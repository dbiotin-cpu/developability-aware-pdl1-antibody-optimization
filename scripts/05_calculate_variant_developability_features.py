#!/usr/bin/env python

from pathlib import Path
import csv
import re
from collections import Counter, defaultdict

from Bio.SeqUtils.ProtParam import ProteinAnalysis


INPUT_VARIANTS_FASTA = Path("data/intermediate/proposed_variants.fasta")
INPUT_CDRS = Path("data/intermediate/cdr_annotation.csv")
INPUT_VARIANT_METADATA = Path("data/intermediate/proposed_variants.csv")

OUT_VARIANT_LIABILITIES = Path("data/intermediate/variant_liability_report.csv")
OUT_VARIANT_FEATURES = Path("data/intermediate/variant_developability_features.csv")
OUT_VARIANT_RANKING = Path("data/intermediate/variant_ranking_pre_structure.csv")

PRIMARY_SCHEME = "imgt"

REGION_ORDER = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC_AA = set("AILMFWVY")
PATCH_HYDROPHOBIC_AA = set("VILFWYM")
AROMATIC_AA = set("FWY")
BASIC_AA = set("KRH")
STRONGLY_BASIC_AA = set("KR")
ACIDIC_AA = set("DE")
POLAR_AA = set("STNQ")

PKA = {
    "Cterm": 2.34,
    "Nterm": 9.69,
    "D": 3.86,
    "E": 4.25,
    "C": 8.33,
    "Y": 10.07,
    "H": 6.00,
    "K": 10.50,
    "R": 12.40,
}


def read_fasta(path):
    records = []
    name = None
    seq_lines = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_lines).replace(" ", "").upper()))
                name = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)

        if name is not None:
            records.append((name, "".join(seq_lines).replace(" ", "").upper()))

    return records


def read_csv(path):
    if not path.exists():
        return []

    with open(path, "r") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def validate_sequence(record_name, sequence):
    invalid = sorted(set(sequence) - CANONICAL_AA)

    if invalid:
        raise ValueError(
            f"{record_name} contains non-canonical residues: {invalid}. "
            "Replace placeholders such as X before calculating features."
        )


def fraction(sequence, aa_set):
    if not sequence:
        return 0.0

    return sum(1 for aa in sequence if aa in aa_set) / len(sequence)


def calculate_net_charge(sequence, ph=7.4):
    counts = Counter(sequence)

    positive = 0.0
    negative = 0.0

    positive += 1 / (1 + 10 ** (ph - PKA["Nterm"]))
    negative += 1 / (1 + 10 ** (PKA["Cterm"] - ph))

    for aa in ["K", "R", "H"]:
        positive += counts[aa] * (1 / (1 + 10 ** (ph - PKA[aa])))

    for aa in ["D", "E", "C", "Y"]:
        negative += counts[aa] * (1 / (1 + 10 ** (PKA[aa] - ph)))

    return positive - negative


def infer_parent_record_name(wt_variant_id):
    if wt_variant_id.endswith("_WT"):
        return wt_variant_id[:-3]
    return wt_variant_id


def get_region_lengths_from_cdr_annotation(cdr_rows, parent_record, scheme=PRIMARY_SCHEME):
    rows = [
        row for row in cdr_rows
        if row["record"] == parent_record and row["scheme"] == scheme
    ]

    region_to_seq = {row["region"]: row["sequence"] for row in rows}

    region_lengths = []
    for region in REGION_ORDER:
        seq = region_to_seq.get(region, "")
        region_lengths.append((region, len(seq), seq))

    return region_lengths


def build_position_region_map_from_lengths(region_lengths, sequence_length):
    """
    Proposed variants are same-length point mutants.
    Therefore we can use WT FR/CDR region lengths to map variant positions.
    """
    region_map = {}
    pos = 1

    for region, length, _seq in region_lengths:
        for _ in range(length):
            region_map[pos] = region
            pos += 1

    mapped_len = len(region_map)

    if mapped_len != sequence_length:
        print(
            f"WARNING: Region map length {mapped_len} does not match sequence length {sequence_length}. "
            "Unmapped residues will be labeled unknown."
        )

    return region_map


def region_for_span(start, end, region_map):
    regions = []
    for pos in range(start, end + 1):
        region = region_map.get(pos, "unknown")
        if region not in regions:
            regions.append(region)

    return ";".join(regions)


def is_cdr_region(region_text):
    return any(region.startswith("CDR") for region in region_text.split(";"))


def add_liability(rows, record, liability_type, motif, start, end,
                  region_map, severity, rationale):
    region = region_for_span(start, end, region_map)
    in_cdr = is_cdr_region(region)

    rows.append(
        {
            "variant_id": record,
            "scheme": PRIMARY_SCHEME,
            "liability_type": liability_type,
            "motif": motif,
            "start": start,
            "end": end,
            "region": region,
            "in_cdr": in_cdr,
            "severity": severity,
            "rationale": rationale,
        }
    )


def scan_liabilities(record, sequence, region_map):
    rows = []

    # N-linked glycosylation: N-X-S/T, where X is not Pro.
    for match in re.finditer(r"(?=(N[^P][ST]))", sequence):
        motif = match.group(1)
        start = match.start(1) + 1
        end = start + len(motif) - 1
        region = region_for_span(start, end, region_map)
        severity = "high" if is_cdr_region(region) else "medium"

        add_liability(
            rows,
            record,
            "N-linked glycosylation motif",
            motif,
            start,
            end,
            region_map,
            severity,
            "N-X-S/T motif may introduce unwanted glycosylation risk.",
        )

    # Asn deamidation motifs.
    for match in re.finditer(r"(?=(N[GSTH]))", sequence):
        motif = match.group(1)
        start = match.start(1) + 1
        end = start + len(motif) - 1
        region = region_for_span(start, end, region_map)
        severity = "medium" if is_cdr_region(region) else "low"

        add_liability(
            rows,
            record,
            "Asn deamidation motif",
            motif,
            start,
            end,
            region_map,
            severity,
            "NG, NS, NT, and NH motifs may be deamidation-prone depending on structure.",
        )

    # Asp isomerization motifs.
    for match in re.finditer(r"(?=(D[GST]))", sequence):
        motif = match.group(1)
        start = match.start(1) + 1
        end = start + len(motif) - 1
        region = region_for_span(start, end, region_map)
        severity = "medium" if is_cdr_region(region) else "low"

        add_liability(
            rows,
            record,
            "Asp isomerization motif",
            motif,
            start,
            end,
            region_map,
            severity,
            "DG, DS, and DT motifs may carry Asp isomerization risk.",
        )

    # Oxidation-prone Met/Trp.
    for i, aa in enumerate(sequence, start=1):
        if aa in {"M", "W"}:
            region = region_for_span(i, i, region_map)
            severity = "medium" if is_cdr_region(region) else "low"

            add_liability(
                rows,
                record,
                "Oxidation-prone residue",
                aa,
                i,
                i,
                region_map,
                severity,
                "Met and Trp can be oxidation-prone, especially if solvent-exposed.",
            )

    # Cysteines for review.
    for i, aa in enumerate(sequence, start=1):
        if aa == "C":
            add_liability(
                rows,
                record,
                "Cysteine for disulfide/free-Cys review",
                aa,
                i,
                i,
                region_map,
                "review",
                "Cysteines should be checked for canonical disulfide pairing or potential free-Cys liability.",
            )

    # Polybasic patches.
    basic_window = 5
    basic_threshold = 3

    for i in range(0, len(sequence) - basic_window + 1):
        window = sequence[i:i + basic_window]
        basic_count = sum(1 for aa in window if aa in STRONGLY_BASIC_AA)

        if basic_count >= basic_threshold:
            start = i + 1
            end = i + basic_window
            region = region_for_span(start, end, region_map)
            severity = "medium" if is_cdr_region(region) else "low"

            add_liability(
                rows,
                record,
                "Polybasic patch",
                window,
                start,
                end,
                region_map,
                severity,
                "Local Lys/Arg enrichment may increase nonspecific electrostatic interactions.",
            )

    # Hydrophobic patches.
    hydro_window = 7
    hydro_threshold = 5

    for i in range(0, len(sequence) - hydro_window + 1):
        window = sequence[i:i + hydro_window]
        hydro_count = sum(1 for aa in window if aa in PATCH_HYDROPHOBIC_AA)

        if hydro_count >= hydro_threshold:
            start = i + 1
            end = i + hydro_window
            region = region_for_span(start, end, region_map)
            severity = "high" if is_cdr_region(region) else "medium"

            add_liability(
                rows,
                record,
                "Hydrophobic patch",
                window,
                start,
                end,
                region_map,
                severity,
                "Local hydrophobic enrichment may increase aggregation, polyspecificity, or solubility risk.",
            )

    return rows


def get_region_sequences(sequence, region_map):
    region_to_seq = defaultdict(list)

    for i, aa in enumerate(sequence, start=1):
        region = region_map.get(i, "unknown")
        region_to_seq[region].append(aa)

    return {region: "".join(aas) for region, aas in region_to_seq.items()}


def count_liabilities(liability_rows):
    severity_counts = Counter(row["severity"] for row in liability_rows)
    type_counts = Counter(row["liability_type"] for row in liability_rows)
    cdr_rows = [row for row in liability_rows if str(row["in_cdr"]) == "True"]
    cdr_type_counts = Counter(row["liability_type"] for row in cdr_rows)

    return {
        "total_liability_flags": len(liability_rows),
        "cdr_liability_flags": len(cdr_rows),
        "high_liability_flags": severity_counts.get("high", 0),
        "medium_liability_flags": severity_counts.get("medium", 0),
        "low_liability_flags": severity_counts.get("low", 0),
        "review_liability_flags": severity_counts.get("review", 0),
        "n_glycosylation_flags": type_counts.get("N-linked glycosylation motif", 0),
        "deamidation_flags": type_counts.get("Asn deamidation motif", 0),
        "asp_isomerization_flags": type_counts.get("Asp isomerization motif", 0),
        "oxidation_flags": type_counts.get("Oxidation-prone residue", 0),
        "cysteine_review_flags": type_counts.get("Cysteine for disulfide/free-Cys review", 0),
        "polybasic_patch_flags": type_counts.get("Polybasic patch", 0),
        "hydrophobic_patch_flags": type_counts.get("Hydrophobic patch", 0),
        "cdr_n_glycosylation_flags": cdr_type_counts.get("N-linked glycosylation motif", 0),
        "cdr_deamidation_flags": cdr_type_counts.get("Asn deamidation motif", 0),
        "cdr_oxidation_flags": cdr_type_counts.get("Oxidation-prone residue", 0),
        "cdr_hydrophobic_patch_flags": cdr_type_counts.get("Hydrophobic patch", 0),
    }


def calculate_simple_risk_score(features):
    score = 0.0

    score += 3.0 * int(features["high_liability_flags"])
    score += 2.0 * int(features["medium_liability_flags"])
    score += 1.0 * int(features["low_liability_flags"])
    score += 0.5 * int(features["review_liability_flags"])

    score += 2.0 * int(features["cdr_liability_flags"])

    pI = float(features["pI"])
    if pI >= 9.5 or pI <= 5.5:
        score += 2.0
    elif pI >= 9.0 or pI <= 6.0:
        score += 1.0

    net_charge_abs = abs(float(features["net_charge_pH_7_4"]))
    if net_charge_abs >= 8:
        score += 2.0
    elif net_charge_abs >= 5:
        score += 1.0

    cdr_hydro = float(features["cdr_hydrophobic_fraction"])
    if cdr_hydro >= 0.45:
        score += 3.0
    elif cdr_hydro >= 0.35:
        score += 1.5

    cdr_aromatic = float(features["cdr_aromatic_fraction"])
    if cdr_aromatic >= 0.35:
        score += 2.0
    elif cdr_aromatic >= 0.25:
        score += 1.0

    return round(score, 3)


def calculate_variant_features(variant_id, sequence, region_map, liability_rows):
    validate_sequence(variant_id, sequence)

    analysis = ProteinAnalysis(sequence)

    region_sequences = get_region_sequences(sequence, region_map)

    cdr1 = region_sequences.get("CDR1", "")
    cdr2 = region_sequences.get("CDR2", "")
    cdr3 = region_sequences.get("CDR3", "")
    all_cdrs = cdr1 + cdr2 + cdr3

    fr_sequence = (
        region_sequences.get("FR1", "")
        + region_sequences.get("FR2", "")
        + region_sequences.get("FR3", "")
        + region_sequences.get("FR4", "")
    )

    liability_counts = count_liabilities(liability_rows)

    features = {
        "variant_id": variant_id,
        "scheme": PRIMARY_SCHEME,
        "sequence_length": len(sequence),
        "molecular_weight_Da": round(analysis.molecular_weight(), 2),
        "pI": round(analysis.isoelectric_point(), 3),
        "net_charge_pH_7_4": round(calculate_net_charge(sequence, ph=7.4), 3),
        "whole_sequence_hydrophobic_fraction": round(fraction(sequence, HYDROPHOBIC_AA), 3),
        "whole_sequence_aromatic_fraction": round(fraction(sequence, AROMATIC_AA), 3),
        "whole_sequence_basic_fraction": round(fraction(sequence, BASIC_AA), 3),
        "whole_sequence_strongly_basic_fraction": round(fraction(sequence, STRONGLY_BASIC_AA), 3),
        "whole_sequence_acidic_fraction": round(fraction(sequence, ACIDIC_AA), 3),
        "whole_sequence_polar_fraction": round(fraction(sequence, POLAR_AA), 3),
        "cdr_total_length": len(all_cdrs),
        "cdr1_length": len(cdr1),
        "cdr2_length": len(cdr2),
        "cdr3_length": len(cdr3),
        "cdr_hydrophobic_fraction": round(fraction(all_cdrs, HYDROPHOBIC_AA), 3),
        "cdr_aromatic_fraction": round(fraction(all_cdrs, AROMATIC_AA), 3),
        "cdr_basic_fraction": round(fraction(all_cdrs, BASIC_AA), 3),
        "cdr_acidic_fraction": round(fraction(all_cdrs, ACIDIC_AA), 3),
        "framework_hydrophobic_fraction": round(fraction(fr_sequence, HYDROPHOBIC_AA), 3),
        "framework_aromatic_fraction": round(fraction(fr_sequence, AROMATIC_AA), 3),
        "cysteine_count": sequence.count("C"),
        "methionine_count": sequence.count("M"),
        "tryptophan_count": sequence.count("W"),
    }

    features.update(liability_counts)
    features["simple_developability_risk_score"] = calculate_simple_risk_score(features)

    return features


def load_variant_metadata(path):
    rows = read_csv(path)
    metadata = {}

    for row in rows:
        metadata[row["variant_id"]] = row

    return metadata


def calculate_pre_structure_ranking(feature_rows, metadata):
    wt_rows = [row for row in feature_rows if row["variant_id"].endswith("_WT")]

    if not wt_rows:
        wt = feature_rows[0]
        print(f"WARNING: No _WT record found. Using {wt['variant_id']} as WT.")
    else:
        wt = wt_rows[0]

    ranking_rows = []

    wt_score = float(wt["simple_developability_risk_score"])
    wt_liabilities = int(wt["total_liability_flags"])
    wt_cdr_liabilities = int(wt["cdr_liability_flags"])
    wt_hydro = float(wt["cdr_hydrophobic_fraction"])
    wt_aromatic = float(wt["cdr_aromatic_fraction"])
    wt_pI = float(wt["pI"])
    wt_charge = float(wt["net_charge_pH_7_4"])

    for row in feature_rows:
        variant_id = row["variant_id"]

        variant_score = float(row["simple_developability_risk_score"])
        variant_liabilities = int(row["total_liability_flags"])
        variant_cdr_liabilities = int(row["cdr_liability_flags"])
        variant_hydro = float(row["cdr_hydrophobic_fraction"])
        variant_aromatic = float(row["cdr_aromatic_fraction"])
        variant_pI = float(row["pI"])
        variant_charge = float(row["net_charge_pH_7_4"])

        delta_risk_score = round(variant_score - wt_score, 3)
        liability_reduction = wt_liabilities - variant_liabilities
        cdr_liability_reduction = wt_cdr_liabilities - variant_cdr_liabilities
        delta_cdr_hydrophobic_fraction = round(variant_hydro - wt_hydro, 3)
        delta_cdr_aromatic_fraction = round(variant_aromatic - wt_aromatic, 3)
        delta_pI = round(variant_pI - wt_pI, 3)
        delta_net_charge = round(variant_charge - wt_charge, 3)

        # Higher is better.
        optimization_score = 0.0

        # The biggest signal: risk score decrease.
        optimization_score += -1.0 * delta_risk_score

        # Direct liability reduction.
        optimization_score += 2.0 * liability_reduction
        optimization_score += 2.0 * cdr_liability_reduction

        # Penalize variants that create more liabilities.
        if liability_reduction < 0:
            optimization_score -= 3.0

        if cdr_liability_reduction < 0:
            optimization_score -= 4.0

        # Reward reduced CDR hydrophobicity/aromaticity.
        if delta_cdr_hydrophobic_fraction < 0:
            optimization_score += 1.0
        elif delta_cdr_hydrophobic_fraction > 0:
            optimization_score -= 1.0

        if delta_cdr_aromatic_fraction < 0:
            optimization_score += 0.5
        elif delta_cdr_aromatic_fraction > 0:
            optimization_score -= 0.5

        # Penalize extreme pI.
        if variant_pI >= 9.5 or variant_pI <= 5.5:
            optimization_score -= 2.0

        # Penalize very large charge magnitude.
        if abs(variant_charge) >= 8:
            optimization_score -= 2.0

        # WT should not be ranked as an optimized candidate.
        if variant_id.endswith("_WT"):
            decision = "Reference"
            optimization_score = 0.0
        else:
            if optimization_score >= 6:
                decision = "Advance to structural re-evaluation"
            elif optimization_score >= 2:
                decision = "Consider for structural re-evaluation"
            elif optimization_score > 0:
                decision = "Low-priority backup"
            else:
                decision = "Do not prioritize before structure review"

        meta = metadata.get(variant_id, {})

        ranking_rows.append(
            {
                "variant_id": variant_id,
                "variant_type": meta.get("variant_type", "WT" if variant_id.endswith("_WT") else ""),
                "mutations": meta.get("mutations", ""),
                "num_mutations": meta.get("num_mutations", "0" if variant_id.endswith("_WT") else ""),
                "risk_score": variant_score,
                "delta_risk_score_vs_WT": delta_risk_score,
                "total_liability_flags": variant_liabilities,
                "liability_reduction_vs_WT": liability_reduction,
                "cdr_liability_flags": variant_cdr_liabilities,
                "cdr_liability_reduction_vs_WT": cdr_liability_reduction,
                "pI": variant_pI,
                "delta_pI_vs_WT": delta_pI,
                "net_charge_pH_7_4": variant_charge,
                "delta_net_charge_vs_WT": delta_net_charge,
                "cdr_hydrophobic_fraction": variant_hydro,
                "delta_cdr_hydrophobic_fraction_vs_WT": delta_cdr_hydrophobic_fraction,
                "cdr_aromatic_fraction": variant_aromatic,
                "delta_cdr_aromatic_fraction_vs_WT": delta_cdr_aromatic_fraction,
                "optimization_score_pre_structure": round(optimization_score, 3),
                "pre_structure_decision": decision,
            }
        )

    ranking_rows = sorted(
        ranking_rows,
        key=lambda x: (
            x["pre_structure_decision"] == "Reference",
            -float(x["optimization_score_pre_structure"]),
            float(x["risk_score"]),
        ),
    )

    return ranking_rows


def main():
    records = read_fasta(INPUT_VARIANTS_FASTA)

    if not records:
        raise ValueError(f"No sequences found in {INPUT_VARIANTS_FASTA}")

    cdr_rows = read_csv(INPUT_CDRS)
    metadata = load_variant_metadata(INPUT_VARIANT_METADATA)

    wt_candidates = [(name, seq) for name, seq in records if name.endswith("_WT")]

    if wt_candidates:
        wt_name, wt_seq = wt_candidates[0]
    else:
        wt_name, wt_seq = records[0]
        print(f"WARNING: No _WT sequence found. Using first FASTA record as WT: {wt_name}")

    parent_record = infer_parent_record_name(wt_name)

    print(f"WT/reference sequence: {wt_name}")
    print(f"Parent record for CDR annotation: {parent_record}")

    region_lengths = get_region_lengths_from_cdr_annotation(
        cdr_rows=cdr_rows,
        parent_record=parent_record,
        scheme=PRIMARY_SCHEME,
    )

    print("\nRegion lengths from WT annotation:")
    for region, length, seq in region_lengths:
        print(f"  {region}: {length} aa")

    all_liability_rows = []
    all_feature_rows = []

    for variant_id, sequence in records:
        validate_sequence(variant_id, sequence)

        print(f"\nProcessing variant: {variant_id}")
        print(f"Sequence length: {len(sequence)} aa")

        region_map = build_position_region_map_from_lengths(
            region_lengths=region_lengths,
            sequence_length=len(sequence),
        )

        liability_rows = scan_liabilities(
            record=variant_id,
            sequence=sequence,
            region_map=region_map,
        )

        all_liability_rows.extend(liability_rows)

        features = calculate_variant_features(
            variant_id=variant_id,
            sequence=sequence,
            region_map=region_map,
            liability_rows=liability_rows,
        )

        all_feature_rows.append(features)

        print(f"  Risk score: {features['simple_developability_risk_score']}")
        print(f"  Total liabilities: {features['total_liability_flags']}")
        print(f"  CDR liabilities: {features['cdr_liability_flags']}")
        print(f"  pI: {features['pI']}")
        print(f"  Net charge pH 7.4: {features['net_charge_pH_7_4']}")

    liability_fields = [
        "variant_id",
        "scheme",
        "liability_type",
        "motif",
        "start",
        "end",
        "region",
        "in_cdr",
        "severity",
        "rationale",
    ]

    feature_fields = [
        "variant_id",
        "scheme",
        "sequence_length",
        "molecular_weight_Da",
        "pI",
        "net_charge_pH_7_4",
        "whole_sequence_hydrophobic_fraction",
        "whole_sequence_aromatic_fraction",
        "whole_sequence_basic_fraction",
        "whole_sequence_strongly_basic_fraction",
        "whole_sequence_acidic_fraction",
        "whole_sequence_polar_fraction",
        "cdr_total_length",
        "cdr1_length",
        "cdr2_length",
        "cdr3_length",
        "cdr_hydrophobic_fraction",
        "cdr_aromatic_fraction",
        "cdr_basic_fraction",
        "cdr_acidic_fraction",
        "framework_hydrophobic_fraction",
        "framework_aromatic_fraction",
        "cysteine_count",
        "methionine_count",
        "tryptophan_count",
        "total_liability_flags",
        "cdr_liability_flags",
        "high_liability_flags",
        "medium_liability_flags",
        "low_liability_flags",
        "review_liability_flags",
        "n_glycosylation_flags",
        "deamidation_flags",
        "asp_isomerization_flags",
        "oxidation_flags",
        "cysteine_review_flags",
        "polybasic_patch_flags",
        "hydrophobic_patch_flags",
        "cdr_n_glycosylation_flags",
        "cdr_deamidation_flags",
        "cdr_oxidation_flags",
        "cdr_hydrophobic_patch_flags",
        "simple_developability_risk_score",
    ]

    ranking_fields = [
        "variant_id",
        "variant_type",
        "mutations",
        "num_mutations",
        "risk_score",
        "delta_risk_score_vs_WT",
        "total_liability_flags",
        "liability_reduction_vs_WT",
        "cdr_liability_flags",
        "cdr_liability_reduction_vs_WT",
        "pI",
        "delta_pI_vs_WT",
        "net_charge_pH_7_4",
        "delta_net_charge_vs_WT",
        "cdr_hydrophobic_fraction",
        "delta_cdr_hydrophobic_fraction_vs_WT",
        "cdr_aromatic_fraction",
        "delta_cdr_aromatic_fraction_vs_WT",
        "optimization_score_pre_structure",
        "pre_structure_decision",
    ]

    ranking_rows = calculate_pre_structure_ranking(
        feature_rows=all_feature_rows,
        metadata=metadata,
    )

    write_csv(OUT_VARIANT_LIABILITIES, all_liability_rows, liability_fields)
    write_csv(OUT_VARIANT_FEATURES, all_feature_rows, feature_fields)
    write_csv(OUT_VARIANT_RANKING, ranking_rows, ranking_fields)

    print("\nSaved:")
    print(f"  {OUT_VARIANT_LIABILITIES}")
    print(f"  {OUT_VARIANT_FEATURES}")
    print(f"  {OUT_VARIANT_RANKING}")

    print("\nTop pre-structure variants:")
    for row in ranking_rows[:10]:
        print(
            f"  {row['variant_id']} | "
            f"score={row['optimization_score_pre_structure']} | "
            f"risk={row['risk_score']} | "
            f"decision={row['pre_structure_decision']}"
        )


if __name__ == "__main__":
    main()
