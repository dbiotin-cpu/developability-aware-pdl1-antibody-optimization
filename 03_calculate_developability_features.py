from pathlib import Path
import csv
from collections import Counter, defaultdict

from Bio.SeqUtils.ProtParam import ProteinAnalysis


INPUT_FASTA = Path("data/input/final_pdl1_candidate.fasta")
INPUT_CDRS = Path("data/intermediate/cdr_annotation.csv")
INPUT_LIABILITIES = Path("data/intermediate/liability_report.csv")

OUT_FEATURES = Path("data/intermediate/developability_features.csv")
OUT_REGION_FEATURES = Path("data/intermediate/region_developability_features.csv")

PRIMARY_SCHEME = "imgt"

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC_AA = set("AILMFWVY")
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
            "Replace placeholders such as X before calculating developability features."
        )


def fraction(sequence, aa_set):
    if not sequence:
        return 0.0

    return sum(1 for aa in sequence if aa in aa_set) / len(sequence)


def count_residues(sequence, aa_set):
    return sum(1 for aa in sequence if aa in aa_set)


def calculate_net_charge(sequence, ph=7.4):
    """
    Approximate net charge using common amino-acid pKa values.

    Positive groups:
    N-terminus, K, R, H

    Negative groups:
    C-terminus, D, E, C, Y
    """
    counts = Counter(sequence)

    positive = 0.0
    negative = 0.0

    # N-terminus
    positive += 1 / (1 + 10 ** (ph - PKA["Nterm"]))

    # C-terminus
    negative += 1 / (1 + 10 ** (PKA["Cterm"] - ph))

    # Side-chain positive groups
    for aa in ["K", "R", "H"]:
        positive += counts[aa] * (1 / (1 + 10 ** (ph - PKA[aa])))

    # Side-chain negative groups
    for aa in ["D", "E", "C", "Y"]:
        negative += counts[aa] * (1 / (1 + 10 ** (PKA[aa] - ph)))

    return positive - negative


def get_regions_for_record(cdr_rows, record_name, scheme=PRIMARY_SCHEME):
    rows = [
        row for row in cdr_rows
        if row["record"] == record_name and row["scheme"] == scheme
    ]

    region_to_seq = {row["region"]: row["sequence"] for row in rows}

    return {
        "FR1": region_to_seq.get("FR1", ""),
        "CDR1": region_to_seq.get("CDR1", ""),
        "FR2": region_to_seq.get("FR2", ""),
        "CDR2": region_to_seq.get("CDR2", ""),
        "FR3": region_to_seq.get("FR3", ""),
        "CDR3": region_to_seq.get("CDR3", ""),
        "FR4": region_to_seq.get("FR4", ""),
    }


def count_liabilities_for_record(liability_rows, record_name):
    rows = [row for row in liability_rows if row["record"] == record_name]

    severity_counts = Counter(row["severity"] for row in rows)
    type_counts = Counter(row["liability_type"] for row in rows)

    cdr_rows = [row for row in rows if row["in_cdr"] == "True"]
    cdr_type_counts = Counter(row["liability_type"] for row in cdr_rows)

    return {
        "total_liability_flags": len(rows),
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
    """
    Simple transparent score for portfolio triage.

    Higher score = higher predicted developability risk.

    This is not a validated therapeutic developability model.
    It is a lightweight ranking heuristic for comparing candidate variants.
    """
    score = 0.0

    score += 3.0 * int(features["high_liability_flags"])
    score += 2.0 * int(features["medium_liability_flags"])
    score += 1.0 * int(features["low_liability_flags"])

    score += 2.0 * int(features["cdr_liability_flags"])

    # pI penalty: very high or very low pI can complicate developability/purification.
    pI = float(features["pI"])
    if pI >= 9.5 or pI <= 5.5:
        score += 2.0
    elif pI >= 9.0 or pI <= 6.0:
        score += 1.0

    # Charge penalty.
    net_charge = abs(float(features["net_charge_pH_7_4"]))
    if net_charge >= 8:
        score += 2.0
    elif net_charge >= 5:
        score += 1.0

    # CDR hydrophobicity penalty.
    cdr_hydro = float(features["cdr_hydrophobic_fraction"])
    if cdr_hydro >= 0.45:
        score += 3.0
    elif cdr_hydro >= 0.35:
        score += 1.5

    # CDR aromatic penalty.
    cdr_aromatic = float(features["cdr_aromatic_fraction"])
    if cdr_aromatic >= 0.35:
        score += 2.0
    elif cdr_aromatic >= 0.25:
        score += 1.0

    return round(score, 3)


def calculate_sequence_features(record_name, sequence, regions, liability_rows):
    validate_sequence(record_name, sequence)

    analysis = ProteinAnalysis(sequence)

    cdr1 = regions.get("CDR1", "")
    cdr2 = regions.get("CDR2", "")
    cdr3 = regions.get("CDR3", "")
    all_cdrs = cdr1 + cdr2 + cdr3

    fr_sequence = (
        regions.get("FR1", "")
        + regions.get("FR2", "")
        + regions.get("FR3", "")
        + regions.get("FR4", "")
    )

    liability_features = count_liabilities_for_record(liability_rows, record_name)

    features = {
        "record": record_name,
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

    features.update(liability_features)

    features["simple_developability_risk_score"] = calculate_simple_risk_score(features)

    return features


def calculate_region_features(record_name, regions, liability_rows):
    rows = []

    for region, region_seq in regions.items():
        if not region_seq:
            continue

        region_liability_rows = [
            row for row in liability_rows
            if row["record"] == record_name and region in row["region"].split(";")
        ]

        severity_counts = Counter(row["severity"] for row in region_liability_rows)

        rows.append(
            {
                "record": record_name,
                "scheme": PRIMARY_SCHEME,
                "region": region,
                "sequence": region_seq,
                "length": len(region_seq),
                "hydrophobic_fraction": round(fraction(region_seq, HYDROPHOBIC_AA), 3),
                "aromatic_fraction": round(fraction(region_seq, AROMATIC_AA), 3),
                "basic_fraction": round(fraction(region_seq, BASIC_AA), 3),
                "acidic_fraction": round(fraction(region_seq, ACIDIC_AA), 3),
                "polar_fraction": round(fraction(region_seq, POLAR_AA), 3),
                "liability_count": len(region_liability_rows),
                "high_liability_count": severity_counts.get("high", 0),
                "medium_liability_count": severity_counts.get("medium", 0),
                "low_liability_count": severity_counts.get("low", 0),
                "review_liability_count": severity_counts.get("review", 0),
            }
        )

    return rows


def main():
    records = read_fasta(INPUT_FASTA)
    cdr_rows = read_csv(INPUT_CDRS)
    liability_rows = read_csv(INPUT_LIABILITIES)

    all_feature_rows = []
    all_region_rows = []

    for record_name, sequence in records:
        print(f"\nCalculating developability features for: {record_name}")
        print(f"Sequence length: {len(sequence)} aa")

        regions = get_regions_for_record(
            cdr_rows=cdr_rows,
            record_name=record_name,
            scheme=PRIMARY_SCHEME,
        )

        cdr_seq = regions["CDR1"] + regions["CDR2"] + regions["CDR3"]

        print(f"CDR1: {regions['CDR1']}")
        print(f"CDR2: {regions['CDR2']}")
        print(f"CDR3: {regions['CDR3']}")
        print(f"Total CDR length: {len(cdr_seq)} aa")

        features = calculate_sequence_features(
            record_name=record_name,
            sequence=sequence,
            regions=regions,
            liability_rows=liability_rows,
        )

        region_rows = calculate_region_features(
            record_name=record_name,
            regions=regions,
            liability_rows=liability_rows,
        )

        all_feature_rows.append(features)
        all_region_rows.extend(region_rows)

        print(f"pI: {features['pI']}")
        print(f"Net charge at pH 7.4: {features['net_charge_pH_7_4']}")
        print(f"CDR hydrophobic fraction: {features['cdr_hydrophobic_fraction']}")
        print(f"CDR aromatic fraction: {features['cdr_aromatic_fraction']}")
        print(f"Total liability flags: {features['total_liability_flags']}")
        print(f"CDR liability flags: {features['cdr_liability_flags']}")
        print(f"Simple risk score: {features['simple_developability_risk_score']}")

    feature_fields = [
        "record",
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

    region_fields = [
        "record",
        "scheme",
        "region",
        "sequence",
        "length",
        "hydrophobic_fraction",
        "aromatic_fraction",
        "basic_fraction",
        "acidic_fraction",
        "polar_fraction",
        "liability_count",
        "high_liability_count",
        "medium_liability_count",
        "low_liability_count",
        "review_liability_count",
    ]

    write_csv(OUT_FEATURES, all_feature_rows, feature_fields)
    write_csv(OUT_REGION_FEATURES, all_region_rows, region_fields)

    print("\nSaved:")
    print(f"  {OUT_FEATURES}")
    print(f"  {OUT_REGION_FEATURES}")


if __name__ == "__main__":
    main()