#!/usr/bin/env python

from pathlib import Path
import csv
import re
from collections import Counter, defaultdict


INPUT_FASTA = Path("data/input/final_pdl1_candidate.fasta")
INPUT_CDRS = Path("data/intermediate/cdr_annotation.csv")

OUT_REPORT = Path("data/intermediate/liability_report.csv")
OUT_SUMMARY = Path("data/intermediate/liability_summary.csv")

PRIMARY_SCHEME = "imgt"

HYDROPHOBIC_AA = set("VILFWYM")
BASIC_AA = set("KR")


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


def read_cdr_annotation(path):
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def build_region_map(record_name, sequence, cdr_rows, scheme=PRIMARY_SCHEME):
    """
    Builds a 1-indexed residue-to-region map from cdr_annotation.csv.

    Expected regions:
    FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4
    """
    region_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

    rows = [
        row for row in cdr_rows
        if row["record"] == record_name and row["scheme"] == scheme
    ]

    region_to_seq = {row["region"]: row["sequence"] for row in rows}

    ordered_segments = []
    for region in region_order:
        seq = region_to_seq.get(region, "")
        if seq:
            ordered_segments.append((region, seq))

    reconstructed = "".join(seq for _, seq in ordered_segments)

    region_map = {}

    # Best case: concatenated FR/CDR sequence exactly matches original sequence.
    if reconstructed == sequence:
        pos = 1
        for region, segment in ordered_segments:
            for aa in segment:
                region_map[pos] = region
                pos += 1
        return region_map

    # Fallback: locate each segment sequentially in the full sequence.
    cursor = 0
    for region, segment in ordered_segments:
        found = sequence.find(segment, cursor)
        if found == -1:
            continue

        start = found + 1
        end = found + len(segment)

        for pos in range(start, end + 1):
            region_map[pos] = region

        cursor = found + len(segment)

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


def add_liability(rows, record, scheme, liability_type, motif, start, end,
                  sequence, region_map, severity, rationale, suggested_action):
    region = region_for_span(start, end, region_map)
    in_cdr = is_cdr_region(region)

    rows.append(
        {
            "record": record,
            "scheme": scheme,
            "liability_type": liability_type,
            "motif": motif,
            "start": start,
            "end": end,
            "region": region,
            "in_cdr": in_cdr,
            "severity": severity,
            "rationale": rationale,
            "suggested_action": suggested_action,
        }
    )


def scan_regex_liabilities(record, sequence, region_map, scheme=PRIMARY_SCHEME):
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
            scheme,
            "N-linked glycosylation motif",
            motif,
            start,
            end,
            sequence,
            region_map,
            severity,
            "N-X-S/T motifs can introduce unwanted glycosylation risk in expressed antibodies.",
            "Consider N->Q, S/T->A, or conservative local redesign if outside key binding contacts.",
        )

    # Deamidation-prone Asn motifs.
    for match in re.finditer(r"(?=(N[GSTH]))", sequence):
        motif = match.group(1)
        start = match.start(1) + 1
        end = start + len(motif) - 1
        region = region_for_span(start, end, region_map)
        severity = "medium" if is_cdr_region(region) else "low"

        add_liability(
            rows,
            record,
            scheme,
            "Asn deamidation motif",
            motif,
            start,
            end,
            sequence,
            region_map,
            severity,
            "NG, NS, NT, and NH motifs can be deamidation-prone depending on structure and solvent exposure.",
            "Review structural exposure; consider N->Q or local conservative substitutions if not binding-critical.",
        )

    # Asp isomerization-prone motifs.
    for match in re.finditer(r"(?=(D[GST]))", sequence):
        motif = match.group(1)
        start = match.start(1) + 1
        end = start + len(motif) - 1
        region = region_for_span(start, end, region_map)
        severity = "medium" if is_cdr_region(region) else "low"

        add_liability(
            rows,
            record,
            scheme,
            "Asp isomerization motif",
            motif,
            start,
            end,
            sequence,
            region_map,
            severity,
            "DG, DS, and DT motifs may carry Asp isomerization risk depending on local flexibility.",
            "Review structural context; avoid changing binding-critical residues without re-evaluation.",
        )

    # Oxidation-prone Met and Trp.
    for i, aa in enumerate(sequence, start=1):
        if aa in {"M", "W"}:
            region = region_for_span(i, i, region_map)
            severity = "medium" if is_cdr_region(region) else "low"

            add_liability(
                rows,
                record,
                scheme,
                "Oxidation-prone residue",
                aa,
                i,
                i,
                sequence,
                region_map,
                severity,
                "Met and Trp can be oxidation-prone, especially when solvent-exposed or located in CDRs.",
                "Review solvent exposure; consider conservative substitution only if not important for binding.",
            )

    # Cysteines for review.
    # This does not prove a cysteine is free; it flags cysteines for structural/disulfide review.
    for i, aa in enumerate(sequence, start=1):
        if aa == "C":
            add_liability(
                rows,
                record,
                scheme,
                "Cysteine for disulfide/free-Cys review",
                aa,
                i,
                i,
                sequence,
                region_map,
                "review",
                "Cysteines should be checked for expected antibody disulfide pairing or potential free-Cys liability.",
                "Confirm whether this cysteine is canonical/disulfide-paired; avoid unpaired cysteines in final candidates.",
            )

    return rows


def scan_patch_liabilities(record, sequence, region_map, scheme=PRIMARY_SCHEME):
    rows = []

    # Polybasic patches: sliding window with many Lys/Arg residues.
    basic_window = 5
    basic_threshold = 3

    for i in range(0, len(sequence) - basic_window + 1):
        window = sequence[i:i + basic_window]
        basic_count = sum(1 for aa in window if aa in BASIC_AA)

        if basic_count >= basic_threshold:
            start = i + 1
            end = i + basic_window
            region = region_for_span(start, end, region_map)
            severity = "medium" if is_cdr_region(region) else "low"

            add_liability(
                rows,
                record,
                scheme,
                "Polybasic patch",
                window,
                start,
                end,
                sequence,
                region_map,
                severity,
                "Local Lys/Arg enrichment may increase nonspecific electrostatic interactions.",
                "Consider reducing exposed basic clustering if outside the binding interface.",
            )

    # Hydrophobic patches: sliding window with many hydrophobic residues.
    hydro_window = 7
    hydro_threshold = 5

    for i in range(0, len(sequence) - hydro_window + 1):
        window = sequence[i:i + hydro_window]
        hydro_count = sum(1 for aa in window if aa in HYDROPHOBIC_AA)

        if hydro_count >= hydro_threshold:
            start = i + 1
            end = i + hydro_window
            region = region_for_span(start, end, region_map)
            severity = "high" if is_cdr_region(region) else "medium"

            add_liability(
                rows,
                record,
                scheme,
                "Hydrophobic patch",
                window,
                start,
                end,
                sequence,
                region_map,
                severity,
                "Local hydrophobic enrichment may increase aggregation, polyspecificity, or poor solubility risk.",
                "Review surface exposure; consider conservative polar substitutions if not involved in binding.",
            )

    return rows


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize_liabilities(rows):
    summary = []

    count_by_type = Counter(row["liability_type"] for row in rows)
    count_by_type_cdr = Counter(
        row["liability_type"] for row in rows if str(row["in_cdr"]) == "True"
    )
    count_by_severity = Counter(row["severity"] for row in rows)

    for liability_type, count in sorted(count_by_type.items()):
        summary.append(
            {
                "category": "liability_type",
                "name": liability_type,
                "total_count": count,
                "cdr_count": count_by_type_cdr.get(liability_type, 0),
            }
        )

    for severity, count in sorted(count_by_severity.items()):
        summary.append(
            {
                "category": "severity",
                "name": severity,
                "total_count": count,
                "cdr_count": "",
            }
        )

    return summary


def main():
    records = read_fasta(INPUT_FASTA)
    cdr_rows = read_cdr_annotation(INPUT_CDRS)

    all_rows = []

    for record_name, sequence in records:
        print(f"\nScanning liabilities for: {record_name}")
        print(f"Sequence length: {len(sequence)} aa")

        region_map = build_region_map(
            record_name=record_name,
            sequence=sequence,
            cdr_rows=cdr_rows,
            scheme=PRIMARY_SCHEME,
        )

        mapped_positions = len(region_map)
        print(f"Mapped positions using {PRIMARY_SCHEME.upper()}: {mapped_positions}/{len(sequence)}")

        regex_rows = scan_regex_liabilities(
            record=record_name,
            sequence=sequence,
            region_map=region_map,
            scheme=PRIMARY_SCHEME,
        )

        patch_rows = scan_patch_liabilities(
            record=record_name,
            sequence=sequence,
            region_map=region_map,
            scheme=PRIMARY_SCHEME,
        )

        rows = regex_rows + patch_rows
        all_rows.extend(rows)

        print(f"Total liability flags: {len(rows)}")

        type_counts = Counter(row["liability_type"] for row in rows)
        for liability_type, count in type_counts.items():
            print(f"  {liability_type}: {count}")

    report_fields = [
        "record",
        "scheme",
        "liability_type",
        "motif",
        "start",
        "end",
        "region",
        "in_cdr",
        "severity",
        "rationale",
        "suggested_action",
    ]

    summary_fields = [
        "category",
        "name",
        "total_count",
        "cdr_count",
    ]

    summary_rows = summarize_liabilities(all_rows)

    write_csv(OUT_REPORT, all_rows, report_fields)
    write_csv(OUT_SUMMARY, summary_rows, summary_fields)

    print("\nSaved:")
    print(f"  {OUT_REPORT}")
    print(f"  {OUT_SUMMARY}")


if __name__ == "__main__":
    main()