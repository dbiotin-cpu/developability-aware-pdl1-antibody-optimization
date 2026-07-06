#!/usr/bin/env python

from pathlib import Path
import csv
from abnumber import Chain


INPUT_FASTA = Path("data/input/final_pdl1_candidate.fasta")
OUT_NUMBERING = Path("data/intermediate/numbered_sequences.csv")
OUT_CDRS = Path("data/intermediate/cdr_annotation.csv")


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


def safe_getattr(obj, attr):
    return getattr(obj, attr, "")


def get_region_from_position(pos):
    """
    AbNumber Position objects may expose region information depending on version.
    This fallback keeps the numbering table usable even if region is unavailable.
    """
    for attr in ["region", "cdr_definition_position", "region_name"]:
        if hasattr(pos, attr):
            value = getattr(pos, attr)
            if value:
                return str(value)

    if hasattr(pos, "get_region"):
        try:
            return str(pos.get_region())
        except Exception:
            pass

    return ""


def annotate_chain(record_name, seq, scheme):
    chain = Chain(seq, scheme=scheme)

    numbering_rows = []
    for pos, aa in chain:
        numbering_rows.append(
            {
                "record": record_name,
                "scheme": scheme,
                "position": str(pos),
                "aa": aa,
                "region_from_position": get_region_from_position(pos),
            }
        )

    cdr_rows = [
        {
            "record": record_name,
            "scheme": scheme,
            "region": "FR1",
            "sequence": safe_getattr(chain, "fr1_seq"),
            "length": len(safe_getattr(chain, "fr1_seq")),
        },
        {
            "record": record_name,
            "scheme": scheme,
            "region": "CDR1",
            "sequence": safe_getattr(chain, "cdr1_seq"),
            "length": len(safe_getattr(chain, "cdr1_seq")),
        },
        {
            "record": record_name,
            "scheme": scheme,
            "region": "FR2",
            "sequence": safe_getattr(chain, "fr2_seq"),
            "length": len(safe_getattr(chain, "fr2_seq")),
        },
        {
            "record": record_name,
            "scheme": scheme,
            "region": "CDR2",
            "sequence": safe_getattr(chain, "cdr2_seq"),
            "length": len(safe_getattr(chain, "cdr2_seq")),
        },
        {
            "record": record_name,
            "scheme": scheme,
            "region": "FR3",
            "sequence": safe_getattr(chain, "fr3_seq"),
            "length": len(safe_getattr(chain, "fr3_seq")),
        },
        {
            "record": record_name,
            "scheme": scheme,
            "region": "CDR3",
            "sequence": safe_getattr(chain, "cdr3_seq"),
            "length": len(safe_getattr(chain, "cdr3_seq")),
        },
        {
            "record": record_name,
            "scheme": scheme,
            "region": "FR4",
            "sequence": safe_getattr(chain, "fr4_seq"),
            "length": len(safe_getattr(chain, "fr4_seq")),
        },
    ]

    return numbering_rows, cdr_rows


def main():
    records = read_fasta(INPUT_FASTA)

    if not records:
        raise ValueError(f"No FASTA records found in {INPUT_FASTA}")

    all_numbering_rows = []
    all_cdr_rows = []

    schemes = ["imgt", "kabat", "chothia"]

    for record_name, seq in records:
        print(f"\nProcessing: {record_name}")
        print(f"Sequence length: {len(seq)} aa")

        for scheme in schemes:
            try:
                numbering_rows, cdr_rows = annotate_chain(record_name, seq, scheme)
                all_numbering_rows.extend(numbering_rows)
                all_cdr_rows.extend(cdr_rows)

                print(f"  {scheme}: success")
                for row in cdr_rows:
                    if row["region"].startswith("CDR"):
                        print(
                            f"    {row['region']}: {row['sequence']} "
                            f"({row['length']} aa)"
                        )

            except Exception as e:
                print(f"  {scheme}: failed")
                print(f"    Error: {e}")

    OUT_NUMBERING.parent.mkdir(parents=True, exist_ok=True)
    OUT_CDRS.parent.mkdir(parents=True, exist_ok=True)

    with open(OUT_NUMBERING, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "record",
                "scheme",
                "position",
                "aa",
                "region_from_position",
            ],
        )
        writer.writeheader()
        writer.writerows(all_numbering_rows)

    with open(OUT_CDRS, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "record",
                "scheme",
                "region",
                "sequence",
                "length",
            ],
        )
        writer.writeheader()
        writer.writerows(all_cdr_rows)

    print("\nSaved:")
    print(f"  {OUT_NUMBERING}")
    print(f"  {OUT_CDRS}")


if __name__ == "__main__":
    main()