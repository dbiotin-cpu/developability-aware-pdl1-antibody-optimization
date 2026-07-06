#!/usr/bin/env python

from pathlib import Path
import csv
import re


INPUT_VARIANTS_FASTA = Path("data/intermediate/proposed_variants.fasta")
INPUT_RANKING = Path("data/intermediate/variant_ranking_pre_structure.csv")
INPUT_PDL1_FASTA = Path("data/input/pdl1_target.fasta")

OUT_DIR = Path("data/structure_eval")
OUT_RF2_DIR = OUT_DIR / "rf2_inputs"
OUT_BOLTZ_DIR = OUT_DIR / "boltz_inputs"
OUT_MPNN_DIR = OUT_DIR / "proteinmpnn_inputs"

OUT_TOP_FASTA = OUT_DIR / "top_variants.fasta"
OUT_TOP_METADATA = OUT_DIR / "top_variants_metadata.csv"
OUT_STRUCTURE_PLAN = Path("docs/structure_evaluation_plan.md")

TOP_N = 5


def read_fasta(path):
    records = []
    name = None
    seq_lines = []

    if not path.exists():
        return records

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
        raise FileNotFoundError(f"Missing required file: {path}")

    with open(path, "r") as f:
        return list(csv.DictReader(f))


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sanitize_filename(name):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return name.strip("_")


def select_top_variants(ranking_rows, top_n=TOP_N):
    non_wt = [
        row for row in ranking_rows
        if not row["variant_id"].endswith("_WT")
    ]

    preferred = [
        row for row in non_wt
        if row.get("pre_structure_decision", "").startswith("Advance")
        or row.get("pre_structure_decision", "").startswith("Consider")
    ]

    if preferred:
        candidates = preferred
    else:
        candidates = non_wt

    for row in candidates:
        row["_optimization_score"] = float(row.get("optimization_score_pre_structure", 0))
        row["_risk_score"] = float(row.get("risk_score", 999))

    candidates = sorted(
        candidates,
        key=lambda row: (-row["_optimization_score"], row["_risk_score"]),
    )

    return candidates[:top_n]


def write_rf2_complex_fasta(path, variant_id, binder_seq, antigen_name, antigen_seq):
    """
    This writes a generic multi-chain FASTA.

    Depending on your RF2/RoseTTAFold2 setup, you may need to adapt the separator
    or input format. This file is intended as a clean starting point.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        f.write(f">{variant_id}_binder_chain_H\n")
        f.write(binder_seq + "\n")
        f.write(f">{antigen_name}_chain_A\n")
        f.write(antigen_seq + "\n")


def write_boltz_yaml(path, variant_id, binder_seq, antigen_seq=None):
    """
    Boltz-style YAML input.

    Chain A = PD-L1 target
    Chain H = antibody/binder variant

    If antigen_seq is unavailable, a binder-only YAML is written.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        f.write("version: 1\n")
        f.write("sequences:\n")

        if antigen_seq:
            f.write("  - protein:\n")
            f.write("      id: A\n")
            f.write(f"      sequence: {antigen_seq}\n")

        f.write("  - protein:\n")
        f.write("      id: H\n")
        f.write(f"      sequence: {binder_seq}\n")


def write_proteinmpnn_readme(path, selected_rows):
    path.parent.mkdir(parents=True, exist_ok=True)

    text = """# ProteinMPNN Preparation Notes

ProteinMPNN is most useful here after a structure model is available.

Recommended workflow:

1. Predict or prepare the WT PD-L1/binder complex structure.
2. Predict or prepare structures for the top developability-improved variants.
3. Identify residues that should remain fixed:
   - PD-L1 chain
   - binder residues contacting PD-L1
   - CDR residues that define the predicted paratope
   - canonical framework residues
   - cysteines involved in expected disulfide pairing

4. Allow design only at selected liability positions or neighboring tolerated positions.
5. Use ProteinMPNN fixed-position constraints to avoid redesigning the whole binder.

The current selected variants for structural re-evaluation are:

"""

    for row in selected_rows:
        text += f"- {row['variant_id']}: {row.get('mutations', '')}\n"

    text += """
Important limitation:
This step does not yet run ProteinMPNN. It prepares the ranked sequence set that should be
evaluated structurally first. ProteinMPNN redesign should be applied only after binding-interface
and structural preservation constraints are defined.
"""

    with open(path, "w") as f:
        f.write(text)


def write_structure_plan(path, selected_rows, has_antigen):
    path.parent.mkdir(parents=True, exist_ok=True)

    antigen_status = (
        "PD-L1 target FASTA was found, so complex input files were generated."
        if has_antigen
        else "PD-L1 target FASTA was not found. Add data/input/pdl1_target.fasta before complex modeling."
    )

    text = f"""# Structure-Based Re-evaluation Plan

## Goal

The goal of this step is to evaluate whether developability-improved variants preserve the intended PD-L1 binding mode.

## Input Status

{antigen_status}

## Selected Variants

"""

    for i, row in enumerate(selected_rows, start=1):
        text += f"""### {i}. {row['variant_id']}

- Mutations: {row.get('mutations', '')}
- Variant type: {row.get('variant_type', '')}
- Pre-structure optimization score: {row.get('optimization_score_pre_structure', '')}
- Risk score: {row.get('risk_score', '')}
- Decision: {row.get('pre_structure_decision', '')}

"""

    text += """## Recommended Evaluation Criteria

For each variant, compare against the original WT PD-L1/binder model:

1. Binder fold preservation
   - Does the binder remain antibody-like or domain-like?
   - Are CDR loop regions structurally reasonable?

2. Binding-mode preservation
   - Does the variant remain close to the original PD-L1 epitope?
   - Does the paratope orientation remain similar?

3. Interface geometry
   - Interface RMSD versus WT model
   - Contact residue overlap
   - Loss or gain of major contact residues

4. Developability improvement
   - Reduced liability count
   - Reduced CDR-localized liabilities
   - Improved or neutral pI/net charge
   - No new obvious hydrophobic or polybasic patches

## Suggested Final Ranking Logic

Final score should combine:

- Epitope retention
- Binding-mode preservation
- Structural confidence
- Developability improvement
- Number of mutations

Recommended decision labels:

- Advance
- Backup
- Watch
- Reject

## Important Limitation

Boltz/RF2/ProteinMPNN results should be treated as in silico triage, not experimental validation.
Top candidates would still require expression, purification, SEC, DSF, SPR/BLI, and cell-binding validation.
"""

    with open(path, "w") as f:
        f.write(text)


def main():
    variant_records = read_fasta(INPUT_VARIANTS_FASTA)
    ranking_rows = read_csv(INPUT_RANKING)
    antigen_records = read_fasta(INPUT_PDL1_FASTA)

    if not variant_records:
        raise ValueError(f"No variant sequences found in {INPUT_VARIANTS_FASTA}")

    variant_seq = {name: seq for name, seq in variant_records}

    selected_rows = select_top_variants(ranking_rows, top_n=TOP_N)

    if not selected_rows:
        raise ValueError("No variants selected. Check variant_ranking_pre_structure.csv")

    antigen_name = None
    antigen_seq = None

    if antigen_records:
        antigen_name, antigen_seq = antigen_records[0]
        print(f"Loaded antigen FASTA: {antigen_name}")
        print(f"Antigen length: {len(antigen_seq)} aa")
    else:
        print(f"WARNING: No antigen FASTA found at {INPUT_PDL1_FASTA}")
        print("Complex RF2/Boltz inputs will be limited until PD-L1 sequence is added.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_RF2_DIR.mkdir(parents=True, exist_ok=True)
    OUT_BOLTZ_DIR.mkdir(parents=True, exist_ok=True)
    OUT_MPNN_DIR.mkdir(parents=True, exist_ok=True)

    top_fasta_records = []
    metadata_rows = []

    print("\nSelected variants for structure re-evaluation:")

    for row in selected_rows:
        variant_id = row["variant_id"]

        if variant_id not in variant_seq:
            print(f"WARNING: {variant_id} not found in proposed_variants.fasta. Skipping.")
            continue

        seq = variant_seq[variant_id]
        safe_id = sanitize_filename(variant_id)

        print(
            f"  {variant_id} | score={row.get('optimization_score_pre_structure', '')} | "
            f"risk={row.get('risk_score', '')}"
        )

        top_fasta_records.append((variant_id, seq))

        metadata_rows.append(
            {
                "variant_id": variant_id,
                "variant_type": row.get("variant_type", ""),
                "mutations": row.get("mutations", ""),
                "num_mutations": row.get("num_mutations", ""),
                "risk_score": row.get("risk_score", ""),
                "optimization_score_pre_structure": row.get("optimization_score_pre_structure", ""),
                "pre_structure_decision": row.get("pre_structure_decision", ""),
                "binder_sequence_length": len(seq),
                "rf2_input_file": str(OUT_RF2_DIR / f"{safe_id}_complex.fasta") if antigen_seq else "",
                "boltz_input_file": str(OUT_BOLTZ_DIR / f"{safe_id}.yaml"),
            }
        )

        # Binder-only Boltz input is still useful for checking fold if antigen sequence is absent.
        write_boltz_yaml(
            path=OUT_BOLTZ_DIR / f"{safe_id}.yaml",
            variant_id=variant_id,
            binder_seq=seq,
            antigen_seq=antigen_seq,
        )

        if antigen_seq:
            write_rf2_complex_fasta(
                path=OUT_RF2_DIR / f"{safe_id}_complex.fasta",
                variant_id=variant_id,
                binder_seq=seq,
                antigen_name=antigen_name,
                antigen_seq=antigen_seq,
            )

    write_fasta(OUT_TOP_FASTA, top_fasta_records)

    metadata_fields = [
        "variant_id",
        "variant_type",
        "mutations",
        "num_mutations",
        "risk_score",
        "optimization_score_pre_structure",
        "pre_structure_decision",
        "binder_sequence_length",
        "rf2_input_file",
        "boltz_input_file",
    ]

    write_csv(OUT_TOP_METADATA, metadata_rows, metadata_fields)

    write_proteinmpnn_readme(
        path=OUT_MPNN_DIR / "README_proteinmpnn_preparation.md",
        selected_rows=selected_rows,
    )

    write_structure_plan(
        path=OUT_STRUCTURE_PLAN,
        selected_rows=selected_rows,
        has_antigen=bool(antigen_seq),
    )

    if not antigen_seq:
        template_path = Path("data/input/pdl1_target.fasta.template")
        with open(template_path, "w") as f:
            f.write(">PDL1_target\n")
            f.write("PASTE_PDL1_SEQUENCE_HERE\n")

    print("\nSaved:")
    print(f"  {OUT_TOP_FASTA}")
    print(f"  {OUT_TOP_METADATA}")
    print(f"  {OUT_BOLTZ_DIR}")
    print(f"  {OUT_RF2_DIR}")
    print(f"  {OUT_MPNN_DIR / 'README_proteinmpnn_preparation.md'}")
    print(f"  {OUT_STRUCTURE_PLAN}")

    if not antigen_seq:
        print("\nReminder:")
        print("  Add PD-L1 sequence to data/input/pdl1_target.fasta before complex modeling.")
        print("  Template created at data/input/pdl1_target.fasta.template")


if __name__ == "__main__":
    main()
