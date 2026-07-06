from pathlib import Path
import csv
from collections import defaultdict


INPUT_FASTA = Path("data/input/final_pdl1_candidate.fasta")
INPUT_LIABILITIES = Path("data/intermediate/liability_report.csv")
INPUT_FEATURES = Path("data/intermediate/developability_features.csv")

OUT_MUTATIONS = Path("data/intermediate/proposed_mutations.csv")
OUT_VARIANTS = Path("data/intermediate/proposed_variants.csv")
OUT_FASTA = Path("data/intermediate/proposed_variants.fasta")


MAX_SINGLE_VARIANTS = 10

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")


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
    with open(path, "r") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def validate_sequence(record_name, sequence):
    invalid = sorted(set(sequence) - CANONICAL_AA)
    if invalid:
        raise ValueError(
            f"{record_name} contains non-canonical residues: {invalid}. "
            "Replace placeholders such as X before proposing variants."
        )


def is_cdr(row):
    return row.get("in_cdr", "") == "True"


def severity_weight(severity):
    weights = {
        "high": 4,
        "medium": 3,
        "low": 2,
        "review": 1,
    }
    return weights.get(severity, 0)


def region_weight(row):
    """
    Framework mutations are safer for first-pass automated proposal.
    CDR mutations are not forbidden, but they require later structure/contact review.
    """
    if is_cdr(row):
        return 1
    return 3


def make_mutation(sequence, position, new_aa):
    """
    position is 1-indexed.
    """
    old_aa = sequence[position - 1]

    if old_aa == new_aa:
        return None

    mutant = sequence[:position - 1] + new_aa + sequence[position:]
    mutation_label = f"{old_aa}{position}{new_aa}"

    return mutation_label, mutant


def choose_n_glycosylation_mutations(sequence, row):
    """
    N-X-S/T motif.
    Conservative first-pass option:
    - N->Q removes glycosylation consensus and preserves side-chain amide character.
    """
    start = int(row["start"])
    motif = row["motif"]

    proposals = []

    if len(motif) == 3 and sequence[start - 1] == "N":
        proposals.append(
            {
                "position": start,
                "new_aa": "Q",
                "mutation_strategy": "N->Q to remove N-linked glycosylation motif while preserving amide-like character",
                "risk_note": "May affect binding if the motif is in a CDR or antigen-contact residue.",
            }
        )

    return proposals


def choose_deamidation_mutations(sequence, row):
    """
    N[GSTH] motifs.
    Conservative first-pass option:
    - N->Q often reduces deamidation risk while preserving amide-like character.
    """
    start = int(row["start"])
    motif = row["motif"]

    proposals = []

    if len(motif) >= 2 and sequence[start - 1] == "N":
        proposals.append(
            {
                "position": start,
                "new_aa": "Q",
                "mutation_strategy": "N->Q to reduce Asn deamidation motif risk",
                "risk_note": "Requires structure review if the Asn is in CDR or antigen-contacting region.",
            }
        )

    return proposals


def choose_asp_isomerization_mutations(sequence, row):
    """
    D[GST] motifs.
    Conservative first-pass option:
    - D->E can reduce Asp isomerization tendency while preserving negative charge.
    """
    start = int(row["start"])
    motif = row["motif"]

    proposals = []

    if len(motif) >= 2 and sequence[start - 1] == "D":
        proposals.append(
            {
                "position": start,
                "new_aa": "E",
                "mutation_strategy": "D->E to reduce Asp isomerization motif risk while preserving negative charge",
                "risk_note": "Requires structure review if the Asp is involved in binding or loop conformation.",
            }
        )

    return proposals


def choose_oxidation_mutations(sequence, row):
    """
    M/W oxidation.
    Conservative policy:
    - M outside CDR: M->L
    - W outside CDR: W->F only as review-level option
    - CDR M/W: do not automatically mutate unless later structure review supports it
    """
    start = int(row["start"])
    aa = sequence[start - 1]

    proposals = []

    if is_cdr(row):
        return proposals

    if aa == "M":
        proposals.append(
            {
                "position": start,
                "new_aa": "L",
                "mutation_strategy": "M->L to reduce oxidation risk outside CDR",
                "risk_note": "Check structural packing; Met may contribute to local stability.",
            }
        )

    elif aa == "W":
        proposals.append(
            {
                "position": start,
                "new_aa": "F",
                "mutation_strategy": "W->F to reduce oxidation risk outside CDR",
                "risk_note": "Tryptophan can be structurally important; treat as lower-priority review mutation.",
            }
        )

    return proposals


def choose_polybasic_patch_mutations(sequence, row):
    """
    Reduce local Lys/Arg enrichment.
    Conservative first-pass option:
    - K->Q or R->Q for one residue in a non-CDR polybasic patch.
    """
    if is_cdr(row):
        return []

    start = int(row["start"])
    end = int(row["end"])

    proposals = []

    for pos in range(start, end + 1):
        aa = sequence[pos - 1]
        if aa in {"K", "R"}:
            proposals.append(
                {
                    "position": pos,
                    "new_aa": "Q",
                    "mutation_strategy": f"{aa}->Q to reduce local polybasic patch outside CDR",
                    "risk_note": "Check whether this residue contributes to structural stability or antigen-independent solubility.",
                }
            )
            break

    return proposals


def choose_hydrophobic_patch_mutations(sequence, row):
    """
    Sequence-only hydrophobic patches are hard to interpret without solvent exposure.
    Therefore this script does not automatically mutate hydrophobic patches.
    Hydrophobic residues can be core-stabilizing, especially in frameworks.

    Later structure-based analysis should decide whether a hydrophobic patch is surface-exposed.
    """
    return []


def propose_from_liability(sequence, row):
    liability_type = row["liability_type"]

    if liability_type == "N-linked glycosylation motif":
        return choose_n_glycosylation_mutations(sequence, row)

    if liability_type == "Asn deamidation motif":
        return choose_deamidation_mutations(sequence, row)

    if liability_type == "Asp isomerization motif":
        return choose_asp_isomerization_mutations(sequence, row)

    if liability_type == "Oxidation-prone residue":
        return choose_oxidation_mutations(sequence, row)

    if liability_type == "Polybasic patch":
        return choose_polybasic_patch_mutations(sequence, row)

    if liability_type == "Hydrophobic patch":
        return choose_hydrophobic_patch_mutations(sequence, row)

    # Cysteines are review-only at this stage.
    return []


def rank_proposed_mutation(row):
    """
    Higher is better.
    Prioritize:
    - high/medium severity
    - non-CDR mutations
    - specific sequence liabilities over broad patch flags
    """
    score = 0

    score += severity_weight(row["source_severity"])
    score += region_weight(row)

    if row["source_liability_type"] in {
        "N-linked glycosylation motif",
        "Asn deamidation motif",
        "Asp isomerization motif",
    }:
        score += 3

    if row["source_liability_type"] == "Oxidation-prone residue":
        score += 1

    if row["source_liability_type"] == "Polybasic patch":
        score += 1

    # CDR mutations are not forbidden, but lower priority for automated proposal.
    if row["in_cdr"] == "True":
        score -= 2

    return score


def apply_mutations(sequence, mutation_rows):
    mutant = sequence
    applied = []

    # Apply from N to C; positions refer to original sequence.
    for row in sorted(mutation_rows, key=lambda x: int(x["position"])):
        pos = int(row["position"])
        old_aa = row["old_aa"]
        new_aa = row["new_aa"]

        current_aa = mutant[pos - 1]

        if current_aa != old_aa:
            raise ValueError(
                f"Mutation mismatch at position {pos}: expected {old_aa}, found {current_aa}"
            )

        mutant = mutant[:pos - 1] + new_aa + mutant[pos:]
        applied.append(f"{old_aa}{pos}{new_aa}")

    return mutant, applied


def generate_variants(record_name, sequence, mutation_rows):
    variants = []

    # Single mutants
    ranked = sorted(
        mutation_rows,
        key=lambda row: (-int(row["priority_score"]), int(row["position"]), row["new_aa"]),
    )

    single_rows = ranked[:MAX_SINGLE_VARIANTS]

    for i, row in enumerate(single_rows, start=1):
        mutant_seq, applied = apply_mutations(sequence, [row])
        variant_id = f"{record_name}_V{i:02d}_{applied[0]}"

        variants.append(
            {
                "variant_id": variant_id,
                "parent_record": record_name,
                "variant_type": "single_mutation",
                "mutations": ";".join(applied),
                "num_mutations": len(applied),
                "rationale": row["mutation_strategy"],
                "risk_note": row["risk_note"],
                "sequence": mutant_seq,
            }
        )

    # Combined conservative variant:
    # only non-CDR mutations, no duplicate positions, highest priority first.
    combined_rows = []
    used_positions = set()

    for row in ranked:
        pos = int(row["position"])

        if row["in_cdr"] == "True":
            continue

        if pos in used_positions:
            continue

        combined_rows.append(row)
        used_positions.add(pos)

        if len(combined_rows) >= 5:
            break

    if len(combined_rows) >= 2:
        mutant_seq, applied = apply_mutations(sequence, combined_rows)

        variant_id = f"{record_name}_V99_combined_conservative"

        rationale = (
            "Combined conservative variant containing multiple non-CDR liability-reducing mutations. "
            "This variant is intended for structural re-evaluation, not direct experimental nomination."
        )

        risk_note = (
            "Because multiple mutations are combined, this variant requires RF2/Boltz/ProteinMPNN "
            "or equivalent structural review before prioritization."
        )

        variants.append(
            {
                "variant_id": variant_id,
                "parent_record": record_name,
                "variant_type": "combined_conservative",
                "mutations": ";".join(applied),
                "num_mutations": len(applied),
                "rationale": rationale,
                "risk_note": risk_note,
                "sequence": mutant_seq,
            }
        )

    return variants


def main():
    records = read_fasta(INPUT_FASTA)
    liability_rows = read_csv(INPUT_LIABILITIES)

    all_mutation_rows = []
    all_variant_rows = []
    fasta_records = []

    for record_name, sequence in records:
        validate_sequence(record_name, sequence)

        print(f"\nProposing rational variants for: {record_name}")
        print(f"Sequence length: {len(sequence)} aa")

        record_liabilities = [
            row for row in liability_rows
            if row["record"] == record_name
        ]

        proposed_by_key = {}

        for liability in record_liabilities:
            proposals = propose_from_liability(sequence, liability)

            for proposal in proposals:
                pos = int(proposal["position"])
                old_aa = sequence[pos - 1]
                new_aa = proposal["new_aa"]

                mutation_result = make_mutation(sequence, pos, new_aa)

                if mutation_result is None:
                    continue

                mutation_label, _ = mutation_result

                key = (record_name, pos, new_aa)

                row = {
                    "record": record_name,
                    "mutation": mutation_label,
                    "position": pos,
                    "old_aa": old_aa,
                    "new_aa": new_aa,
                    "region": liability["region"],
                    "in_cdr": liability["in_cdr"],
                    "source_liability_type": liability["liability_type"],
                    "source_motif": liability["motif"],
                    "source_start": liability["start"],
                    "source_end": liability["end"],
                    "source_severity": liability["severity"],
                    "mutation_strategy": proposal["mutation_strategy"],
                    "risk_note": proposal["risk_note"],
                    "priority_score": 0,
                }

                row["priority_score"] = rank_proposed_mutation(row)

                # Deduplicate same mutation if it addresses multiple overlapping liabilities.
                if key not in proposed_by_key:
                    proposed_by_key[key] = row
                else:
                    existing = proposed_by_key[key]
                    existing["source_liability_type"] += f";{liability['liability_type']}"
                    existing["source_motif"] += f";{liability['motif']}"

        mutation_rows = list(proposed_by_key.values())

        mutation_rows = sorted(
            mutation_rows,
            key=lambda row: (-int(row["priority_score"]), int(row["position"]), row["new_aa"]),
        )

        all_mutation_rows.extend(mutation_rows)

        variants = generate_variants(record_name, sequence, mutation_rows)

        all_variant_rows.extend(variants)

        fasta_records.append((record_name + "_WT", sequence))
        for variant in variants:
            fasta_records.append((variant["variant_id"], variant["sequence"]))

        print(f"Candidate liability flags: {len(record_liabilities)}")
        print(f"Proposed mutations: {len(mutation_rows)}")
        print(f"Generated variants: {len(variants)}")

        if mutation_rows:
            print("\nTop proposed mutations:")
            for row in mutation_rows[:10]:
                print(
                    f"  {row['mutation']} | {row['region']} | "
                    f"{row['source_liability_type']} | priority={row['priority_score']}"
                )
        else:
            print("No automatic mutations proposed. Review liability report manually.")

    mutation_fields = [
        "record",
        "mutation",
        "position",
        "old_aa",
        "new_aa",
        "region",
        "in_cdr",
        "source_liability_type",
        "source_motif",
        "source_start",
        "source_end",
        "source_severity",
        "mutation_strategy",
        "risk_note",
        "priority_score",
    ]

    variant_fields = [
        "variant_id",
        "parent_record",
        "variant_type",
        "mutations",
        "num_mutations",
        "rationale",
        "risk_note",
        "sequence",
    ]

    write_csv(OUT_MUTATIONS, all_mutation_rows, mutation_fields)
    write_csv(OUT_VARIANTS, all_variant_rows, variant_fields)
    write_fasta(OUT_FASTA, fasta_records)

    print("\nSaved:")
    print(f"  {OUT_MUTATIONS}")
    print(f"  {OUT_VARIANTS}")
    print(f"  {OUT_FASTA}")


if __name__ == "__main__":
    main()
