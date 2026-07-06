from pathlib import Path
import csv
import math


PRE_STRUCTURE_RANKING = Path("data/intermediate/variant_ranking_pre_structure.csv")
BOLTZ_SUMMARY = Path("data/structure_eval/boltz_structure_summary.csv")
EPITOPE_SCORES = Path("data/structure_eval/epitope_retention_scores.csv")

OUT_FINAL_RANKING = Path("data/output/final_ranked_variants.csv")
OUT_FINAL_RECOMMENDATION = Path("data/output/final_recommendation.md")


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")

    with open(path, "r") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value, default=None):
    if value is None or value == "":
        return default

    try:
        return float(value)
    except ValueError:
        return default


def safe_int(value, default=0):
    if value is None or value == "":
        return default

    try:
        return int(float(value))
    except ValueError:
        return default


def index_by_variant(rows):
    return {row["variant_id"]: row for row in rows}


def normalize_positive(value, max_value):
    if value is None:
        return 0.0

    if max_value is None or max_value <= 0:
        return 0.0

    return max(0.0, min(1.0, value / max_value))


def mutation_burden_component(num_mutations):
    """
    Fewer mutations are preferred.
    0 mutations = 1.0
    6 or more mutations = 0.0
    """
    num_mutations = safe_int(num_mutations, default=0)
    return round(max(0.0, 1.0 - min(num_mutations, 6) / 6.0), 3)


def structure_component(epitope_row):
    """
    Prefer explicit epitope_retention_score if available.
    If unavailable, fall back to partial metrics.
    """
    score = safe_float(epitope_row.get("epitope_retention_score", ""), default=None)

    if score is not None:
        return max(0.0, min(1.0, score))

    epitope = safe_float(epitope_row.get("epitope_contact_retention", ""), default=None)
    binder = safe_float(epitope_row.get("binder_contact_retention", ""), default=None)
    jaccard = safe_float(epitope_row.get("contact_pair_jaccard_vs_reference", ""), default=None)

    available = [x for x in [epitope, binder, jaccard] if x is not None]

    if available:
        return round(sum(available) / len(available), 3)

    return 0.0


def interface_component(epitope_row):
    """
    Contact-pair Jaccard is a strong interface-preservation metric.
    If missing, use epitope retention.
    """
    jaccard = safe_float(epitope_row.get("contact_pair_jaccard_vs_reference", ""), default=None)

    if jaccard is not None:
        return max(0.0, min(1.0, jaccard))

    epitope = safe_float(epitope_row.get("epitope_contact_retention", ""), default=None)

    if epitope is not None:
        return max(0.0, min(1.0, epitope))

    return 0.0


def rmsd_penalty(epitope_row):
    rmsd = safe_float(epitope_row.get("binder_ca_rmsd_after_target_alignment_A", ""), default=None)

    if rmsd is None:
        return 0.0

    if rmsd <= 3:
        return 0.0
    elif rmsd <= 6:
        return 0.05
    elif rmsd <= 10:
        return 0.10
    else:
        return 0.20


def decide_final(row):
    final_score = safe_float(row["final_integrated_score"], default=0.0)
    structure_score = safe_float(row["structure_retention_component"], default=0.0)
    dev_score = safe_float(row["developability_component"], default=0.0)
    epitope_decision = row.get("structure_retention_decision", "")

    if row.get("is_reference") == "True":
        return "Reference parent"

    if "Low retention" in epitope_decision:
        return "Reject: binding mode not retained"

    if structure_score >= 0.75 and dev_score >= 0.40 and final_score >= 0.65:
        return "Advance"

    if structure_score >= 0.60 and final_score >= 0.50:
        return "Backup"

    if structure_score >= 0.40 and final_score >= 0.35:
        return "Watch"

    return "Reject"


def make_markdown_report(rows):
    ranked = [row for row in rows if row["final_decision"] != "Reference parent"]

    advanced = [row for row in ranked if row["final_decision"] == "Advance"]
    backups = [row for row in ranked if row["final_decision"] == "Backup"]
    watches = [row for row in ranked if row["final_decision"] == "Watch"]

    lines = []

    lines.append("# Final Variant Recommendation")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(
        "Developability-optimized variants were ranked by combining pre-structure "
        "developability improvement with Boltz-based structure/interface retention metrics."
    )
    lines.append("")
    lines.append("The final integrated score used four components:")
    lines.append("")
    lines.append("- Developability improvement relative to the parent candidate")
    lines.append("- Epitope/contact retention relative to the WT/reference complex")
    lines.append("- Interface contact preservation")
    lines.append("- Mutation burden")
    lines.append("")
    lines.append("Boltz/RF2-style structural scores are treated as in silico triage, not experimental validation.")
    lines.append("")

    if advanced:
        lines.append("## Recommended Variants to Advance")
        lines.append("")
        for row in advanced:
            lines.append(f"### {row['final_rank']}. {row['variant_id']}")
            lines.append("")
            lines.append(f"- Mutations: {row['mutations']}")
            lines.append(f"- Final integrated score: {row['final_integrated_score']}")
            lines.append(f"- Developability component: {row['developability_component']}")
            lines.append(f"- Structure retention component: {row['structure_retention_component']}")
            lines.append(f"- Interface component: {row['interface_component']}")
            lines.append(f"- Mutation burden component: {row['mutation_burden_component']}")
            lines.append(f"- Final decision: {row['final_decision']}")
            lines.append("")
    else:
        lines.append("## Recommended Variants to Advance")
        lines.append("")
        lines.append("No variant met the strict Advance threshold. Review Backup candidates for possible experimental follow-up.")
        lines.append("")

    if backups:
        lines.append("## Backup Variants")
        lines.append("")
        for row in backups[:5]:
            lines.append(
                f"- **{row['variant_id']}**: {row['mutations']} "
                f"(score={row['final_integrated_score']})"
            )
        lines.append("")

    if watches:
        lines.append("## Watch List")
        lines.append("")
        for row in watches[:5]:
            lines.append(
                f"- **{row['variant_id']}**: {row['mutations']} "
                f"(score={row['final_integrated_score']})"
            )
        lines.append("")

    lines.append("## Recommended Experimental Validation")
    lines.append("")
    lines.append("Top-ranked variants should be tested experimentally by:")
    lines.append("")
    lines.append("- Transient mammalian expression")
    lines.append("- Protein A or affinity purification")
    lines.append("- SEC purity and aggregation analysis")
    lines.append("- DSF or nanoDSF thermal stability")
    lines.append("- SPR/BLI binding to PD-L1")
    lines.append("- Cell-surface PD-L1 binding by flow cytometry")
    lines.append("- Specificity or polyspecificity assessment")
    lines.append("")

    return "\n".join(lines)


def main():
    pre_rows = read_csv(PRE_STRUCTURE_RANKING)
    boltz_rows = read_csv(BOLTZ_SUMMARY)
    epitope_rows = read_csv(EPITOPE_SCORES)

    boltz_by_variant = index_by_variant(boltz_rows)
    epitope_by_variant = index_by_variant(epitope_rows)

    non_ref_pre_rows = [
        row for row in pre_rows
        if not row["variant_id"].endswith("_WT")
    ]

    max_opt = max(
        [
            safe_float(row.get("optimization_score_pre_structure", ""), default=0.0)
            for row in non_ref_pre_rows
        ]
        or [1.0]
    )

    final_rows = []

    for pre in pre_rows:
        variant_id = pre["variant_id"]

        boltz = boltz_by_variant.get(variant_id, {})
        epitope = epitope_by_variant.get(variant_id, {})

        is_reference = variant_id.endswith("_WT")

        pre_opt = safe_float(pre.get("optimization_score_pre_structure", ""), default=0.0)

        developability_component = normalize_positive(pre_opt, max_opt)
        struct_component = structure_component(epitope)
        int_component = interface_component(epitope)
        mut_component = mutation_burden_component(pre.get("num_mutations", "0"))

        penalty = rmsd_penalty(epitope)

        if is_reference:
            final_score = 0.0
        else:
            final_score = (
                0.40 * developability_component
                + 0.35 * struct_component
                + 0.15 * int_component
                + 0.10 * mut_component
                - penalty
            )

        final_score = round(max(0.0, min(1.0, final_score)), 3)

        row = {
            "variant_id": variant_id,
            "is_reference": str(is_reference),
            "variant_type": pre.get("variant_type", ""),
            "mutations": pre.get("mutations", ""),
            "num_mutations": pre.get("num_mutations", "0"),
            "pre_structure_decision": pre.get("pre_structure_decision", ""),
            "pre_structure_risk_score": pre.get("risk_score", ""),
            "delta_risk_score_vs_WT": pre.get("delta_risk_score_vs_WT", ""),
            "total_liability_flags": pre.get("total_liability_flags", ""),
            "liability_reduction_vs_WT": pre.get("liability_reduction_vs_WT", ""),
            "cdr_liability_flags": pre.get("cdr_liability_flags", ""),
            "cdr_liability_reduction_vs_WT": pre.get("cdr_liability_reduction_vs_WT", ""),
            "pI": pre.get("pI", ""),
            "net_charge_pH_7_4": pre.get("net_charge_pH_7_4", ""),
            "cdr_hydrophobic_fraction": pre.get("cdr_hydrophobic_fraction", ""),
            "cdr_aromatic_fraction": pre.get("cdr_aromatic_fraction", ""),
            "optimization_score_pre_structure": pre.get("optimization_score_pre_structure", ""),
            "developability_component": round(developability_component, 3),
            "structure_retention_component": round(struct_component, 3),
            "interface_component": round(int_component, 3),
            "mutation_burden_component": round(mut_component, 3),
            "rmsd_penalty": round(penalty, 3),
            "epitope_contact_retention": epitope.get("epitope_contact_retention", ""),
            "binder_contact_retention": epitope.get("binder_contact_retention", ""),
            "contact_pair_jaccard_vs_reference": epitope.get("contact_pair_jaccard_vs_reference", ""),
            "interface_contact_count_ratio_vs_reference": epitope.get("interface_contact_count_ratio_vs_reference", ""),
            "binder_ca_rmsd_after_target_alignment_A": epitope.get("binder_ca_rmsd_after_target_alignment_A", ""),
            "epitope_retention_score": epitope.get("epitope_retention_score", ""),
            "structure_retention_decision": epitope.get("structure_retention_decision", ""),
            "boltz_structure_file": boltz.get("representative_structure_copied", ""),
            "available_boltz_metrics": boltz.get("available_boltz_metrics", ""),
            "final_integrated_score": final_score,
            "final_decision": "",
            "final_rank": "",
        }

        row["final_decision"] = decide_final(row)

        final_rows.append(row)

    # Sort: reference last, then by final score descending.
    final_rows = sorted(
        final_rows,
        key=lambda row: (
            row["is_reference"] == "True",
            -safe_float(row["final_integrated_score"], default=0.0),
            safe_float(row["pre_structure_risk_score"], default=999.0),
        ),
    )

    rank = 1
    for row in final_rows:
        if row["is_reference"] == "True":
            row["final_rank"] = "Reference"
        else:
            row["final_rank"] = rank
            rank += 1

    fieldnames = [
        "final_rank",
        "variant_id",
        "is_reference",
        "final_decision",
        "final_integrated_score",
        "developability_component",
        "structure_retention_component",
        "interface_component",
        "mutation_burden_component",
        "rmsd_penalty",
        "variant_type",
        "mutations",
        "num_mutations",
        "pre_structure_decision",
        "pre_structure_risk_score",
        "delta_risk_score_vs_WT",
        "total_liability_flags",
        "liability_reduction_vs_WT",
        "cdr_liability_flags",
        "cdr_liability_reduction_vs_WT",
        "pI",
        "net_charge_pH_7_4",
        "cdr_hydrophobic_fraction",
        "cdr_aromatic_fraction",
        "optimization_score_pre_structure",
        "epitope_contact_retention",
        "binder_contact_retention",
        "contact_pair_jaccard_vs_reference",
        "interface_contact_count_ratio_vs_reference",
        "binder_ca_rmsd_after_target_alignment_A",
        "epitope_retention_score",
        "structure_retention_decision",
        "boltz_structure_file",
        "available_boltz_metrics",
    ]

    write_csv(OUT_FINAL_RANKING, final_rows, fieldnames)

    report = make_markdown_report(final_rows)

    OUT_FINAL_RECOMMENDATION.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_FINAL_RECOMMENDATION, "w") as f:
        f.write(report)

    print("Saved:")
    print(f"  {OUT_FINAL_RANKING}")
    print(f"  {OUT_FINAL_RECOMMENDATION}")

    print("\nTop final ranked variants:")
    for row in final_rows[:10]:
        print(
            f"  {row['final_rank']} | {row['variant_id']} | "
            f"score={row['final_integrated_score']} | "
            f"decision={row['final_decision']}"
        )


if __name__ == "__main__":
    main()
