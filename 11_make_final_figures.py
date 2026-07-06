from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


INPUT_FINAL_RANKING = Path("data/output/final_ranked_variants.csv")

FIGURE_DIR = Path("figures")
OUT_README_SUMMARY = Path("data/output/final_ranked_variants_readme_summary.csv")

OUT_FINAL_SCORE = FIGURE_DIR / "final_integrated_score_by_variant.png"
OUT_SCORE_COMPONENTS = FIGURE_DIR / "final_score_components_top_variants.png"
OUT_DEV_VS_STRUCTURE = FIGURE_DIR / "developability_vs_structure_retention.png"
OUT_DECISION_COUNTS = FIGURE_DIR / "final_decision_counts.png"


TOP_N = 10


def shorten_label(label, max_len=36):
    label = str(label)

    if len(label) <= max_len:
        return label

    return label[:max_len - 3] + "..."


def load_final_ranking(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run scripts/10_final_variant_ranking.py first."
        )

    df = pd.read_csv(path)

    numeric_cols = [
        "final_integrated_score",
        "developability_component",
        "structure_retention_component",
        "interface_component",
        "mutation_burden_component",
        "rmsd_penalty",
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
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "is_reference" in df.columns:
        df["is_reference"] = df["is_reference"].astype(str)
    else:
        df["is_reference"] = "False"

    return df


def get_non_reference(df):
    return df[df["is_reference"] != "True"].copy()


def make_final_score_plot(df):
    plot_df = get_non_reference(df)

    plot_df = plot_df.sort_values(
        ["final_integrated_score", "pre_structure_risk_score"],
        ascending=[False, True],
    ).head(TOP_N)

    labels = [shorten_label(x) for x in plot_df["variant_id"]]

    plt.figure(figsize=(10, max(4, 0.45 * len(plot_df))))
    plt.barh(labels, plot_df["final_integrated_score"])
    plt.xlabel("Final integrated score, higher is better")
    plt.title("Final Integrated Score by Variant")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(OUT_FINAL_SCORE, dpi=300)
    plt.close()


def make_score_component_plot(df):
    plot_df = get_non_reference(df)

    plot_df = plot_df.sort_values(
        ["final_integrated_score", "pre_structure_risk_score"],
        ascending=[False, True],
    ).head(TOP_N)

    labels = [shorten_label(x) for x in plot_df["variant_id"]]

    component_cols = [
        "developability_component",
        "structure_retention_component",
        "interface_component",
        "mutation_burden_component",
    ]

    component_labels = [
        "Developability",
        "Structure retention",
        "Interface preservation",
        "Mutation burden",
    ]

    plt.figure(figsize=(11, max(4, 0.48 * len(plot_df))))

    left = pd.Series([0.0] * len(plot_df), index=plot_df.index)

    for col, label in zip(component_cols, component_labels):
        values = plot_df[col].fillna(0.0)
        plt.barh(labels, values, left=left, label=label)
        left = left + values

    plt.xlabel("Component value before final weighting")
    plt.title("Final Score Components for Top Variants")
    plt.legend(frameon=False)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(OUT_SCORE_COMPONENTS, dpi=300)
    plt.close()


def make_developability_vs_structure_plot(df):
    plot_df = get_non_reference(df)

    plt.figure(figsize=(8, 6))
    plt.scatter(
        plot_df["developability_component"],
        plot_df["structure_retention_component"],
    )

    for _, row in plot_df.iterrows():
        label = shorten_label(row["variant_id"], max_len=18)
        x = row["developability_component"]
        y = row["structure_retention_component"]

        if pd.notna(x) and pd.notna(y):
            plt.annotate(label, (x, y), fontsize=8, xytext=(4, 4), textcoords="offset points")

    plt.xlabel("Developability improvement component")
    plt.ylabel("Structure retention component")
    plt.title("Developability Improvement vs Structure Retention")
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.tight_layout()
    plt.savefig(OUT_DEV_VS_STRUCTURE, dpi=300)
    plt.close()


def make_decision_count_plot(df):
    plot_df = get_non_reference(df)

    counts = (
        plot_df["final_decision"]
        .fillna("Unknown")
        .value_counts()
        .reset_index()
    )

    counts.columns = ["final_decision", "count"]

    plt.figure(figsize=(9, 5))
    plt.bar(counts["final_decision"], counts["count"])
    plt.xlabel("Final decision")
    plt.ylabel("Number of variants")
    plt.title("Final Variant Decision Counts")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(OUT_DECISION_COUNTS, dpi=300)
    plt.close()


def save_readme_summary(df):
    plot_df = get_non_reference(df)

    plot_df = plot_df.sort_values(
        ["final_integrated_score", "pre_structure_risk_score"],
        ascending=[False, True],
    )

    keep_cols = [
        "final_rank",
        "variant_id",
        "final_decision",
        "final_integrated_score",
        "mutations",
        "num_mutations",
        "pre_structure_risk_score",
        "liability_reduction_vs_WT",
        "cdr_liability_reduction_vs_WT",
        "epitope_contact_retention",
        "contact_pair_jaccard_vs_reference",
        "binder_ca_rmsd_after_target_alignment_A",
        "epitope_retention_score",
    ]

    keep_cols = [col for col in keep_cols if col in plot_df.columns]

    OUT_README_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    plot_df[keep_cols].to_csv(OUT_README_SUMMARY, index=False)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    df = load_final_ranking(INPUT_FINAL_RANKING)

    print(f"Loaded final ranking: {INPUT_FINAL_RANKING}")
    print(f"Rows: {len(df)}")

    make_final_score_plot(df)
    make_score_component_plot(df)
    make_developability_vs_structure_plot(df)
    make_decision_count_plot(df)
    save_readme_summary(df)

    print("\nSaved figures:")
    print(f"  {OUT_FINAL_SCORE}")
    print(f"  {OUT_SCORE_COMPONENTS}")
    print(f"  {OUT_DEV_VS_STRUCTURE}")
    print(f"  {OUT_DECISION_COUNTS}")

    print("\nSaved README summary table:")
    print(f"  {OUT_README_SUMMARY}")


if __name__ == "__main__":
    main()