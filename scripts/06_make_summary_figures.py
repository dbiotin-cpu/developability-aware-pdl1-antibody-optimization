#!/usr/bin/env python

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


RANKING_FILE = Path("data/intermediate/variant_ranking_pre_structure.csv")
FEATURE_FILE = Path("data/intermediate/variant_developability_features.csv")
LIABILITY_FILE = Path("data/intermediate/variant_liability_report.csv")

FIGURE_DIR = Path("figures")

OUT_RISK_SCORE = FIGURE_DIR / "pre_structure_risk_score_by_variant.png"
OUT_OPT_SCORE = FIGURE_DIR / "pre_structure_optimization_score_by_variant.png"
OUT_LIABILITY_COUNT = FIGURE_DIR / "liability_count_by_variant.png"
OUT_CDR_HYDRO = FIGURE_DIR / "cdr_hydrophobicity_by_variant.png"
OUT_LIABILITY_TYPES = FIGURE_DIR / "liability_type_summary.png"

OUT_TOP_TABLE = Path("data/output/top_pre_structure_variants.csv")


def shorten_label(label, max_len=32):
    label = str(label)

    if len(label) <= max_len:
        return label

    return label[:max_len - 3] + "..."


def save_bar_plot(df, x_col, y_col, title, ylabel, output_path, top_n=12, ascending=True):
    plot_df = df.copy()

    if "variant_id" in plot_df.columns:
        wt_rows = plot_df[plot_df["variant_id"].str.endswith("_WT", na=False)]
        non_wt_rows = plot_df[~plot_df["variant_id"].str.endswith("_WT", na=False)]

        non_wt_rows = non_wt_rows.sort_values(y_col, ascending=ascending).head(top_n)

        if not wt_rows.empty:
            plot_df = pd.concat([wt_rows, non_wt_rows], axis=0)
        else:
            plot_df = non_wt_rows
    else:
        plot_df = plot_df.sort_values(y_col, ascending=ascending).head(top_n)

    labels = [shorten_label(x) for x in plot_df[x_col]]

    plt.figure(figsize=(10, max(4, 0.45 * len(plot_df))))
    plt.barh(labels, plot_df[y_col])
    plt.xlabel(ylabel)
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def make_risk_score_plot(ranking):
    df = ranking.copy()
    df["risk_score"] = pd.to_numeric(df["risk_score"], errors="coerce")

    save_bar_plot(
        df=df,
        x_col="variant_id",
        y_col="risk_score",
        title="Pre-Structure Developability Risk Score",
        ylabel="Risk score, lower is better",
        output_path=OUT_RISK_SCORE,
        top_n=12,
        ascending=True,
    )


def make_optimization_score_plot(ranking):
    df = ranking.copy()
    df["optimization_score_pre_structure"] = pd.to_numeric(
        df["optimization_score_pre_structure"], errors="coerce"
    )

    # For optimization score, higher is better.
    save_bar_plot(
        df=df,
        x_col="variant_id",
        y_col="optimization_score_pre_structure",
        title="Pre-Structure Optimization Score",
        ylabel="Optimization score, higher is better",
        output_path=OUT_OPT_SCORE,
        top_n=12,
        ascending=False,
    )


def make_liability_count_plot(ranking):
    df = ranking.copy()
    df["total_liability_flags"] = pd.to_numeric(
        df["total_liability_flags"], errors="coerce"
    )

    save_bar_plot(
        df=df,
        x_col="variant_id",
        y_col="total_liability_flags",
        title="Total Sequence Liability Flags by Variant",
        ylabel="Liability flag count, lower is better",
        output_path=OUT_LIABILITY_COUNT,
        top_n=12,
        ascending=True,
    )


def make_cdr_hydrophobicity_plot(ranking):
    df = ranking.copy()
    df["cdr_hydrophobic_fraction"] = pd.to_numeric(
        df["cdr_hydrophobic_fraction"], errors="coerce"
    )

    save_bar_plot(
        df=df,
        x_col="variant_id",
        y_col="cdr_hydrophobic_fraction",
        title="CDR Hydrophobic Fraction by Variant",
        ylabel="CDR hydrophobic fraction",
        output_path=OUT_CDR_HYDRO,
        top_n=12,
        ascending=True,
    )


def make_liability_type_summary(liabilities):
    if liabilities.empty:
        print("No liability rows found. Skipping liability type summary figure.")
        return

    summary = (
        liabilities.groupby("liability_type")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=True)
    )

    labels = [shorten_label(x, max_len=45) for x in summary["liability_type"]]

    plt.figure(figsize=(10, max(4, 0.5 * len(summary))))
    plt.barh(labels, summary["count"])
    plt.xlabel("Count across WT and proposed variants")
    plt.title("Liability Type Summary")
    plt.tight_layout()
    plt.savefig(OUT_LIABILITY_TYPES, dpi=300)
    plt.close()


def save_top_variant_table(ranking):
    df = ranking.copy()

    df["optimization_score_pre_structure"] = pd.to_numeric(
        df["optimization_score_pre_structure"], errors="coerce"
    )
    df["risk_score"] = pd.to_numeric(df["risk_score"], errors="coerce")

    top = df[~df["variant_id"].str.endswith("_WT", na=False)].copy()

    top = top.sort_values(
        ["optimization_score_pre_structure", "risk_score"],
        ascending=[False, True],
    )

    keep_cols = [
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
        "net_charge_pH_7_4",
        "cdr_hydrophobic_fraction",
        "cdr_aromatic_fraction",
        "optimization_score_pre_structure",
        "pre_structure_decision",
    ]

    keep_cols = [col for col in keep_cols if col in top.columns]

    OUT_TOP_TABLE.parent.mkdir(parents=True, exist_ok=True)
    top[keep_cols].to_csv(OUT_TOP_TABLE, index=False)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    if not RANKING_FILE.exists():
        raise FileNotFoundError(
            f"Missing {RANKING_FILE}. Run scripts/05_calculate_variant_developability_features.py first."
        )

    ranking = pd.read_csv(RANKING_FILE)

    if LIABILITY_FILE.exists():
        liabilities = pd.read_csv(LIABILITY_FILE)
    else:
        liabilities = pd.DataFrame()

    print(f"Loaded ranking table: {RANKING_FILE}")
    print(f"Number of ranked entries: {len(ranking)}")

    make_risk_score_plot(ranking)
    make_optimization_score_plot(ranking)
    make_liability_count_plot(ranking)
    make_cdr_hydrophobicity_plot(ranking)
    make_liability_type_summary(liabilities)
    save_top_variant_table(ranking)

    print("\nSaved figures:")
    print(f"  {OUT_RISK_SCORE}")
    print(f"  {OUT_OPT_SCORE}")
    print(f"  {OUT_LIABILITY_COUNT}")
    print(f"  {OUT_CDR_HYDRO}")
    print(f"  {OUT_LIABILITY_TYPES}")

    print("\nSaved summary table:")
    print(f"  {OUT_TOP_TABLE}")


if __name__ == "__main__":
    main()
