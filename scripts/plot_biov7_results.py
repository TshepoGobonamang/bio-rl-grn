"""Publication-style descriptive figures for the locked BIO-V7 outputs."""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

LABELS = {
    "expression_only": "Expression only", "fixed_sources": "Fixed sources",
    "global_source": "Global source gate", "no_graph_state": "No graph state",
    "no_parent_state": "No parent state", "no_action_state": "No action state",
    "gcev": "GCEV",
}
ORDER = ["expression_only", "fixed_sources", "global_source", "no_graph_state",
         "no_parent_state", "no_action_state", "gcev"]
NAVY, TEAL, ORANGE, GRAY = "#16324F", "#197C80", "#D47A35", "#698096"

def build(results: Path, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(results / "summary.csv")
    contrasts = pd.read_csv(results / "contrasts.csv")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": NAVY, "text.color": NAVY,
                         "xtick.color": NAVY, "ytick.color": NAVY,
                         "savefig.facecolor": "white", "pdf.fonttype": 42})
    primary = summary.loc[summary.primary_case].set_index("method").loc[ORDER]
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    y = np.arange(7)
    values = primary.macro_nmse_64.to_numpy()
    colors = [TEAL if m == "gcev" else GRAY for m in ORDER]
    ax.scatter(values, y, s=78, c=colors, zorder=3)
    ax.axvline(values[0], color=GRAY, linestyle="--", linewidth=1.0, alpha=.7)
    for value, yy, method in zip(values, y, ORDER):
        ax.annotate(f"{value:.6f}", (value, yy), xytext=(8, 0),
                    textcoords="offset points", va="center",
                    fontweight="bold" if method == "gcev" else "normal", fontsize=10)
    ax.set_yticks(y, [LABELS[m] for m in ORDER])
    ax.invert_yaxis()
    ax.set_ylim(6.6, -.6)
    ax.set_xlim(.4009, .40355)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.set_xlabel("Mean normalized squared prediction error (lower is better)")
    ax.grid(axis="x", color="#DFE5EB", linewidth=.6)
    fig.suptitle("BIO-V7 | Primary real K562 comparison", x=.05, ha="left",
                 fontsize=16, fontweight="bold", color=NAVY)
    fig.text(.05, .90, "Original log1p expression + original DNase processing; 131 internal-score cells",
             fontsize=10)
    fig.text(.05, .025, "Same 64 targets and 192-edge graph budget per method. Point estimates; no regulatory truth evaluated.",
             fontsize=9, color=GRAY)
    fig.subplots_adjust(left=.24, right=.97, top=.82, bottom=.16)
    for ext in ("png", "pdf"):
        fig.savefig(output / f"F01_BIOV7_Primary_Model_Comparison.{ext}", dpi=220)
    plt.close(fig)

    primary_contrasts = ["gcev_minus_no_graph_state", "gcev_minus_fixed_sources",
                         "gcev_minus_expression_only"]
    # Bracket indexing is deliberate: DataFrame.transform is a method.
    c = contrasts.loc[(contrasts["transform"] == "log1p_original") &
                      (contrasts["variant"] == "primary_mixed")].set_index("contrast").loc[primary_contrasts]
    sensitivity = contrasts.loc[contrasts.contrast == "gcev_minus_no_graph_state"].copy()
    condition_order = [(t,v) for t in ("log1p_original","log1p_cp10k")
                       for v in ("primary_mixed","matched_v1","matched_v2","motif_only")]
    sensitivity = sensitivity.set_index(["transform","variant"]).loc[condition_order]
    fig, axes = plt.subplots(1,2,figsize=(12.6,5.8),gridspec_kw={"width_ratios":[1,1.45]})
    for ax in axes:
        ax.axvline(0, color=NAVY, linewidth=1)
        ax.grid(axis="x", color="#E2E7EC", linewidth=.6)
        ax.set_xlabel("GCEV minus comparator NMSE")
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax=axes[0]
    for i,row in enumerate(c.itertuples()):
        ax.errorbar(row.delta_macro_nmse,i,
                    xerr=np.array([[row.delta_macro_nmse-row.cell_bootstrap_ci95_low],
                                   [row.cell_bootstrap_ci95_high-row.delta_macro_nmse]]),
                    fmt="o",color=TEAL,capsize=4,markersize=6,lw=1.5)
    ax.set_yticks(range(3),["No graph state","Fixed sources","Expression only"])
    ax.set_ylim(2.6,-.6)
    ax.set_title("A  Primary case: three paired contrasts",loc="left",fontsize=11,fontweight="bold",pad=14)
    ax=axes[1]
    for i,row in enumerate(sensitivity.itertuples()):
        color=TEAL if row.Index[0]=="log1p_original" else ORANGE
        ax.errorbar(row.delta_macro_nmse,i,
                    xerr=np.array([[row.delta_macro_nmse-row.cell_bootstrap_ci95_low],
                                   [row.cell_bootstrap_ci95_high-row.delta_macro_nmse]]),
                    fmt="o",color=color,capsize=3,markersize=5,lw=1.4)
    labels=["Original | mixed","Original | matched v1","Original | matched v2","Original | motif only",
            "CP10k | mixed","CP10k | matched v1","CP10k | matched v2","CP10k | motif only"]
    ax.set_yticks(range(8),labels,fontsize=9)
    ax.set_ylim(7.6,-.6)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.axhline(3.5,color="#D6DEE7",lw=.7)
    ax.set_title("B  GCEV versus no graph state: all cases",loc="left",fontsize=11,fontweight="bold",pad=14)
    fig.suptitle("BIO-V7 | Effect size and processing sensitivity",x=.04,ha="left",fontsize=16,fontweight="bold")
    fig.text(.04,.89,"Negative values favor GCEV. Every primary-case interval includes zero.",fontsize=10)
    fig.text(.04,.05,"95% paired cell-bootstrap intervals (10,000 resamples); conditional on this fitted cohort and these genes.",
             fontsize=9,color=GRAY)
    fig.text(.04,.018,"Processing variants reuse the same cells and biosamples. They are sensitivity analyses, not biological replications.",
             fontsize=9,color=GRAY)
    fig.subplots_adjust(left=.14,right=.98,wspace=.85,top=.79,bottom=.20)
    for ext in ("png","pdf"):
        fig.savefig(output / f"F02_BIOV7_Paired_Effects_and_Sensitivity.{ext}",dpi=220)
    plt.close(fig)

if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--results",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    build(a.results,a.output)

