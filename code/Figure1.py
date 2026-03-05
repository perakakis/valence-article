"""
Figure 1: Convergence of P2N as a robust predictor across methods.

Three panels:
  A) LASSO selection frequency — summary bars + outcome-group dots
  B) R² — three-model grouped bars (Solo, over P2N, over P2N+N2P) + group dots
  C) Relative importance (LMG) — summary bars + outcome-group dots

Dots represent outcome-group averages (weighted by sample size):
  Depression (DASS-D, PHQ-9), Anxiety (DASS-A, GAD-7),
  Life Satisfaction (SWLS, SWLS*), and individual outcomes (AAQ, BRS, FS).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =============================================================================
# Configuration
# =============================================================================
predictor_order = ["M_PosA", "M_NegA", "SD_PosA", "SD_NegA", "N2P", "P2N"]
predictor_labels = ["mPA", "mNA", "sdPA", "sdNA", "N2P", "P2N"]
n_pred = len(predictor_order)

# Outcome groups: which outcomes belong to which group
# Grouped categories get averaged; singletons stay as-is
outcome_to_group = {
    "DASSd": "Depression",
    "PHQ":   "Depression",
    "DASSa": "Anxiety",
    "GAD":   "Anxiety",
    "SWL":   "Life Sat.",
    "SWL1":  "Life Sat.",
    "AAQ":   "AAQ",
    "BRS":   "BRS",
    "FS":    "FS",
}

# Visual encoding per group: (color, marker, size)
group_style = {
    "Depression": ("#c0392b", "o",  32),   # red circles
    "Anxiety":    ("#e67e22", "D",  28),   # orange diamonds
    "Life Sat.":  ("#27ae60", "^",  32),   # green triangles
    "AAQ":        ("#2980b9", "s",  26),   # blue squares
    "BRS":        ("#8e44ad", "P",  32),   # purple plus
    "FS":         ("#16a085", "X",  30),   # teal X
}

legend_order = ["Depression", "Anxiety", "Life Sat.", "AAQ", "BRS", "FS"]

# R² model colors (matching MATLAB / LinearRegression.R)
R2_COLORS = {
    "solo": "#0072BD",
    "over_bench": "#D95319",
    "over_ext": "#EDB120",
}

# =============================================================================
# Load data
# =============================================================================
lasso = pd.read_csv("results/LassoSelectionFreqIndividual.csv")
r2_avg_csv = pd.read_csv("results/R2Avg.csv")
r2_solo_csv = pd.read_csv("results/R2.csv")
relimp_raw = pd.read_csv("results/RelativeImportance.csv")

# =============================================================================
# Helper: compute group-level averages (weighted by sample size)
# =============================================================================
def compute_group_dots(per_outcome_df, predictor_col, value_col, weight_col=None):
    """
    Given per-outcome data, compute weighted averages within each outcome group.
    Returns dict: {group_name: {predictor: value}}
    """
    per_outcome_df = per_outcome_df.copy()
    per_outcome_df["_group"] = per_outcome_df["outcome"].map(outcome_to_group)
    per_outcome_df = per_outcome_df.dropna(subset=["_group"])

    result = {}
    for grp in per_outcome_df["_group"].unique():
        grp_data = per_outcome_df[per_outcome_df["_group"] == grp]
        pred_vals = {}
        for pred in predictor_order:
            sub = grp_data[grp_data[predictor_col] == pred]
            if len(sub) == 0:
                continue
            if weight_col and weight_col in sub.columns:
                pred_vals[pred] = np.average(sub[value_col], weights=sub[weight_col])
            else:
                pred_vals[pred] = sub[value_col].mean()
        result[grp] = pred_vals
    return result


# =============================================================================
# Panel A: LASSO selection frequency
# =============================================================================
# Overall weighted average (bar heights)
lasso_summary = (
    lasso.groupby("predictor")
    .apply(lambda g: pd.Series({
        "mean": np.average(g["weighted_selection"], weights=g["n_total"]),
        "se": np.sqrt(np.average(
            (g["weighted_selection"] - np.average(g["weighted_selection"], weights=g["n_total"]))**2,
            weights=g["n_total"]) / len(g)) if len(g) > 1 else 0
    }))
    .reindex(predictor_order)
)

# Group dots
lasso_group_dots = compute_group_dots(
    lasso.rename(columns={"weighted_selection": "val"}),
    "predictor", "val", "n_total"
)

# =============================================================================
# Panel B: R² three-model comparison
# =============================================================================
# Bars: overall "Well-being" from R2Avg.csv
r2_wb = r2_avg_csv[r2_avg_csv["outcome_type"] == "Well-being"].copy()
r2_solo_mod = r2_wb[r2_wb["model"] == "solo"].set_index("predictor")
r2_ob_mod = r2_wb[r2_wb["model"] == "over_bench"].set_index("predictor")
r2_oe_mod = r2_wb[r2_wb["model"] == "over_ext"].set_index("predictor")

# Group dots: Solo R² per outcome, averaged across datasets, then within groups
r2_solo_long = r2_solo_csv.melt(
    id_vars=["dataset", "outcome", "ind"],
    value_vars=["M_PosA", "M_NegA", "SD_PosA", "SD_NegA", "P2N", "N2P"],
    var_name="predictor", value_name="r2"
)
# First: weighted average across datasets per outcome (weight = 1 per dataset here)
r2_per_outcome = r2_solo_long.groupby(["outcome", "predictor"])["r2"].mean().reset_index()
r2_group_dots = compute_group_dots(r2_per_outcome, "predictor", "r2")

# =============================================================================
# Panel C: Relative Importance (LMG)
# =============================================================================
# Per outcome (weighted across datasets)
relimp_per_outcome = (
    relimp_raw.groupby(["outcome", "predictor"])
    .apply(lambda g: pd.Series({
        "lmg": np.average(g["lmg_value"], weights=g["n_obs"]),
        "n_total": g["n_obs"].sum(),
    }))
    .reset_index()
)

# Overall weighted average (bar heights)
relimp_summary = (
    relimp_raw.groupby("predictor")
    .apply(lambda g: pd.Series({
        "mean": np.average(g["lmg_value"], weights=g["n_obs"]),
        "se": np.sqrt(np.average(
            (g["lmg_value"] - np.average(g["lmg_value"], weights=g["n_obs"]))**2,
            weights=g["n_obs"]) / len(g)) if len(g) > 1 else 0
    }))
    .reindex(predictor_order)
)

# Group dots
relimp_group_dots = compute_group_dots(
    relimp_per_outcome, "predictor", "lmg", "n_total"
)

# =============================================================================
# Helper: plot group dots on an axis
# =============================================================================
np.random.seed(42)

def plot_group_dots(ax, group_dots, x_offset=0.0, spread=0.28):
    """Plot one dot per outcome group, jittered symmetrically."""
    groups_present = [g for g in legend_order if g in group_dots]
    n_groups = len(groups_present)
    for pred_idx, pred in enumerate(predictor_order):
        for g_idx, grp in enumerate(groups_present):
            if pred not in group_dots[grp]:
                continue
            val = group_dots[grp][pred]
            color, marker, ms = group_style[grp]
            # Spread dots evenly across the bar width
            offset = (g_idx - (n_groups - 1) / 2) * (spread / max(n_groups - 1, 1))
            ax.scatter(
                pred_idx + x_offset + offset, val,
                s=ms, alpha=0.75, color=color, marker=marker,
                edgecolors="white", linewidths=0.4, zorder=6
            )

# =============================================================================
# Figure
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=False)
plt.subplots_adjust(wspace=0.32, bottom=0.18)

x = np.arange(n_pred)
bar_width = 0.55

def get_bar_colors(predictors):
    return ["#2c3e50" if p == "P2N" else "#b0b8c1" for p in predictors]

# =============================================================================
# Panel A: LASSO Selection Frequency
# =============================================================================
ax = axes[0]
means_a = lasso_summary["mean"].values
ses_a = lasso_summary["se"].values

ax.bar(x, means_a, width=bar_width, color=get_bar_colors(predictor_order),
       edgecolor="white", linewidth=0.5, zorder=3)
ax.errorbar(x, means_a, yerr=ses_a, fmt="none", ecolor="#333333",
            capsize=3, capthick=1, linewidth=1, zorder=4)
plot_group_dots(ax, lasso_group_dots)

ax.set_ylabel("Selection Frequency", fontsize=10)
ax.set_ylim(0, 1.12)
ax.set_title("A   LASSO Selection Frequency", fontsize=11, fontweight="bold",
             loc="left", pad=8)

# =============================================================================
# Panel B: R² (three-model grouped bars + group dots on Solo)
# =============================================================================
ax = axes[1]
bar_w = 0.22
offsets = [-bar_w, 0, bar_w]
model_configs = [
    ("solo",       r2_solo_mod, "weighted_R2",      "weighted_R2_se",      predictor_order),
    ("over_bench", r2_ob_mod,   "weighted_R2_over", "weighted_R2_over_se",
     [p for p in predictor_order if p != "P2N"]),
    ("over_ext",   r2_oe_mod,   "weighted_R2_over", "weighted_R2_over_se",
     [p for p in predictor_order if p not in ("P2N", "N2P")]),
]
model_labels = ["Solo", "Over P2N", "Over P2N & N2P"]

for m_idx, (model_key, data_src, val_col, se_col, valid_preds) in enumerate(model_configs):
    vals, errs, positions = [], [], []
    for p_idx, pred in enumerate(predictor_order):
        if pred in valid_preds and pred in data_src.index:
            vals.append(data_src.loc[pred, val_col])
            errs.append(data_src.loc[pred, se_col])
            positions.append(p_idx + offsets[m_idx])

    ax.bar(positions, vals, width=bar_w, color=R2_COLORS[model_key],
           edgecolor="white", linewidth=0.5, label=model_labels[m_idx], zorder=3)
    ax.errorbar(positions, vals, yerr=errs, fmt="none", ecolor="#333333",
                capsize=2, capthick=0.8, linewidth=0.8, zorder=4)

# Group dots overlaid on the Solo bar positions
plot_group_dots(ax, r2_group_dots, x_offset=offsets[0], spread=0.32)

ax.set_ylabel("R\u00b2", fontsize=10)
ax.set_title("B   Variance Explained (R\u00b2)", fontsize=11, fontweight="bold",
             loc="left", pad=8)
ax.legend(fontsize=7, loc="upper left", framealpha=0.92,
          handlelength=1.2, handletextpad=0.4)

# =============================================================================
# Panel C: Relative Importance (LMG)
# =============================================================================
ax = axes[2]
means_c = relimp_summary["mean"].values
ses_c = relimp_summary["se"].values

ax.bar(x, means_c, width=bar_width, color=get_bar_colors(predictor_order),
       edgecolor="white", linewidth=0.5, zorder=3)
ax.errorbar(x, means_c, yerr=ses_c, fmt="none", ecolor="#333333",
            capsize=3, capthick=1, linewidth=1, zorder=4)
plot_group_dots(ax, relimp_group_dots)

ax.set_ylabel("Relative Importance (LMG)", fontsize=10)
ax.set_title("C   Relative Importance", fontsize=11, fontweight="bold",
             loc="left", pad=8)

# =============================================================================
# Shared formatting
# =============================================================================
for ax in axes:
    ax.set_xticks(x)
    ax.set_xticklabels(predictor_labels, fontsize=8.5, rotation=45, ha="right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.set_xlim(-0.5, n_pred - 0.5)
    ax.yaxis.grid(True, alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)

# =============================================================================
# Shared outcome-group legend (below figure)
# =============================================================================
legend_elements = []
for grp in legend_order:
    color, marker, ms = group_style[grp]
    legend_elements.append(
        Line2D([0], [0], marker=marker, color="w", markerfacecolor=color,
               markeredgecolor="white", markeredgewidth=0.4,
               markersize=6.5, label=grp, linestyle="None")
    )

fig.legend(
    handles=legend_elements, ncol=len(legend_order),
    loc="lower center", bbox_to_anchor=(0.5, -0.01),
    fontsize=8, frameon=True, framealpha=0.9,
    handletextpad=0.3, columnspacing=1.2,
    title="Outcome groups", title_fontsize=8.5
)

plt.savefig("plots/Figure1.pdf", dpi=300, bbox_inches="tight", pad_inches=0.35)
plt.savefig("plots/Figure1.png", dpi=300, bbox_inches="tight", pad_inches=0.35)
print("Figure saved to plots/Figure1.pdf and plots/Figure1.png")
plt.close()
