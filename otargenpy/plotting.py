"""Visualization functions for otargenpy query outputs.

All functions accept a pandas DataFrame (output of a query function)
and return a matplotlib Figure.
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# plot_adverse_events
# ---------------------------------------------------------------------------

def plot_adverse_events(df: pd.DataFrame, top_n: int = 20):
    """Lollipop chart of adverse events ranked by log-likelihood ratio.

    Args:
        df: Output of :func:`~otargenpy.drug.adverse_events_query`.
        top_n: Max events to show.

    Returns:
        matplotlib Figure.
    """
    df = df.nlargest(top_n, "logLR").sort_values("logLR")
    crit = df["criticalValue"].iloc[0]

    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.35)))
    ax.hlines(y=df["name"], xmin=0, xmax=df["logLR"], color="grey", linewidth=0.6)
    scatter = ax.scatter(df["logLR"], df["name"], c=df["logLR"],
                         cmap="RdYlBu_r", s=40, zorder=3)
    ax.axvline(x=crit, color="firebrick", linestyle="--", linewidth=0.8,
               label=f"Critical value ({crit:.2f})")
    ax.set_xlabel("Log-likelihood ratio (logLR)")
    ax.set_title(f"Adverse Events (top {len(df)})")
    ax.legend(loc="lower right", fontsize=8)
    fig.colorbar(scatter, ax=ax, label="logLR", shrink=0.6)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# plot_interactions
# ---------------------------------------------------------------------------

def plot_interactions(df: pd.DataFrame, top_n: int = 20):
    """Circular network graph of protein interaction partners.

    Args:
        df: Output of :func:`~otargenpy.gene.interactions_query`.
        top_n: Max interactions to show.

    Returns:
        matplotlib Figure.
    """
    df = df.nlargest(top_n, "score")
    nodes = list(dict.fromkeys(
        list(df["targetA.approvedSymbol"]) + list(df["targetB.approvedSymbol"])
    ))
    n = len(nodes)
    angles = [2 * math.pi * i / n for i in range(n)]
    pos = {node: (math.cos(a), math.sin(a)) for node, a in zip(nodes, angles)}

    fig, ax = plt.subplots(figsize=(8, 8))
    # Edges
    for _, row in df.iterrows():
        a = pos[row["targetA.approvedSymbol"]]
        b = pos[row["targetB.approvedSymbol"]]
        ax.plot([a[0], b[0]], [a[1], b[1]], color="grey",
                alpha=0.4, linewidth=row["score"] * 3)
    # Nodes
    xs = [pos[n][0] for n in nodes]
    ys = [pos[n][1] for n in nodes]
    ax.scatter(xs, ys, s=120, color="steelblue", zorder=3)
    for node in nodes:
        x, y = pos[node]
        ax.annotate(node, (x * 1.12, y * 1.12), ha="center", va="center", fontsize=8)
    ax.set_aspect("equal")
    ax.set_title("Protein Interaction Network")
    ax.axis("off")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# plot_l2g
# ---------------------------------------------------------------------------

def plot_l2g(df: pd.DataFrame, top_n: int = 15):
    """Horizontal bar chart of locus-to-gene prediction scores.

    Args:
        df: Output of :func:`~otargenpy.genetics.locus2gene_query`.
        top_n: Max genes to show.

    Returns:
        matplotlib Figure.
    """
    df = df.nlargest(top_n, "score").sort_values("score")

    fig, ax = plt.subplots(figsize=(7, max(3, len(df) * 0.4)))
    colors = plt.cm.Blues(np.linspace(0.3, 1.0, len(df)))
    ax.barh(df["target.approvedSymbol"], df["score"], color=colors, height=0.6)
    ax.set_xlabel("L2G Score")
    ax.set_title("Locus-to-Gene (L2G) Predictions")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# plot_colocalisation
# ---------------------------------------------------------------------------

def plot_colocalisation(df: pd.DataFrame, h4_threshold: float = 0.8):
    """Scatter plot of H4 posterior vs number of colocalising variants.

    Args:
        df: Output of :func:`~otargenpy.genetics.gwas_colocalisation`.
        h4_threshold: Dashed threshold line for H4.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(9, 6))
    scatter = ax.scatter(df["numberColocalisingVariants"], df["h4"],
                         c=df["h4"], cmap="plasma", s=50, alpha=0.8, zorder=3)
    ax.axhline(y=h4_threshold, color="firebrick", linestyle="--", linewidth=0.7,
               label=f"H4 threshold ({h4_threshold})")
    # Label points (truncate long trait names)
    for _, row in df.iterrows():
        trait = row.get("study.traitReported", "")
        label = (trait[:32] + "...") if len(str(trait)) > 35 else str(trait)
        ax.annotate(label, (row["numberColocalisingVariants"], row["h4"]),
                    fontsize=6, alpha=0.7, xytext=(4, 4),
                    textcoords="offset points")
    ax.set_xlabel("Number of colocalising variants")
    ax.set_ylabel("H4 posterior")
    ax.set_title("GWAS Colocalisation")
    ax.legend(fontsize=8)
    fig.colorbar(scatter, ax=ax, label="H4", shrink=0.6)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# plot_indications
# ---------------------------------------------------------------------------

# Clinical stage ordering and display labels
_STAGE_MAP = {
    "IND": "IND",
    "EARLY_PHASE_1": "Early Phase 1",
    "PHASE_1": "Phase 1",
    "PHASE_1_2": "Phase 1/2",
    "PHASE_2": "Phase 2",
    "PHASE_2_3": "Phase 2/3",
    "PHASE_3": "Phase 3",
    "APPROVAL": "Approved",
}


def plot_indications(df: pd.DataFrame, top_n: int = 10):
    """Faceted bar chart of drug indications grouped by clinical stage.

    Args:
        df: Output of :func:`~otargenpy.drug.indications_query`.
        top_n: Max diseases per stage panel.

    Returns:
        matplotlib Figure.
    """
    df = df.copy()
    df["stage_label"] = df["maxClinicalStage"].map(_STAGE_MAP).fillna(df["maxClinicalStage"])
    # Order from stage map
    ordered = [v for v in _STAGE_MAP.values() if v in df["stage_label"].values]
    if not ordered:
        ordered = df["stage_label"].unique().tolist()

    # Color ramp yellow -> green
    n_stages = len(ordered)
    cmap = plt.cm.YlGn(np.linspace(0.3, 0.9, n_stages))
    color_map = dict(zip(ordered, cmap))

    ncols = 2
    nrows = math.ceil(n_stages / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, max(3, nrows * 3)))
    axes = axes.flatten() if n_stages > 1 else [axes]

    for i, stage in enumerate(ordered):
        ax = axes[i]
        sub = df[df["stage_label"] == stage].head(top_n).sort_values("disease.name")
        ax.barh(sub["disease.name"], [1] * len(sub), color=color_map[stage], height=0.6)
        ax.set_title(stage, fontweight="bold", fontsize=10)
        ax.set_xlim(0, 1.1)
        ax.get_xaxis().set_visible(False)
        ax.tick_params(axis="y", labelsize=8)

    # Hide unused panels
    for j in range(n_stages, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Drug Indications by Clinical Stage", fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig
