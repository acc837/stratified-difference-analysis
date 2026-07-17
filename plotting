from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from analysis_core import Genie3Result, StratificationResult


def plot_stratification(result: StratificationResult):
    frame = pd.DataFrame(
        {
            "sample": result.expression.index,
            "expression": result.expression.values,
            "group": result.groups.reindex(result.expression.index).values,
        }
    ).sort_values("expression")

    color_map = {"Low": "#3B82F6", "High": "#DC2626", "Excluded": "#C7CDD4"}
    colors = frame["group"].map(color_map)

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.bar(np.arange(len(frame)), frame["expression"], color=colors, width=0.9)
    ax.axhline(result.low_cutoff, linestyle="--", linewidth=1, color="#3B82F6")
    ax.axhline(result.high_cutoff, linestyle="--", linewidth=1, color="#DC2626")
    ax.set_xlabel("Samples ordered by target-gene expression")
    ax.set_ylabel(f"{result.target_gene} normalized expression")
    ax.set_title(f"{result.target_gene}: high/low sample stratification")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def plot_pair_scatter(expression_sg: pd.DataFrame, gene_a: str, gene_b: str):
    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    ax.scatter(expression_sg[gene_a], expression_sg[gene_b], s=38, alpha=0.8)
    ax.set_xlabel(f"{gene_a} expression used by GENIE3")
    ax.set_ylabel(f"{gene_b} expression used by GENIE3")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def plot_target_network(result: Genie3Result, target_gene: str, top_n_each: int = 8):
    edges = result.edges
    selected = pd.concat(
        [
            edges.loc[edges["regulator"] == target_gene].head(top_n_each),
            edges.loc[edges["target"] == target_gene].head(top_n_each),
        ],
        ignore_index=True,
    ).drop_duplicates(["regulator", "target"])

    graph = nx.DiGraph()
    for row in selected.itertuples(index=False):
        graph.add_edge(row.regulator, row.target, weight=row.weight)

    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    if graph.number_of_nodes() == 0:
        ax.text(0.5, 0.5, "No target-related edges", ha="center", va="center")
        ax.axis("off")
        return fig

    pos = nx.spring_layout(graph, seed=2026, k=1.25)
    node_sizes = [1500 if node == target_gene else 850 for node in graph.nodes]
    node_colors = ["#1D4ED8" if node == target_gene else "#E8EEF6" for node in graph.nodes]
    widths = [0.8 + 5.0 * graph[u][v]["weight"] for u, v in graph.edges]

    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors="#1F2937",
        linewidths=0.8,
        ax=ax,
    )
    # Draw labels separately to allow a contrasting label color for the target node.
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors="#1F2937",
        linewidths=0.8,
        ax=ax,
    )
    nx.draw_networkx_edges(
        graph,
        pos,
        width=widths,
        edge_color="#667085",
        arrows=True,
        arrowsize=15,
        connectionstyle="arc3,rad=0.05",
        ax=ax,
    )
    for node, (x, y) in pos.items():
        ax.text(
            x,
            y,
            node,
            ha="center",
            va="center",
            fontsize=9,
            color="white" if node == target_gene else "black",
            fontweight="bold" if node == target_gene else "normal",
        )
    ax.set_title(f"Top predicted edges connected with {target_gene}")
    ax.axis("off")
    fig.tight_layout()
    return fig
