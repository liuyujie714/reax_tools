"""Transfer-flow graph construction and plotting."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors
from matplotlib.patches import FancyArrowPatch, Rectangle
import networkx as nx
import pandas as pd

from .chemistry import formula_to_subscript


def build_transfer_graph(flow: pd.DataFrame, label: str = "formula") -> nx.DiGraph:
    graph = nx.DiGraph()
    for _, row in flow.iterrows():
        source = str(row["source_id"])
        target = str(row["target_id"])
        source_label = row["source_label"] if label == "formula" else source
        target_label = row["target_label"] if label == "formula" else target
        graph.add_node(source, label=formula_to_subscript(source_label))
        graph.add_node(target, label=formula_to_subscript(target_label))
        graph.add_edge(
            source,
            target,
            count=int(row["count"]),
            atom_transfer=int(row["atom_transfer"]),
        )
    return graph


def _layout_graph(graph: nx.DiGraph, layout: str) -> dict:
    layout = layout.lower()
    undirected = graph.to_undirected()
    n_nodes = max(1, graph.number_of_nodes())
    if layout in {"spring", "fdp"}:
        return nx.spring_layout(graph, seed=42, k=max(1.1, math.sqrt(n_nodes) / 2.4), iterations=220)
    if layout in {"kamada", "kk"}:
        pos = nx.kamada_kawai_layout(undirected)
        return _expand_layout(pos, max(1.0, math.sqrt(n_nodes) / 4.0))
    if layout in {"circular", "circle"}:
        return nx.circular_layout(graph)
    if layout == "shell":
        return nx.shell_layout(graph)
    if layout == "spectral":
        return nx.spectral_layout(undirected)
    if layout in {"layered", "dot"}:
        return _layered_layout(graph)
    return nx.kamada_kawai_layout(undirected)


def _layered_layout(graph: nx.DiGraph) -> dict:
    condensed = nx.condensation(graph)
    layers = {}
    for node in nx.topological_sort(condensed):
        preds = list(condensed.predecessors(node))
        layers[node] = 0 if not preds else max(layers[pred] + 1 for pred in preds)
    component_members = condensed.graph.get("mapping", {})
    by_layer: dict[int, list[str]] = {}
    for original, component in component_members.items():
        by_layer.setdefault(layers[component], []).append(original)
    pos = {}
    for x, nodes in by_layer.items():
        nodes = sorted(nodes)
        offset = (len(nodes) - 1) / 2
        for i, node in enumerate(nodes):
            pos[node] = (float(x) * 3.2, float(offset - i) * 1.25)
    return pos or nx.spring_layout(graph, seed=42)


def _expand_layout(pos: dict, factor: float) -> dict:
    return {node: (coords[0] * factor, coords[1] * factor) for node, coords in pos.items()}


def _auto_figsize(graph: nx.DiGraph, pos: dict, config: dict) -> tuple[float, float]:
    base_width = float(config.get("figure_width", 9.0))
    base_height = float(config.get("figure_height", 7.0))
    n_nodes = max(1, graph.number_of_nodes())
    xs = [coords[0] for coords in pos.values()]
    ys = [coords[1] for coords in pos.values()]
    span_x = max(xs) - min(xs) if xs else 1.0
    span_y = max(ys) - min(ys) if ys else 1.0
    width = max(base_width, min(18.0, 2.0 + span_x * 1.35, 6.0 + math.sqrt(n_nodes) * 1.15))
    height = max(base_height, min(18.0, 2.0 + span_y * 0.45, 5.0 + math.sqrt(n_nodes) * 0.95))
    return width, height


def _scale(values: list[float], low: float, high: float) -> list[float]:
    if not values:
        return []
    min_v = min(values)
    max_v = max(values)
    if max_v == min_v:
        return [(low + high) / 2 for _ in values]
    return [low + (value - min_v) * (high - low) / (max_v - min_v) for value in values]


def _node_box_size(pos: dict, config: dict) -> tuple[float, float]:
    xs = [coords[0] for coords in pos.values()]
    ys = [coords[1] for coords in pos.values()]
    span_x = max(xs) - min(xs) if xs else 1.0
    span_y = max(ys) - min(ys) if ys else 1.0
    height = float(config.get("node_box_height", 0.62))
    if span_y <= 3.0:
        height = max(height, span_y * 0.08)
    width = height * float(config.get("node_box_aspect", 2.0))
    if span_x <= 3.0:
        width = max(width, span_x * 0.12)
    return width, height


def _edge_boundary_points(source_pos, target_pos, box_w: float, box_h: float) -> tuple[tuple[float, float], tuple[float, float]]:
    sx, sy = source_pos
    tx, ty = target_pos
    dx = tx - sx
    dy = ty - sy
    if dx == 0 and dy == 0:
        return (sx, sy), (tx, ty)

    def boundary(cx, cy, sign):
        candidates = []
        if dx != 0:
            t = (box_w / 2) / abs(dx)
            candidates.append(t)
        if dy != 0:
            t = (box_h / 2) / abs(dy)
            candidates.append(t)
        t = min(candidates) if candidates else 0
        pad = 0.10
        return cx + sign * dx * (t + pad * t), cy + sign * dy * (t + pad * t)

    return boundary(sx, sy, 1), boundary(tx, ty, -1)


def _draw_top_reactions(ax, graph: nx.DiGraph, top_n: int) -> None:
    ax.axis("off")
    edges = sorted(graph.edges(data=True), key=lambda item: int(item[2]["count"]), reverse=True)[:top_n]
    if not edges:
        return
    ax.text(0.0, 0.98, f"Top {len(edges)} transfers", fontsize=12, fontweight="bold", va="top")
    for idx, (source, target, data) in enumerate(edges, start=1):
        src = graph.nodes[source].get("label", source)
        tgt = graph.nodes[target].get("label", target)
        text = f"{idx}. {src} -> {tgt}: {int(data['count'])}"
        x = 0.0 if idx <= math.ceil(top_n / 2) else 0.52
        y = 0.82 - ((idx - 1) % math.ceil(top_n / 2)) * 0.16
        ax.text(x, y, text, fontsize=12, va="top")


def draw_transfer_graph(
    graph: nx.DiGraph,
    output: str | Path,
    *,
    layout: str = "kamada",
    config: dict | None = None,
    use_colors: bool = True,
) -> None:
    if graph.number_of_nodes() == 0:
        raise ValueError("Cannot draw an empty transfer graph")
    config = config or {}
    pos = _layout_graph(graph, layout)
    labels = nx.get_node_attributes(graph, "label")
    box_w, box_h = _node_box_size(pos, config)

    edge_counts = [float(graph.edges[e]["count"]) for e in graph.edges]
    widths = _scale(
        edge_counts,
        float(config.get("edge_width_min", 0.9)),
        float(config.get("edge_width_max", 2.7)),
    )
    node_strength = {
        node: sum(float(data["count"]) for _, _, data in graph.in_edges(node, data=True))
        + sum(float(data["count"]) for _, _, data in graph.out_edges(node, data=True))
        for node in graph.nodes
    }
    node_size_high = float(config.get("node_size_max", 1800))
    if graph.number_of_nodes() > 25:
        node_size_high *= max(0.45, 25 / graph.number_of_nodes())
    node_size_low = min(float(config.get("node_size_min", 900)), node_size_high * 0.65)
    node_sizes = _scale(
        [node_strength[node] for node in graph.nodes],
        node_size_low,
        node_size_high,
    )

    cmap = plt.get_cmap(str(config.get("color_map", "YlGnBu")))
    if use_colors:
        norm = mcolors.Normalize(vmin=min(edge_counts), vmax=max(edge_counts)) if edge_counts else mcolors.Normalize(0, 1)
        edge_colors = [cmap(0.28 + 0.62 * norm(value)) for value in edge_counts]
        node_norm = _scale([node_strength[node] for node in graph.nodes], 0.22, 0.78)
        node_colors = [cmap(value) for value in node_norm]
    else:
        norm = mcolors.Normalize(vmin=min(edge_counts), vmax=max(edge_counts)) if edge_counts else mcolors.Normalize(0, 1)
        edge_colors = ["#555555" for _ in edge_counts]
        node_colors = ["#eeeeee" for _ in graph.nodes]

    show_top = int(config.get("top_reactions", 0) or 0)
    if show_top > 0:
        fig = plt.figure(figsize=(
            _auto_figsize(graph, pos, config)[0],
            _auto_figsize(graph, pos, config)[1] + float(config.get("top_reactions_height", 1.35)),
        ))
        gs = fig.add_gridspec(2, 1, height_ratios=[_auto_figsize(graph, pos, config)[1], float(config.get("top_reactions_height", 1.35))])
        ax = fig.add_subplot(gs[0, 0])
        list_ax = fig.add_subplot(gs[1, 0])
    else:
        fig, ax = plt.subplots(figsize=_auto_figsize(graph, pos, config))
        list_ax = None

    for (node, (x, y)), color in zip(pos.items(), node_colors):
        rect = Rectangle(
            (x - box_w / 2, y - box_h / 2),
            box_w,
            box_h,
            facecolor=color,
            edgecolor="#333333",
            linewidth=0.8,
            alpha=float(config.get("alpha", 0.78)),
            zorder=2,
        )
        ax.add_patch(rect)
    label_font = 12 if graph.number_of_nodes() <= 35 else 10
    for node, (x, y) in pos.items():
        ax.text(x, y, labels.get(node, node), fontsize=label_font, ha="center", va="center", zorder=3)

    for idx, (source, target, data) in enumerate(graph.edges(data=True)):
        start, end = _edge_boundary_points(pos[source], pos[target], box_w, box_h)
        arrow = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=widths[idx],
            color=edge_colors[idx],
            alpha=0.85,
            connectionstyle="arc3,rad=0.08",
            zorder=1,
        )
        ax.add_patch(arrow)
        if bool(config.get("show_edge_labels", False)):
            mx = (start[0] + end[0]) / 2
            my = (start[1] + end[1]) / 2
            ax.text(mx, my, str(int(data["count"])), fontsize=12, ha="center", va="center", zorder=4)

    if use_colors and bool(config.get("show_colorbar", True)) and edge_counts:
        scalar = cm.ScalarMappable(norm=norm, cmap=cmap)
        scalar.set_array(edge_counts)
        cbar = fig.colorbar(scalar, ax=ax, fraction=0.032, pad=0.02)
        cbar.set_label("Transfer count", fontsize=14)
        cbar.ax.tick_params(labelsize=12)
    if list_ax is not None:
        _draw_top_reactions(list_ax, graph, show_top)
    ax.axis("off")
    xs = [coords[0] for coords in pos.values()]
    ys = [coords[1] for coords in pos.values()]
    ax.set_xlim(min(xs) - box_w * 1.2, max(xs) + box_w * 1.2)
    ax.set_ylim(min(ys) - box_h * 1.5, max(ys) + box_h * 1.5)
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)
