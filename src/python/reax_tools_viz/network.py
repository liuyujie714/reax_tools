"""Transfer-flow graph construction and Graphviz rendering."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
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
            source_label=formula_to_subscript(source_label),
            target_label=formula_to_subscript(target_label),
        )
    return graph


def _scale(value: float, values: list[float], low: float, high: float) -> float:
    if not values:
        return (low + high) / 2.0
    min_v = min(values)
    max_v = max(values)
    if min_v == max_v:
        return (low + high) / 2.0
    return low + (value - min_v) * (high - low) / (max_v - min_v)


def _rgba_hex(color: tuple[float, float, float, float], alpha: float) -> str:
    r, g, b, _ = color
    return "#{:02x}{:02x}{:02x}{:02x}".format(
        round(r * 255),
        round(g * 255),
        round(b * 255),
        round(max(0.0, min(1.0, alpha)) * 255),
    )


def _edge_color(value: float, values: list[float], cfg: dict[str, Any]) -> str:
    cmap = matplotlib.colormaps[str(cfg.get("color_map", "YlGnBu"))]
    low = float(cfg.get("color_low", 0.2))
    high = float(cfg.get("color_high", 1.0))
    fraction = _scale(value, values, low, high)
    return _rgba_hex(cmap(fraction), float(cfg.get("alpha", 0.75)))


def _node_strength(graph: nx.DiGraph, node: str, weight: str) -> float:
    return sum(float(data.get(weight, 0.0)) for _, _, data in graph.in_edges(node, data=True)) + sum(
        float(data.get(weight, 0.0)) for _, _, data in graph.out_edges(node, data=True)
    )


def _graphviz_attrs(graph_cfg: dict[str, Any]) -> dict[str, str]:
    return {str(k): str(v) for k, v in graph_cfg.items() if v is not None and str(v) != ""}


def _node_attrs(graph: nx.DiGraph, node: str, strengths: list[float], config: dict[str, Any], weight: str) -> dict[str, str]:
    node_cfg = dict(config.get("nodes", {}))
    strength = _node_strength(graph, node, weight)
    width = _scale(strength, strengths, float(node_cfg.get("width_min", 0.7)), float(node_cfg.get("width_max", 2.0)))
    height = _scale(strength, strengths, float(node_cfg.get("height_min", 0.35)), float(node_cfg.get("height_max", 0.88)))
    return {
        "label": str(graph.nodes[node].get("label", node)) if bool(node_cfg.get("label_enabled", True)) else "",
        "shape": str(node_cfg.get("shape", "ellipse")),
        "style": str(node_cfg.get("style", "filled")),
        "fixedsize": str(node_cfg.get("fixedsize", "false")),
        "width": f"{width:.3f}",
        "height": f"{height:.3f}",
        "fillcolor": str(node_cfg.get("fillcolor_fixed", "#9ed9c7")),
        "color": str(node_cfg.get("outline_color", "#263238")),
        "penwidth": str(node_cfg.get("outline_width", 1.4)),
        "fontname": str(node_cfg.get("fontname", "Arial")),
        "fontsize": str(node_cfg.get("fontsize", 20)),
        "fontcolor": str(node_cfg.get("fontcolor", "#111111")),
    }


def _edge_attrs(
    source: str,
    target: str,
    data: dict[str, Any],
    values: list[float],
    config: dict[str, Any],
    weight: str,
    labelled_edges: set[tuple[str, str]],
) -> dict[str, str]:
    edge_cfg = dict(config.get("edges", {}))
    label_cfg = dict(config.get("edge_labels", {}))
    value = float(data.get(weight, 0.0))
    attrs = {
        "penwidth": f"{_scale(value, values, float(edge_cfg.get('width_min', 2.0)), float(edge_cfg.get('width_max', 6.0))):.3f}",
        "color": _edge_color(value, values, edge_cfg),
        "arrowsize": f"{_scale(value, values, float(edge_cfg.get('arrowsize_min', 1.0)), float(edge_cfg.get('arrowsize_max', 2.0))):.3f}",
        "arrowhead": str(edge_cfg.get("arrowhead", "vee")),
        "style": str(edge_cfg.get("style", "solid")),
        "constraint": str(edge_cfg.get("constraint", "true")),
        "weight": f"{_scale(value, values, float(edge_cfg.get('weight_graphviz_min', 1.0)), float(edge_cfg.get('weight_graphviz_max', 10.0))):.3f}",
    }
    if bool(label_cfg.get("enabled", False)) and (source, target) in labelled_edges:
        label_value = float(data.get(str(label_cfg.get("field", weight)), value))
        attrs.update(
            label=str(label_cfg.get("format", "{value:.0f}")).format(value=label_value),
            fontname=str(label_cfg.get("fontname", "Arial")),
            fontsize=str(label_cfg.get("fontsize", 14)),
            fontcolor=str(label_cfg.get("fontcolor", "#333333")),
            decorate=str(label_cfg.get("decorate", "false")),
            labeldistance=str(label_cfg.get("labeldistance", "1.6")),
            labelangle=str(label_cfg.get("labelangle", "0")),
        )
    return attrs


def _build_agraph(graph: nx.DiGraph, *, layout: str, config: dict[str, Any], title: str | None) -> Any:
    try:
        import pygraphviz as pgv
    except ImportError as exc:
        raise RuntimeError("pygraphviz is required for ReaxTools transfer diagrams") from exc

    graph_cfg = dict(config.get("graph", {}))
    title_cfg = dict(config.get("title", {}))
    weight = str(config.get("weight", "atom_transfer"))
    agraph = pgv.AGraph(
        strict=bool(graph_cfg.pop("strict", False)),
        directed=bool(graph_cfg.pop("directed", True)),
    )
    agraph.graph_attr.update(_graphviz_attrs(graph_cfg))
    if config.get("dpi"):
        agraph.graph_attr.update(dpi=str(config["dpi"]))
    if title and bool(title_cfg.get("enabled", True)):
        agraph.graph_attr.update(
            label=title,
            labelloc=str(title_cfg.get("labelloc", "t")),
            labeljust=str(title_cfg.get("labeljust", "c")),
            fontname=str(title_cfg.get("fontname", "Arial")),
            fontsize=str(title_cfg.get("fontsize", 24)),
            fontcolor=str(title_cfg.get("fontcolor", "#222222")),
        )

    strengths = [_node_strength(graph, node, weight) for node in graph.nodes]
    for node in graph.nodes:
        agraph.add_node(node, **_node_attrs(graph, node, strengths, config, weight))

    edges = list(graph.edges(data=True))
    values = [float(data.get(weight, 0.0)) for _, _, data in edges]
    label_cfg = dict(config.get("edge_labels", {}))
    label_count = int(label_cfg.get("top_edges_for_labelling", 0))
    labelled_edges = {
        (source, target)
        for source, target, _ in sorted(edges, key=lambda item: float(item[2].get(weight, 0.0)), reverse=True)[:label_count]
    }
    for source, target, data in sorted(edges, key=lambda item: float(item[2].get(weight, 0.0))):
        agraph.add_edge(source, target, **_edge_attrs(source, target, data, values, config, weight, labelled_edges))
    agraph.layout(prog=layout)
    return agraph


def draw_transfer_graph(
    graph: nx.DiGraph,
    output: str | Path,
    *,
    layout: str = "dot",
    config: dict | None = None,
    use_colors: bool = True,
    title: str | None = None,
) -> None:
    if graph.number_of_nodes() == 0:
        raise ValueError("Cannot draw an empty transfer graph")
    config = dict(config or {})
    if not use_colors:
        edge_cfg = dict(config.get("edges", {}))
        edge_cfg["color_map"] = "Greys"
        edge_cfg["color_low"] = 0.35
        edge_cfg["color_high"] = 0.75
        config["edges"] = edge_cfg
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    agraph = _build_agraph(graph, layout=layout, config=config, title=title)
    agraph.draw(str(path), format=path.suffix.lstrip(".") or "png")
