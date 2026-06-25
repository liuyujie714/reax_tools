"""Sankey-style flow diagrams for transfer_flow.csv."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle
import pandas as pd

from .chemistry import formula_to_subscript
from .filters import prepare_transfer_flow


def _labels(flow: pd.DataFrame) -> dict[str, str]:
    labels = {}
    for _, row in flow.iterrows():
        labels[str(row["source_id"])] = formula_to_subscript(row["source_label"])
        labels[str(row["target_id"])] = formula_to_subscript(row["target_label"])
    return labels


def _plain(label: str) -> str:
    return label.translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789"))


def _resolve_node(flow: pd.DataFrame, value: str) -> str | None:
    value = str(value)
    labels = _labels(flow)
    if value in labels:
        return value
    for node, label in labels.items():
        if _plain(label) == value:
            return node
    return None


def _draw_band(ax, start, end, width: float, color, alpha: float = 0.55, zorder: int = 1) -> None:
    x0, y0 = start
    x1, y1 = end
    dx = max(0.4, abs(x1 - x0) * 0.5)
    verts = [(x0, y0), (x0 + dx, y0), (x1 - dx, y1), (x1, y1)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    patch = PathPatch(
        MplPath(verts, codes),
        facecolor="none",
        edgecolor=color,
        lw=width,
        alpha=alpha,
        capstyle="round",
        joinstyle="round",
        zorder=zorder,
    )
    ax.add_patch(patch)


def _scaled_width(values: list[float], min_width: float, max_width: float):
    if not values:
        return lambda _: (min_width + max_width) / 2
    low = min(values)
    high = max(values)
    if low == high:
        return lambda _: (min_width + max_width) / 2
    return lambda value: min_width + (value - low) * (max_width - min_width) / (high - low)


def _draw_box(ax, x: float, y: float, w: float, h: float, text: str, *, face: str, fontsize: float) -> None:
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, facecolor=face, edgecolor="#444", lw=0.8, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, zorder=4)


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", _plain(text))


def _formula_atom_count(label: str) -> int:
    text = _plain(label)
    total = 0
    for _, count in re.findall(r"([A-Z][a-z]*)(\d*)", text):
        total += int(count or 1)
    return max(1, total)


def _weight_column(config: dict) -> str:
    return str(config.get("weight", "atom_transfer"))


def _add_colorbar(fig, ax, cmap, norm, label: str, fontsize: float) -> None:
    scalar = cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    cbar = fig.colorbar(scalar, ax=ax, fraction=0.032, pad=0.02)
    cbar.set_label(label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=max(8, fontsize - 2))


def _node_flux_color(cmap, norm, value: float):
    return cmap(0.10 + 0.58 * norm(value))


def _flow_cmap(config: dict):
    name = str(config.get("color_map", "fresh"))
    if name == "fresh":
        return LinearSegmentedColormap.from_list(
            "reax_fresh",
            ["#f1f5a6", "#b8e3b3", "#67c8c4", "#5a8fd6", "#7f68c7"],
        )
    return plt.get_cmap(name)


def _allocate_heights(values: list[float], total_height: float, min_height: float) -> list[float]:
    if not values:
        return []
    total = sum(values)
    if total <= 0:
        return [total_height / len(values)] * len(values)
    heights = [total_height * value / total for value in values]
    if min_height <= 0 or len(values) * min_height > total_height * 0.82:
        return heights

    fixed = set()
    while True:
        changed = False
        for idx, height in enumerate(heights):
            if idx not in fixed and height < min_height:
                fixed.add(idx)
                changed = True
        if not changed:
            break
        remaining_height = total_height - len(fixed) * min_height
        remaining_value = sum(value for idx, value in enumerate(values) if idx not in fixed)
        if remaining_height <= 0 or remaining_value <= 0:
            break
        for idx, value in enumerate(values):
            heights[idx] = min_height if idx in fixed else remaining_height * value / remaining_value
    return heights


def _edge_color(cmap, norm, value: float):
    return cmap(0.14 + 0.50 * norm(value))


def _draw_layered_sankey(
    nodes: list[dict],
    edges: list[tuple[str, str, float]],
    output: str | Path,
    *,
    config: dict,
    metadata: dict | None = None,
) -> bool:
    if not nodes:
        return False
    metadata = metadata or {}
    node_by_id = {str(node["id"]): node for node in nodes}
    used_ids = {str(node_id) for edge in edges for node_id in edge[:2]}
    used_ids.update(str(node["id"]) for node in nodes if float(node.get("flux", 0)) > 0)
    nodes = [node for node in nodes if str(node["id"]) in used_ids]
    if not nodes:
        return False

    by_layer = defaultdict(list)
    for node in nodes:
        by_layer[int(node["level"])].append(node)
    layers = sorted(by_layer)
    for layer in layers:
        by_layer[layer].sort(key=lambda item: float(item.get("flux", 0)), reverse=True)

    layer_totals = {
        layer: sum(float(node.get("flux", 0)) for node in group)
        for layer, group in by_layer.items()
    }
    max_total = max(layer_totals.values()) if layer_totals else 1.0
    if max_total <= 0:
        max_total = 1.0

    fig, ax = plt.subplots(figsize=(float(config.get("figure_width", 8.0)), float(config.get("figure_height", 4.8))))
    ax.axis("off")

    cmap = _flow_cmap(config)
    values = [float(weight) for _, _, weight in edges] + [float(node.get("flux", 0)) for node in nodes]
    values = [value for value in values if value > 0] or [1.0]
    norm = mcolors.Normalize(vmin=min(values), vmax=max(values))
    width_for = _scaled_width(
        [float(weight) for _, _, weight in edges] or values,
        float(config.get("edge_width_min", 2.0)),
        float(config.get("edge_width_max", 16.0)),
    )

    bar_width = float(config.get("bar_width", config.get("node_box_width", 0.84)))
    max_bar_height = float(config.get("max_bar_height", 3.8))
    min_bar_height = float(config.get("min_bar_height", 0.65))
    segment_gap = float(config.get("segment_gap", 0.035))
    min_segment_height = float(config.get("min_segment_height", 0.18))
    layer_spacing = float(config.get("layer_spacing", 2.25))
    font = float(config.get("font_size", 10))
    label_min_height = float(config.get("label_min_height", 0.23))
    inside_label_max_chars = int(config.get("inside_label_max_chars", 5))
    outside_label_gap = float(config.get("outside_label_gap", 0.22))
    alpha = float(config.get("alpha", 0.58))

    positions = {}
    for idx, layer in enumerate(layers):
        group = by_layer[layer]
        raw_height = max_bar_height * layer_totals[layer] / max_total
        total_height = max(min_bar_height, raw_height) if group else 0.0
        total_gap = segment_gap * max(0, len(group) - 1)
        usable_height = max(0.1, total_height - total_gap)
        heights = _allocate_heights(
            [float(node.get("flux", 0)) for node in group],
            usable_height,
            min_segment_height,
        )
        top = total_height / 2
        x = idx * layer_spacing
        for node, height in zip(group, heights):
            y_top = top
            y_bottom = top - height
            y_center = (y_top + y_bottom) / 2
            node_id = str(node["id"])
            positions[node_id] = {
                "x": x,
                "y": y_center,
                "height": height,
                "top": y_top,
                "bottom": y_bottom,
            }
            top = y_bottom - segment_gap

    incoming_edges = defaultdict(list)
    outgoing_edges = defaultdict(list)
    for source, target, weight in edges:
        source = str(source)
        target = str(target)
        outgoing_edges[source].append((target, float(weight)))
        incoming_edges[target].append((source, float(weight)))

    def _ports(grouped_edges):
        ports = {}
        for node_id, items in grouped_edges.items():
            if node_id not in positions:
                continue
            items = sorted(items, key=lambda item: item[1], reverse=True)
            total = sum(weight for _, weight in items)
            pos = positions[node_id]
            if total <= 0:
                for other, _ in items:
                    ports[(node_id, other)] = pos["y"]
                continue
            cursor = pos["top"]
            for other, weight in items:
                height = pos["height"] * weight / total
                ports[(node_id, other)] = cursor - height / 2
                cursor -= height
        return ports

    source_ports = _ports(outgoing_edges)
    target_ports = _ports(incoming_edges)

    for source, target, weight in edges:
        source = str(source)
        target = str(target)
        if source not in positions or target not in positions:
            continue
        s = positions[source]
        t = positions[target]
        _draw_band(
            ax,
            (s["x"] + bar_width / 2, source_ports.get((source, target), s["y"])),
            (t["x"] - bar_width / 2, target_ports.get((target, source), t["y"])),
            width_for(float(weight)),
            _edge_color(cmap, norm, float(weight)),
            alpha=alpha,
            zorder=1,
        )

    outside_labels = []
    for node in nodes:
        node_id = str(node["id"])
        if node_id not in positions:
            continue
        pos = positions[node_id]
        label = str(node["label"])
        face = _node_flux_color(cmap, norm, float(node.get("flux", 0)))
        ax.add_patch(
            Rectangle(
                (pos["x"] - bar_width / 2, pos["bottom"]),
                bar_width,
                pos["height"],
                facecolor=face,
                edgecolor="#444",
                lw=0.65,
                zorder=3,
            )
        )
        if pos["height"] >= label_min_height and len(_plain(label)) <= inside_label_max_chars:
            ax.text(
                pos["x"],
                pos["y"],
                label,
                ha="center",
                va="center",
                fontsize=font,
                zorder=4,
            )
        else:
            outside_labels.append((pos["x"], pos["y"], pos["top"], pos["bottom"], label))

    by_label_column = defaultdict(list)
    for item in outside_labels:
        by_label_column[item[0]].append(item)
    for x, items in by_label_column.items():
        items = sorted(items, key=lambda item: item[1], reverse=True)
        adjusted = []
        last_y = None
        for item in items:
            y = item[1] if last_y is None else min(item[1], last_y - outside_label_gap)
            adjusted.append((item, y))
            last_y = y
        overflow = min((y for _, y in adjusted), default=0) - min((item[3] for item, _ in adjusted), default=0)
        if overflow < 0:
            adjusted = [(item, y - overflow) for item, y in adjusted]
        for (x0, y0, _top, _bottom, label), y in adjusted:
            label_x = x0 + bar_width / 2 + 0.10
            ax.plot([x0 + bar_width / 2, label_x - 0.03], [y0, y], color="#666", lw=0.45, zorder=4)
            ax.text(label_x, y, label, ha="left", va="center", fontsize=font, zorder=5)

    layer_labels = metadata.get("layer_labels") or {}
    for idx, layer in enumerate(layers):
        x = idx * layer_spacing
        title = layer_labels.get(layer, f"Stage {idx + 1}")
        ax.text(x, max_bar_height / 2 + 0.38, title, ha="center", va="bottom", fontsize=font, fontweight="bold")

    xs = [pos["x"] for pos in positions.values()]
    ys = [pos["top"] for pos in positions.values()] + [pos["bottom"] for pos in positions.values()]
    if bool(config.get("show_colorbar", True)):
        _add_colorbar(fig, ax, cmap, norm, "Atom transfer", font)
    note = metadata.get("note")
    if note:
        ax.text(min(xs) - 0.2, min(ys) - 0.42, note, fontsize=max(7, font - 2), ha="left", va="top")
    ax.set_xlim(min(xs) - layer_spacing * 0.55, max(xs) + layer_spacing * 0.85)
    ax.set_ylim(min(ys) - 0.7, max(max(ys), max_bar_height / 2) + 0.7)
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return True


def draw_focus_sankey(
    raw_flow: pd.DataFrame,
    center: str,
    output: str | Path,
    *,
    config: dict | None = None,
    include_self_loops: bool = False,
) -> bool:
    config = config or {}
    flow = prepare_transfer_flow(raw_flow, include_self_loops=include_self_loops)
    weight_col = _weight_column(config)
    center_id = _resolve_node(flow, center)
    if center_id is None:
        return False
    labels = _labels(flow)
    incoming = flow[flow["target_id"].astype(str) == center_id].copy()
    outgoing = flow[flow["source_id"].astype(str) == center_id].copy()
    incoming = incoming.sort_values(weight_col, ascending=False).head(int(config.get("max_in", 8)))
    outgoing = outgoing.sort_values(weight_col, ascending=False).head(int(config.get("max_out", 8)))
    if incoming.empty and outgoing.empty:
        return False

    nodes = []
    edges = []
    center_draw_id = f"focus:{center_id}"
    for _, row in incoming.iterrows():
        source = str(row["source_id"])
        source_draw_id = f"in:{source}"
        value = float(row[weight_col])
        nodes.append({"id": source_draw_id, "label": labels.get(source, source), "level": 0, "flux": value})
        edges.append((source_draw_id, center_draw_id, value))
    center_flux = max(incoming[weight_col].astype(float).sum(), outgoing[weight_col].astype(float).sum())
    nodes.append({"id": center_draw_id, "label": labels.get(center_id, center_id), "level": 1, "flux": center_flux})
    for _, row in outgoing.iterrows():
        target = str(row["target_id"])
        target_draw_id = f"out:{target}"
        value = float(row[weight_col])
        nodes.append({"id": target_draw_id, "label": labels.get(target, target), "level": 2, "flux": value})
        edges.append((center_draw_id, target_draw_id, value))

    metadata = {"layer_labels": {0: "Inflow", 1: "Focus", 2: "Outflow"}}
    return _draw_layered_sankey(nodes, edges, output, config=config, metadata=metadata)


def draw_focus_sankey_set(
    raw_flow: pd.DataFrame,
    output_dir: str | Path,
    *,
    centers: list[str] | None = None,
    config: dict | None = None,
    include_self_loops: bool = False,
) -> list[Path]:
    config = config or {}
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    prepared = prepare_transfer_flow(raw_flow, include_self_loops=include_self_loops)
    centers = centers or role_focus_centers(prepared, int(config.get("centers", 5)), config=config)
    labels = _labels(prepared)
    written = []
    for center in centers:
        center_id = _resolve_node(prepared, center) or center
        label = labels.get(center_id, str(center))
        path = out / f"focus_{_safe_name(label)}.png"
        if draw_focus_sankey(prepared, center, path, config=config, include_self_loops=True):
            written.append(path)
    return written


def role_focus_centers(flow: pd.DataFrame, n: int, *, config: dict | None = None) -> list[str]:
    config = config or {}
    weight_col = _weight_column(config)
    scores = Counter()
    for _, row in flow.iterrows():
        weight = float(row[weight_col])
        scores[str(row["source_id"])] += weight
        scores[str(row["target_id"])] += weight
    return [node for node, _ in scores.most_common(n)]


def _role_breaks(config: dict) -> list[float]:
    return [float(value) for value in config.get("role_breaks", [-0.65, -0.25, 0.25, 0.65])]


def _role_level(role: float, breaks: list[float]) -> int:
    for idx, value in enumerate(breaks):
        if role <= value:
            return idx
    return len(breaks)


def role_flow_model(raw_flow: pd.DataFrame, *, config: dict | None = None, include_self_loops: bool = False):
    config = config or {}
    flow = prepare_transfer_flow(raw_flow, include_self_loops=include_self_loops)
    if flow.empty:
        return [], [], {}
    weight_col = _weight_column(config)
    labels = _labels(flow)
    incoming = Counter()
    outgoing = Counter()
    all_nodes = set()
    for _, row in flow.iterrows():
        source = str(row["source_id"])
        target = str(row["target_id"])
        weight = float(row[weight_col])
        outgoing[source] += weight
        incoming[target] += weight
        all_nodes.update([source, target])

    breaks = _role_breaks(config)
    stats = {}
    by_level = defaultdict(list)
    for node in all_nodes:
        win = incoming[node]
        wout = outgoing[node]
        total = win + wout
        if total <= 0:
            continue
        role = (win - wout) / total
        level = _role_level(role, breaks)
        flux = max(win, wout)
        item = {"id": node, "label": labels.get(node, node), "level": level, "role": role, "flux": flux}
        stats[node] = item
        by_level[level].append(item)

    max_species = int(config.get("max_species_per_level", 6))
    terminal_keep = int(config.get("keep_terminal_species", 4))
    forced = set()
    if terminal_keep > 0:
        size_penalty = float(config.get("terminal_size_penalty", 0.65))

        def product_score(item):
            node = item["id"]
            net = incoming[node] - outgoing[node]
            if net <= 0:
                return -1.0
            size = _formula_atom_count(item["label"])
            return net * max(float(item["role"]), 0.05) / (size ** size_penalty)

        def reactant_score(item):
            node = item["id"]
            net = outgoing[node] - incoming[node]
            if net <= 0:
                return -1.0
            size = _formula_atom_count(item["label"])
            return net * max(-float(item["role"]), 0.05) / (size ** size_penalty)

        products = sorted(
            stats.values(),
            key=lambda item: (product_score(item), incoming[item["id"]] - outgoing[item["id"]]),
            reverse=True,
        )
        reactants = sorted(
            stats.values(),
            key=lambda item: (reactant_score(item), outgoing[item["id"]] - incoming[item["id"]]),
            reverse=True,
        )
        forced.update(item["id"] for item in products[:terminal_keep] if incoming[item["id"]] > outgoing[item["id"]])
        forced.update(item["id"] for item in reactants[:terminal_keep] if outgoing[item["id"]] > incoming[item["id"]])

    display = {}
    nodes = []
    other_flux = Counter()
    other_ids = {}
    for level, items in by_level.items():
        items = sorted(items, key=lambda item: item["flux"], reverse=True)
        selected_ids = {item["id"] for item in items[:max_species]}
        selected_ids.update(item["id"] for item in items if item["id"] in forced)
        selected = [item for item in items if item["id"] in selected_ids]
        for item in selected:
            display[item["id"]] = item["id"]
            nodes.append(item)
        rest = [item for item in items if item["id"] not in selected_ids]
        if rest:
            other_id = f"__other_level_{level}"
            other_ids[level] = other_id
            flux = sum(item["flux"] for item in rest)
            other_flux[level] = flux
            nodes.append({"id": other_id, "label": f"Other L{level + 1}", "level": level, "role": None, "flux": flux})
            for item in rest:
                display[item["id"]] = other_id

    edge_weights = Counter()
    discarded = Counter()
    for _, row in flow.iterrows():
        source = str(row["source_id"])
        target = str(row["target_id"])
        if source not in stats or target not in stats:
            continue
        source_level = stats[source]["level"]
        target_level = stats[target]["level"]
        weight = float(row[weight_col])
        if target_level <= source_level:
            discarded["same_or_backward"] += weight
            continue
        ds = display.get(source)
        dt = display.get(target)
        if ds is None or dt is None or ds == dt:
            discarded["collapsed"] += weight
            continue
        edge_weights[(ds, dt)] += weight

    max_links = int(config.get("max_links", 24))
    edge_items = edge_weights.most_common()
    selected_edges = {(s, t) for (s, t), _ in edge_items[:max_links]}
    anchor_links = int(config.get("keep_anchor_links", 4))
    if anchor_links > 0 and forced:
        for node_id in forced:
            connected = [((s, t), w) for (s, t), w in edge_items if s == node_id or t == node_id]
            for (s, t), _ in connected[:anchor_links]:
                selected_edges.add((s, t))
    edges = [(s, t, edge_weights[(s, t)]) for s, t in edge_weights if (s, t) in selected_edges]
    edges.sort(key=lambda item: item[2], reverse=True)
    used = {node for edge in edges for node in edge[:2]}
    nodes = [node for node in nodes if node["id"] in used]
    return nodes, edges, {"discarded": dict(discarded), "weight": weight_col, "forced": list(forced)}


def draw_role_flow(
    raw_flow: pd.DataFrame,
    output: str | Path,
    *,
    config: dict | None = None,
    include_self_loops: bool = False,
) -> bool:
    config = config or {}
    nodes, edges, metadata = role_flow_model(raw_flow, config=config, include_self_loops=include_self_loops)
    if not nodes:
        return False
    discarded = metadata.get("discarded", {})
    metadata["note"] = None
    if bool(config.get("show_discarded_note", True)):
        metadata["note"] = f"omitted back/same-layer atom transfer: {int(discarded.get('same_or_backward', 0))}"
    return _draw_layered_sankey(nodes, edges, output, config=config, metadata=metadata)
