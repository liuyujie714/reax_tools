"""Optional molecule drawing from molecules.json graph records."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx


ELEMENT_COLORS = {
    "H": "#f7f7f7",
    "C": "#555555",
    "N": "#4c78a8",
    "O": "#e45756",
    "S": "#f2cf5b",
    "P": "#b279a2",
}


def draw_molecule_record(record: dict, output: str | Path) -> None:
    example = record.get("example", {})
    atoms = example.get("atoms", [])
    bonds = example.get("bonds", [])
    if not atoms:
        raise ValueError(f"molecule {record.get('id')} has no example atoms")

    graph = nx.Graph()
    for atom in atoms:
        atom_id = str(atom["id"])
        element = str(atom["element"])
        graph.add_node(atom_id, element=element, label=f"{element}{atom_id}")
    for bond in bonds:
        graph.add_edge(str(bond["a"]), str(bond["b"]), order=int(bond.get("order", 1)))

    pos = nx.spring_layout(graph, seed=7)
    colors = [ELEMENT_COLORS.get(graph.nodes[n]["element"], "#cccccc") for n in graph.nodes]
    labels = {n: graph.nodes[n]["element"] for n in graph.nodes}

    plt.figure(figsize=(4, 4))
    nx.draw_networkx_edges(graph, pos, width=1.5, edge_color="#555")
    nx.draw_networkx_nodes(graph, pos, node_color=colors, edgecolors="#222", node_size=500)
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=10)
    plt.title(str(record.get("formula", record.get("id", ""))))
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def draw_selected_molecules(molecules: list[dict], output_dir: str | Path, species: list[str] | None = None) -> list[Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    wanted = set(species or [])
    written: list[Path] = []
    for record in molecules:
        if wanted and str(record.get("id")) not in wanted and str(record.get("formula")) not in wanted:
            continue
        filename = f"{record.get('formula', record.get('id'))}_{record.get('id')}.png"
        path = out / filename
        draw_molecule_record(record, path)
        written.append(path)
    return written
