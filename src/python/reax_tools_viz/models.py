"""Data models for ReaxTools Python processing.

These classes describe loaded C++ outputs. They intentionally contain no file
parsing, plotting, or filtering logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MoleculeRecord:
    id: str
    formula: str
    atom_counts: dict[str, int]


@dataclass(frozen=True)
class TransferEdge:
    source_id: str
    target_id: str
    source_label: str
    target_label: str
    count: int
    atom_transfer: int
    self_loop: bool


@dataclass(frozen=True)
class ReactionEvent:
    event_id: int
    frame: int
    reactant_ids: tuple[str, ...]
    product_ids: tuple[str, ...]
    tracked_reactant_ids: tuple[str, ...]
    tracked_product_ids: tuple[str, ...]
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    atom_transfer: int
    atom_conserved: bool


@dataclass(frozen=True)
class RunBundle:
    root: Path
    manifest: dict
    molecules: dict[str, MoleculeRecord]

