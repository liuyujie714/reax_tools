"""Manifest-driven readers for ReaxTools raw output directories."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .models import MoleculeRecord, RunBundle


REQUIRED_FILES = {
    "reax_tools_manifest.json",
    "molecules.json",
    "species_count.csv",
    "bond_count.csv",
    "atom_bonded_num_count.csv",
    "ring_count.csv",
    "reaction_events.csv",
    "reaction_event_pairs.csv",
    "transfer_flow.csv",
}


def load_bundle(output_dir: str | Path, validate: bool = True) -> RunBundle:
    root = Path(output_dir)
    manifest_path = root / "reax_tools_manifest.json"
    manifest = json.loads(manifest_path.read_text())

    molecules_json = json.loads((root / "molecules.json").read_text())
    molecules = {
        str(item["id"]): MoleculeRecord(
            id=str(item["id"]),
            formula=str(item["formula"]),
            atom_counts={str(k): int(v) for k, v in item.get("atom_counts", {}).items()},
        )
        for item in molecules_json.get("molecules", [])
    }

    bundle = RunBundle(root=root, manifest=manifest, molecules=molecules)
    if validate:
        validate_bundle(bundle)
    return bundle


def validate_bundle(bundle: RunBundle) -> None:
    missing = [name for name in REQUIRED_FILES if not (bundle.root / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing ReaxTools output files: {', '.join(missing)}")

    molecule_ids = set(bundle.molecules)

    species = pd.read_csv(bundle.root / "species_count.csv", dtype={"molecule_id": str})
    unknown = set(species["molecule_id"]) - molecule_ids
    if unknown:
        raise ValueError(f"species_count.csv references unknown molecule ids: {sorted(unknown)[:10]}")

    flow = pd.read_csv(bundle.root / "transfer_flow.csv", dtype={"source_id": str, "target_id": str})
    unknown = (set(flow["source_id"]) | set(flow["target_id"])) - molecule_ids
    if unknown:
        raise ValueError(f"transfer_flow.csv references unknown molecule ids: {sorted(unknown)[:10]}")

    events = pd.read_csv(bundle.root / "reaction_events.csv", dtype=str)
    for column in [
        "reactant_hashes",
        "product_hashes",
        "reactant_formulas",
        "product_formulas",
    ]:
        if column not in events.columns:
            raise ValueError(f"reaction_events.csv missing {column}")

    unknown_event_ids: set[str] = set()
    for column in ["reactant_hashes", "product_hashes"]:
        for value in events[column].dropna():
            unknown_event_ids.update(part for part in value.split("+") if part and part not in molecule_ids)
    if unknown_event_ids:
        raise ValueError(f"reaction_events.csv references unknown molecule ids: {sorted(unknown_event_ids)[:10]}")


def read_count_table(bundle: RunBundle, name: str) -> pd.DataFrame:
    return pd.read_csv(bundle.root / name)


def read_transfer_flow(bundle: RunBundle) -> pd.DataFrame:
    return pd.read_csv(bundle.root / "transfer_flow.csv", dtype={"source_id": str, "target_id": str})


def read_reaction_events(bundle: RunBundle) -> pd.DataFrame:
    return pd.read_csv(bundle.root / "reaction_events.csv", dtype=str)
