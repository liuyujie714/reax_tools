#!/usr/bin/env python3
"""Basic regression checks for ReaxTools raw C++ outputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


REQUIRED_FILES = [
    "species_count.csv",
    "bond_count.csv",
    "atom_bonded_num_count.csv",
    "ring_count.csv",
    "reaction_events.csv",
    "reaction_event_pairs.csv",
    "transfer_flow.csv",
    "molecules.json",
    "reax_tools.log",
    "reax_tools_manifest.json",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fp:
        return list(csv.DictReader(fp))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--expect-events", type=int)
    parser.add_argument("--expect-transfer-edges", type=int)
    parser.add_argument("--expect-self-loop-count", type=int)
    args = parser.parse_args()

    out = args.output_dir
    missing = [name for name in REQUIRED_FILES if not (out / name).exists()]
    if missing:
        raise SystemExit(f"missing required files: {', '.join(missing)}")

    manifest = json.loads((out / "reax_tools_manifest.json").read_text())
    molecules = json.loads((out / "molecules.json").read_text())
    events = read_csv(out / "reaction_events.csv")
    transfer = read_csv(out / "transfer_flow.csv")

    manifest_files = set(manifest.get("files", []))
    for name in REQUIRED_FILES:
        if name not in manifest_files and name != "reax_tools_manifest.json":
            raise SystemExit(f"manifest does not list {name}")

    for count_file in ["bond_count.csv", "atom_bonded_num_count.csv", "ring_count.csv"]:
        header = (out / count_file).read_text().splitlines()[0].split(",")
        if not header or header[0] != "frame":
            raise SystemExit(f"{count_file} must start with frame column")

    species_header = (out / "species_count.csv").read_text().splitlines()[0].split(",")
    if species_header != ["frame", "molecule_id", "formula", "count"]:
        raise SystemExit("species_count.csv must use frame,molecule_id,formula,count schema")

    molecule_ids = {str(molecule["id"]) for molecule in molecules.get("molecules", [])}
    species_rows = read_csv(out / "species_count.csv")
    unknown_species_ids = {row["molecule_id"] for row in species_rows if row["molecule_id"] not in molecule_ids}
    if unknown_species_ids:
        raise SystemExit(f"species_count.csv references unknown molecule ids: {sorted(unknown_species_ids)[:10]}")

    unknown_transfer_ids = {
        row[key]
        for row in transfer
        for key in ["source_id", "target_id"]
        if row[key] not in molecule_ids
    }
    if unknown_transfer_ids:
        raise SystemExit(f"transfer_flow.csv references unknown molecule ids: {sorted(unknown_transfer_ids)[:10]}")

    for required in ["reactant_ids", "product_ids", "tracked_reactant_ids", "tracked_product_ids"]:
        if events and required not in events[0]:
            raise SystemExit(f"reaction_events.csv missing {required}")

    unknown_event_ids = set()
    for row in events:
        for key in ["reactant_ids", "product_ids"]:
            for molecule_id in row[key].split("+"):
                if molecule_id and molecule_id not in molecule_ids:
                    unknown_event_ids.add(molecule_id)
    if unknown_event_ids:
        raise SystemExit(f"reaction_events.csv references unknown molecule ids: {sorted(unknown_event_ids)[:10]}")

    bad_conservation = [row["event_id"] for row in events if row.get("atom_conserved") != "1"]
    if bad_conservation:
        raise SystemExit(f"atom conservation failed for events: {bad_conservation[:10]}")

    if args.expect_events is not None and len(events) != args.expect_events:
        raise SystemExit(f"event count mismatch: got {len(events)}, expected {args.expect_events}")

    if args.expect_transfer_edges is not None and len(transfer) != args.expect_transfer_edges:
        raise SystemExit(f"transfer edge count mismatch: got {len(transfer)}, expected {args.expect_transfer_edges}")

    self_loop_sum = sum(int(row["count"]) for row in transfer if row.get("self_loop") == "1")
    if args.expect_self_loop_count is not None and self_loop_sum != args.expect_self_loop_count:
        raise SystemExit(f"self-loop count mismatch: got {self_loop_sum}, expected {args.expect_self_loop_count}")

    print("ReaxTools output check passed")
    print(f"  events: {len(events)}")
    print(f"  transfer edges: {len(transfer)}")
    print(f"  self-loop count: {self_loop_sum}")
    print(f"  molecules: {len(molecules.get('molecules', []))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
