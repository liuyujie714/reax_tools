"""Export helpers for filtered ReaxTools Python products."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def export_transfer_flow_dot(flow: pd.DataFrame, output: str | Path) -> None:
    with Path(output).open("w") as fp:
        fp.write("digraph TransferFlow {\n")
        fp.write("  rankdir=LR;\n")
        for _, row in flow.iterrows():
            fp.write(f'  "{row["source_id"]}" [label="{row["source_label"]}"];\n')
            fp.write(f'  "{row["target_id"]}" [label="{row["target_label"]}"];\n')
            fp.write(
                f'  "{row["source_id"]}" -> "{row["target_id"]}" '
                f'[label="{int(row["count"])}"];\n'
            )
        fp.write("}\n")

