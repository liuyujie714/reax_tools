"""Plot template loading and merging."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


DEFAULT_TEMPLATE = Path(__file__).with_name("default_plot_template.yaml")


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_plot_template(output_dir: str | Path | None = None, template_file: str | Path | None = None) -> dict[str, Any]:
    config = yaml.safe_load(DEFAULT_TEMPLATE.read_text()) or {}
    if output_dir is not None:
        for name in ["reax_tools_plot.yaml", "reax_tools_template.yaml"]:
            local = Path(output_dir) / name
            if local.exists():
                config = _deep_merge(config, yaml.safe_load(local.read_text()) or {})
                break
    if template_file is not None:
        config = _deep_merge(config, yaml.safe_load(Path(template_file).read_text()) or {})
    return config
