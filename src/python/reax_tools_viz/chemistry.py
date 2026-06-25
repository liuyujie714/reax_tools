"""Small chemistry-facing formatting helpers."""

from __future__ import annotations

import re


SUBSCRIPT_TRANS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def formula_to_subscript(formula: object) -> str:
    """Render only digits in formula-like labels as unicode subscripts."""
    return str(formula).translate(SUBSCRIPT_TRANS)


def reaction_to_subscript(label: object) -> str:
    return " + ".join(formula_to_subscript(part) for part in str(label).split("+"))


def natural_count_key(label: object) -> tuple:
    text = str(label)
    parts = re.split(r"(\d+)", text)
    return tuple(int(part) if part.isdigit() else part for part in parts)


def parse_atom_degree(label: str) -> tuple[str, int] | None:
    match = re.fullmatch(r"([A-Z][a-z]?)-(\d+)", str(label))
    if not match:
        return None
    return match.group(1), int(match.group(2))


def parse_bond(label: str) -> tuple[str, str] | None:
    match = re.fullmatch(r"([A-Z][a-z]?)-([A-Z][a-z]?)", str(label))
    if not match:
        return None
    return match.group(1), match.group(2)
