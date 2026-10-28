#!/usr/bin/env python3

import json
import os
import sys
from typing import Dict, List, Tuple

MOL2_HEADER = """@<TRIPOS>MOLECULE
{name}
{num_atoms} {num_bonds}
SMALL
NO_CHARGES


@<TRIPOS>ATOM
"""


def generate_mol2_atoms(atom_ids: List[int], coords: List[Dict]) -> Tuple[str, Dict[int, int]]:
    """Return the ATOM section plus a map of old atom id -> new sequential id."""
    coord_map = {atom['id']: atom for atom in coords}
    id_map = {old: new for new, old in enumerate(atom_ids, 1)}

    atom_lines = []
    for old_id, new_id in id_map.items():
        atom = coord_map.get(old_id)
        if atom is None:
            continue
        element = atom['element']
        atom_lines.append(
            f"{new_id:>6} {element}{new_id:<6} "
            f"{atom['x']:>10.6f} {atom['y']:>10.6f} {atom['z']:>10.6f} "
            f"{element:<4} 1 UNL 0.0000")
    return "\n".join(atom_lines), id_map


def generate_mol2_bonds(bonds: List[Dict], id_map: Dict[int, int]) -> Tuple[str, int]:
    """Return the BOND section (bonds between molecule atoms only) and its count."""
    kept = [b for b in bonds if b['a'] in id_map and b['b'] in id_map]
    lines = [f"{i:>6} {id_map[b['a']]:>6} {id_map[b['b']]:>6} {b['order']}"
             for i, b in enumerate(kept, 1)]
    return "\n".join(lines), len(kept)


def generate_mol2(atom_ids: List[int], coords: List[Dict], bonds: List[Dict],
                  name: str) -> Tuple[str, int]:
    atoms_str, id_map = generate_mol2_atoms(atom_ids, coords)
    bonds_str, bond_count = generate_mol2_bonds(bonds, id_map)

    mol2 = MOL2_HEADER.format(name=name, num_atoms=len(atom_ids), num_bonds=bond_count) + atoms_str
    if bonds_str:
        mol2 += f"\n\n@<TRIPOS>BOND\n{bonds_str}"
    return mol2 + "\n", bond_count


def process_type(data: Dict, mol_type: str, output_dir: str) -> int:
    """Write .mol2 files for every reactant/product molecule and return the count."""
    prefix = 'reactant' if mol_type == 'reactant' else 'product'
    label = '反应物' if prefix == 'reactant' else '产物'
    molecules = data.get(prefix + 's', [])
    coords = data.get(prefix + '_coords', [])
    bonds = data.get(prefix + '_bonds', [])

    if not molecules:
        print(f"  [!] 未找到{label}数据")
        return 0

    for mol in molecules:
        atom_ids = mol.get('atom_ids', [])
        name = f"{prefix}_{mol.get('formula', 'UNKNOWN')}_{mol.get('tracked_molecule_id', 0)}"
        print(f"  生成分子: {name}")
        print(f"    原子ID: {atom_ids}")

        mol2, bond_count = generate_mol2(atom_ids, coords, bonds, name)
        output_file = os.path.join(output_dir, f"{name}.mol2")
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(mol2)
        print(f"  [OK] 已生成: {output_file}")
        print(f"    原子数: {len(atom_ids)}, 键数: {bond_count}")
    return len(molecules)


def json_to_mol2(json_file_path: str, output_dir: str = "."):
    os.makedirs(output_dir, exist_ok=True)

    try:
        with open(json_file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"错误: 文件 '{json_file_path}' 未找到")
        return
    except json.JSONDecodeError as e:
        print(f"错误: JSON解析失败 - {e}")
        return

    missing = [f for f in ['reactants', 'products', 'reactant_coords', 'product_coords']
               if f not in data]
    if missing:
        print(f"错误: JSON缺少必需字段 '{missing[0]}'")
        return

    print(f"开始处理JSON文件: {json_file_path}")
    print("=" * 50)

    print("处理反应物...")
    n_reactants = process_type(data, 'reactant', output_dir)

    print("\n处理产物...")
    n_products = process_type(data, 'product', output_dir)

    print("\n" + "=" * 50)
    print(f"完成! 共生成 {n_reactants} 个反应物和 {n_products} 个产物的MOL2文件")
    print(f"文件保存在: {os.path.abspath(output_dir)}")


def main():
    if len(sys.argv) < 2:
        print("用法: python reaction_to_mol2.py <json文件路径> [输出目录]")
        print("示例: python reaction_to_mol2.py reaction.json")
        print("示例: python reaction_to_mol2.py reaction.json ./mol2_output")
        return

    json_to_mol2(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else ".")


if __name__ == "__main__":
    main()
