#!/usr/bin/env python3

import json
import sys
import os
from typing import Dict, List, Any, Set

def generate_mol2_header(molecule_name: str, num_atoms: int, num_bonds: int) -> str:
    header = f"""@<TRIPOS>MOLECULE
{molecule_name}
{num_atoms} {num_bonds}
SMALL
NO_CHARGES


@<TRIPOS>ATOM
"""
    return header

def generate_mol2_atoms(atom_ids: List[int], coords: List[Dict]) -> tuple:
    atom_lines = []
    coord_map = {atom['id']: atom for atom in coords}
    
    id_map = {}
    for new_id, old_id in enumerate(atom_ids, 1):
        id_map[old_id] = new_id
        
        if old_id in coord_map:
            atom = coord_map[old_id]
            element = atom['element']
            x, y, z = atom['x'], atom['y'], atom['z']
            
            atom_name = f"{element}{new_id}"
            atom_lines.append(f"{new_id:>6} {atom_name:<6} {x:>10.6f} {y:>10.6f} {z:>10.6f} {element:<4} 1 UNL 0.0000")
    
    return "\n".join(atom_lines), id_map

def generate_mol2_bonds(atom_ids: List[int], bonds: List[Dict], id_map: Dict[int, int]) -> str:
    if not bonds:
        return ""
    
    atom_set = set(atom_ids)
    
    bond_lines = []
    bond_counter = 1
    
    for bond in bonds:
        a1_old = bond['a']
        a2_old = bond['b']
        order = bond['order']
        
        if a1_old in atom_set and a2_old in atom_set:
            a1_new = id_map.get(a1_old)
            a2_new = id_map.get(a2_old)
            
            if a1_new is not None and a2_new is not None:
                bond_lines.append(f"{bond_counter:>6} {a1_new:>6} {a2_new:>6} {order}")
                bond_counter += 1
    
    return "\n".join(bond_lines)

def generate_mol2_for_molecule(atom_ids: List[int], coords: List[Dict], 
                               bonds: List[Dict], name: str) -> str:
    atoms_str, id_map = generate_mol2_atoms(atom_ids, coords)
    
    bonds_str = generate_mol2_bonds(atom_ids, bonds, id_map)
    bond_count = len(bonds_str.split('\n')) if bonds_str else 0
    
    final_mol2 = f"""@<TRIPOS>MOLECULE
{name}
{len(atom_ids)} {bond_count}
SMALL
NO_CHARGES


@<TRIPOS>ATOM
{atoms_str}"""

    if bonds_str:
        final_mol2 += f"""

@<TRIPOS>BOND
{bonds_str}"""
    
    final_mol2 += "\n"
    
    return final_mol2

def extract_molecule_from_json(data: Dict, mol_type: str) -> Dict:
    if mol_type == 'reactant':
        molecules = data.get('reactants', [])
        coords = data.get('reactant_coords', [])
        bonds = data.get('reactant_bonds', [])
    else:
        molecules = data.get('products', [])
        coords = data.get('product_coords', [])
        bonds = data.get('product_bonds', [])
    
    if not molecules:
        return None
    
    result = {
        'molecules': [],
        'bonds': bonds,
        'coords': coords
    }
    
    for mol in molecules:
        atom_ids = mol.get('atom_ids', [])
        formula = mol.get('formula', 'UNKNOWN')
        mol_hash = mol.get('molecule_hash', '')
        mol_id = mol.get('tracked_molecule_id', 0)
        
        result['molecules'].append({
            'atom_ids': atom_ids,
            'formula': formula,
            'hash': mol_hash,
            'id': mol_id
        })
    
    return result

def json_to_mol2(json_file_path: str, output_dir: str = "."):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        with open(json_file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"错误: 文件 '{json_file_path}' 未找到")
        return
    except json.JSONDecodeError as e:
        print(f"错误: JSON解析失败 - {e}")
        return
    
    required_fields = ['reactants', 'products', 'reactant_coords', 'product_coords']
    for field in required_fields:
        if field not in data:
            print(f"错误: JSON缺少必需字段 '{field}'")
            return
    
    print(f"开始处理JSON文件: {json_file_path}")
    print("=" * 50)
    
    print("处理反应物...")
    reactant_data = extract_molecule_from_json(data, 'reactant')
    if reactant_data:
        for mol in reactant_data['molecules']:
            formula = mol['formula']
            mol_id = mol['id']
            atom_ids = mol['atom_ids']
            
            name = f"reactant_{formula}_{mol_id}"
            print(f"  生成分子: {name}")
            print(f"    原子ID: {atom_ids}")
            
            mol2_content = generate_mol2_for_molecule(
                atom_ids, 
                reactant_data['coords'], 
                reactant_data['bonds'],
                name
            )
            
            output_file = os.path.join(output_dir, f"{name}.mol2")
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(mol2_content)
            print(f"  ✓ 已生成: {output_file}")
            
            atom_set = set(atom_ids)
            bond_count = sum(1 for bond in reactant_data['bonds'] 
                           if bond['a'] in atom_set and bond['b'] in atom_set)
            print(f"    原子数: {len(atom_ids)}, 键数: {bond_count}")
    else:
        print("  ⚠ 未找到反应物数据")
    
    print("\n处理产物...")
    product_data = extract_molecule_from_json(data, 'product')
    if product_data:
        for mol in product_data['molecules']:
            formula = mol['formula']
            mol_id = mol['id']
            atom_ids = mol['atom_ids']
            
            name = f"product_{formula}_{mol_id}"
            print(f"  生成分子: {name}")
            print(f"    原子ID: {atom_ids}")
            
            mol2_content = generate_mol2_for_molecule(
                atom_ids,
                product_data['coords'],
                product_data['bonds'],
                name
            )
            
            output_file = os.path.join(output_dir, f"{name}.mol2")
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(mol2_content)
            print(f"  ✓ 已生成: {output_file}")
            
            atom_set = set(atom_ids)
            bond_count = sum(1 for bond in product_data['bonds'] 
                           if bond['a'] in atom_set and bond['b'] in atom_set)
            print(f"    原子数: {len(atom_ids)}, 键数: {bond_count}")
    else:
        print("  ⚠ 未找到产物数据")
    
    print("\n" + "=" * 50)
    print(f"完成! 共生成 {len(reactant_data['molecules']) if reactant_data else 0} 个反应物和 "
          f"{len(product_data['molecules']) if product_data else 0} 个产物的MOL2文件")
    print(f"文件保存在: {os.path.abspath(output_dir)}")

def main():
    if len(sys.argv) < 2:
        print("用法: python reaction_to_mol2.py <json文件路径> [输出目录]")
        print("示例: python reaction_to_mol2.py reaction.json")
        print("示例: python reaction_to_mol2.py reaction.json ./mol2_output")
        return
    
    json_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    
    json_to_mol2(json_file, output_dir)

if __name__ == "__main__":
    main()
