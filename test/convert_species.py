#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import os


def generate_frame_matrix(csv_file_path: str, output_file: str = "frame_matrix.csv"):
    try:
        data = np.genfromtxt(
            csv_file_path, 
            delimiter=',', 
            skip_header=1, 
            dtype=[('frame', int), ('molecule_id', int), ('formula', 'U100'), ('count', int)],
            encoding='utf-8'
        )
    except FileNotFoundError:
        print(f"Error: File '{csv_file_path}' not found")
        return
    except Exception as e:
        print(f"Error: Failed to read CSV - {e}")
        return

    frames = np.unique(data['frame'])
    formulas = np.unique(data['formula'])
    formula_to_idx = {f: i for i, f in enumerate(formulas)}
    
    matrix = np.zeros((len(frames), len(formulas)), dtype=int)
    
    for row in data:
        frame_idx = np.where(frames == row['frame'])[0][0]
        formula_idx = formula_to_idx[row['formula']]
        matrix[frame_idx, formula_idx] = row['count']
    
    try:
        with open(output_file, 'w', encoding='utf-8', newline='\n') as f:
            f.write(','.join(formulas) + '\n')
            for row in matrix:
                f.write(','.join(map(str, row)) + '\n')
        print(f"[OK] Generated: {output_file}")
        print(f"    Frames: {len(frames)}, Formulas: {len(formulas)}")
        print(f"    Output: {os.path.abspath(output_file)}")
    except Exception as e:
        print(f"Error: Failed to write file - {e}")
        return


def main():
    parser = argparse.ArgumentParser(
        description='Generate frame-species matrix from new species_count.csv'
    )
    parser.add_argument(
        'csv_file',
        help='Input CSV file path (species_count.csv)'
    )
    parser.add_argument(
        '-o', '--output',
        default='frame_matrix.csv',
        help='Output matrix file path (default: frame_matrix.csv)'
    )
    
    args = parser.parse_args()
    
    print(f"Processing: {args.csv_file}")
    print("=" * 50)
    generate_frame_matrix(args.csv_file, args.output)
    print("=" * 50)


if __name__ == "__main__":
    main()
