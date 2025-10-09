#!/usr/bin/env python3

"""
Copy the top-ranked AlphaFold model for each cluster into a foldtree workspace
Usage:
    python prepare_foldtree_structs.py --models Structural_analysis/AlphaFold_models --output foldtree/structs
"""

import argparse
import json
import shutil
from pathlib import Path


def find_best_model_dir(cluster_dir: Path) -> Path:
    ranking_file = cluster_dir / "ranking_debug.json"
    
    data = json.loads(ranking_file.read_text())
    order = data.get("order")


    best_index = 0
    ranked_candidate = cluster_dir / f"ranked_{best_index}.pdb"
    if ranked_candidate.exists():
        return ranked_candidate

    relaxed_candidate = cluster_dir / f"relaxed_model_{best_index + 1}.pdb"
    if relaxed_candidate.exists():
        return relaxed_candidate

    pdb_files = list(cluster_dir.glob("*.pdb"))
    
    return pdb_files[0]


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare structures for foldtree and foldseek")
    parser.add_argument("--models", required=True, help="Directory with AlphaFold model folders")
    parser.add_argument("--output", required=True, help="Destination structs directory")
    args = parser.parse_args()

    models_dir = Path(args.models)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    for existing in output_dir.glob("*.pdb"):
        existing.unlink()

    for cluster_dir in sorted(models_dir.iterdir()):
        model_path = find_best_model_dir(cluster_dir)
        

        destination = output_dir / f"{cluster_dir.name}.pdb"
        shutil.copyfile(model_path, destination)


if __name__ == "__main__":
    main()
