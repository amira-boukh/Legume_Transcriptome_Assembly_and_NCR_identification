#!/bin/bash

#Structural alignment workflow for each Foldseek supercluster using MUSTANG

mkdir -p Structural_analysis/Structural_alignments

for cluster_dir in Structural_analysis/Foldseek_superclusters/*/; do
  cluster_name=$(basename "$cluster_dir")
  alignment_dir="Structural_analysis/Structural_alignments/${cluster_name}"
  mkdir -p "$alignment_dir"

  mapfile -t pdb_files < <(ls "$cluster_dir"/*.pdb 2>/dev/null)
  if [ "${#pdb_files[@]}" -eq 0 ]; then
    continue
  fi

  representative=$(ls -S "$cluster_dir"/*.pdb | head -n 1)
  ordered_pdbs=("$representative")
  for pdb in "${pdb_files[@]}"; do
    if [ "$pdb" != "$representative" ]; then
      ordered_pdbs+=("$pdb")
    fi
  done

  MUSTANG -f "${ordered_pdbs[@]}" -o "$alignment_dir"

  aligned_pdb="$alignment_dir/aligned.pdb"
  if [ ! -f "$aligned_pdb" ]; then
    echo "No aligned.pdb produced for $cluster_name" >&2
    continue
  fi

  pymol_script="${alignment_dir}/${cluster_name}_alignment.pml"
  {
    echo "reinitialize"
    echo "bg_color white"
    echo "load ${aligned_pdb}, ${cluster_name}"
    echo "split_states ${cluster_name}, ${cluster_name}"
    echo "delete ${cluster_name}"
    echo "set_name ${cluster_name}_0001, ${cluster_name}_representative"
    echo "color red, ${cluster_name}_representative"
    echo "show cartoon, ${cluster_name}_*"
    echo "cartoon putty, ${cluster_name}_*"
    echo "spectrum b, ${cluster_name}_*"
    echo "orient ${cluster_name}_*"
    echo "ray 2000,1500"
    echo "png ${cluster_name}_alignment.png, dpi=300"
    echo "quit"
  } > "$pymol_script"

done

#USalign was also tested