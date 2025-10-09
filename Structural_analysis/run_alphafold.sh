#!/bin/bash

#Run AlphaFold2 for each NCR cluster until one sequence reaches mean pLDDT equal or higher than 70
#Update the AlphaFold paths before use


mkdir -p Structural_analysis/AlphaFold_models

for cluster_dir in Structural_analysis/Cluster_fastas/*/; do
  cluster_name=$(basename "$cluster_dir")
  msa_file="Structural_analysis/Cluster_msas/${cluster_name}.a3m"

  rm -rf "Structural_analysis/AlphaFold_models/${cluster_name}" \
         Structural_analysis/AlphaFold_models/${cluster_name}_attempt_*

  attempt=1

  for sequence_fasta in "${cluster_dir}"*.fasta; do
    out_dir="Structural_analysis/AlphaFold_models/${cluster_name}_attempt_${attempt}"
    rm -rf "$out_dir"
    mkdir -p "$out_dir"

    python /path/to/alphafold/run_alphafold.py \
      --data_dir=/path/to/alphafold/data \
      --output_dir="$out_dir" \
      --fasta_paths="$sequence_fasta" \
      --msa_paths="$msa_file" \
      --model_preset=monomer \
      --db_preset=full_dbs \
      --max_template_date=2023-12-31 \
      --uniref90_database_path=/path/to/databases/uniref90/uniref90.fasta \
      --mgnify_database_path=/path/to/databases/mgnify/mgy_clusters_2018_12.fa \
      --uniref30_database_path=/path/to/databases/uniref30/uniref30_2103_db.fasta \
      --bfd_database_path=/path/to/databases/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
      --pdb70_database_path=/path/to/databases/pdb70/pdb70 \
      --template_mmcif_dir=/path/to/databases/pdb_mmcif/mmcif_files \
      --obsolete_pdbs_path=/path/to/databases/pdb_mmcif/obsolete.dat \
      --pdb_seqres_database_path=/path/to/databases/pdb_seqres/pdb_seqres.txt \
      --use_gpu_relax=true

    result=$(python Structural_analysis/check_plddt.py "$out_dir/ranking_debug.json")
    if [ "$result" = "True" ]; then
      rm -rf "Structural_analysis/AlphaFold_models/${cluster_name}"
      mv "$out_dir" "Structural_analysis/AlphaFold_models/${cluster_name}"
      break
    fi

    rm -rf "$out_dir"
    attempt=$((attempt + 1))
  done

done
