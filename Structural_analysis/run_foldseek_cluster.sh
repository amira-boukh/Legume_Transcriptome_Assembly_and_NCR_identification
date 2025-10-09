#!/bin/bash

#Foldseek clustering of NCR AlphaFold structures using direct commands
#Input: Structural_analysis/AlphaFold_models (directory with AlphaFold PDB models)
#Output: Structural_analysis/Foldseek_results (directory with clustering results)

#Install Foldseek:
conda create -n foldseek -c conda-forge -c bioconda foldseek
conda activate foldseek

#Or download precompiled binaries
#Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xvzf foldseek-linux-avx2.tar.gz
export PATH=$(pwd)/foldseek/bin/:$PATH

#Prepare structures for Foldseek
python3 Structural_analysis/prepare_foldtree_structs.py \
  --models Structural_analysis/AlphaFold_models \
  --output Structural_analysis/foldtree_workspace/structs

mkdir -p Structural_analysis/Foldseek_results
mkdir -p Structural_analysis/foldseek_tmp

#Create database and run clustering
foldseek createdb \
  Structural_analysis/foldtree_workspace/structs \
  Structural_analysis/Foldseek_results/ncr_structures_db

#Run clustering with Foldseek cluster (Foldseek easy-cluster was also tested)
#Different parameters were tested for more or less stringent clustering 
foldseek cluster \
  Structural_analysis/Foldseek_results/ncr_structures_db \
  Structural_analysis/Foldseek_results/ncr_structures_cluster \
  Structural_analysis/foldseek_tmp \
  --tmscore-threshold 0.5 \
  --alignment-type 1 \
  --cov-mode 0 \
  -c 0.8

#Create TSV file with cluster representatives and members
foldseek createtsv \
  Structural_analysis/Foldseek_results/ncr_structures_db \
  Structural_analysis/Foldseek_results/ncr_structures_db \
  Structural_analysis/Foldseek_results/ncr_structures_cluster \
  Structural_analysis/Foldseek_results/ncr_structures_clusters.tsv

rm -rf Structural_analysis/Foldseek_superclusters
mkdir -p Structural_analysis/Foldseek_superclusters

#Copy structures into folders according to superclusters
while read -r representative member; do
  cluster_folder="Structural_analysis/Foldseek_superclusters/${representative}"
  mkdir -p "$cluster_folder"
  for structure_id in "$representative" "$member"; do
    source_pdb="Structural_analysis/foldtree_workspace/structs/${structure_id}.pdb"
    if [ -f "$source_pdb" ]; then
      cp "$source_pdb" "$cluster_folder/"
    fi
  done
done < Structural_analysis/Foldseek_results/ncr_structures_clusters.tsv

rm -rf Structural_analysis/foldseek_tmp

#This command was used to align two structures and visualize the results like in figure 4B
#Example
foldseek easy-search Foldseek_results/SC312/NCR247.pdb Foldseek_results/SC37/Indigo_0000.pdb result_NCR247_Indigo.html foldseek_tmp --format-mode 3

#Command used to compute TM score inside superclusters
foldseek easy-search Foldseek_results/SC312/ Foldseek_results/SC312/ result_SC312.tsv foldseek_tmp --alignment-type 1 

#Command used to compute TM score between different superclusters
foldseek easy-search Foldseek_results/SC312/ Foldseek_results/SC37/ result_SC312_vs_SC37.tsv foldseek_tmp --alignment-type 1
