#!/bin/bash

#Run the DessimozLab fold_tree Snakemake workflow on NCR structures.
#Install dependencies from the github repository: https://github.com/DessimozLab/fold_tree/

#Update the config file with the correct paths
python3 Structural_analysis/update_foldtree_config.py \
  Structural_analysis/fold_tree/workflow/config/config_vars.yaml

#Prepare the input structures (normally done before when running Foldseek)
python3 Structural_analysis/prepare_foldtree_structs.py \
  --models Structural_analysis/AlphaFold_models \
  --output Structural_analysis/foldtree_workspace/structs

mkdir -p Structural_analysis/foldtree_workspace
mkdir -p Structural_analysis/Newick_files

touch Structural_analysis/foldtree_workspace/identifiers.txt

snakemake --cores 4 --use-conda -s Structural_analysis/fold_tree_repo/workflow/fold_tree \
  --config folder="$(cd Structural_analysis/foldtree_workspace && pwd)" \
  filter=False custom_structs=True

cp Structural_analysis/foldtree_workspace/foldtree_struct_tree.PP.nwk.rooted.final \
   Structural_analysis/Newick_files/foldtree_struct_tree.PP.nwk.rooted.final

#The above analysis was for all the structures
#Do the same for each supercluster also using the folders created at the end of foldseek analysis
#I used iTOL for the annotation of the trees: https://itol.embl.de/
#The annotation files are in Structural_analysis/iTol_annotation_files