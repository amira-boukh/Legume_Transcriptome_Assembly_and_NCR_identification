# Structure-based phylogenetic analysis of NCR peptides

This repository contains the command-line workflows and helper scripts that was used. in the analyses described in the manuscript *Structure-based phylogenetic analysis reveals multiple events of convergent evolution of cysteine-rich antimicrobial peptides in legume-rhizobium symbiosis* (https://doi.org/10.1101/2025.09.09.675119). Each subdirectory contains one stage of the study, from raw RNA-seq processing to structural comparisons.

## The analysis 

1. **Transcriptome_assembly/** 
Bash workflows for cleaning root and nodule RNA-seq reads (fastp, Rcorrector, Trim Galore), removing bacterial contamination, running Trinity de novo assembly, merging assemblies with CD-HIT, clustering with Corset and SuperTranscript, annotating coding sequences with TransDecoder and BLAST/GFAP, and building Salmon indices/quantifications for expression analysis.
2. **DEGs/** 
Expression analysis.  
   - `DEGs_indigo.R` keeps the original code style: Use the salmon quantifications and run  DESeq2 and edgeR analysis and generate plots.  
   - `salmon_merge.py` and `tpm_plot_indigo.R` merge the salmon outputs, calculate TPM and plot box plot visualizations. 
3. **NCR_detection/** 
The SPADA-based pipeline for discovering NCR peptides, filtering on signal peptides/cysteine patterns, nodule expression and validating motifs.
4. **Orthology_clustering/** 
All-vs-all BLAST, OrthAgogue, and MCL clustering scripts for organising known NCR peptides into orthologous groups.
5. **Structural_analysis/** 
AlphaFold modelling, Foldseek clustering, and FoldTree/structural alignment workflows. (`run_alphafold.sh`, `check_plddt.py`, `prepare_foldtree_structs.py`, `update_foldtree_config.py`, `run_foldseek_cluster.sh`, `run_foldseek_foldtree.sh`, `run_structural_alignments.sh`).

## Getting started

The pipelines assume a Unix-like environment with standard bioinformatics tools on your `PATH`. A non-exhaustive list of software used throughout the repository includes:

- Read processing and assembly: `fastqc`, `fastp`, `rcorrector`, `TrimGalore`, `bowtie2`, `samtools`, `Trinity`, `cd-hit-est`, `corset`, `lace`.
- Annotation and peptide discovery: `SignalP`, `TransDecoder`, `BLAST+`, `GFAP`, `SPADA`, `hmmsearch`.
- Structure and clustering: `Foldseek`, `FoldTree`, `fastme`, `quicktree`, `AlphaFold2` 
- Orthology and clustering: `orthAgogue`, `mcl`.
- Expression analysis: `salmon`, `tximport`, `DESeq2`, `edgeR`, `ggplot2`, `pheatmap`.

Some steps were executed on an HPC cluster. Adjust thread counts and memory flags to suit your hardware.

## Reproducing the main analyses

1. **Assemble transcriptomes**  
   Run `Transcriptome_assembly/Transcriptome_assembly.sh`, editing the file paths to your FASTQ inputs. The script is heavily commented so you can enable or skip individual stages (quality control, host read removal, Trinity assembly, clustering, and annotation).

2. **Quantify expression and call DEGs**  
   - Quantify each sample with Salmon (mapping-based mode) against the merged Trinity transcriptome. The core commands are in `Transcriptome_assembly/Transcriptome_assembly.sh`. 

   Adjust the sample IDs/paths if your filenames differ.  
   - Merge/preview TPMs with `DEGs/salmon_merge.py` (e.g. `python DEGs/salmon_merge.py -i Transcriptome_assembly/salmon/quants -o DEGs -p indigo`) to generate replicate summaries along with the `Quants.csv` and `geneID.csv` inputs expected by the R workflow.
   - Edit the paths at the top of `DEGs/DEGs_indigo.R` and run it to generate DESeq2/edgeR DEG tables along with MA, heatmap, volcano, and marker gene plots tailored to Indigofera.

3. **Identify NCR peptides**  
   Use the SPADA wrapper scripts in `NCR_detection/` to screen assemblies and genomes with IRLC, Dalbergioid, and CRP profiles. Post-processing scripts filter candidates by signal peptides, motif composition, length, and nodule-specific expression.

4. **Cluster NCRs across species**  
   Execute `Orthology_clustering/KnownNCRs_classification.sh` to normalise protein identifiers, build full proteomes, compute BLASTP similarities, and cluster sequences with OrthAgogue + MCL. Outputs include cluster FASTA files and ID mappings used later for structural work.

5. **Explore structural convergence**  
   - Edit the hard-coded paths inside `Structural_analysis/run_alphafold.sh`, then run it to iterate through per-sequence FASTA files and keep the first AlphaFold2 model reaching mean pLDDT of 70 or higher.  
   - Execute `Structural_analysis/run_foldseek_cluster.sh` to generate Foldseek structure clusters (and supercluster-specific PDB folders), then use `Structural_analysis/run_structural_alignments.sh` to run MUSTANG on each supercluster and produce PyMOL sausage scripts/images centred on the largest representative.  
   - Launch `Structural_analysis/run_foldseek_foldtree.sh` to clone the DessimozLab `fold_tree` Snakemake workflow, patch its configuration, prepare structures, and run the FoldTree pipeline (Snakemake + Foldseek + `fastme` + `quicktree`) to obtain rooted structural trees saved in `Structural_analysis/Newick_files/ncr_structures_foldtree.nwk`.

## Tips for adapting the workflows

- Scripts rely on relative paths, so run them from within their respective directories or update the variables accordingly.
- When working with new species, start by adjusting the species lists and file naming conventions inside the shell scripts.
- The differential expression pipeline accepts additional covariates via the `--design` option, allowing you to include batch or treatment effects beyond the nodule vs root contrast.

- Manuscript figures can be regenerated with the helper scripts under `Scripts_figures/`. For example, `Rscript Scripts_figures/Figure4/figure4b_violin.R` exports the simplified violin version of Figure 4b using the same count table as the original heatmap.

Please reach out if you have any questions or want to contribute. 
