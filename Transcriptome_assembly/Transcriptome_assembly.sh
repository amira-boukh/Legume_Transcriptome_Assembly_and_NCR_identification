#!/bin/bash

#Nodule and root transcriptome assembly of Indigofera argentea

#Check the quality of our data
mkdir Indigo_fastqc
fastqc fastq_after_trimming/* -o Indigo_fastqc -t 16

#Clean the data with fastp 
#Run the fastp command for each sample (3 nodule samples and 3 roots samples)
mkdir fastq_after_fastp
fastp -i fastq_after_trimming/IndigoNod4_S4_all_R1_001_cutadapt.fastq.gz -o fastq_after_fastp/out_IndigoNod4_S4_all_R1_001_cutadapt.fastq.gz -I fastq_after_trimming/IndigoNod4_S4_all_R2_001_cutadapt.fastq.gz -O fastq_after_fastp/out_IndigoNod4_S4_all_R2_001_cutadapt.fastq.gz -w 16 -V -h -r --cut_right -W 4 -M 15 -c --detect_adapter_for_pe

#Check the quality of preprocessed data
mkdir fastqc_after_fastp
fastqc fastq_after_fastp/* -o fastqc_after_fastp -t 16

#Run Rcorrector to find the errors 
#Nodules
rcorrector -t 48 -1 fastq_after_fastp/out_IndigoNod4_S4_all_R1_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoNod5_S5_all_R1_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoNod7_S6_all_R1_001_cutadapt.fastq.gz \
-2 fastq_after_fastp/out_IndigoNod4_S4_all_R2_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoNod5_S5_all_R2_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoNod7_S6_all_R2_001_cutadapt.fastq.gz \
-od Rcorrector
#Roots
rcorrector -t 48 -1 fastq_after_fastp/out_IndigoRoots4_S1_all_R1_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoRoots5_S5_all_R1_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoRoots7_S6_all_R1_001_cutadapt.fastq.gz \
-2 fastq_after_fastp/out_IndigoRoots4_S1_all_R1_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoRoots5_S5_all_R1_001_cutadapt.fastq.gz,fastq_after_fastp/out_IndigoRoots7_S6_all_R1_001_cutadapt.fastq.gz \
-od Rcorrector

#Run to fix the errors for each sample
#First download the FilterUncorrectabledPEfastq.py 
python2 FilterUncorrectabledPEfastq.py -1 IndigoNod4_S4_all_R1_001_cutadapt.cor.fq.gz -2 IndigoNod4_S4_all_R2_001_cutadapt.cor.fq.gz -s ../FilterUncorrectabledPEfastq.py

#Run trimgalore for each sample 
perl run_trim_galore.pl --paired --fastqc -j 12 --phred33 --output_dir TrimGalore --length 25 -q 20 --stringency 1 -e 0.1 unfixrm_IndigoNod4_S4_all_R1_001_cutadapt.cor.fq.gz unfixrm_IndigoNod4_S4_all_R2_001_cutadapt.cor.fq.gz

#Assess the quality of the preprocessed reads 
mkdir fastqc_preprocessed
fastqc TrimGalore/* -o fastqc_preprocessed -t 16

#Remove possible bacterial contamination from I. argentea root and nodule reads
#Do the same for each sample also
bowtie2-build Bradyrhizobium_elkanii_SA281.fasta Brady_SA281_index 
bowtie2 -x Brady_SA281_index -1 TrimGalore/unfixrm_IndigoNod4_S4_all_R1_001_cutadapt.cor_val_1.fq.gz -2 TrimGalore/unfixrm_IndigoNod4_S4_all_R2_001_cutadapt.cor_val_2.fq.gz --threads 8 -S Indigo_nodS4_brady.sam
samtools view -bS Indigo_nodS4_brady.sam > Indigo_nodS4_brady.bam
samtools view -b -f 12 -F 256 Indigo_nodS4_brady.bam > Indigo_nodS4_brady_unmapped.bam
samtools fastq -1 Indigo_nodS4_HostRemoved_R1.fastq -2 Indigo_nodS4_HostRemoved_R2.fastq -n Indigo_nodS4_brady_unmapped.bam

#Run nodule and root transcriptome assembly 
Trinity --seqType fq --left  Indigo_nodS4_HostRemoved_R1.fastq,Indigo_nodS5_HostRemoved_R1.fastq,Indigo_nodS6_HostRemoved_R1.fastq \
--right  Indigo_nodS4_HostRemoved_R2.fastq,Indigo_nodS5_HostRemoved_R2.fastq,Indigo_nodS6_HostRemoved_R2.fastq \
--trimmomatic --quality_trimming_params "SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:25" --max_memory 100G --CPU 24
#Same for roots

#Check the alignment rate of assemblies with bowtie2
bowtie2 -p 10 -q --no-unal -k 20 -x Indigo_Nods/Trinity.fasta -1 Indigo_nodS4_HostRemoved_R1.fastq,Indigo_nodS5_HostRemoved_R1.fastq,Indigo_nodS6_HostRemoved_R1.fastq -2 Indigo_nodS4_HostRemoved_R2.fastq,Indigo_nodS5_HostRemoved_R2.fastq,Indigo_nodS6_HostRemoved_R2.fastq 2>align_stats_nods.txt| samtools view -@10 -Sb -o bowtie2_nods.bam 
samtools sort -@10 bowtie2_nods.bam -o bowtie2_nods_sorted.bam
samtools index bowtie2_nods_sorted.bam
#Do the same for roots
#run Trinity_stats.pl to check basic stats

#Merge assemblies (root and nodule) with CD-HIT to remove redundancies
cat Indigo_Nods/Trinity.fasta Indigo_Roots/Trinity.fasta > merged_Indigo_transcripts.fasta
#Run CD-HIT-EST to cluster the sequences at 95% identity
cd-hit-est -i merged_Indigo_transcripts.fasta -o Indigo_merged_transcripts_cdhit95.fasta -c 0.95 -T 16 
bowtie2 -p 10 -q --no-unal -k 20 -x Indigo_merged_transcripts_cdhit95.fasta -1 Indigo_nodS4_HostRemoved_R1.fastq,Indigo_nodS5_HostRemoved_R1.fastq,Indigo_nodS6_HostRemoved_R1.fastq,Indigo_RootsS4_HostRemoved_R1.fastq,Indigo_RootsS5_HostRemoved_R1.fastq,Indigo_RootsS6_HostRemoved_R1.fastq -2 Indigo_nodS4_HostRemoved_R2.fastq,Indigo_nodS5_HostRemoved_R2.fastq,Indigo_nodS6_HostRemoved_R2.fastq,Indigo_RootsS4_HostRemoved_R2.fastq,Indigo_RootsS5_HostRemoved_R2.fastq,Indigo_RootsS6_HostRemoved_R2.fastq 2>align_stats_nods.txt| samtools view -@10 -Sb -o bowtie2_merged.bam 
samtools sort -@10 bowtie2_merged.bam -o bowtie2_merged_sorted.bam
samtools index bowtie2_merged_sorted.bam

#Cluster transcripts into gene groups using Corset
corset -D 999999999 bowtie2_merged_sorted.bam
corset -D 999999999 bowtie2_nods_sorted.bam
corset -D 999999999 bowtie2_roots_sorted.bam

#Run Lace to produce supertranscript (clusters.txt is the output from Corset)
Lace --cores 12 -a --outputDir merged/new_lacer Indigo_merged_transcripts_cdhit95.fasta merged/clusters.txt
Lace --cores 12 -a --outputDir Indigo_Nods/new_lacer Indigo_Nods/Trinity.fasta Indigo_Nods/clusters.txt
Lace --cores 12 -a --outputDir Indigo_Roots/new_lacer Indigo_Roots/Trinity.fasta Indigo_Roots/clusters.txt

#Run STAR to check the alignment rate of supertranscripts
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir star_index_merged --genomeFastaFiles merged/new_lacer/SuperDuper.fasta --sjdbGTFfile merged/new_lacer/SuperDuper.gff --sjdbGTFtagExonParentTranscript gene_id --limitGenomeGenerateRAM 61000000000 --genomeSAindexNbases 12
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir star_index_nods --genomeFastaFiles Indigo_Nods/new_lacer/SuperDuper.fasta --sjdbGTFfile Indigo_Nods/new_lacer/SuperDuper.gff --sjdbGTFtagExonParentTranscript gene_id --limitGenomeGenerateRAM 61000000000 --genomeSAindexNbases 12
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir star_index_roots --genomeFastaFiles Indigo_Roots/new_lacer/SuperDuper.fasta --sjdbGTFfile Indigo_Roots/new_lacer/SuperDuper.gff --sjdbGTFtagExonParentTranscript gene_id --limitGenomeGenerateRAM 61000000000 --genomeSAindexNbases 12

mkdir -p merged/star

#STAR two-pass alignment as separate runs
#First pass: discover splice junctions only (no BAM output)
STAR --runMode alignReads --runThreadN 48 --genomeDir star_index_merged \
--outSAMtype None --alignIntronMax 6000 --alignIntronMin 60 \
--outSAMunmapped Within --outFileNamePrefix merged/star/first_pass_ \
--readFilesIn \
Indigo_nodS4_HostRemoved_R1.fastq,Indigo_nodS5_HostRemoved_R1.fastq,Indigo_nodS6_HostRemoved_R1.fastq,Indigo_RootsS4_HostRemoved_R1.fastq,Indigo_RootsS5_HostRemoved_R1.fastq,Indigo_RootsS6_HostRemoved_R1.fastq \
Indigo_nodS4_HostRemoved_R2.fastq,Indigo_nodS5_HostRemoved_R2.fastq,Indigo_nodS6_HostRemoved_R2.fastq,Indigo_RootsS4_HostRemoved_R2.fastq,Indigo_RootsS5_HostRemoved_R2.fastq,Indigo_RootsS6_HostRemoved_R2.fastq

#Second pass: provide junctions from first pass and output sorted BAM
STAR --runMode alignReads --runThreadN 48 --genomeDir star_index_merged \
--sjdbFileChrStartEnd merged/star/first_pass_SJ.out.tab \
--outSAMtype BAM SortedByCoordinate --alignIntronMax 6000 --alignIntronMin 60 \
--outSAMunmapped Within --outFileNamePrefix merged/star/second_pass_ \
--readFilesIn \
Indigo_nodS4_HostRemoved_R1.fastq,Indigo_nodS5_HostRemoved_R1.fastq,Indigo_nodS6_HostRemoved_R1.fastq,Indigo_RootsS4_HostRemoved_R1.fastq,Indigo_RootsS5_HostRemoved_R1.fastq,Indigo_RootsS6_HostRemoved_R1.fastq \
Indigo_nodS4_HostRemoved_R2.fastq,Indigo_nodS5_HostRemoved_R2.fastq,Indigo_nodS6_HostRemoved_R2.fastq,Indigo_RootsS4_HostRemoved_R2.fastq,Indigo_RootsS5_HostRemoved_R2.fastq,Indigo_RootsS6_HostRemoved_R2.fastq

#Do the same for nodules and roots separately

#Assess completeness with BUSCO 
#Using both Fabales and Viridiplantae lineage datasets (transcriptome mode)
mkdir -p busco

#Trinity (merged) vs Fabales
busco -i Indigo_merged_transcripts_cdhit95.fasta -l fabales_odb10 \
-o merged_trinity_fabales -m transcriptome --cpu 48 --out_path busco

#Trinity (merged) vs Viridiplantae
busco -i Indigo_merged_transcripts_cdhit95.fasta -l viridiplantae_odb10 \
-o merged_trinity_viridiplantae -m transcriptome --cpu 48 --out_path busco

#SuperTranscripts (merged) vs Fabales
busco -i merged/new_lacer/SuperDuper.fasta -l fabales_odb10 \
-o merged_super_fabales -m transcriptome --cpu 48 --out_path busco

#SuperTranscripts (merged) vs Viridiplantae
busco -i merged/new_lacer/SuperDuper.fasta -l viridiplantae_odb10 \
-o merged_super_viridiplantae -m transcriptome --cpu 48 --out_path busco

#Generate combined BUSCO summary plot across all runs (Figure S2 in the paper)
generate_plot.py -wd busco

#Do the same for nodules and roots separately

#Predict coding regions (ORFs) with TransDecoder and annotate
#Requires: TransDecoder v5.7.1, BLAST+ v2.12.0+, GFAP installed
#Inputs expected:
#UniProt (release-2023_05) protein FASTA: uniprot_sprot_2023_05.fasta
#Known NCR peptides FASTA: known_NCR_peptides.fasta

mkdir -p annotation/db annotation/transdecoder annotation/blast annotation/gfap

#Build combined UniProt+NCR protein database for BLASTP
cat uniprot_sprot_2023_05.fasta known_NCR_peptides.fasta > annotation/db/uniprot_sprot_2023_05_plus_NCR.fasta
makeblastdb -in annotation/db/uniprot_sprot_2023_05_plus_NCR.fasta -dbtype prot -parse_seqids -out annotation/db/uniprot_sprot_2023_05_plus_NCR

#TransDecoder on merged Trinity transcripts
TransDecoder.LongOrfs -t Indigo_merged_transcripts_cdhit95.fasta -O annotation/transdecoder/merged_trinity
blastp -query annotation/transdecoder/merged_trinity/longest_orfs.pep \
  -db annotation/db/uniprot_sprot_2023_05_plus_NCR \
  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 24 \
  -out annotation/blast/merged_trinity_vs_uniprotNCR.blastp.outfmt6
TransDecoder.Predict -t Indigo_merged_transcripts_cdhit95.fasta -O annotation/transdecoder/merged_trinity \
  --retain_blastp_hits annotation/blast/merged_trinity_vs_uniprotNCR.blastp.outfmt6

##SuperTranscripts coding region prediction following Trinity util-based flow
#Convert GFF to GFF3
gtf_to_alignment_gff3.pl merged/new_lacer/SuperDuper.gff > merged/new_lacer/SuperDuper.gff3
#Extract cDNA sequences from supertranscripts
gtf_genome_to_cdna_fasta.pl merged/new_lacer/SuperDuper.gff merged/new_lacer/SuperDuper.fasta > merged/new_lacer/transcript.fasta
#TransDecoder.LongOrfs
TransDecoder.LongOrfs -t merged/new_lacer/transcript.fasta -O annotation/transdecoder/merged_super
#BLASTP against UniProt+NCR
blastp -query annotation/transdecoder/merged_super/longest_orfs.pep \
  -db annotation/db/uniprot_sprot_2023_05_plus_NCR \
  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 24 \
  -out annotation/blast/merged_super_vs_uniprotNCR.blastp.outfmt6
#Final prediction with single_best_only and retain BLAST hits
TransDecoder.Predict --single_best_only -t merged/new_lacer/transcript.fasta -O annotation/transdecoder/merged_super \
  --retain_blastp_hits annotation/blast/merged_super_vs_uniprotNCR.blastp.outfmt6
#Map ORFs back to supertranscripts (GFF3 with ORFs)
cdna_alignment_orf_to_genome_orf.pl \
  annotation/transdecoder/merged_super/transcript.fasta.transdecoder.gff3 \
  merged/new_lacer/SuperDuper.gff3 \
  merged/new_lacer/transcript.fasta \
  > annotation/transdecoder/merged_super/supertranscripts.wOrfs.gff3

#GFAP functional annotation using Glycine max as reference
#Adjust the GFAP executable and species key to your installation if needed
TRINITY_PEP=$(ls annotation/transdecoder/merged_trinity/*.transdecoder.pep | head -n 1)
SUPER_PEP=$(ls annotation/transdecoder/merged_super/*.transdecoder.pep | head -n 1)

#If GFAP is installed as a Python script on PATH
GFAP.py -i "$TRINITY_PEP" -o annotation/gfap/merged_trinity -s Glycine_max -t 24
GFAP.py -i "$SUPER_PEP"   -o annotation/gfap/merged_super   -s Glycine_max -t 24


#Nodule transcriptome assembly of Lupinus species was performed with the same method
#Refer to the Supplementary dataset S1 to find the public RNAseq data and use the same method
#Only nodule RNAseq datasets were used for the assembly in Lupinus species
#The root RNAseq data were used in downstream analysis to filter the nodule-specific NCR peptides 
