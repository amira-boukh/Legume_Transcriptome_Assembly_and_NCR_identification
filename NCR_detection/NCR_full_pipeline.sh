#!/usr/bin/env bash

#peptide and CDS sequences per NCR cluster
#Use the per-cluster FASTAs produced in Orthology_clustering/KnownNCRs_classification.sh under results/IRLC and results/DAL

#IRLC: build species CDS and extract CDS per cluster, run SignalP notm per cluster, align and build HMMs per cluster
mkdir -p results/IRLC/cds results/IRLC/clusters_cds results/IRLC/alignments results/IRLC/profiles_mature results/IRLC/signalp results/IRLC/mature results/IRLC/clusters_prot
gffread -g data/Medicago_truncatula/genome.fa -x results/IRLC/cds/Medicago_truncatula.cds.fa -y results/IRLC/cds/Medicago_truncatula.prot.fa data/Medicago_truncatula/genes.gtf
gffread -g data/Medicago_sativa/genome.fa      -x results/IRLC/cds/Medicago_sativa.cds.fa      -y results/IRLC/cds/Medicago_sativa.prot.fa      data/Medicago_sativa/genes.gtf
gffread -g data/Pisum_sativum/genome.fa        -x results/IRLC/cds/Pisum_sativum.cds.fa        -y results/IRLC/cds/Pisum_sativum.prot.fa        data/Pisum_sativum/genes.gtf
gffread -g data/Cicer_arietinum/genome.fa      -x results/IRLC/cds/Cicer_arietinum.cds.fa      -y results/IRLC/cds/Cicer_arietinum.prot.fa      data/Cicer_arietinum/genes.gtf
cat results/IRLC/cds/*.cds.fa > results/IRLC/all_species.cds.fa
for cf in results/IRLC/clusters_prot/*.fasta; do bn=$(basename "$cf" .fasta); grep "^>" "$cf" | sed 's/^>//; s/[\t ].*$//' > results/IRLC/clusters_cds/${bn}.ids; seqkit grep -f results/IRLC/clusters_cds/${bn}.ids results/IRLC/all_species.cds.fa > results/IRLC/clusters_cds/${bn}.fna; signalp -s notm -t euk -f summary -n results/IRLC/signalp/${bn}.summary -m results/IRLC/mature/${bn}.mature.faa -v "$cf" > results/IRLC/signalp/${bn}.short; macse -prog alignSequences -seq results/IRLC/clusters_cds/${bn}.fna -out_NT results/IRLC/alignments/${bn}.cds.macse.fna -out_AA results/IRLC/alignments/${bn}.aa.macse.faa; hmmbuild results/IRLC/profiles_mature/${bn}.hmm results/IRLC/alignments/${bn}.aa.macse.faa; done
cat results/IRLC/profiles_mature/*.hmm > profiles/IRLC.profiles.hmm
hmmpress profiles/IRLC.profiles.hmm
perl spada_soft/spada/build_profile.pl -i results/IRLC/alignments -o results/SPADA/IRLC_profiles

#Dalbergioids: build species CDS and extract CDS per cluster, run SignalP notm per cluster, align and build HMMs per cluster
mkdir -p results/DAL/cds results/DAL/clusters_cds results/DAL/alignments results/DAL/profiles_mature results/DAL/signalp results/DAL/mature results/DAL/clusters_prot
gffread -g data/Aeschynomene_evenia/genome.fa   -x results/DAL/cds/Aeschynomene_evenia.cds.fa   -y results/DAL/cds/Aeschynomene_evenia.prot.fa   data/Aeschynomene_evenia/genes.gtf
gffread -g data/Arachis_hypogaea/genome.fa      -x results/DAL/cds/Arachis_hypogaea.cds.fa      -y results/DAL/cds/Arachis_hypogaea.prot.fa      data/Arachis_hypogaea/genes.gtf
cat results/DAL/cds/*.cds.fa > results/DAL/all_species.cds.fa
for cf in results/DAL/clusters_prot/*.fasta; do bn=$(basename "$cf" .fasta); grep "^>" "$cf" | sed 's/^>//; s/[\t ].*$//' > results/DAL/clusters_cds/${bn}.ids; seqkit grep -f results/DAL/clusters_cds/${bn}.ids results/DAL/all_species.cds.fa > results/DAL/clusters_cds/${bn}.fna; signalp -s notm -t euk -f summary -n results/DAL/signalp/${bn}.summary -m results/DAL/mature/${bn}.mature.faa -v "$cf" > results/DAL/signalp/${bn}.short; macse -prog alignSequences -seq results/DAL/clusters_cds/${bn}.fna -out_NT results/DAL/alignments/${bn}.cds.macse.fna -out_AA results/DAL/alignments/${bn}.aa.macse.faa; hmmbuild results/DAL/profiles_mature/${bn}.hmm results/DAL/alignments/${bn}.aa.macse.faa; done
cat results/DAL/profiles_mature/*.hmm > profiles/DAL.profiles.hmm
hmmpress profiles/DAL.profiles.hmm
perl spada_soft/spada/build_profile.pl -i results/DAL/alignments -o results/SPADA/DAL_profiles

#CRP: prepare CRP profiles. If alignments available build with build_profile else use provided hmm.crp
mkdir -p results/CRP/alignments results/SPADA/CRP_profiles
perl spada_soft/spada/build_profile.pl -i results/CRP/alignments -o results/SPADA/CRP_profiles
cat spada_soft/spada/hmm.crp > profiles/CRP.profiles.hmm
hmmpress profiles/CRP.profiles.hmm

#SPADA pipeline
#Genome sanity check (do it for each genome/transcriptome to analyze)
perl spada_soft/spada/seq.check.pl -i data/Medicago_truncatula/genome.fa -o data/Medicago_truncatula/genome.checked.fa

#Profiles already built per clade above. Use profiles/IRLC.profiles.hmm, profiles/DAL.profiles.hmm, profiles/CRP.profiles.hmm

#Run SPADA three times per genome/transcriptome (IRLC, Dalbergioids, CRP)
perl spada_soft/spada/spada.pl -c spada_soft/spada/cfg.txt -d results/SPADA/Medtr_IRLC --fas data/Medicago_truncatula/genome.checked.fa --hmm profiles/IRLC.profiles.hmm --org Mtruncatula --sp --thread 8
perl spada_soft/spada/spada.pl -c spada_soft/spada/cfg.txt -d results/SPADA/Medtr_DAL  --fas data/Medicago_truncatula/genome.checked.fa --hmm profiles/DAL.profiles.hmm  --org Mtruncatula --sp --thread 8
perl spada_soft/spada/spada.pl -c spada_soft/spada/cfg.txt -d results/SPADA/Medtr_CRP  --fas data/Medicago_truncatula/genome.checked.fa --hmm spada_soft/spada/hmm.crp         --org Mtruncatula --sp --thread 8

#Repeat the three commands above for each target species/transcriptome by editing the paths

#Merge SPADA results from the three runs and extract FASTA; remove duplicates
cat results/SPADA/Medtr_IRLC/61_final.tbl results/SPADA/Medtr_DAL/61_final.tbl results/SPADA/Medtr_CRP/61_final.tbl > results/SPADA/Medtr_merged.tbl
cp results/SPADA/Medtr_merged.tbl 61_final.tbl
python3 NCR_detection/extract_seq.py
mv seq.fasta results/SPADA/Medtr_merged.seq.fasta
seqkit rmdup -s results/SPADA/Medtr_merged.seq.fasta > results/SPADA/Medtr_merged.noDup.fasta

#Post-SPADA filtering
#Self-BLAST to collapse highly similar sequences and remove reciprocal duplicates
makeblastdb -in results/SPADA/Medtr_merged.noDup.fasta -dbtype prot -parse_seqids
blastp -query results/SPADA/Medtr_merged.noDup.fasta -db results/SPADA/Medtr_merged.noDup.fasta -culling_limit 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > results/SPADA/Medtr_merged.self_blast.tsv
python3 NCR_detection/NCR_check.py -i results/SPADA/Medtr_merged.self_blast.tsv | awk '{print $2}' > results/SPADA/Medtr_merged.rm_ids.txt
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_merged.rm_ids.txt results/SPADA/Medtr_merged.noDup.fasta results/SPADA/Medtr_merged.noDup_noHS.fasta

#Signal peptide and mature peptide extraction (SignalP 4.1g notm)
signalp -s notm -t euk -f summary -n results/SPADA/Medtr.signalp41.summary -m results/SPADA/Medtr.mature.fasta -v results/SPADA/Medtr_merged.noDup_noHS.fasta > results/SPADA/Medtr.signalp41.short

#Length < 100 aa and >= 4 Cys in mature peptides
python3 NCR_detection/NCR_final_check.py -i results/SPADA/Medtr.mature.fasta > results/SPADA/Medtr_mature.final_check.txt
grep "crp\|cluster" results/SPADA/Medtr_mature.final_check.txt > results/SPADA/Medtr_mature.rm_ids.txt
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_mature.rm_ids.txt results/SPADA/Medtr.mature.fasta results/SPADA/Medtr.mature.checked.fasta
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_mature.rm_ids.txt results/SPADA/Medtr_merged.noDup_noHS.fasta results/SPADA/Medtr.full.checked.fasta

#RNA-seq support: SRA download, mapping (STAR), and counting (htseq-count)
#Download SRA and convert to FASTQ 
fastq-dump --split-files SRR******  #replace ***** with actual SRA run ID (repeat for each run)

#STAR genome index (example genome and GTF)
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir results/STAR_index/Medtr --genomeFastaFiles data/Medicago_truncatula/genome.fa --sjdbGTFfile data/Medicago_truncatula/genes.gtf

#STAR mapping (paired-end example)
STAR --runThreadN 8 --genomeDir results/STAR_index/Medtr --readFilesIn SRR******_1.fastq SRR******_2.fastq --outFileNamePrefix results/STAR_map/Medtr/ --outSAMtype BAM SortedByCoordinate

#Prepare GTF filtered to retained IDs
grep "^>" results/SPADA/Medtr.full.checked.fasta | sed 's/>//' > results/SPADA/Medtr_kept_ids.txt
awk -F '"' 'NR==FNR { ids[$1]=1; next } ids[$2]' results/SPADA/Medtr_kept_ids.txt data/Medicago_truncatula/genes.gtf > results/SPADA/Medtr_kept.gtf

#Count reads on transcripts (htseq-count)
htseq-count --format=bam --order=pos --stranded=yes --idattr=transcript_id --mode=union results/STAR_map/Medtr/Aligned.sortedByCoord.out.bam results/SPADA/Medtr_kept.gtf > results/SPADA/Medtr_counts_htseq.txt

#Remove unexpressed candidates and keep expressed NCRs
awk -F '\t' '$2 == 0' results/SPADA/Medtr_counts_htseq.txt | awk '{print $1}' > results/SPADA/Medtr_counts_zero_ids.txt
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_counts_zero_ids.txt results/SPADA/Medtr.full.checked.fasta results/SPADA/Medtr.full.checked_RNA.fasta
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_counts_zero_ids.txt results/SPADA/Medtr.mature.checked.fasta results/SPADA/Medtr.mature.checked_RNA.fasta

#Motif-based classification
#For IRLC species: exclude defensin-like (>=8 Cys) from annotation set
python3 NCR_detection/NCR_final_check2.py -i results/SPADA/Medtr.mature.checked_RNA.fasta > results/SPADA/Medtr_mature.motif_check.txt
grep "crp\|cluster" results/SPADA/Medtr_mature.motif_check.txt > results/SPADA/Medtr_mature.motif_rm_ids.txt
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_mature.motif_rm_ids.txt results/SPADA/Medtr.mature.checked_RNA.fasta results/SPADA/Medtr.mature.checked_RNA_motif.fasta
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_mature.motif_rm_ids.txt results/SPADA/Medtr.full.checked_RNA.fasta results/SPADA/Medtr.full.checked_RNA_motif.fasta

#Optional motif reporting (NCR_motif.py)
python3 NCR_detection/NCR_motif.py -i results/SPADA/Medtr.mature.checked_RNA.fasta

#Compute pI values (R package pIR) on a table of IDs and sequences (edit pI.R as needed)
Rscript NCR_detection/pI.R

#Classify CRP-found sequences via hmmsearch against our cluster profiles and merge
#Prepare CRP sequences (unclassified) and search with hmmsearch per cluster
awk '/^>/{print $0; next} !/^>/{print}' results/SPADA/Medtr.mature.checked_RNA.fasta > results/SPADA/Medtr_mature_clean.fasta
mkdir -p results/ALL/profiles_mature results/ALL/hmmsearch_by_cluster
cat results/IRLC/profiles_mature/*.hmm results/DAL/profiles_mature/*.hmm > results/ALL/profiles_mature/all_clusters.hmm
hmmpress results/ALL/profiles_mature/all_clusters.hmm
for hmm in results/IRLC/profiles_mature/*.hmm results/DAL/profiles_mature/*.hmm; do bn=$(basename "$hmm" .hmm); hmmsearch --cpu 8 --tblout results/ALL/hmmsearch_by_cluster/${bn}.tbl "$hmm" results/SPADA/Medtr_mature_clean.fasta; done
cat results/ALL/hmmsearch_by_cluster/*.tbl | grep -v '^#' > results/ALL/hmmsearch_by_cluster/all.tbl

#Keep best hit per sequence across clusters and annotate FASTA headers with best cluster
perl -ane '$key=$F[0];$eval=$F[4];$cluster=$F[2];if(!exists $data{$key}||$eval<$data{$key}[1]){$data{$key}=[$F[1],$eval,$cluster]} END {foreach $key (sort keys %data){print "$key $data{$key}[0] $data{$key}[2]\n";}}' results/ALL/hmmsearch_by_cluster/all.tbl | sed 's/ -//' > results/SPADA/Medtr_CRPs_best_hits.txt
awk 'BEGIN{while(getline<"results/SPADA/Medtr_CRPs_best_hits.txt"){a[$1]=$2}} /^>/{print $0" "a[substr($1,2)]}!/^>/{print}' < results/SPADA/Medtr_mature_clean.fasta | sed 's/\s/\|/' > results/SPADA/Medtr_mature_crps_classified.fasta
perl -ne 'if (/^>/) { ($id) = $_ =~ /^(>\S+)/; $found = 0; if ($id =~ /Cluster/) { $found = 1; print; } } elsif ($found) { print; }' results/SPADA/Medtr_mature_crps_classified.fasta > results/SPADA/Medtr_mature_crps_classified_clusters.fasta

#Split classified CRP sequences per cluster for downstream per-cluster handling
mkdir -p results/SPADA/Medtr_CRPs_classified_by_cluster
seqkit split --by-id --id-regexp "(Cluster_[0-9]+)" results/SPADA/Medtr_mature_crps_classified_clusters.fasta -O results/SPADA/Medtr_CRPs_classified_by_cluster

#Mark sequences with no cluster hit as CladeSpecific
grep '^>' results/SPADA/Medtr_mature_clean.fasta | sed 's/^>//; s/[\t ].*$//' | sort -u > results/SPADA/Medtr_all_ids.txt
cut -d' ' -f1 results/SPADA/Medtr_CRPs_best_hits.txt | sort -u > results/SPADA/Medtr_hit_ids.txt
comm -23 results/SPADA/Medtr_all_ids.txt results/SPADA/Medtr_hit_ids.txt > results/SPADA/Medtr_nohit_ids.txt
awk 'NR==FNR{no[$1]=1; next} /^>/{h=$0; sub(/^>/,"",h); id=h; if(id in no){print ">"$0"|CladeSpecific"} else {print}} !/^>/{print}' results/SPADA/Medtr_nohit_ids.txt results/SPADA/Medtr_mature_clean.fasta > results/SPADA/Medtr_mature_clade_specific_annot.fasta

#Filter clade-specific sequences by NCR criteria and keep as singletons
python3 NCR_detection/NCR_final_check.py -i results/SPADA/Medtr_mature_clade_specific_annot.fasta > results/SPADA/Medtr_clade_specific_check.txt
grep -v "please check" results/SPADA/Medtr_clade_specific_check.txt | sed '/^$/d' > results/SPADA/Medtr_clade_specific_fail_ids.txt
python3 NCR_detection/Exclude_ids.py results/SPADA/Medtr_clade_specific_fail_ids.txt results/SPADA/Medtr_mature_clade_specific_annot.fasta results/SPADA/Medtr_clade_specific_pass.fasta
mkdir -p results/FINAL/singletons/Medtr
seqkit split --by-id results/SPADA/Medtr_clade_specific_pass.fasta -O results/FINAL/singletons/Medtr

#Merge SPADA-classified and hmmsearch-classified NCRs as the final set
cat results/SPADA/Medtr.mature.checked_RNA_motif.fasta results/SPADA/Medtr_mature_crps_classified_clusters.fasta > results/SPADA/Medtr_final_NCRs.fasta

#Special case: I. argentea (re-run SPADA with expanded profiles and species-specific HMM)
perl spada_soft/spada/spada.pl -c spada_soft/spada/cfg.txt -d results/SPADA/Iarg_expanded --fas data/Ionopsidium_argenteum/genome.checked.fa --hmm profiles/expanded_clusters.hmm --org Iargentea --sp --thread 8
perl spada_soft/spada/spada.pl -c spada_soft/spada/cfg.txt -d results/SPADA/Iarg_specific --fas data/Ionopsidium_argenteum/genome.checked.fa --hmm profiles/Iargentea_5peptides.hmm --org Iargentea --sp --thread 8
cat results/SPADA/Iarg_expanded/61_final.tbl results/SPADA/Iarg_specific/61_final.tbl > results/SPADA/Iarg_merged.tbl
cp results/SPADA/Iarg_merged.tbl 61_final.tbl
python3 NCR_detection/extract_seq.py
mv seq.fasta results/SPADA/Iarg_merged.seq.fasta
signalp -s notm -t euk -f summary -n results/SPADA/Iarg.signalp41.summary -m results/SPADA/Iarg.mature.fasta -v results/SPADA/Iarg_merged.seq.fasta > results/SPADA/Iarg.signalp41.short
python3 NCR_detection/NCR_final_check.py -i results/SPADA/Iarg.mature.fasta > results/SPADA/Iarg_mature.final_check.txt
grep "NCR sequence too long\|invalid number of Cys" results/SPADA/Iarg_mature.final_check.txt | awk '{print $NF}' > results/SPADA/Iarg_mature.rm_ids.txt
python3 NCR_detection/Exclude_ids.py results/SPADA/Iarg_mature.rm_ids.txt results/SPADA/Iarg.mature.fasta results/SPADA/Iarg.mature.flex.fasta

#Extract CDS for retained NCRs (per-species example using genome+GTF)
gffread -g data/Medicago_truncatula/genome.fa -x results/SPADA/Medtr.cds.all.fa -y results/SPADA/Medtr.prot.all.fa data/Medicago_truncatula/genes.gtf
grep "^>" results/SPADA/Medtr.full.checked_RNA_motif.fasta | sed 's/>//' > results/SPADA/Medtr_final_ids.txt
seqkit grep -f results/SPADA/Medtr_final_ids.txt results/SPADA/Medtr.cds.all.fa > results/SPADA/Medtr_final_NCRs.cds.fna

#Merge expanded clusters across species using merge_clusters script
bash NCR_detection/Merge_clusters.sh
