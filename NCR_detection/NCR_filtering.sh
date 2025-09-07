#Analysis Spada 

#Check genome
perl spada_soft/spada/seq.check.pl -i genome.fa -o genome.checked.fa

#Spada with CRP clusters 
/home/amira.boukherissa/spada_soft/spada/spada.pl -c /home/amira.boukherissa/spada_soft/spada/cfg.txt -d /home/amira.boukherissa/spada_soft/spada/spada.Vfaba.CRP --fas Vfaba.checked.fa --hmm /home/amira.boukherissa/spada_soft/spada/hmm.crp --org Mtruncatula --sp --thread 12

#Spada with IRLC NCR clusters : 


#Spada with Dalbergioids NCR clusters : 


#Analysis after SPADA 

#1- Merge spada all sequences 
#extract sequences from output files
#run extract_seq.py script in the spada output folder. It will generate the seq.fasta file contains all sequences from spada output
#then use https://arn.ugr.es/srnatoolbox/helper/removedup/ to remove any repeated sequences


#2- filtred : 
#A - blast : 
blastp -query *_spada_high_confidence.fasta -db *_spada_high_confidence.fasta -culling_limit 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > *_blast.txt
python3 NCR_check.py -i *_blast.txt | awk '{print $2}' > *_blast_ids.txt
python3 Exclude_ids.py *_blast_ids.txt *_spada_high_confidence.fasta *_spada_high_confidence_noDup.fasta
#B- signalP : 
/home/amira.boukherissa/spada_soft/signalp-4.1/signalp -s notm -t euk -f summary -n *_spada_HC_noDup_signal.gff -m *_spada_HC_noDup_mature.fasta -v *_spada_high_confidence_noDup.fasta > *_signal_short.out
#C- final check (mature length < 100, Cys mature>=4): 
python3 NCR_final_check.py -i *_spada_HC_noDup_mature.fasta > final_check_mature.txt
grep "crp\|cluster" final_check_mature.txt > final_check_mature_ids.txt
python3 Exclude_ids.py final_check_mature_ids.txt *_spada_HC_noDup_mature.fasta *_spada_HC_noDup_mature_checked.fasta
python3 Exclude_ids.py final_check_mature_ids.txt *_spada_high_confidence_noDup.fasta *_spada_high_confidence_noDup_checked.fasta


#3- RNAseq expression : 
cat LupAl_noDup_checked.fasta | sed 's/.*\=/>/' | sed 's/\s.*//' > LupAl_noDup_checked_clean.fasta #Keep only one ID (delete =) 
grep "^>" *_spada_high_confidence_noDup_checked.fasta | sed 's/>//' > ids_kept.txt
awk -F'"' 'NR==FNR { ids[$1]=1;next } ids[$2]' ids_kept.txt ../../../Medtr/61_final.gtf > 61_final_filtred.gtf

#Use htseqcount or featurecounts
htseq-count --format=bam --order=pos --stranded=yes --idattr=transcript_id --mode=union /store/EQUIPES/PBI/Amira_Boukherissa/Genomes/Genomes_spada/*/star_aligned_pass2Aligned.sortedByCoord.out.bam /store/EQUIPES/PBI/Amira_Boukherissa/results_spada/*/new_analysis/61_final_filtred.gtf > /store/EQUIPES/PBI/Amira_Boukherissa/results_spada/*/new_analysis/*_counts_htseq.txt

featureCounts -T 8 -s 0 -p -t transcript -g transcript_id -a /store/EQUIPES/PBI/Amira_Boukherissa/results_spada/pfam/*/new_analysis/61_final_filtred.gtf -o /store/EQUIPES/PBI/Amira_Boukherissa/results_spada/pfam/*/new_analysis/*.featurecounts.txt /store/EQUIPES/PBI/Amira_Boukherissa/Genomes/Genomes_spada/*/star_align2_Aligned.sortedByCoord.out.bam

#if using htseq_count (Supprimer les lignes à la fin)
awk -F'\t' '$2 == 0' *_counts_htseq.txt | awk '{print $1}'  > *_counts_Zero.txt
#if using featurecounts
awk -F'\t' '$7 == 0' *.featurecounts.txt | awk '{print $1}'  > *_counts_Zero.txt
 
python3 ../../pfam/Exclude_ids.py *_counts_Zero.txt ../../pfam/*/new_analysis/*_spada_high_confidence_noDup_checked.fasta ../../pfam/*/new_analysis/*_spada_high_confidence_noDup_checked_RNA.fasta
python3 ../../pfam/Exclude_ids.py *_counts_Zero.txt ../../pfam/*/new_analysis/*_spada_HC_noDup_mature_checked.fasta ../../pfam/*/new_analysis/*_spada_HC_noDup_mature_checked_RNA.fasta

#4- Motif check for IRLC (Cys 4 or 6) : 
python3 ../../NCR_final_check2.py -i *_spada_HC_noDup_mature_checked_RNA.fasta > motif_check_mature.txt
grep "crp\|cluster" motif_check_mature.txt > motif_check_mature_ids.txt
python3 ../../Exclude_ids.py motif_check_mature_ids.txt *_spada_HC_noDup_mature_checked_RNA.fasta *_spada_HC_noDup_mature_checked_RNA_motif.fasta
python3 ../../Exclude_ids.py final_check_mature_ids.txt *_spada_high_confidence_noDup_checked_RNA.fasta *_spada_high_confidence_noDup_checked_RNA_motif.fasta

#Search new NCRs 


#Search NCRs to add 



#If studied plant with known NCR peptides, continue the analysis with "AddNCR", else if the NCR peptides are not known yet continue with NCR_RNA or NCR_RNA_motif

#5-Split sequences according to clusters: 
seqkit split --by-id --id-regexp "(cluster[0-9]+)" meloff_spada_noDup_mature_checked_motif.fasta 

#Remove crp files : 
rm *crp* 

#Rename cluster files : 
for file in *.fasta ; do mv $file ${file//TriPr_spada_noDup_mature_checked_motif.id_cluster/Cluster_} ; done 

#Extract crp sequences : 
cd TriPr_spada_noDup_mature_checked_motif.fasta.split
cat *.fasta > ../All_clusters.fasta
cd .. 
grep "^>" All_clusters.fasta | sed 's/>//' | sed 's/\s.*//' > clusters_ids.txt 
python3 ../../Exclude_ids.py clusters_ids.txt TriPr_spada_noDup_mature_checked_RNA_motif.fasta TriPr_mature_crps.fasta
cat TriPr_mature_crps.fasta | sed 's/\s.*//' > TriPr_mature_crp.fasta #Delete the description after the fasta IDs 


#Hmmsearch crp sequences in our clusters 
cd PBI/Amira_Boukherissa/Genomes/MCL_IRLC/new_analysis
mkdir profiles_mature 
cd mature_alignments
for file in *.fasta ; do hmmbuild ../profiles_mature/$file.hmm $file ; done
hmmsearch --cpu 24 --tblout TriPr_CRPs_hmmsearch.out all_clusters.hmm /home/amira.boukherissa/PBI/Amira_Boukherissa/results_spada/pfam/TriPr/all_spada_analysis/TriPr_mature_crps.fasta 

#Extract best hit for each sequence of Hmmsearch
perl -ane '$key = $F[0]; $eval = $F[4]; $cluster = $F[2]; if (!exists $data{$key} || $eval < $data{$key}[1]) {$data{$key} = [$F[1],$eval,$cluster]} END {foreach $key (sort keys %data) {print "$key $data{$key}[0] $data{$key}[2]\n";}}' TriPr_CRPs_hmmsearch.out | sed 's/ -//' > only_best_hits.txt

#Extract fasta sequences of the best hits
awk 'BEGIN{while(getline<"TriPr_CRPs_best_hits.txt"){a[$1]=$2}} /^>/{print $0" "a[substr($1,2)]}!/^>/{print}' < TriPr_mature_crp.fasta | sed 's/\s/\|/' > TriPr_mature_crps_classified.fasta

#Keep only classified sequences
perl -ne 'if (/^>/) { ($id) = $_ =~ /^>(\S+)/; $found = 0; if ($id =~ /Cluster/) { $found = 1; print; } } elsif ($found) { print; }' TriPr_mature_crps_classified.fasta > TriPr_mature_crps_classified_clusters.fasta

#Split sequences according to clusters : 
cat Meloff_mature_crps_classified_clusters.fasta | sed 's/Cluster_/Cluster/' > Meloff_mature_crps_classified_cluster.fasta #delete "_"
seqkit split --by-id --id-regexp "(Cluster[0-9]+)" TriPr_mature_crps_classified_cluster.fasta
#Rename cluster files 
for file in *.fasta ; do mv $file ${file//TriPr_mature_crps_classified_clusters.id_Cluster/Cluster_} ; done


#Merge the same clusters : 
mkdir merged_clusters 

#Change the directories names in this script
../../Merge_clusters.sh #Merge clusters between clusters classified with spada and those classified with hmmsearch (copy the files if there is no intersection) 







