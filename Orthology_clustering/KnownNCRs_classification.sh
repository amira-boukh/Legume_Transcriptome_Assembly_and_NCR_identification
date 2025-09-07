#!/usr/bin/env bash

#DATA_DIR/<SpeciesName>/proteome.faa
#DATA_DIR/<SpeciesName>/genome.fa
#DATA_DIR/<SpeciesName>/known_ncrs.faa
#Where <SpeciesName> is one of:
#Medicago_truncatula
#Medicago_sativa
#Pisum_sativum
#Cicer_arietinum
#Aeschynomene_evenia
#Arachis_hypogaea
#
#Usage: ./KnownNCRs_classification.sh DATA_DIR OUT_DIR [threads]

DATA_DIR="${1:-data}"
OUT_ROOT="${2:-results}"
THREADS="${3:-8}"
EVALUE="1e-5"
MCL_I="1.5"

IRLC=("Medicago_truncatula" "Medicago_sativa" "Pisum_sativum" "Cicer_arietinum")
DAL=("Aeschynomene_evenia" "Arachis_hypogaea")

for SET in IRLC DAL COMBINED; do
  if [ "$SET" = "IRLC" ]; then SPECIES=("${IRLC[@]}"); fi
  if [ "$SET" = "DAL" ]; then SPECIES=("${DAL[@]}"); fi
  if [ "$SET" = "COMBINED" ]; then SPECIES=("${IRLC[@]}" "${DAL[@]}"); fi

  OUTDIR="$OUT_ROOT/$SET"
  mkdir -p "$OUTDIR"
  > "$OUTDIR/all_augmented_proteomes.list"
  > "$OUTDIR/all_known_NCR_ids.txt"

  for species in "${SPECIES[@]}"; do
    sdir="$OUTDIR/$species"
    mkdir -p "$sdir"
    proteome_faa="$DATA_DIR/$species/proteome.faa"
    genome_fa="$DATA_DIR/$species/genome.fa"
    known_ncr_faa="$DATA_DIR/$species/known_ncrs.faa"

    #Normalize headers to Species|ID
    awk -v S="$species" 'BEGIN{OFS=""} /^>/ {sub(/^>/,""); print ">",S,"|",$0; next} {print}' "$proteome_faa" > "$sdir/${species}.proteome.tag.faa"
    awk -v S="$species" 'BEGIN{OFS=""} /^>/ {sub(/^>/,""); print ">",S,"|",$0; next} {print}' "$known_ncr_faa" > "$sdir/${species}.knownNCRs.tag.faa"

    #Make BLAST DBs
    makeblastdb -in "$sdir/${species}.proteome.tag.faa" -dbtype prot -parse_seqids 
    makeblastdb -in "$genome_fa" -dbtype nucl -parse_seqids 

    #BLASTP known NCRs vs proteome
    blastp -query "$sdir/${species}.knownNCRs.tag.faa" -db "$sdir/${species}.proteome.tag.faa" -out "$sdir/${species}.ncr_vs_proteome.blastp.tsv" -outfmt 6 -evalue "$EVALUE" -num_threads "$THREADS"

    #Keep top hit per query and collect proteome NCR IDs
    sort -k1,1 -k12,12nr "$sdir/${species}.ncr_vs_proteome.blastp.tsv" | awk '!(q[$1]++)' > "$sdir/${species}.ncr_vs_proteome.top.tsv"
    cut -f2 "$sdir/${species}.ncr_vs_proteome.top.tsv" | sort -u > "$sdir/${species}.proteome_ncr_ids.txt"

    #Queries without proteome hit
    grep '^>' "$sdir/${species}.knownNCRs.tag.faa" | sed 's/^>//; s/[\t ].*$//' | sort -u > "$sdir/${species}.ncr_queries.txt"
    cut -f1 "$sdir/${species}.ncr_vs_proteome.top.tsv" | sort -u > "$sdir/${species}.ncr_with_hit.txt"
    comm -23 "$sdir/${species}.ncr_queries.txt" "$sdir/${species}.ncr_with_hit.txt" > "$sdir/${species}.ncr_without_hit.txt"

    #Start proteome with known NCR hits
    cp "$sdir/${species}.proteome.tag.faa" "$sdir/${species}.proteome.aug.faa"

    #If some NCRs missing in proteome, search in genome via TBLASTN and append supported peptides
    if [ -s "$sdir/${species}.ncr_without_hit.txt" ]; then
      awk 'NR==FNR{a[$1]=1; next} /^>/{h=$0; sub(/^>/,"",h); split(h,b,/\s+/); keep=a[b[1]]} {if(/^>/){if(keep) print ">"$0; else next} else {if(keep) print}}' "$sdir/${species}.ncr_without_hit.txt" "$sdir/${species}.knownNCRs.tag.faa" > "$sdir/${species}.missing_queries.faa"
      tblastn -query "$sdir/${species}.missing_queries.faa" -db "$genome_fa" -out "$sdir/${species}.missing_vs_genome.tblastn.tsv" -outfmt 6 -evalue "$EVALUE" -num_threads "$THREADS"
      if [ -s "$sdir/${species}.missing_vs_genome.tblastn.tsv" ]; then
        cut -f1 "$sdir/${species}.missing_vs_genome.tblastn.tsv" | sort -u > "$sdir/${species}.ncr_tblastn_supported.txt"
        awk 'NR==FNR{a[$1]=1; next} /^>/{h=$0; sub(/^>/,"",h); split(h,b,/\s+/); keep=a[b[1]]} {if(/^>/){if(keep) print ">"$0; else next} else {if(keep) print}}' "$sdir/${species}.ncr_tblastn_supported.txt" "$sdir/${species}.knownNCRs.tag.faa" > "$sdir/${species}.ncr_to_append.faa"
        cat "$sdir/${species}.ncr_to_append.faa" >> "$sdir/${species}.proteome.aug.faa"
        grep '^>' "$sdir/${species}.ncr_to_append.faa" | sed 's/^>//; s/[\t ].*$//' >> "$sdir/${species}.proteome_ncr_ids.txt"
        sort -u "$sdir/${species}.proteome_ncr_ids.txt" -o "$sdir/${species}.proteome_ncr_ids.txt"
      fi
    fi

    echo "$sdir/${species}.proteome.aug.faa" >> "$OUTDIR/all_augmented_proteomes.list"
    cat "$sdir/${species}.proteome_ncr_ids.txt" >> "$OUTDIR/all_known_NCR_ids.txt"
  done

  sort -u "$OUTDIR/all_known_NCR_ids.txt" -o "$OUTDIR/all_known_NCR_ids.txt"

  #Concatenate all augmented proteomes for this set
  ALL_FASTA="$OUTDIR/all_species.proteins.aug.faa"
  > "$ALL_FASTA"
  while read -r f; do cat "$f" >> "$ALL_FASTA"; done < "$OUTDIR/all_augmented_proteomes.list"

  #All-vs-all BLASTP
  makeblastdb -in "$ALL_FASTA" -dbtype prot -parse_seqids >/dev/null
  blastp -query "$ALL_FASTA" -db "$ALL_FASTA" -out "$OUTDIR/all_vs_all.blastp.tsv" -outfmt 6 -evalue "$EVALUE" -num_threads "$THREADS"
  sort -k1,1 -k12,12nr "$OUTDIR/all_vs_all.blastp.tsv" > "$OUTDIR/all_vs_all.sorted.tsv"

  #orthAgogue and MCL
  mkdir -p "$OUTDIR/orthagogue"
  orthAgogue --input "$OUTDIR/all_vs_all.sorted.tsv" -s '|' -c 12 -u -O "$OUTDIR/orthagogue"
  ABC_FILE=$(ls -1 "$OUTDIR/orthagogue"/*.abc 2>/dev/null | head -n1)
  mcl "$ABC_FILE" --abc -I "$MCL_I" -o "$OUTDIR/mcl_clusters.txt"

  #Classify clusters
  awk -v OFS="\t" 'ARGIND==1{n[$1]=1; next} {cid++; m=split($0,a,/\s+/); nc=0; for(i=1;i<=m;i++){if(a[i] in n) nc++} no=m-nc; if(nc>=2 && no==0){print cid, m, nc, "NCR", $0 > ncr} else if(nc>=1 && no>=1){print cid, m, nc, "NCR-mixed", $0 > mix}}' "$OUTDIR/all_known_NCR_ids.txt" "$OUTDIR/mcl_clusters.txt" ncr="$OUTDIR/clusters_ncr.txt" mix="$OUTDIR/clusters_ncr_mixed.txt"

  #All clusters containing at least 2 NCRs (pure or mixed)
  cat "$OUTDIR/clusters_ncr.txt" "$OUTDIR/clusters_ncr_mixed.txt" 2>/dev/null | awk '$3>=2' > "$OUTDIR/clusters_ge2NCR.txt"

  #Extract NCR IDs from those clusters and the corresponding FASTA
  awk 'ARGIND==1{n[$1]=1; next} {for(i=5;i<=NF;i++){if($i in n) print $i}}' "$OUTDIR/all_known_NCR_ids.txt" "$OUTDIR/clusters_ge2NCR.txt" | sort -u > "$OUTDIR/ge2NCR_ncr.ids"
  awk 'NR==FNR{ids[$1]=1; next} /^>/{h=$0; sub(/^>/,"",h); split(h,a,/\s+/); keep=ids[a[1]]} {if(/^>/){if(keep) print ">"$0}else{if(keep) print}}' "$OUTDIR/ge2NCR_ncr.ids" "$ALL_FASTA" > "$OUTDIR/ge2NCR_ncr.faa"

  #Per-cluster FASTAs: build ID->Cluster map and split the NCR peptides into one file per cluster
  awk 'ARGIND==1{n[$1]=1; next} {cid=$1; for(i=5;i<=NF;i++){id=$i; if(id in n) print id"\tCluster_"cid}}' "$OUTDIR/all_known_NCR_ids.txt" "$OUTDIR/clusters_ge2NCR.txt" > "$OUTDIR/id2cluster.tsv"
  awk 'BEGIN{while((getline<ARGV[1])>0){m[$1]=$2} ARGV[1]=""} /^>/{id=$1; sub(/^>/,"",id); split(id,a,/\s+/); print $0"|"m[a[1]]; next} {print}' "$OUTDIR/id2cluster.tsv" "$OUTDIR/ge2NCR_ncr.faa" > "$OUTDIR/ge2NCR_ncr.cluster_annot.faa"
  mkdir -p "$OUTDIR/clusters_prot"
  seqkit split --by-id --id-regexp "(Cluster_[0-9]+)" "$OUTDIR/ge2NCR_ncr.cluster_annot.faa" -O "$OUTDIR/clusters_prot"
  for f in "$OUTDIR"/clusters_prot/*.fasta; do mv "$f" "${f/*id_Cluster_/Cluster_}"; done

done
