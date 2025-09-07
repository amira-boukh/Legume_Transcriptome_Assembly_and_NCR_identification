mkdir -p results/FINAL/merged_clusters/IRLC
mkdir -p results/FINAL/merged_clusters/DAL

# Merge IRLC base clusters with newly classified clusters
for file1 in results/IRLC/clusters_prot/*; do
    base1=$(basename "$file1")
    if [ -f "results/FINAL/merged_clusters/IRLC/$base1" ]; then
        cat "$file1" >> "results/FINAL/merged_clusters/IRLC/$base1"
    else
        cp "$file1" "results/FINAL/merged_clusters/IRLC/$base1"
    fi
done
for file2 in results/SPADA/*_CRPs_classified_by_cluster/*.fasta; do
    base2=$(basename "$file2")
    if [ -f "results/FINAL/merged_clusters/IRLC/$base2" ]; then
        cat "$file2" >> "results/FINAL/merged_clusters/IRLC/$base2"
    else
        cp "$file2" "results/FINAL/merged_clusters/IRLC/$base2"
    fi
done

# Merge DAL base clusters with newly classified clusters
for file1 in results/DAL/clusters_prot/*; do
    base1=$(basename "$file1")
    if [ -f "results/FINAL/merged_clusters/DAL/$base1" ]; then
        cat "$file1" >> "results/FINAL/merged_clusters/DAL/$base1"
    else
        cp "$file1" "results/FINAL/merged_clusters/DAL/$base1"
    fi
done
for file2 in results/SPADA/*_CRPs_classified_by_cluster/*.fasta; do
    base2=$(basename "$file2")
    if [ -f "results/FINAL/merged_clusters/DAL/$base2" ]; then
        cat "$file2" >> "results/FINAL/merged_clusters/DAL/$base2"
    else
        cp "$file2" "results/FINAL/merged_clusters/DAL/$base2"
    fi
done

