import os
import csv
from Bio import SeqIO

def determine_cysteine_motif(sequence):
    count_c = sequence.count('C')
    if count_c in [4, 5]:
        return "4 cysteines"
    elif count_c in [6, 7]:
        return "6 cysteines"
    elif count_c >= 8:
        return "8 cysteines"
    else:
        return "Less than 4 cysteines"

# Define the root directory containing the superclusters
root_dir = "/media/amira.boukherissa/Extreme Pro/CoEvolSym/NCR_peptides/By_supercluster"

# Initialize an empty list to store the table data
table_data = []

# Traverse through each supercluster folder
for supercluster in os.listdir(root_dir):
    supercluster_path = os.path.join(root_dir, supercluster)
    if os.path.isdir(supercluster_path):
        # Traverse through each cluster file in the supercluster
        for cluster_file in os.listdir(supercluster_path):
            if cluster_file.endswith(".fasta") or cluster_file.endswith(".fa"):
                cluster_path = os.path.join(supercluster_path, cluster_file)
                cluster_name = os.path.splitext(cluster_file)[0]
                
                # Parse the FASTA file
                for record in SeqIO.parse(cluster_path, "fasta"):
                    sequence_id = record.id
                    sequence = str(record.seq)
                    species_name = sequence_id.split('|')[0]  # Assuming species name is the first part of the sequence ID
                    cysteine_motif = determine_cysteine_motif(sequence)
                    
                    # Append the data to the table_data list
                    table_data.append([sequence_id, cluster_name, supercluster, species_name, cysteine_motif])

# Write the table data to a CSV file
output_csv = "output_table.csv"
with open(output_csv, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    # Write the header
    csvwriter.writerow(['Sequence_ID', 'Cluster_Name', 'Supercluster_Name', 'Species_Name', 'Cysteine_Motif'])
    # Write the rows
    csvwriter.writerows(table_data)

print(f"Table saved to {output_csv}")
