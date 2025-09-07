from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='final check of length and number of Cys residues of NCR peptide sequences')
parser.add_argument('-i',required=True,help='NCR peptide sequences')
args = vars(parser.parse_args())



def check(infile):
    for record in SeqIO.parse(infile,"fasta"):
        id = record.id
        sequence = str(record.seq)
        length = len(sequence)
        n_cys = sequence.count("C")
        if length > 100:
            print("NCR sequence too long please check")
            print(id + "\n")

        if n_cys < 4 :
            print("invalid number of Cys residues please check")
            print(id + "\n")


if __name__ == '__main__':
    alfalfa_NCR_file = args['i']
    check(alfalfa_NCR_file)
    
