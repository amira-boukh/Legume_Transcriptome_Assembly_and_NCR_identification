import argparse

parser = argparse.ArgumentParser(description='Check blast result of NCR peptides and output high similar sequences ids')
parser.add_argument('-i',required=True,help='input blast table')

args = vars(parser.parse_args())


def check(file):
    with open(file,"r") as f:
        for l in f:
            info = l.strip("\n").split("\t")
            id1 = info[0]
            id2 = info[1]
            if id1 != id2:
                print(id1 + " " + id2)

#check("/home/rui-huang/Documents/RNA_seq_doc/pfam/alfalfa_NCR_blast.txt")
check(args['i'])