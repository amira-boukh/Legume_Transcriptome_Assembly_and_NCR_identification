import os


with open("./61_final.tbl","r") as f:
    w = open("./seq.fasta","w")
    for l in f:
        info = l.strip("\n").split("\t")
        id = info[0]
        seq = info[10]
    
        w.write(">" + id + "\n")
        w.write(seq + "\n")
    w.close()