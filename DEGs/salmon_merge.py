import os
from numpy import mean
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Merge salmon results for replicates')
parser.add_argument('-i',required=True,help='folder of salmon outputs')
parser.add_argument('-o',required=True,help='output folder')
parser.add_argument('-p',required=True,help='plant name (alfalfa or clover)')
args = vars(parser.parse_args())

def merge_salmon(in_folder,out,plant):
    replicate_folders = ["Lupinus_lut1","Lupinus_lut2","Lupinus_lut3","Lupinus_mar1","Lupinus_mar2","Lupinus_mar3"]
    out_file = open(out + "/" + plant + "_triplicate_salmon_tpm.txt","w")
    if plant == "alfalfa":
        out_file.write("Geneid\t{}\t{}\t{}\t{}\t{}\t{}\tAverage_count\tAverage_tpm\n".format(replicate_folders[0] + "_count",replicate_folders[1] + "_count",replicate_folders[2] + "_count",replicate_folders[0] + "_tpm",replicate_folders[1] + "_tpm",replicate_folders[2] + "_tpm"))
    elif plant == "clover":
        out_file.write("Geneid\t{}\t{}\t{}\t{}\t{}\t{}\tAverage_count\tAverage_tpm\n".format(replicate_folders[3] + "_count",replicate_folders[4] + "_count",replicate_folders[5] + "_count",replicate_folders[3] + "_tpm",replicate_folders[4] + "_tpm",replicate_folders[5] + "_tpm"))

    total_id_list = []
    dict_1 = {}
    dict_2 = {}
    dict_3 = {}
    for f in replicate_folders:
        if plant in f:
            with open(salmon_folder + "/" + f + "/" + "quant.sf") as file:
                lines = file.readlines()
                lines = lines[1::]
                for l in lines:
                    info = l.strip("\n").split("\t")
                    total_id_list.append(info[0])
                    if f == replicate_folders[0]:
                        dict_1.update({info[0]:info})
                    if f == replicate_folders[1]:
                        dict_2.update({info[0]:info})
                    if f == replicate_folders[2]:
                        dict_3.update({info[0]:info})
                    if f == replicate_folders[3]:
                        dict_1.update({info[0]:info})
                    if f == replicate_folders[4]:
                        dict_2.update({info[0]:info})
                    if f == replicate_folders[5]:
                        dict_3.update({info[0]:info})
    total_id_list = list(set(total_id_list))
    for id in total_id_list:
        try:
            tpm_1 = [dict_1[id][4],dict_1[id][3]]
        except:
            tpm_1 = [0,0]
        try:
            tpm_2 = [dict_2[id][4],dict_2[id][3]]
        except:
            tpm_2 = [0,0]
        try:
            tpm_3 = [dict_3[id][4],dict_3[id][3]]
        except:
            tpm_3 = [0,0]
        count_list = [float(tpm_1[0]),float(tpm_2[0]),float(tpm_3[0])]
        tpm_list = [float(tpm_1[1]),float(tpm_2[1]),float(tpm_3[1])]
        count_average = mean(count_list)
        tpm_average = mean(tpm_list)
        out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(id,count_list[0],count_list[1],count_list[2],tpm_list[0],tpm_list[1],tpm_list[2],count_average,tpm_average))
    
    length_list = []
    for id in total_id_list:
        try:
            length = dict_1[id][1]
        except:
            try:
                length = dict_2[id][1]
            except:
                try:
                    length = dict_3[id][1]
                except:
                    length = 0
        length_list.append(length)
    df = pd.DataFrame(total_id_list,columns = [plant + "Geneid"])
    df["Length"] = length_list
    df.to_csv(folder + "/" + plant + "_length.txt", index= False)


if __name__ == '__main__':
    salmon_folder = args["i"]
    folder = args["o"]
    plant = args["p"]
    
    merge_salmon(salmon_folder,folder,plant)



    
