import os
from Bio import SeqIO
import multiprocessing

input_dir="/user1/scl1/yanzeli/aptmp/Paper_pipeline/Sources/"

def Calc_avglen(input_path):
    input_seqr_arr=list(SeqIO.parse(input_path,"fastq"))
    total_len=sum(len(seqr) for seqr in input_seqr_arr)
    avglen=total_len/len(input_seqr_arr)
    return avglen

def Calc_PP_avglen(PP):
    PP_avglen = 0
    for suffix in suffix_arr:
        input_path = "{}{}_{}.fastq".format(input_dir, PP, suffix)
        PP_avglen += Calc_avglen(input_path) / 2
    return PP,PP_avglen
PP_arr=sorted(set(name.split("_")[0] for name in os.listdir(input_dir)))
suffix_arr=sorted(set(name.split(".")[0].split("_")[1] for name in os.listdir(input_dir)))
pool=multiprocessing.Pool(20)
res_gen=pool.imap_unordered(Calc_PP_avglen,PP_arr)
for PP,PP_avglen in res_gen:
    print(PP,int(PP_avglen))