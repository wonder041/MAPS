import sys
import os
from Bio import SeqIO
from Bio import pairwise2

merged_path,r1_path,r2_path,Output_path=sys.argv[1:]

def Trim_merged(merged_seqr):
    r1_seqr,r2_seqr=id_r1r2_dic[merged_seqr.id]    
    index1=len(merged_seqr)-len(r1_seqr)
    correct_rate1=sum(1 for trimmed_char,r1_char in zip(str(merged_seqr.seq[index1:]),str(r1_seqr.seq)) if trimmed_char==r1_char or r1_char =="N")/len(r1_seqr)
    if correct_rate1 < 0.9:
        print("correct_rate1",correct_rate1)
        if str(r1_seqr.seq).count("N") >5:
            return
        print(pairwise2.align.globalxx(merged_seqr.seq[index1:],r1_seqr.seq)[0])
        raise
    index2=len(r2_seqr)
    correct_rate2=sum(1 for trimmed_char,r2_char in zip(str(merged_seqr.seq[:index2]),str(r2_seqr.seq.reverse_complement())) if trimmed_char==r2_char or r2_char=="N")/len(r2_seqr)
    if correct_rate2 < 0.9:
        print("correct_rate2",correct_rate2)
        if str(r2_seqr.seq).count("N") >5:
            return
        print(pairwise2.align.globalxx(merged_seqr.seq[:index2],r2_seqr.seq.reverse_complement()))
        raise
    return(merged_seqr[index1:index2])

    
r1r2_gen=zip(SeqIO.parse(r1_path,"fastq"),SeqIO.parse(r2_path,"fastq"))
id_r1r2_dic=dict((r1_seqr.id,(r1_seqr,r2_seqr)) for r1_seqr,r2_seqr in r1r2_gen)    
Output_seqr_gen=(Trim_merged(merged_seqr) for merged_seqr in SeqIO.parse(merged_path,"fastq") if str(merged_seqr.seq).count("N") <1)
SeqIO.write((seqr for seqr in Output_seqr_gen if seqr),Output_path,"fastq")
