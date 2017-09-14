import os
import pandas as pd
import multiprocessing

resouce_dir="/aptmp/yanzeli/Megaviridae/Source_FL/"
# outputs_dir="/lustre1/aptmp/ideas2/yanzeli/Megaviridae/All_170905_2/"
outputs_dir="/lustre1/aptmp/ideas2/yanzeli/Megaviridae/All_170907_0/"
barcode_arr=sorted(set(file_name.split("_")[0] for file_name in os.listdir(resouce_dir)))
    
def Count_fastq(input_path):
    try:
        input_file=open(input_path,"r")
        total=sum(1 for line in input_file if line=="+\n")
        return total
    except:
        return 0
    
def Count_fasta(input_path):
    try:
        input_file=open(input_path,"r")
        total=sum(1 for line in input_file if line.startswith(">"))
        return total
    except:
        return 0
    
def Count_faa(input_path):
    try:
        input_file=open(input_path,"r")
        total=sum(int(line.split("S")[2].split("F")[0]) for line in input_file if line.startswith(">"))
        return total
    except:
        return 0

    
# barcode_arr_dic={}
# for barcode in barcode_arr:
def Stat_barcode(barcode):
    raw=Count_fastq(resouce_dir+barcode+"_R1.fastq")
    hq=Count_fastq(outputs_dir+"1_A5G40/"+barcode+"_R1.fastq")
    wp=Count_fastq(outputs_dir+"2_CUTADAPT_G40/"+barcode+"_R1.fastq")
    merged=Count_fasta(outputs_dir+"3_MERGE/"+barcode+".fna")
    blast=Count_faa(outputs_dir+"6_BLASTP/"+barcode+".faa")
    pplacer=Count_faa(outputs_dir+"7_ALIGNMENT/"+barcode+"_Pplacer.faa")
    return (barcode,[raw,hq,wp,merged,blast,pplacer])
    # barcode_arr_dic[barcode]=[raw,hq,wp,merged,blast,pplacer]

pool = multiprocessing.Pool(20)
str_gen = pool.imap_unordered(Stat_barcode, barcode_arr)    
barcode_arr_dic=dict(str_gen)    
    
    
df=pd.DataFrame(barcode_arr_dic).T
df.to_csv("All_170907_0.stat")

    

    
    
    

