import sys
from Bio import SeqIO

#revieve argvs from pipeline
primer_path,mixture_path,barcode=sys.argv[1:]

#read primers
PP_seq_dic=dict((seqr.id,str(seqr.seq)) for seqr in SeqIO.parse(primer_path,"fasta"))

#read mixtures
with open(mixture_path,"r") as mixture_file:
    PP_MP_mix_arr_dic=dict((line.strip().split("\t")[0],line.strip().split("\t")[1:]) for line in mixture_file if not line.startswith("#"))

#convert PP->MP to MP->PP
MP_mix_PP_dic_dic={}
for nMP,MP in enumerate(("MP5","MP10","MP20")):
    mix_set=set(MP_mix_arr[nMP] for MP_mix_arr in PP_MP_mix_arr_dic.values())
    MP_mix_PP_dic_dic[MP]=dict((mix,sorted([PP for PP,MP_mix_arr in PP_MP_mix_arr_dic.items() if MP_mix_arr[nMP]==mix]))for mix in mix_set)

#build PP list from reading barcode
if "PP" in barcode:
    PP_arr=["PP"+str(int(barcode.split("PP")[1].split("c")[0]))]
elif "MP" in barcode:
    MP=barcode.split("-")[-2]
    mix=barcode.split("-")[-1]
    PP_arr=MP_mix_PP_dic_dic[MP][mix]
else:
    PP_arr=sorted(list(set(key.split("_")[0] for key in PP_seq_dic.keys())))

#add mark by orientation
orientation_mark_dic={"_L":"g","_R":"G"}
primer_cmd="".join(" -"+orientation_mark_dic[ori]+" '^"+PP_seq_dic[PP+ori]+"'" for ori in sorted(orientation_mark_dic.keys()) for PP in PP_arr)

print(primer_cmd)
