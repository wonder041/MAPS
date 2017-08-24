from Bio import SeqIO
import decode_clstr

# import Make_OTU_table
# import sys

# input_path="/lustre1/aptmp/yanzeli/Megaviridae/Outputs_170310/9_CDHIT/A0_90.faa"
# output_path="/lustre1/aptmp/yanzeli/Megaviridae/Outputs_170310/9_CDHIT/A0_90_1000.faa"

# rep_psdic_dic = Make_OTU_table.Get_dic(input_path+".clstr")
# rep_size_tup_arr = sorted([(rep,sum(psdic.values())) for rep,psdic in rep_psdic_dic.items()],key=lambda tup:tup[1],reverse=True)
# rep_picked_arr=[rep for (rep,size) in rep_size_tup_arr if size > rep_size_tup_arr[1000][1]]

# output_seqr_gen=(seqr for seqr in SeqIO.parse(input_path,"fasta") if seqr.id in rep_picked_arr)
# SeqIO.write(output_seqr_gen,output_path,"fasta")

input_path="/aptmp/yanzeli/Paper_pipeline/Outputs_170812/8_POSTPROCESS/90.faa"
OTU_barcode_table=decode_clstr.Calc_table(input_path+".clstr")
loc_arr=OTU_barcode_table.sum(1).sort_values().iloc[::-1].index.tolist()
tar_set=set(loc_arr[:1000])


def Rename(nseqr,seqr):
    seqr.id=seqr.name=seqr.description=seqr.id+"U"+str(OTU_barcode_table.iloc[nseqr].sum())
    return seqr

subset_gen=(Rename(nseqr,seqr) for nseqr,seqr in enumerate(SeqIO.parse(input_path,"fasta")) if nseqr in tar_set)
SeqIO.write(subset_gen,"".join(input_path.split(".")[:-1])+".sub."+input_path.split(".")[-1],"fasta")
    