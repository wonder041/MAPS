from Bio import SeqIO
import decode_clstr

input_path="/aptmp/yanzeli/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/S090.faa"
OTU_barcode_table=decode_clstr.Calc_table(input_path+".clstr")
loc_arr=OTU_barcode_table.sum(1).sort_values().iloc[::-1].index.tolist()
tar_set=set(loc_arr[:1000])
