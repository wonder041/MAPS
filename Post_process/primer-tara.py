import pandas as pd
csv=pd.read_csv("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/Reftree/MEGAPRIMER_v1_PolB_primer_pairs_annotations.csv")

with open("PP_pos","r") as PP_pos_file:
    PP_pos_lst=list(set(line.strip().split("c")[0] for line in PP_pos_file))



PP_tar_dic=dict(zip(csv.iloc[:,0].tolist(), csv.iloc[:,-2].tolist()))
PP_out_dic=dict(zip(csv.iloc[:,0].tolist(), csv.iloc[:,-1].tolist()))
PP_out_dic={PP:out_lst for PP,out_lst in PP_out_dic.items() if type(out_lst)==type("")}




PP_neg_lst=sorted(PP for PP in PP_tar_dic if PP not in PP_pos_lst)


seq_lst=[]

for PP in PP_neg_lst:
    seq_lst.extend(PP_tar_dic[PP].split(","))
    seq_lst.extend(PP_out_dic.get(PP,"").split(","))

tara_lst=sorted(set(seq for seq in seq_lst if seq and "MIMI" not in seq))


collections.Counter(tara.split("_")[1] for tara in tara_lst)
#Counter({'DCM': 73, 'MES': 5, 'SUR': 214})

collections.Counter(map(int,(tara.split("_")[0] for tara in tara_lst)))
# Counter({4: 6,
#          7: 1,
#          9: 6,
#          23: 1,
#          25: 6,
#          30: 14,
#          32: 1,
#          34: 2,
#          36: 6,
#          38: 4,
#          41: 2,
#          42: 1,
#          45: 1,
#          48: 3,
#          58: 5,
#          64: 2,
#          65: 8,
#          66: 7,
#          67: 38,
#          68: 6,
#          70: 6,
#          72: 2,
#          76: 11,
#          78: 8,
#          82: 4,
#          84: 1,
#          85: 5,
#          93: 6,
#          94: 3,
#          96: 2,
#          98: 1,
#          99: 3,
#          100: 9,
#          102: 9,
#          109: 7,
#          111: 3,
#          112: 4,
#          122: 10,
#          123: 5,
#          124: 10,
#          125: 13,
#          128: 2,
#          132: 4,
#          133: 16,
#          137: 6,
#          138: 2,
#          140: 2,
#          142: 2,
#          145: 1,
#          146: 3,
#          148: 1,
#          149: 1,
#          150: 3,
#          152: 7})
