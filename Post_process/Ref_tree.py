from ete3 import Tree,NodeStyle,TreeFace,faces,TreeStyle,TreeStyle,AttrFace,ClusterTree,ProfileFace,RectFace,CircleFace,TextFace
import pandas as pd
import math
import itertools
import sys

t=Tree("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/Reftree/ref.tre")

csv=pd.read_csv("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/Reftree/MEGAPRIMER_v1_PolB_primer_pairs_annotations.csv")
PP_tar_dic=dict(zip(csv.iloc[:,0].tolist(), csv.iloc[:,-2].tolist()))
PP_out_dic=dict(zip(csv.iloc[:,0].tolist(), csv.iloc[:,-1].tolist()))
PP_out_dic={PP:out_lst for PP,out_lst in PP_out_dic.items() if type(out_lst)==type("")}



PP_tara_dic={}
with open("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/ref_PP_cov/tara_10.primersearch","r") as primersearch_file:
    for line in primersearch_file:
        if line.startswith("Primer name"):
            PP=line.strip().split()[-1]
            continue
        if line.startswith("\tSequence"):
            tara=line.strip().split()[-1]
            PP_tara_dic[PP]=PP_tara_dic.get(PP,[])+[tara]

tara_PP_dic={}
for PP,tara_lst in PP_tara_dic.items():
    for tara in tara_lst:
        tara_PP_dic[tara]=tara_PP_dic.get(tara,[])+[PP]


















# tar_PP_dic={tar:PP for PP,tar_lst in PP_tar_dic.items() for tar in tar_lst.split(",")}
# out_PP_dic={out:PP for PP,out_lst in PP_out_dic.items() for out in out_lst.split(",")}

# seq_PP_dic={**tar_PP_dic,**out_PP_dic}

# mimi_PP_dic={}
# for PP,mimi_lst in PP_out_dic.items():
#     for mimi in mimi_lst.split(","):
#         try:
#             mimi_PP_dic[mimi].append(PP)
#         except KeyError:
#             mimi_PP_dic[mimi]=[PP]






line_width=1

ns_com = NodeStyle(size=0,hz_line_width=line_width,vt_line_width=line_width)

ts_sub_c=TreeStyle()
ts_sub_c.mode="c"

ts_sub_c.scale=100
ts_sub_c.draw_guiding_lines= False
ts_sub_c.show_leaf_name = False
ts_sub_c.show_scale = False
ts_sub_c.allow_face_overlap = True
ts_sub_c.legend_position=4
ts_sub_c.extra_branch_line_color="White"
# ts_sub_c.complete_branch_lines_when_necessary=False



for n in t.traverse():
   n.set_style(ns_com)


for leaf in t:
    PP_lst=tara_PP_dic.get(leaf.name,[])
    length_lst=[int(PP[2:]) for PP in PP_lst]
    for pos in range(83):
        # color="#FF"+str(hex(255-pos*2-80)).split("x")[-1].upper().zfill(2)+"00" if pos==length else "White"
        color="#FF"+str(hex(255-pos*2-80)).split("x")[-1].upper().zfill(2)+"00"
        l=20 if (pos in length_lst) else 1
        # l=20
        leaf.add_face(RectFace(l,13,color,color), pos, position="aligned")
    if "MIMI" in leaf.name:
        leaf.add_face(TextFace(leaf.name,fsize=20,fgcolor="Black"),83,position="aligned")





# t.show(tree_style=ts_sub_c)
t.render("/user1/scl1/yanzeli/Exchange/Ref_tree_180402.png",w=6400,units="px",tree_style=ts_sub_c)




sys.exit(0)

# PP_lst=tara_PP_dic.get(leaf.name.replace("MEGA_",""),["PP62"])
