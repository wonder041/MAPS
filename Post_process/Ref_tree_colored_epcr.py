from ete3 import Tree,NodeStyle,TreeFace,faces,TreeStyle,TreeStyle,AttrFace,ClusterTree,ProfileFace,RectFace,CircleFace,TextFace
import pandas as pd
import math
import itertools
import matplotlib.pyplot as plt

colors_f=plt.get_cmap("tab10").colors

convert_color=lambda color:"#{:X}{:X}{:X}".format(*map(lambda x:int(256*x),color))

colors=list(map(convert_color,colors_f))

t=Tree("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/Reftree/MEGA21_TARA906_full.tre")


PP_tara_dic={}
with open("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/ref_PP_cov/tara_10.primersearch","r") as primersearch_file:
    for line in primersearch_file:
        if line.startswith("Primer name"):
            PP=line.strip().split()[-1]
            continue
        if line.startswith("\tSequence"):
            tara=line.strip().split()[1]
            PP_tara_dic[PP]=PP_tara_dic.get(PP,[])+[tara]

tara_PP_dic={}
for PP,tara_lst in PP_tara_dic.items():
    for tara in tara_lst:
        tara_PP_dic[tara]=tara_PP_dic.get(tara,[])+[PP]


PP_mega_dic={}
with open("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/ref_PP_cov/MEGA21_10.primersearch","r") as primersearch_file:
    for line in primersearch_file:
        if line.startswith("Primer name"):
            PP=line.strip().split()[-1]
            continue
        if line.startswith("\tSequence"):
            mega=line.strip().split()[1]
            PP_mega_dic[PP]=PP_mega_dic.get(PP,[])+[mega]

mega_PP_dic={}
for PP,mega_lst in PP_mega_dic.items():
    for mega in mega_lst:
        mega_PP_dic[mega]=mega_PP_dic.get(mega,[])+[PP]




seq_PP_dic={**tara_PP_dic,**mega_PP_dic}













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




def get_depth(name):
    if name.isnumeric():
        return(int(name))
    else:
        return(0)

def cal_depth(node):
    if node.is_leaf():
        return
    [cal_depth(child) for child in node.get_children()]
    node.name = node.name+str(max(get_depth(child.name)
                                  for child in node.get_children())+1)


cal_depth(t)
layer_sum = int(t.name)+1

# total length = 3
intervention = 3/layer_sum

# appoint length


def set_dist(node):
    if node.is_leaf():
        node.dist = intervention*int(node.up.name)
        return

    [set_dist(child) for child in node.get_children()]

    if int(node.name) == layer_sum-1:
        node.dist = intervention
        return
    node.dist = intervention*(int(node.up.name)-int(node.name))


set_dist(t)






# n_gen_circle_name=0
# def gen_circle_name():
#     global n_gen_circle_name
#     n_gen_circle_name+=1
#     return chr(n_gen_circle_name+9311)



for n in t.traverse():
   n.set_style(ns_com)


for leaf in t:
    if "MIMI" in leaf.name:
        PP_lst=seq_PP_dic.get(leaf.name ,[]) 
    else:    
        PP_lst=seq_PP_dic.get("MEGA_"+leaf.name ,[]) 
    if not PP_lst:
        print(leaf.name)
    length_lst=[int(PP[2:]) for PP in PP_lst]
    for pos in range(84):
        # color="#FF"+str(hex(255-pos*2-80)).split("x")[-1].upper().zfill(2)+"00" if pos==length else "White"
        # color="#FF"+str(hex(255-pos*2-80)).split("x")[-1].upper().zfill(2)+"00"
        color=colors[pos%len(colors)]
        l=20 if (pos in length_lst) else 1
        # l=20
        leaf.add_face(RectFace(l,int(0.2*pos+2),color,color), pos, position="aligned")
    # if "MIMI" in leaf.name:
        # char=gen_circle_name()
        # leaf.add_face(TextFace(char,fsize=100,fgcolor="Black"),83,position="aligned")
        # leaf.add_face(TextFace(leaf.name.replace("MIMI_",""),fsize=20,fgcolor="Black"),83,position="aligned")





# t.show(tree_style=ts_sub_c)
t.render("/user1/scl1/yanzeli/Exchange/Ref_tree_180612.png",w=6400,units="px",tree_style=ts_sub_c)