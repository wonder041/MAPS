#!/aptmp/yanzeli/miniconda3/envs/full/bin/python


import pandas as pd
import numpy as np
import os

def Calc_table(clstr_path):
    abstract_barcode=lambda line:line.split(">")[1].split("N")[0]
    abstract_size=lambda line:int(line.split("S")[2].split("F")[0].split("L")[0])
    abstract_name=lambda line:line.split(">")[1].split(".")[0]

    series_arr=[]
    barcode_size_series=pd.Series(dtype="int")


    with open(clstr_path,"r") as clstr_file:
        next(clstr_file)
        for line in clstr_file:
            if line.startswith('>Cluster'):
                series_arr+=[barcode_size_series]
                barcode_size_series=pd.Series(dtype="int")
                continue
            barcode=abstract_barcode(line)
            size=abstract_size(line)
            if line.strip().endswith("*"):
                barcode_size_series.name=abstract_name(line)
            try:
                barcode_size_series[barcode]+=size
            except KeyError:
                barcode_size_series[barcode]=size
        else:
            series_arr+=[barcode_size_series]
     
    OTU_barcode_table=pd.DataFrame(series_arr).fillna(0).astype("int")
    try:
        OTU_barcode_table.to_csv(clstr_path+".csv")
    except:
        pass
    return OTU_barcode_table

clstr_path="/aptmp/yanzeli/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/S090.faa.clstr"
OTU_barcode_table=Calc_table(clstr_path)



from ete3 import Tree,NodeStyle,TreeFace,faces,TreeStyle,TreeStyle,AttrFace,ClusterTree,ProfileFace,RectFace,CircleFace,TextFace
import math


for nPP in range(15,60):


    tree = Tree("/aptmp/yanzeli/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/S090.sub.combo.aln.rearranged_rename.tre")
    output_fig_path="/user1/scl1/yanzeli/Exchange/nPP/PP_tree_"+str(nPP)+"_180319.png"

    #########################################
    line_width=6
    stick_width=15
    scale=330
    # scale=600
    log_scale=480
    color_scale=6 
    color_MIMI="Black"  

    ns_MIMI = NodeStyle(size=0,hz_line_color="Red",vt_line_color="Red",hz_line_width=line_width,vt_line_width=line_width)
    ns_com = NodeStyle(size=0,hz_line_width=line_width,vt_line_width=line_width)

    ts_sub_c=TreeStyle()
    ts_sub_c.mode="c"
    ts_sub_c.arc_start=0
    ts_sub_c.arc_span=357
    ts_sub_c.scale=scale
    ts_sub_c.draw_guiding_lines= True
    # ts_sub_c.extra_branch_line_type=0
    # ts_sub_c.extra_branch_line_color="Black"
    ts_sub_c.guiding_lines_type=2
    # ts_sub_c.force_topology=True
    # ts_sub_c.guiding_lines_color=color_query
    ts_sub_c.guiding_lines_color="Grey"
    ts_sub_c.show_leaf_name = False
    ts_sub_c.show_scale = False
    ts_sub_c.allow_face_overlap = True
    ts_sub_c.legend_position=4
    #########################################

    mimi_min=6
    mimi_max=18
    mimi_lab=14
    olpv_min=3
    olpv_max=olpv_min+1

    #########################################


    query_arr=[leaf for leaf in tree if "PP" in leaf.name]
    MIMI_arr=[leaf for leaf in tree if "MIMI" in leaf.name]
    MIMI_ancestor=tree.get_common_ancestor(MIMI_arr)

    PHYCO_arr=[leaf for leaf in tree if "PHYCO" in leaf.name]
    PHYCO_ancestor=tree.get_common_ancestor(PHYCO_arr)
    PHYCO_ancestor.detach()

    #test from MIMI ancestor
    test_node=MIMI_ancestor
    if test_node.up:
        while not [leaf for leaf in test_node.up if ("PP" not in leaf.name) and ("MIMI" not in leaf.name)]:
            test_node=test_node.up

    MIMI_query_reasonable_ancestor=test_node


    #calculate depth
    def get_depth(name):
        if name.isnumeric():
            return(int(name))
        else:
            return(0)

    def cal_depth(node):
        if node.is_leaf():
            return
        [cal_depth(child) for child in node.get_children()]
        node.name = node.name+str(max(get_depth(child.name) for child in node.get_children())+1)

    cal_depth(MIMI_query_reasonable_ancestor)        
    layer_sum=int(MIMI_query_reasonable_ancestor.name)+1

    #total length = 3
    intervention=3/layer_sum

    #appoint length
    def set_dist(node):
        if node.is_leaf():
            node.dist=intervention*int(node.up.name)
            return
            
        [set_dist(child) for child in node.get_children()]
        
        if int(node.name)==layer_sum-1:
            node.dist=intervention
            return
        node.dist=intervention*(int(node.up.name)-int(node.name))

    set_dist(MIMI_query_reasonable_ancestor)

    #gen from left leaf
    def Gen_from_left(node):        
        yield node
        if not node.is_leaf():
            for leaf in Gen_from_left(node.get_children()[0]):
                yield leaf
            for leaf in Gen_from_left(node.get_children()[1]):
                yield leaf

                
    #
    nMIMI=1
  
    for node in Gen_from_left(MIMI_query_reasonable_ancestor):
        if not node.is_leaf():
            node.set_style(ns_com)
            continue
        if "MIMI" in node.name:
            node.set_style(ns_MIMI)
                        
            for j in range(int(6*log_scale//color_scale)):
                color=color_MIMI if j%4 else "White"
                node.add_face(RectFace(color_scale,line_width,color,color), j, position="aligned")
            else:
                j+=1

                
                if nMIMI > mimi_min and nMIMI < mimi_max and nMIMI != mimi_lab:
                    pass
                elif nMIMI == olpv_max:
                    pass
                else:
                    if nMIMI==mimi_lab:
                        face_circle_number=TextFace("〜",fsize=200,fgcolor=color_MIMI)
                    elif nMIMI==olpv_min:
                        face_circle_number=TextFace(chr(nMIMI+9311)+chr(nMIMI+9312),fsize=200,fgcolor=color_MIMI)
                    else:
                        face_circle_number=TextFace(chr(nMIMI+9311),fsize=200,fgcolor=color_MIMI)
                    node.add_face(face_circle_number,j,position="aligned")
            nMIMI+=1
            
        else: 
            node.set_style(ns_com)
            size=int(node.name.split("U")[1])
            if size==1:
                node.detach()
                continue
            # stick_length=log_scale*math.log(size,10)
            stick_length=log_scale*6

            for j in range(int(6*log_scale//color_scale)):
                color="White"
                if OTU_barcode_table.iloc[:,j//int(6*log_scale//color_scale//60)][node.name.split("U")[0]] > 1:
                    if j//int(6*log_scale//color_scale//60)==nPP:
                        # color="#FF"+str(hex(255-j)).split("x")[-1].upper().zfill(2)+"00"
                        color="OrangeRed"

                node.add_face(RectFace(color_scale,stick_width,color,color), j, position="aligned")


    MIMI_query_reasonable_ancestor.render(output_fig_path,w=6400,units="px",tree_style=ts_sub_c)  