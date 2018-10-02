#!/aptmp/yanzeli/miniconda3/envs/full/bin/python

######################################################
import pandas as pd
import numpy as np
import os

clstr_path = "/aptmp/yanzeli/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/S090.faa.clstr"


def Calc_table(clstr_path):
    def abstract_barcode(line): return line.split(">")[1].split("N")[0]

    def abstract_size(line): return int(
        line.split("S")[2].split("F")[0].split("L")[0])

    def abstract_name(line): return line.split(">")[1].split(".")[0]

    series_lst = []
    barcode_size_series = pd.Series(dtype="int")

    with open(clstr_path, "r") as clstr_file:
        next(clstr_file)
        for line in clstr_file:
            if line.startswith('>Cluster'):
                series_lst += [barcode_size_series]
                barcode_size_series = pd.Series(dtype="int")
                continue
            barcode = abstract_barcode(line)
            size = abstract_size(line)
            if line.strip().endswith("*"):
                barcode_size_series.name = abstract_name(line)
            try:
                barcode_size_series[barcode] += size
            except KeyError:
                barcode_size_series[barcode] = size
        else:
            series_lst += [barcode_size_series]

    OTU_barcode_table = pd.DataFrame(series_lst).fillna(0).astype("int")
    try:
        OTU_barcode_table.to_csv(clstr_path+".csv")
    except:
        pass
    return OTU_barcode_table


OTU_barcode_table = Calc_table(clstr_path)

######################################################

import matplotlib.pyplot as plt

colors_f = plt.get_cmap("tab10").colors


def convert_color(color): return "#{:X}{:X}{:X}".format(
    *map(lambda x: int(256*x), color))


colors = list(map(convert_color, colors_f))

######################################################

from ete3 import Tree, NodeStyle, TreeFace, faces, TreeStyle, TreeStyle, AttrFace, ClusterTree, ProfileFace, RectFace, CircleFace, TextFace
import math


tree = Tree("/aptmp/yanzeli/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/S090.sub.combo.aln.rearranged_rename.tre")
# tree = Tree("/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/ref_amplicon/MIMI_S090_sub.tre")

output_fig_path = "/user1/scl1/yanzeli/Exchange/PP_tree_180327.png"

#########################################
line_width = 6
stick_width = 24
scale = 330
# scale=600
log_scale = 480
color_scale = 6
color_MEGA = "Black"

ns_MIMI = NodeStyle(size=0, hz_line_color="Red", vt_line_color="Red",
                    hz_line_width=line_width, vt_line_width=line_width)
ns_com = NodeStyle(size=0, hz_line_width=line_width, vt_line_width=line_width)
ns_MIMI_ancestor = NodeStyle(
    size=0, hz_line_width=line_width, vt_line_width=line_width, bgcolor="DarkSeaGreen")

ts = TreeStyle()
ts.mode = "c"
ts.arc_start = 0
ts.arc_span = 357
ts.scale = scale
ts.draw_guiding_lines = True
# ts_sub_c.extra_branch_line_type=0
# ts_sub_c.extra_branch_line_color="Black"
ts.guiding_lines_type = 2
# ts_sub_c.force_topology=True
# ts_sub_c.guiding_lines_color=color_query
ts.guiding_lines_color = "Grey"
ts.show_leaf_name = False
ts.show_scale = False
ts.allow_face_overlap = True
ts.legend_position = 4
#########################################

MIMI_min = 6
MIMI_max = 18
MIMI_lable = 14
OLPV_min = 3
OLPV_max = OLPV_min+1

#########################################
query_lst = [leaf for leaf in tree if "PP" in leaf.name]
MEGA_lst = [leaf for leaf in tree if "MIMI" in leaf.name]
MEGA_ancestor = tree.get_common_ancestor(MEGA_lst)

PHYCO_lst = [leaf for leaf in tree if "PHYCO" in leaf.name]
PHYCO_ancestor = tree.get_common_ancestor(PHYCO_lst)
PHYCO_ancestor.detach()

# test from MIMI ancestor
test_node = MEGA_ancestor
if test_node.up:
    while not [leaf for leaf in test_node.up if ("PP" not in leaf.name) and ("MIMI" not in leaf.name)]:
        test_node = test_node.up

MEGA_query_ancestor = test_node
#########################################
# calculate depth


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


cal_depth(MEGA_query_ancestor)
layer_sum = int(MEGA_query_ancestor.name)+1

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


set_dist(MEGA_query_ancestor)
#########################################
# gen from left leaf


def Gen_from_left(node):
    yield node
    if not node.is_leaf():
        for leaf in Gen_from_left(node.get_children()[0]):
            yield leaf
        for leaf in Gen_from_left(node.get_children()[1]):
            yield leaf


nMEGA = 1
MIMI_node_lst = []

for node in Gen_from_left(MEGA_query_ancestor):
    if not node.is_leaf():
        node.set_style(ns_com)
        continue
    if "MIMI" in node.name:
        node.set_style(ns_MIMI)

        for j in range(int(6*log_scale//color_scale)):
            color = color_MEGA if j % 4 else "White"
            node.add_face(RectFace(color_scale, line_width,
                                   color, color), j, position="aligned")
        else:
            j += 1

            if nMEGA > MIMI_min and nMEGA < MIMI_max and nMEGA != MIMI_lable:
                pass
            elif nMEGA == OLPV_max:
                pass
            else:
                if nMEGA == MIMI_lable:
                    face_circle_number = TextFace(
                        "〜", fsize=200, fgcolor=color_MEGA)
                elif nMEGA == OLPV_min:
                    face_circle_number = TextFace(
                        chr(nMEGA+9311)+chr(nMEGA+9312), fsize=200, fgcolor=color_MEGA)
                else:
                    face_circle_number = TextFace(
                        chr(nMEGA+9311), fsize=200, fgcolor=color_MEGA)
                node.add_face(face_circle_number, j, position="aligned")

        if 17 > nMEGA > 5:
            MIMI_node_lst.append(node)
        nMEGA += 1

    else:
        node.set_style(ns_com)
        size = int(node.name.split("U")[1])
        if size == 1:
            node.detach()
            continue
        # stick_length=log_scale*math.log(size,10)
        stick_length = log_scale*6

        for j in range(int(6*log_scale//color_scale)):
            color = "White"
            if OTU_barcode_table.iloc[:, j//int(6*log_scale//color_scale//60)][node.name.split("U")[0]] >= 1:
                # if j//int(6*log_scale//color_scale//60)==nPP:
                    # color="#FF"+str(hex(255-j)).split("x")[-1].upper().zfill(2)+"00"
                    # color="OrangeRed"
                color = colors[j //
                               int(6*log_scale//color_scale//60) % len(colors)]
            if j % 8 == 0:
                color = colors[j //
                               int(6*log_scale//color_scale//60) % len(colors)]

            node.add_face(RectFace(color_scale, stick_width,
                                   color, color), j, position="aligned")


MIMI_ancestor = MEGA_query_ancestor.get_common_ancestor(MIMI_node_lst)
# MIMI_ancestor.set_style(ns_MIMI_ancestor)


MEGA_query_ancestor.render(output_fig_path, w=6400, units="px", tree_style=ts)









def calc_pr():
    crov=[leaf for leaf in MEGA_query_ancestor if "CroV" in leaf.name][0]
    knv=[leaf for leaf in MEGA_query_ancestor if "KNV" in leaf.name][0]
    MEGA_query_ancestor.get_common_ancestor(knv,crov)
    crov_an=crov.up.up.up.up.up.up.up
    megamimi_an=MEGA_query_ancestor.get_common_ancestor(knv,crov)
    mid_an=megamimi_an.children[0].children[1]
    mimi_an=megamimi_an.children[0].children[0].children[1]

    megamimi_name_lst=[leaf.name for leaf in crov_an]+[leaf.name for leaf in mid_an]+[leaf.name for leaf in mimi_an]
    megamimi_name_lst=list(filter(lambda name:"S0" in name,megamimi_name_lst))
    all_name_lst=list(filter(lambda name:"S0" in name,(leaf.name for leaf in MEGA_query_ancestor)))


    all_num=sum(OTU_barcode_table.sum(1)[name.split("U")[0]] for name in all_name_lst)
    megamimi_num=sum(OTU_barcode_table.sum(1)[name.split("U")[0]] for name in megamimi_name_lst)
