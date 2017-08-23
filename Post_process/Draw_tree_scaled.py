from ete3 import Tree,NodeStyle,TreeFace,faces,TreeStyle,TreeStyle,AttrFace,ClusterTree,ProfileFace,RectFace,CircleFace,TextFace
import math
import sys
import Make_OTU_table

tree = Tree("/lustre1/aptmp/yanzeli/Megaviridae/Outputs_170310/Tree/test.tre")
rep_psdic_dic = Make_OTU_table.Get_dic("/lustre1/aptmp/yanzeli/Megaviridae/Outputs_170310/9_CDHIT/A0_90.faa.clstr")
output_fig_path="/user1/scl1/yanzeli/Megaviridae/Figures/Tree_1000_scaled.png"
#########################################
line_width=2
stick_width=15   #35
scale=600
log_scale=300   #150
color_scale=6   #4
color_MIMI="Black"  #"#4A19FF"
color_query="#12B217"
color_query_arr=["#12B217","Lime","LimeGreen","Green","DarkGreen","DarkGreen"]



ns_MIMI = NodeStyle(size=0,hz_line_color=color_MIMI,vt_line_color=color_MIMI,hz_line_width=line_width,vt_line_width=line_width)
ns_com = NodeStyle(size=0,hz_line_width=line_width,vt_line_width=line_width)

ts_sub_c=TreeStyle()
ts_sub_c.mode="c"
ts_sub_c.arc_start=-10
ts_sub_c.arc_span=360
ts_sub_c.scale=scale
ts_sub_c.draw_guiding_lines= True
ts_sub_c.guiding_lines_type=0
ts_sub_c.guiding_lines_color="LightGrey"
ts_sub_c.show_leaf_name = False
ts_sub_c.show_scale = False
ts_sub_c.allow_face_overlap = True
ts_sub_c.legend_position=4
#########################################


query_arr=[leaf for leaf in tree if "PP" in leaf.name]
MIMI_arr=[leaf for leaf in tree if "MIMI" in leaf.name]
MIMI_ancestor=tree.get_common_ancestor(MIMI_arr)
test_node=MIMI_ancestor
while not [leaf for leaf in test_node.up if ("PP" not in leaf.name) and ("MIMI" not in leaf.name)]:
    test_node=test_node.up
MIMI_query_reasonable_ancestor=test_node

# def get_depth(name):
    # if name.isnumeric():
        # return(int(name))
    # else:
        # return(0)

# def cal_depth(node):
    # if node.is_leaf():
        # return
    # child0,child1=node.get_children()
    # cal_depth(child0)
    # cal_depth(child1)
    # node.name =node.name+str(max([get_depth(child0.name),get_depth(child1.name)])+1)

# cal_depth(MIMI_query_reasonable_ancestor)        
# layer_sum=int(MIMI_query_reasonable_ancestor.name)+1
# intervention=3/layer_sum

# def set_dist(node):
    # if node.is_leaf():
        # node.dist=intervention*int(node.up.name)
        # return
        
    # child0,child1=node.get_children()
    # set_dist(child0)
    # set_dist(child1)
    
    # if int(node.name)==layer_sum-1:
        # node.dist=intervention
        # return
    # node.dist=intervention*(int(node.up.name)-int(node.name))

# set_dist(MIMI_query_reasonable_ancestor)

def Gen_from_left(node):        
    yield node
    if not node.is_leaf():
        for leaf in Gen_from_left(node.get_children()[0]):
            yield leaf
        for leaf in Gen_from_left(node.get_children()[1]):
            yield leaf

nMIMI=0
for node in Gen_from_left(MIMI_query_reasonable_ancestor):
    if not node.is_leaf():
        node.set_style(ns_com)
        continue
    if "MIMI" in node.name:
        node.set_style(ns_MIMI)        
        circle_number=chr(nMIMI+9312)   # if nMIMI >= 0 else chr(nMIMI+9312+18)
        face0 = TextFace(circle_number,fsize=100,fgcolor=color_MIMI)
        
        # node.dist=0
        # node.add_face(RectFace(node.dist*scale+5*log_scale,line_width,color_MIMI,color_MIMI),0,position="float")
        # node.add_face(face0,1,position="float")
        
        for j in range(int(5*log_scale//color_scale)):
            color=color_MIMI if j%4==0 else "White"
            node.add_face(RectFace(color_scale,line_width,color,color), j, position="aligned")
        # else:
            # j+=1
            # node.add_face(face0,j,position="aligned")
                      
        legned_face = TextFace(circle_number+":"+node.name.split("_")[-1],fsize=120,fgcolor="Black")
        # ts_sub_c.legend.add_face(legned_face,column=nMIMI//5)
        
        
        nMIMI+=1
        
    else: 
        node.set_style(ns_com)
        size=sum(rep_psdic_dic[node.name].values())
        if size==1:
            node.detach()
            continue
        stick_length=log_scale*math.log(size,10)
        for j in range(int(stick_length/color_scale)):
            color="#FF"+str(hex(255-j)).split("x")[-1].upper().zfill(2)+"00" if j%(log_scale//color_scale) else "Grey"
            node.add_face(RectFace(color_scale,stick_width,color,color), j, position="aligned")
    
for j in range(int(6*log_scale//color_scale)+1):
    color="#FF"+str(hex(255-j)).split("x")[-1].upper().zfill(2)+"00" if j%(log_scale//color_scale) else "Black"
    ts_sub_c.legend.add_face(RectFace(color_scale,2*stick_width,color,color),column=j+1)
else:
    ts_sub_c.legend.add_face(RectFace(log_scale,50*stick_width,"White","White"),column=j+2)
    
    circle_number_str="".join(chr(i+9312) for i in range(20))
    # ts_sub_c.legend.add_face(TextFace(circle_number_str,fsize=100,fgcolor=color_MIMI),column=j+1)
    
#################################################################

MIMI_query_reasonable_ancestor.render(output_fig_path,w=12800,units="px",tree_style=ts_sub_c)  