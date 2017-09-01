from ete3 import Tree,NodeStyle,TreeFace,faces,TreeStyle,TreeStyle,AttrFace,ClusterTree,ProfileFace,RectFace,CircleFace,TextFace
import math
# import sys
# import Make_OTU_table

tree = Tree("/aptmp/yanzeli/Paper_pipeline/Outputs_170812/8_POSTPROCESS/90.sub.combo.faa.aln.rearranged.tre")
output_fig_path="/user1/scl1/yanzeli/Megaviridae/Figures/Tree_1000_170824.png"



#########################################
line_width=6
stick_width=15   #35
scale=600
log_scale=300   #150
color_scale=6   #4
color_MIMI="Black"  #"#4A19FF"

ns_MIMI = NodeStyle(size=0,hz_line_color=color_MIMI,vt_line_color=color_MIMI,hz_line_width=line_width,vt_line_width=line_width)
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


query_arr=[leaf for leaf in tree if "PP" in leaf.name]
MIMI_arr=[leaf for leaf in tree if "MIMI" in leaf.name]
MIMI_ancestor=tree.get_common_ancestor(MIMI_arr)

#test from MIMI ancestor
test_node=MIMI_ancestor
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
    child0,child1=node.get_children()
    cal_depth(child0)
    cal_depth(child1)
    node.name =node.name+str(max([get_depth(child0.name),get_depth(child1.name)])+1)

cal_depth(MIMI_query_reasonable_ancestor)        
layer_sum=int(MIMI_query_reasonable_ancestor.name)+1

#total length = 3
intervention=3/layer_sum

#appoint length
def set_dist(node):
    if node.is_leaf():
        node.dist=intervention*int(node.up.name)
        return
        
    child0,child1=node.get_children()
    set_dist(child0)
    set_dist(child1)
    
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
nMIMI=0
MIMI_name_arr=[]
for node in Gen_from_left(MIMI_query_reasonable_ancestor):
    if not node.is_leaf():
        node.set_style(ns_com)
        continue
    if "MIMI" in node.name:
        node.set_style(ns_MIMI)
               
        # node.dist=0
        # node.add_face(RectFace(node.dist*scale+5*log_scale,line_width,color_MIMI,color_MIMI),0,position="float")
        # node.add_face(face0,1,position="float")
        
        for j in range(int(6*log_scale//color_scale)):
            color=color_MIMI if j%4 else "White"
            node.add_face(RectFace(color_scale,line_width,color,color), j, position="aligned")
        else:
            j+=1
            nMIMI_re=nMIMI
        
            circle_number=chr(nMIMI_re+9312)
            
            
            if nMIMI_re > 3 and nMIMI_re < 15 and nMIMI_re != 6:
                pass
            elif nMIMI == 17:
                pass
            else:
                if nMIMI_re==6:
                    face_circle_number=TextFace("〜",fsize=200,fgcolor=color_MIMI)
                elif nMIMI_re==16:
                    face_circle_number=TextFace(chr(nMIMI_re+9312)+chr(nMIMI_re+9313),fsize=200,fgcolor=color_MIMI)
                else:
                    face_circle_number=TextFace(chr(nMIMI_re+9312),fsize=200,fgcolor=color_MIMI)
                node.add_face(face_circle_number,j,position="aligned")
                print(nMIMI_re)
               
            # face0 = TextFace(circle_number,fsize=200,fgcolor=color_MIMI)
            
            # node.add_face(TextFace(circle_number,fsize=200,fgcolor=color_MIMI),j,position="aligned")
            
            # if nMIMI_re not in range(1,12):
                # if nMIMI_re == 16:
                    # node.add_face(TextFace(chr(nMIMI_re+9312)+chr(nMIMI_re+9313),fsize=200,fgcolor=color_MIMI),j,position="aligned")
                # elif nMIMI_re == 17:
                    # pass
                # else:
                    # node.add_face(face0,j,position="aligned")
            # if nMIMI_re == 9:
                # node.add_face(TextFace("~",fsize=300,fgcolor=color_MIMI),j,position="aligned")
            
                
        legned_face = TextFace(circle_number+":"+node.name.split("_")[-1],fsize=120,fgcolor="Black")
        ts_sub_c.legend.add_face(legned_face,column=nMIMI//7)
        MIMI_name_arr.append(node.name.split("_")[-1])
        
        
        nMIMI+=1
        
    else: 
        node.set_style(ns_com)
        size=int(node.name.split("U")[1])
        if size==1:
            node.detach()
            continue
        stick_length=log_scale*math.log(size,10)
        # for j in range(int(stick_length/color_scale)):
        for j in range(int(6*log_scale//color_scale)):
            if j%(log_scale//color_scale)==0:
                color="LightGrey"
                node.add_face(RectFace(color_scale,1.5*stick_width,color,color), j, position="aligned")
            elif j < stick_length/color_scale:
                color="#FF"+str(hex(255-j)).split("x")[-1].upper().zfill(2)+"00"
                node.add_face(RectFace(color_scale,stick_width,color,color), j, position="aligned")
            else:
                color="White"
                node.add_face(RectFace(color_scale,stick_width,color,color), j, position="aligned")

# for nMIMI_re,MIMI_name in enumerate(MIMI_name_arr[6:]+MIMI_name_arr[:6]):
    # ts_sub_c.legend.add_face(TextFace(chr(nMIMI_re+9312)+":"+MIMI_name+" ",fsize=120,fgcolor="Black"),column=nMIMI_re//5)
    # print(chr(nMIMI_re+9312)+":"+MIMI_name)
    
num_arr=["1","10","100","1k","10k","100k"]

ts_sub_c.legend.add_face(RectFace(log_scale,2*stick_width,"White","White"),column=3)                
for j in range(int(6*log_scale//color_scale)+1):
    color="#FF"+str(hex(255-j)).split("x")[-1].upper().zfill(2)+"00" if j%(log_scale//color_scale) else "Black"
    ts_sub_c.legend.add_face(RectFace(color_scale,2*stick_width,color,color),column=j+4)
else:
    for k in range(6):
        ts_sub_c.legend.add_face(RectFace(log_scale,2*stick_width,"White","White"),column=j+5+k)
        ts_sub_c.legend.add_face(TextFace(num_arr[k],fsize=100,fgcolor="Black"),column=j+5+k)
        # ts_sub_c.legend.add_face(TextFace("1E"+str(k),fsize=100,fgcolor="Black"),column=j+6+k)
    
    circle_number_str="".join(chr(i+9312) for i in range(20))
    # ts_sub_c.legend.add_face(TextFace(circle_number_str,fsize=100,fgcolor=color_MIMI),column=j+1)
    
#################################################################

MIMI_query_reasonable_ancestor.render(output_fig_path,w=6400,units="px",tree_style=ts_sub_c)  