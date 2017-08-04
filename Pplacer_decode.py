import sys
import json
import re
from Bio import Phylo
from io import StringIO

jplace_path = sys.argv[1]
jplace_file = open(jplace_path, "r")
jplace_dic = json.load(jplace_file)


tree = Phylo.read(StringIO(jplace_dic["tree"]), "newick")

name_position_pattern = re.compile(r"([-\w]*):[\d\.\-e]+\{(\d+)\}")
name_position_tuple_arr = name_position_pattern.findall(jplace_dic["tree"])
terminal_position_dic = dict((name_position_tuple[1], name_position_tuple[
                             0]) for name_position_tuple in name_position_tuple_arr if name_position_tuple[0] != "")

remove_mark = lambda x: x.replace("{", "").replace("}", "")


def get_nonterminal_position(nonterminal_clade):
    name_arr = [terminal_position_dic[remove_mark(
        clade.name)] for clade in nonterminal_clade.get_terminals()]
    return (remove_mark(nonterminal_clade.name), name_arr)

nonterminal_position_gen = (get_nonterminal_position(
    nonterminal_clade) for nonterminal_clade in tree.get_nonterminals())

name_position_dic = dict((position, [terminal_position_dic[
                         position]]) for position in terminal_position_dic)
name_position_dic.update(nonterminal_position_gen)

name_cate_dic = {}
for placements_dic in jplace_dic["placements"]:
    query_name_arr = [nm[0]for nm in placements_dic["nm"]]
    query_position = str(placements_dic["p"][0][1])
    query_cate_set = set(cate_name.split(
        "_")[0] for cate_name in name_position_dic[query_position])

    if query_cate_set == {"MIMI"}:
        print(*query_name_arr,sep="\n")

    # for query_name in query_name_arr:
       # name_cate_dic[query_name]=query_cate_set

# import Check_OTU_size
# name_size_dic=Check_OTU_size.Get_dic(jplace_path.replace("combo.jplace","clstr"))

# seq_num=0
# mimi_num=0
# for name in name_cate_dic:
    # num=name_size_dic[name]
    # seq_num+=num
    # if name_cate_dic[name]=={"MIMI"}:mimi_num+=num

# print(mimi_num,seq_num,mimi_num/seq_num)
