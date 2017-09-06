import decode_clstr
import numpy as np
import pandas as pd
import itertools
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles

output_path="/user1/scl1/yanzeli/Megaviridae/Figures/venn_170906.png"
ngroups=3
target_tup_arr=list(itertools.chain.from_iterable(itertools.combinations(range(ngroups),npick) for npick in range(1,ngroups+1)))


for PP in (16,45):
    clstr_path="/lustre1/aptmp/ideas2/yanzeli/Megaviridae/All_170905/8_POST/"+str(PP)+"97.fna.clstr"
    raw_df=decode_clstr.Calc_table(clstr_path)
    
    label_arr=raw_df.columns.tolist()
    
    ns_df=raw_df[raw_df.sum(1) > 1]
    
    range_df_dic={"ns":ns_df}
    for ntop in (100,500):

        jug_arr_gen=(raw_df.iloc[:,ncol] > sorted(raw_df.iloc[:,ncol],reverse=True)[ntop] for ncol in range(ngroups))
        jug_arr_res=[any(jug_arr) for jug_arr in zip(*jug_arr_gen)]    
        top_df=raw_df[jug_arr_res]
        range_df_dic["top"+str(ntop)]=top_df

    
    

    for dfrange,query_df in range_df_dic.items():
        value_OTU_arr=[]
        value_read_arr=[]
        for target_tup in target_tup_arr:
            judge_gen=((query_df.iloc[:,ngroup] > 0 if ngroup in target_tup else query_df.iloc[:,ngroup] == 0) for ngroup in range(ngroups))
            judge_res=[all(judge_item) for judge_item in zip(*judge_gen)]
            
            value_OTU=query_df[judge_res].shape[0]
            value_OTU_arr.append(value_OTU)
            
            value_read=query_df[judge_res].sum().sum()
            value_read_arr.append(value_read)
            
        
        plt.clf()
        v = venn3(subsets=value_OTU_arr, set_labels=label_arr)
        for npos,pos in enumerate(("100","010","001","110","101","011","111")):
            if v.get_label_by_id(pos):
                v.get_label_by_id(pos).set_text(str(value_OTU_arr[npos])+" OTUs: "+str(int(value_OTU_arr[npos]/sum(value_OTU_arr)*100))+"%\n"+str(value_read_arr[npos])+" reads: "+str(int(value_read_arr[npos]/sum(value_read_arr)*100))+"%")
            
        plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/venn_170906"+str(PP)+dfrange+".png",bbox_inches='tight',pad_inches=0.2,dpi=1000,orientation="portrait")

    
    



    
