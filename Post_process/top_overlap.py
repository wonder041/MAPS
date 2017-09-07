import decode_clstr
import numpy as np
import pandas as pd
import itertools
from matplotlib import pyplot as plt

for PP in (16,45):
    clstr_path="/lustre1/aptmp/ideas2/yanzeli/Megaviridae/All_170905/8_POST/"+str(PP)+"97.fna.clstr"
    raw_df=decode_clstr.Calc_table(clstr_path)
    
    OTU_arr=[]
    read_arr=[]
    xlim=max(raw_df.astype(bool).sum())
    for top in range(xlim):
        jud_arr_gen=(raw_df.iloc[:,nitem] >= max(1,sorted(raw_df.iloc[:,nitem],reverse=True)[top]) for nitem in range(raw_df.shape[1]))
        jud_res=[any(jud) for jud in zip(*jud_arr_gen)]
        total_OTU=raw_df.iloc[jud_res].shape[0]
        overlap_OTU=raw_df.iloc[jud_res].all(1).sum()
        total_read=raw_df.iloc[jud_res].sum().sum()
        overlap_read=raw_df.iloc[jud_res].iloc[raw_df.iloc[jud_res].all(1).tolist()].sum().sum()
        OTU_arr.append(overlap_OTU/total_OTU*100)
        read_arr.append(overlap_read/total_read*100)
    
    plt.clf()
    plt.plot(range(xlim),OTU_arr, label="OTU")
    plt.plot(range(xlim),read_arr, label="read")
    plt.legend()
    plt.ylabel("Percentage of overlap for all samples")
    plt.xlabel("Number of Top OTU in each sample")
    plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/top_overlap_"+str(PP)+".png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")
