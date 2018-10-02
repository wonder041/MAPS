#!/aptmp/yanzeli/miniconda3/envs/mainenv/bin/python
# import random
# import pandas as pd
# import numpy as np
# import os
import matplotlib.pyplot as plt
# import multiprocessing
# import itertools
import json

def Calc_slope(sample_size_arr,OTU_size_arr):
    dif_sample_size=sample_size_arr[0]-sample_size_arr[1]
    dif_OTU_size=OTU_size_arr[0]-OTU_size_arr[1]
    return dif_OTU_size/dif_sample_size*1000000
    

with open("/aptmp/yanzeli/Paper_pipeline/Scripts/Post_process/identity_point_dic") as dic_file:
    identity_point_dic=json.load(dic_file)


clstr_path="/aptmp/yanzeli/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/S097.fna.clstr"

def Calc_size_arr(clstr_path):
    abstract_size=lambda line:int(line.split("S")[2].split("F")[0].split("L")[0])

    size_arr=[]

    with open(clstr_path,"r") as clstr_file:
        size=0
        next(clstr_file)
        for line in clstr_file:
            if line.startswith('>Cluster'):
                size_arr.append(size)
                size=0
                continue
            size += abstract_size(line)
        else:
            size_arr.append(size)
    return size_arr
    
size_arr=Calc_size_arr(clstr_path)
size_arr=[size for size in size_arr if size > 1]




plt.close("all")
fig=plt.figure(figsize=(7,3),dpi=1200)
n_space=42
fig.suptitle("A"+" "*n_space+"B"+" "*n_space,size="xx-large")

ax = fig.add_subplot(121)
for identity in sorted(identity_point_dic,reverse=True):
    sample_size_arr, OTU_size_arr=identity_point_dic[identity]
    slope=Calc_slope(sample_size_arr,OTU_size_arr)
    lable=format(slope, '0.2f')+"/1M"
    ax.plot([x/1000000 for x in sample_size_arr],OTU_size_arr, label=str(identity)+"%")
    ax.annotate(xy=[sample_size_arr[0]/1000000-2.1,OTU_size_arr[0]-1400],s=lable)

ax.xaxis.label.set_size("large")
ax.yaxis.label.set_size("large")
ax.set_ylabel('Number of OTUs Observed',size="large")
ax.set_xlabel("Number of Reads Sampled (million reads)",size="large")
# plt.legend(title="Percent Identity",bbox_to_anchor=(1.05, 1), loc=2, fontsize="large")
ax.legend(title="Identity",bbox_to_anchor=(0, 1), loc=2)

ax = fig.add_subplot(122)
ax.xaxis.label.set_size("large")
ax.yaxis.label.set_size("large")
ax.set_xlabel('Rank',size="large")
ax.set_ylabel("Abundance",size="large")
ax.semilogy(range(len(size_arr)),sorted(size_arr,reverse=True), label="111")

fig.tight_layout()
fig.subplots_adjust(top=0.8)

plt.savefig("/user1/scl1/yanzeli/Exchange/Figure_2_20180613.png",bbox_inches='tight',pad_inches=0.1,orientation="portrait")
plt.savefig("/user1/scl1/yanzeli/Exchange/Figure_2_20180613.eps",bbox_inches='tight',pad_inches=0.1,orientation="portrait")