import sys
import os
import random
import matplotlib.pyplot as plt

# input_path_arr=["/user1/scl1/yanzeli/Megaviridae/Outputs_presentation/9_CDHIT/97.fna.clstr","/user1/scl1/florian/megaviridae_data/osaka_primermix_feb/Lis_pipeline_OUTPUT/9_CDHIT/Osaka_Bay_97.fna.clstr"]
input_path_arr=["/user1/scl1/yanzeli/Megaviridae/Outputs_presentation/9_CDHIT/97.fna.clstr"]
cate_str_dic={'A3':'Hot Sprint', 'AMP10':'Osaka Bay 3 MP10', 'A2':'Japan Sea', 'AMP20':'Osaka Bay 3 MP20', 'A0':'Osaka Bay Nov', 'AMP5':'Osaka Bay 3 MP5', 'A4':'Mangrove', 'A1':'Osaka Bay Aug'}

OTU_cate_size_dic_arr=[]

for input_path in input_path_arr:
    with open(input_path) as clstr_file:
        next(clstr_file)
        cate_size_dic={}
        for line in clstr_file:
            if line.startswith(">"):
                OTU_cate_size_dic_arr.append(cate_size_dic)
                cate_size_dic={}
                continue
            else:
                cate=line.split(">")[1].split("U")[0]
                size=int(line.split("S")[1].split("L")[0])
                cate_size_dic[cate]=cate_size_dic.get(cate,0)+size            
        else:
            OTU_cate_size_dic_arr.append(cate_size_dic)
cate_arr=list(set(cate for cate_size_dic in OTU_cate_size_dic_arr for cate in cate_size_dic.keys()))
# cate_size_arr_dic=dict((cate,[cate_size_dic[cate] for cate_size_dic in OTU_cate_size_dic_arr if cate_size_dic.get(cate,0) >1]) for cate in cate_arr)
cate_size_arr_dic=dict((cate,[cate_size_dic[cate] for cate_size_dic in OTU_cate_size_dic_arr if cate_size_dic.get(cate,0) >0]) for cate in cate_arr)
# print(cate_size_arr_dic.keys())

cate_arr=sorted([cate for cate in cate_size_arr_dic.keys()])
# [print(cate_str_dic[cate]) for cate in cate_arr]
# [print(len(cate_size_arr_dic[cate])) for cate in cate_arr]
# [print(sum(1 for size in cate_size_arr_dic[cate] if size > 1)) for cate in cate_arr]




# sys.exit(0)
   
def Get_point(size_arr_with_singleton):
    nboot=10
    size_arr=[size for size in size_arr_with_singleton if size >1]
    total_number=sum(size_arr)
    sample_size_port_arr=[j*(0.5**i) for i in range(10) for j in (1,0.875,0.75,0.625)]
    sample_size_port_arr.append(0)
    sample_size_arr=[int(total_number*port) for port in sample_size_port_arr]
    
    OTU_arr=[]
    [OTU_arr.append(OTU) for OTU,size in enumerate(size_arr) for i in range(size)]
    
    OTU_size_arr=[0 for i in range(len(sample_size_arr))]
    
    for i in range(nboot):
        random.shuffle(OTU_arr)
        for sample_number,sample_size in enumerate(sample_size_arr):
            OTU_size_arr[sample_number]+=len(set(OTU_arr[:sample_size]))/nboot 
    
    return sample_size_arr, OTU_size_arr


    
# plt.figure(figsize=(4.5,6))    
# plt.figure(figsize=(6,4.5))    
# for cate,size_arr in cate_size_arr_dic.items():
    # sample_size_arr, OTU_size_arr=Get_point(size_arr)    
    # plt.plot([size/1000000 for size in sample_size_arr],OTU_size_arr, label=cate_str_dic[cate])
    # if cate!="A0":
        # plt.annotate(xy=[sample_size_arr[0]/1000000+0.1,OTU_size_arr[0]],s=cate_str_dic[cate])
    # else:
        # plt.annotate(xy=[sample_size_arr[0]/1000000-1.5,OTU_size_arr[0]-400],s=cate_str_dic[cate])
# plt.ylabel('Number of OTUs Observed')
# plt.xlabel("Number of Reads Sampled (million reads)")
# plt.legend(title="Sample",bbox_to_anchor=(1, 0), loc=4,fontsize=7)
# plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/8samples_rare.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")

# plt.figure(figsize=(6,4.5))    
# for cate,size_arr in cate_size_arr_dic.items():
    # if cate =="A0" or cate == "A1":
        # continue
    # sample_size_arr, OTU_size_arr=Get_point(size_arr)    
    # plt.plot([size/1000 for size in sample_size_arr],OTU_size_arr, label=cate_str_dic[cate])
    # if cate!="A2":
        # plt.annotate(xy=[sample_size_arr[0]/1000+100,OTU_size_arr[0]],s=cate_str_dic[cate])
    # else:
        # plt.annotate(xy=[sample_size_arr[0]/1000-200,OTU_size_arr[0]-200],s=cate_str_dic[cate])
# plt.ylabel('Number of OTUs Observed')
# plt.xlabel("Number of Reads Sampled (thousand reads)")
# plt.legend(title="Sample",bbox_to_anchor=(1, 0), loc=4,fontsize=7)
# plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/6samples_rare.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")

    

plt.figure(figsize=(6,4.5))    
for cate,size_arr in cate_size_arr_dic.items():
    sample_size_arr, OTU_size_arr=Get_point(size_arr)    
    plt.plot([size/1000000 for size in sample_size_arr],OTU_size_arr, label=cate_str_dic[cate])
    if cate!="A0":
        plt.annotate(xy=[sample_size_arr[0]/1000000+0.1,OTU_size_arr[0]],s=cate_str_dic[cate])
    else:
        plt.annotate(xy=[sample_size_arr[0]/1000000-1.5,OTU_size_arr[0]-400],s=cate_str_dic[cate])
plt.ylabel('Number of OTUs Observed')
plt.xlabel("Number of Reads Sampled (million reads)")
plt.legend(title="Sample",bbox_to_anchor=(1, 0), loc=4,fontsize=7)
plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/5samples_rare.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
