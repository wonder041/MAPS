import random
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

def Calculate_point(size_arr):
    #discard singleton OTUs
    size_arr=[size for size in size_arr if size >1]
    total_number=sum(size_arr)
    
    #calculate loop number
    nboot=int(10000000/total_number)+1
    
    #calculate sample size
    sample_size_port_arr=[j*(0.5**i) for i in range(10) for j in [k/16 for k in range(17,9,-1)]]
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
        
def Darw_various_identity():

    input_dir="/aptmp/yanzeli/Paper_pipeline/Outputs_170812/8_POSTPROCESS/"
    plt.figure(figsize=(4.5,6))
    
    identity_point_dic={}
    identity_tup=(99,97,95,90)

    for identity in identity_tup:
        clstr_path=input_dir+"S0_"+str(identity)+".fna.clstr"
        if os.path.exists(clstr_path+".csv"):
            OTU_barcode_table=pd.read_csv(clstr_path+".csv")
        else:
            import decode_clstr
            OTU_barcode_table=decode_clstr.Calc_table(clstr_path)
            

        S0_table=OTU_barcode_table.loc[:,OTU_barcode_table.columns.to_series().apply(lambda x:x.startswith("S0")).tolist()] 
        S0_nosingleton_table= S0_table.loc[lambda x:x.sum(1)>1]
        size_arr=S0_nosingleton_table.sum(1).tolist()
        
        sample_size_arr, OTU_size_arr=Calculate_point(size_arr)
        
        # print(sample_size_arr, OTU_size_arr)
        plt.plot([x/1000000 for x in sample_size_arr],OTU_size_arr, label=str(identity)+"%")
        plt.annotate(xy=[sample_size_arr[0]/1000000-2,OTU_size_arr[0]-700],s=str(identity)+"%:"+str(int(OTU_size_arr[0]+0.5)))

    plt.ylabel('Number of OTUs Observed')
    plt.xlabel("Number of Reads Sampled (million reads)")
    plt.legend(title="Percent Identity",bbox_to_anchor=(1.05, 1), loc=2, fontsize=12)
    plt.legend(title="Identity",bbox_to_anchor=(0, 1), loc=2,fontsize=7)
    plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/Whole_rare4_170824.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")
        

def Darw_various_PP():
    clstr_path="/aptmp/yanzeli/Paper_pipeline/Outputs_170812/8_POSTPROCESS/S0_97.fna.clstr"
    
    if os.path.exists(clstr_path+".csv"):
        OTU_barcode_table=pd.read_csv(clstr_path+".csv")
    else:
        import decode_clstr
        OTU_barcode_table=decode_clstr.Calc_table(clstr_path)
    
    
    PP_points_dic={}
    for table_item in OTU_barcode_table.iteritems():
        barcode=table_item[1].name
        if "PP" not in barcode:
            continue
        PP=barcode.split("-")[-1]
        size_arr=table_item[1][lambda x:x>1].tolist()
        PP_points_dic[PP]=Calculate_point(size_arr)
    
    PP_arr=sorted(PP_points_dic.keys(),key=lambda x:PP_points_dic[x][0][0],reverse=True)
    
    
    fig, axes = plt.subplots(12, 5)
    fig.set_size_inches(21,30)
    for num,PP in enumerate(PP_arr):
        points=PP_points_dic[PP]
        axes[num//5,num%5].plot([x/1000 for x in points[0]],points[1])
        axes[num//5,num%5].locator_params(axis='x', nbins=6)
        axes[num//5,num%5].locator_params(axis='y', nbins=4)
        axes[num//5,num%5].text(1,0,PP,verticalalignment='bottom',horizontalalignment='right',transform=axes[num//5,num%5].transAxes,size=24 if "c" in PP else 48)
    
    num+=1
    axes[num//5,num%5].axis('off')
    axes[num//5,num%5].text(1,0,"x-axis: Number of Reads Sampled \n(thousand reads)\ny-axis: Number of OTUs Observed\n(singleton OTUs are excluded)",verticalalignment='bottom',horizontalalignment='right',transform=axes[num//5,num%5].transAxes,size=12)
    fig.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/PP_rare_170824.png",bbox_inches="tight",dpi=600)

    
def Darw_various_samples():
    clstr_path="/aptmp/yanzeli/Megaviridae/Outputs_170821/8_POSTPROCESS/97.fna.clstr"
    sample_description_dic={'S3':'Hot Spring', 'S5':'Osaka Bay Nov-2', 'S6':'Osaka Bay Nov-3', 'S2':'Japan Sea', 'S0':'Osaka Bay Nov', 'S4':'Mangrove', 'S1':'Osaka Bay Aug'}
    
    if os.path.exists(clstr_path+".csv"):
        OTU_barcode_table=pd.read_csv(clstr_path+".csv")
    else:
        import decode_clstr
        OTU_barcode_table=decode_clstr.Calc_table(clstr_path)

    sample_arr=sorted(set(barcode.split("-")[0] for barcode in OTU_barcode_table.columns))[:5]
    
    sample_points_dic={}
    for sample in sample_arr:
        sample_table=OTU_barcode_table.loc[:,OTU_barcode_table.columns.to_series().apply(lambda x:x.startswith(sample)).tolist()]
        sample_nosingleton_table= sample_table.loc[lambda x:x.sum(1)>1]
        
        size_arr=sample_nosingleton_table.sum(1).tolist()
        
        sample_points_dic[sample]=Calculate_point(size_arr)
        
    plt.clf()
    
    for sample,(sample_size_arr, OTU_size_arr) in sample_points_dic.items():
        # sample_size_arr, OTU_size_arr=Calculate_point(size_arr)
        
        
        # print(sample_size_arr, OTU_size_arr)
        plt.plot([x/1000000 for x in sample_size_arr],OTU_size_arr, label=sample_description_dic[sample])
        # print([sample_size_arr[0]/1000000-2,OTU_size_arr[0]-700])
        if sample=="S0":
            plt.annotate(xy=[sample_size_arr[0]/1000000-2.5,OTU_size_arr[0]-400],s=sample_description_dic[sample])
            plt.annotate(xy=[sample_size_arr[0]/1000000-2.5,OTU_size_arr[0]-650],s="reads:"+format(sample_size_arr[0], ',d')+" OTUs:"+format(int(OTU_size_arr[0]), ',d'),fontsize=6)
        elif sample=="S3":
            plt.annotate(xy=[sample_size_arr[0]/1000000+0.2,OTU_size_arr[0]-150],s=sample_description_dic[sample])
            plt.annotate(xy=[sample_size_arr[0]/1000000+0.2,OTU_size_arr[0]-400],s="reads:"+format(sample_size_arr[0], ',d')+" OTUs:"+format(int(OTU_size_arr[0]), ',d'),fontsize=6)
        else:
            plt.annotate(xy=[sample_size_arr[0]/1000000+0.2,OTU_size_arr[0]],s=sample_description_dic[sample])
            plt.annotate(xy=[sample_size_arr[0]/1000000+0.2,OTU_size_arr[0]-250],s="reads:"+format(sample_size_arr[0], ',d')+" OTUs:"+format(int(OTU_size_arr[0]), ',d'),fontsize=6)

    plt.ylabel('Number of OTUs Observed')
    plt.xlabel("Number of Reads Sampled (million reads)")
    plt.legend(title="Sampling point",bbox_to_anchor=(1, 0), loc=4, fontsize=10)
    plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/5samples_170824.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")


        
    
if __name__ == "__main__":
    # Darw_various_identity()
    # Darw_various_PP()
    Darw_various_samples()
    
    


