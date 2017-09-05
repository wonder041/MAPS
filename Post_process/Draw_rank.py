import matplotlib.pyplot as plt
import math

clstr_path="/aptmp/yanzeli/Paper_pipeline/Outputs_170812/8_POSTPROCESS/S0_97.fna.clstr"

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

print("\t".join(str(size) for size in size_arr))  

size_arr=[math.log10(size) for size in size_arr if size > 1]
plt.clf()    
plt.plot(range(len(size_arr)),sorted(size_arr,reverse=True), label="111")
# plt.ylabel('Number of OTUs Observed')
# plt.xlabel("Number of Reads Sampled (million reads)")
# plt.legend(title="Percent Identity",bbox_to_anchor=(1.05, 1), loc=2, fontsize=12)
# plt.legend(title="Identity",bbox_to_anchor=(0, 1), loc=2,fontsize=7)
# plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/Rank_170904.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")