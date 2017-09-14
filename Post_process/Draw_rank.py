import matplotlib.pyplot as plt
import math

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

# print("\t".join(str(size) for size in size_arr))  

size_arr=[size for size in size_arr if size > 1]
plt.clf()    
plt.xlabel('Rank')
plt.ylabel("Abundance")
plt.semilogy(range(len(size_arr)),sorted(size_arr,reverse=True), label="111")
# plt.plot(range(len(size_arr)),sorted(size_arr,reverse=True), label="111")
plt.savefig("/user1/scl1/yanzeli/Megaviridae/Figures/Rank_170912.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")