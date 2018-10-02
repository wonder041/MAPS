import decode_clstr
import collections
# from operator import itemgetter
import matplotlib.pyplot as plt
import sys

clstr_path="/user1/scl1/yanzeli/aptmp/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/PP_overlap/S0_du.fna.clstr"
OTU_barcode_table=decode_clstr.Calc_table(clstr_path)




bar_counter=collections.Counter((OTU_barcode_table > 0).sum(1))
x1,y1=zip(*sorted(bar_counter.items()))
x2,y2=(OTU_barcode_table >0).sum(1),OTU_barcode_table.sum(1)

plt.close("all")
fig=plt.figure(figsize=(6,2.5),dpi=1200)
n_space=35
fig.suptitle("A"+" "*n_space+"B"+" "*n_space,size="xx-large")


ax = fig.add_subplot(121)
ax.set_yscale('log')
ax.xaxis.label.set_size("large")
ax.yaxis.label.set_size("large")
ax.bar(x1,y1,color="Black")
ax.set_xlabel('Number of primer pairs',size="large")
ax.set_ylabel('Number of genotypes',size="large")

ax = fig.add_subplot(122)
ax.set_yscale('log')
ax.xaxis.label.set_size("large")
ax.yaxis.label.set_size("large")
ax.plot(x2,y2,".",color="Black")
ax.set_xlabel('Number of primer pairs',size="large")
ax.set_ylabel('Number of reads',size="large")

fig.tight_layout()
fig.subplots_adjust(top=0.8)
# plt.figure(figsize=(6,2.5),dpi=1200)
plt.savefig("/user1/scl1/yanzeli/Exchange/PP_overlap_180409.eps",bbox_inches='tight',pad_inches=0.1,orientation="portrait")
plt.savefig("/user1/scl1/yanzeli/Exchange/PP_overlap_180409.png",bbox_inches='tight',pad_inches=0.1,orientation="portrait")

sys.exit(0)


# plt.xticks(fontsize=5)
# plt.yticks(fontsize=5)

# plt.ylim(ymin=0.5)
# PP_nread_series=OTU_barcode_table.groupby((OTU_barcode_table > 0).sum(1)).sum().sum(1)
# plt.bar(PP_nread_series.index,PP_nread_series,color="Black")

# ax.subplot(122)


# plt.clf()
# plt.yscale('log')
# plt.plot(PP_nread_series.index,PP_nread_series,"o")
# plt.plot(x,y,".")
# plt.savefig("/user1/scl1/yanzeli/Exchange/PP_overlap_point_180309.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")

# plt.clf()
# plt.yscale('log')
# all_data=[a[1].sum(1).tolist() for a in OTU_barcode_table.groupby((OTU_barcode_table > 0).sum(1))]
# plt.boxplot(all_data,sym="r.",labels=PP_nread_series.index)
# plt.savefig("/user1/scl1/yanzeli/Exchange/PP_overlap_box_180313.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")

def draw_pie():
    PP_nread_series=OTU_barcode_table.groupby((OTU_barcode_table > 0).sum(1)).sum().sum(1)

    plt.clf()
    plt.pie(PP_nread_series,radius=1.2,labels=PP_nread_series.index, colors=plt.get_cmap("Set3").colors)
    plt.pie(PP_nread_series.groupby(PP_nread_series.index//5).sum(), radius=1, autopct='%1.1f%%',colors=plt.get_cmap("Set3").colors)
    # plt.pie(PP_nread_series.groupby(PP_nread_series.index//5).sum(),labels=["1-5","6-10","10-15","16-20","20-25","26-30","30-35","36+"], radius=1, autopct='%1.1f%%',colors=plt.get_cmap("Set3").colors)
    plt.savefig("/user1/scl1/yanzeli/Exchange/PP_overlap_pie_180314.png",bbox_inches='tight',pad_inches=0.1,dpi=1000,orientation="portrait")
