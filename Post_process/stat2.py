import os
import pandas as pd
import multiprocessing
import decode_clstr
import itertools

target_dir="/aptmp/yanzeli/Paper_pipeline/Outputs_ALL_170907_0/Post_pipeline/"
identity_tup=(99,97,95,90)

def Calc_clstr(identity): 
    clstr_path=target_dir+"S0"+str(identity)+".fna.clstr"
    OTU_barcode_table=decode_clstr.Calc_table(clstr_path)
    
    if identity==97:
        PP_arr_gen=((PP,(OTU_barcode_table[PP] == 1).sum(),(OTU_barcode_table[PP] > 1).sum()) for PP in sorted(OTU_barcode_table.columns))
        df=pd.DataFrame(PP_arr_gen)
        df.to_csv("All_170907_0_PP_OTU.stat")
    nOTU=OTU_barcode_table.shape[0]
    nsin=(OTU_barcode_table.sum(1) == 1).sum() 
    nnosin=(OTU_barcode_table.sum(1) > 1).sum() 
    nnosinread=OTU_barcode_table.iloc[(OTU_barcode_table.sum(1) > 1).tolist()].sum().sum()
    
    return (identity,nOTU,nsin,nnosin,nnosinread)

pool = multiprocessing.Pool(4)
identity_arr_arr = list(pool.imap_unordered(Calc_clstr, identity_tup))

df=pd.DataFrame(identity_arr_arr).T
df.to_csv("All_170907_0_OTU.stat")

    

    
    
    

