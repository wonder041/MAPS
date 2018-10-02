import pandas as pd
import numpy as np
import os

def Calc_table(clstr_path,reload=False):
    try:
        return pd.read_csv(clstr_path+".csv",index_col=0)
    except:
        pass


    abstract_barcode=lambda line:line.split(">")[1].split("N")[0]
    abstract_size=lambda line:int(line.split("S")[2].split("F")[0].split("L")[0])
    abstract_name=lambda line:line.split(">")[1].split(".")[0]

    series_arr=[]
    barcode_size_series=pd.Series(dtype="int")


    with open(clstr_path,"r") as clstr_file:
        next(clstr_file)
        for line in clstr_file:
            if line.startswith('>Cluster'):
                series_arr+=[barcode_size_series]
                barcode_size_series=pd.Series(dtype="int")
                continue
            barcode=abstract_barcode(line)
            size=abstract_size(line)
            if line.strip().endswith("*"):
                barcode_size_series.name=abstract_name(line)
            try:
                barcode_size_series[barcode]+=size
            except KeyError:
                barcode_size_series[barcode]=size
        else:
            series_arr+=[barcode_size_series]
     
    OTU_barcode_table=pd.DataFrame(series_arr).fillna(0).astype("int")
    try:
        OTU_barcode_table.to_csv(clstr_path+".csv")
    except:
        pass
    return OTU_barcode_table


if __name__ == "__main__":
    import sys
    clstr_path=sys.argv[1]
    OTU_barcode_table=Calc_table(clstr_path)
    # OTU_barcode_table=pd.read_csv(clstr_path+".csv")
