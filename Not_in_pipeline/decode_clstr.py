import sys
import pandas as pd
import numpy as np

clstr_path=sys.argv[1]


abstract_barcode=lambda line:line.split(">")[1].split("N")[0]
abstract_size=lambda line:int(line.split("S")[2].split("F")[0].split("L")[0])

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
        try:
            barcode_size_series[barcode]+=size
        except "KeyError":
            barcode_size_series[barcode]=size
 
OTU_barcode_table=pd.DataFrame(series_arr).fillna(0).astype("int")
OTU_barcode_table.to_csv(clstr_path+".csv")

OTU_barcode_table=pd.read_csv(clstr_path+".csv")