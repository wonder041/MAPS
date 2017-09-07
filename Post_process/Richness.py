from skbio.diversity import *
import pandas as pd

df=pd.read_csv("/aptmp/yanzeli/Paper_pipeline/Outputs_170812/8_POSTPROCESS/S0_97.fna.clstr.csv")
df1=df.iloc[:,1:]
print(alpha_diversity('michaelis_menten_fit',df1.sum(1)))

