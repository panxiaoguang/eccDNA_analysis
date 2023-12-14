import pandas as pd
from functools import reduce
from sys import argv
def readIN(upfile,mainfile,downFile,sample,outdir):
    df1=pd.read_csv(upfile,sep="\t")
    df1.columns=["ecc","upstream"]
    df2=pd.read_csv(mainfile,sep="\t")
    df2.columns=["ecc","self"]
    df3=pd.read_csv(downFile,sep="\t")
    df3.columns=["ecc","downstream"]
    fin=reduce(lambda x,y:pd.merge(x,y,on="ecc",how="outer"), [df1,df2,df3])
    fin.to_csv(f"{outdir}/{sample}.gcContents.txt",sep="\t",index=False)

readIN(argv[1],argv[2],argv[3],argv[4],argv[5])
