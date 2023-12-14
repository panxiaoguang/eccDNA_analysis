import pandas as pd
import numpy as np
import glob
from sys import argv
## filter parameters: split >=2 circlescore >=200 sci>=0.33 eci >=0.33 sd<mean coveragecon <=0.1


def filt(file, sample, outdir, split, circlescore):
    df=pd.read_table(file)
    df["length"] = df.end- df.start
    df = df.loc[(df["soft-clipped"] >= split) & (df.score >= circlescore) & (df["start_ratio"] >=0.33) & (df["end_ratio"] >=0.33) & (df.continuity <=0.2) &(df["std"] < df["mean"])]
    df['ft']=np.where(df.length<=2000, 100,df.discordants)
    df = df.loc[df["ft"]>0]
    df = df.drop("ft",axis=1)
    df.to_csv(f'{outdir}/{sample}_circle_site.filter.tsv',sep="\t",index=None)

### 
### argv[1] = input file
### argv[2] = sample name
### argv[3] = output directory
### argv[4] = split
### argv[5] = circlescore

filt(argv[1],argv[2],argv[3],int(argv[4]),int(argv[5]))