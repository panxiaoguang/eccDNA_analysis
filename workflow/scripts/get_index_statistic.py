import pysam as ps
import pandas as pd
from sys import argv

chroms=["chr"+str(x) for x in range(1,23)]
chroms.append("chrM")
chroms.append("chrX")
chroms.append("chrY")

def get_mappings(bamFile,sample,outdir,chr):
    samfile=ps.AlignmentFile(bamFile,"rb")
    df=pd.DataFrame(samfile.get_index_statistics()).set_index("contig")
    total_reads=sum(df.total)
    fin=df.loc[chr]
    fin["mappedRatio"]=fin["mapped"]/total_reads*100
    fin.to_csv(f"{outdir}/{sample}.mapping.stats.tsv",sep="\t")

## argv[1] is the bam file
## argv[2] is the sample name
## argv[3] is the output directory

get_mappings(argv[1],argv[2],argv[3],chroms)