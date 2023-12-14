import pysam as ps
from sys import argv

def get_mappings(bamFile,sample):
    samfile=ps.AlignmentFile(bamFile,'rb')
    totals=sum([x.mapped for x in samfile.get_index_statistics() if x.contig !="chrM"])
    return sample, totals

def main(bamFile,sample,outdir):
    with open(f"{outdir}/{sample}_total_mappings.txt","w+") as out:
        sample,dt = get_mappings(bamFile,sample)
        out.write(f"{sample}\t{dt}\n")
    
## argv[1] is the bam file
## argv[2] is the sample name
## argv[3] is the output directory

main(argv[1],argv[2],argv[3])