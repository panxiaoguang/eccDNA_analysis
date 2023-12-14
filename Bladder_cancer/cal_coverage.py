import pysam as ps
import numpy as np

def readBED(bed):
    bedlst=[]
    with open(bed,"r") as f:
        for line in f:
            chrom,start,end = line.split("\t")
            end=int(end.strip("\n"))
            start=int(start)
            bedlst.append([chrom,start,end])
    return bedlst
'''
def cal_coverage(sam,contig, start, end,totals):
    cov = sam.count_coverage(contig, start, end, quality_threshold = 0)
    summarized_cov = np.array([cov[0], cov[1], cov[2], cov[3]]).sum(axis=0)
    mean_cov = np.mean(summarized_cov)
    lgth = (end-start)/1000
    RPKM = mean_cov/(totals*lgth)
    return RPKM
'''

def cal_coverage(sam,contig, start, end,totals):
    reads = sam.count(contig,start,end)
    lgth = (end-start)/1000
    RPKM = reads/(totals*lgth)
    return RPKM

def comp_cov(sample):
    print("sample\tregion\tWGST\tWGSN\tcircleT\tcircleN")
    bedlst=readBED(f"/zfsqd1/ST_LBI/PROJECT/panxiaoguang/tmp/ecDNA/{sample}.merge.ecDNA.bed")
    case=ps.AlignmentFile(f"WGS/alignData/{sample}T.cs.rmdup.sort.bam","rb")
    control=ps.AlignmentFile(f"WGS/alignData/{sample}N.cs.rmdup.sort.bam","rb")
    case2=ps.AlignmentFile(f"circleSeq/alignData/sorted_cBca_{sample}T_circle.bam","rb")
    control2=ps.AlignmentFile(f"circleSeq/alignData/sorted_cBca_{sample}N_circle.bam","rb")
    total_mappings1 = np.sum([x.mapped for x in case.get_index_statistics()])/1000000
    total_mappings2 = np.sum([x.mapped for x in control.get_index_statistics()])/1000000
    total_mappings3 = np.sum([x.mapped for x in case2.get_index_statistics()])/1000000
    total_mappings4 = np.sum([x.mapped for x in control2.get_index_statistics()])/1000000
    for region in bedlst:
        cs,ct,cs2,ct2 = cal_coverage(case,region[0],region[1],region[2],total_mappings1),cal_coverage(control,region[0],region[1],region[2],total_mappings2),cal_coverage(case2,region[0],region[1],region[2],total_mappings3),cal_coverage(control2,region[0],region[1],region[2],total_mappings4)
        print(sample,"\t",f"{region[0]}:{region[1]}-{region[2]}","\t",cs,"\t",ct,"\t",cs2,"\t",ct2,sep="")

samples=["13","14","15","16","17","1",
        "21","23","24","25","28","29",
        "2","31","34","36","37","40",
        "41","46","47","50","51","5",
        "60","63","65","67","69","70",
        "71","72","74","75","76","77",
        "79","80","84","85","86","87",
        "8","9"]

for sm in samples:
    comp_cov(sm)
