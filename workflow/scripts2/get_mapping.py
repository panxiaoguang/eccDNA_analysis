import pysam as ps


def getMapping(ipt,out,sp):
   samfile=ps.AlignmentFile(ipt,"rb")
   mapped_chr = sum([x.mapped for x in samfile.get_index_statistics() if x.contig !="chrM"])
   mapped_MT = [x.mapped for x in samfile.get_index_statistics() if x.contig =="chrM"][0]
   unmapped = samfile.unmapped
   with open(out,"w") as f:
      f.write(sp+"\t"+str(mapped_chr)+"\t"+str(mapped_MT)+"\t"+str(unmapped)+"\n")
   
getMapping(snakemake.input[0],snakemake.output[0],snakemake.wildcards.sample)