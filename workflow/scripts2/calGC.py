from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqUtils import GC

refs = Fasta(snakemake.input['ref'])
with open(snakemake.input['ecc'],'r') as f , open(snakemake.output[0],'w') as o:
    for line in f:
        if line.startswith("chrom"):
            continue
        else:
            line = line.split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            seq = refs[chrom][start:end].seq
            gcs = GC(Seq(seq))
            o.write(str(gcs)+"\n")

