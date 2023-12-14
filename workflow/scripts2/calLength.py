
with open(snakemake.input[0],'r') as f, open(snakemake.output[0],'w') as o:
    for line in f:
        if line.startswith("chrom"):
            continue
        else:
            line = line.split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            length = end - start
            o.write(str(length)+"\n")