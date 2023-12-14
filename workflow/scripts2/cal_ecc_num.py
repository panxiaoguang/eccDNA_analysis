
i = 0
with open(snakemake.input[0],'r') as f:
    for line in f:
        i+=1

with open(snakemake.output[0],'w') as f:
    f.write(str(i))