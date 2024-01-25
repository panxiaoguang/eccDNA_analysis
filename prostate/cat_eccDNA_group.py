import os

def change_sufix(sample):
    filename = f"filter_bed/{sample}.ecc.bed"
    cmd_prefix = "awk -v OFS=\"\\t\" \'{sub(\"_[0-9]+$\", \"\", $4); print }\'"
    os.system(cmd_prefix + " " + filename + " > " + f"overlaps/{sample}.ecc.bed")

def merge_bed(samples,group):
    cmd = "cat "
    for sample in samples:
        cmd += f"overlaps/{sample}.ecc.bed "
    cmd += f"| sort -k1,1 -k2,2n  > overlaps/{group}.all_ecc.bed"
    os.system(cmd)

def clean():
    os.system("rm -rf overlaps/*.ecc.bed")

def readSamples(group):
    filename = f"{group}_sampels.tsv"
    samples = []
    with open(filename) as f:
        for line in f:
            if line.startswith("samples"):
                continue
            else:
                sample = line.strip("\n")
                samples.append(sample)
    return samples

def main():
    Case_samples = readSamples("Case")
    Control_samples = readSamples("Control")
    for sample in Case_samples:
        change_sufix(sample)
    for sample in Control_samples:
        change_sufix(sample)
    merge_bed(Case_samples,"Case")
    merge_bed(Control_samples,"Control")
    clean()

if __name__ == "__main__":
    main()