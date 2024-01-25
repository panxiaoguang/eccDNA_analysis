import os

def change_sufix(sample):
    filename = f"filter_bed/{sample}.ecc.bed"
    cmd_prefix = "awk -v OFS=\"\\t\" \'{sub(\"_[0-9]+$\", \"\", $4); print }\'"
    os.system(cmd_prefix + " " + filename + " > " + f"overlaps/{sample}.ecc.bed")

def merge_bed(samples,group):
    cmd = "cat "
    for sample in samples:
        cmd += f"overlaps/{sample}.ecc.bed "
    cmd += f"| sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o count_distinct > overlaps/{group}.merged_ecc.bed"
    os.system(cmd)

def add_group_info(group):
    cmd_prefix = "awk -v OFS=\"\\t\" \'{print $0, \"" + group + "\"}\' "
    os.system(cmd_prefix + " " + f"overlaps/{group}.merged_ecc.bed > overlaps/{group}.merged_ecc.tmp.bed")

def find_cluster():
    cmd = "cat overlaps/Case.merged_ecc.tmp.bed overlaps/Control.merged_ecc.tmp.bed |sort -k1,1 -k2,2n | bedtools cluster -i - > overlaps/Finally.out.bed"
    os.system(cmd)

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
    ## check if the path overlaps exists
    if not os.path.exists("overlaps"):
        os.mkdir("overlaps")
    Case_samples = readSamples("Case")
    Control_samples = readSamples("Control")
    ## step1:
    print("step1: change the sufix of the sample name")
    for sample in Case_samples:
        change_sufix(sample)
    for sample in Control_samples:
        change_sufix(sample)
    ## step2:
    print("step2: merge the bed files of each group")
    merge_bed(Case_samples,"Case")
    merge_bed(Control_samples,"Control")
    ## step3:
    print("step3: add group information")
    add_group_info("Case")
    add_group_info("Control")
    ## step4:
    print("step4: find the clusters")
    find_cluster()

if __name__ == "__main__":
    main()