import os

folder_names = ["dbs","filter_bed","filterData","mappingState","start_anno","GCs","element_anno"]

for folder in folder_names:
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder {folder} created.")
    else:
        print(f"Folder {folder} check passed.")

db_names = ["dbs/hg38.chromo.size","dbs/hg38.coding.bed","dbs/hg38.fa","dbs/hg38.fa.fai","dbs/hg38.genetic_elements.merge.bed"]

if os.path.exists(db_names[0]):
    print(f"File {db_names[0]} check passed.")
else:
    print(f"File {db_names[0]} not found, you can download from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes")

if os.path.exists(db_names[1]):
    print(f"File {db_names[1]} check passed.")
else:
    print(f"File {db_names[1]} not found, you can download from biomart (https://www.ensembl.org/biomart/)")

if os.path.exists(db_names[2]) and os.path.exists(db_names[3]):
    print(f"File {db_names[2]} check passed.")
elif os.path.exists(db_names[2]) and not os.path.exists(db_names[3]):
    print(f"Genome {db_names[3]} index not found, you can build index by samtools faidx {db_names[2]}")
else:
    print(f"File {db_names[2]} not found, you can download from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")

if os.path.exists(db_names[4]):
    print(f"File {db_names[4]} check passed.")
else:
    print(f"File {db_names[4]} not found, you can download from UCSC TableBrowser.")


