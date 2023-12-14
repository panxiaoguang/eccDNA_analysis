library(maftools)
counts = maftools::gtMarkers(t_bam = "Kidney_WGS/alignData/11T.cs.rmdup.sort.bam",
                             n_bam = "Kidney_WGS/alignData/11N.cs.rmdup.sort.bam",
                             loci = "11N.germline.het.txt",
                             fa = "/dellfsqd2/ST_LBI/USER/panxiaoguang/app/data_repo/GRCh38/hg38full.fa")
                            
ascat.bc = maftools::prepAscat(t_counts = "11T.cs.rmdup.sort_nucleotide_counts.tsv",
                               n_counts = "11N.cs.rmdup.sort_nucleotide_counts.tsv",
                               sample_name = "11")
library(ASCAT)

ascat.bc = ascat.loadData(Tumor_LogR_file = "11.tumour.logR.txt", Tumor_BAF_file = "11.tumour.BAF.txt", Germline_LogR_file = "11.normal.logR.txt", Germline_BAF_file = "11.normal.BAF.txt", gender = 'XX', genomeVersion = "hg38") # isTargetedSeq=T for targeted sequencing data
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_example.txt", replictimingfile = "RT_example.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc) # penalty=25 for targeted sequencing data
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, write_segments = T, gamma = 1) # gamma=1 for HTS data
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')