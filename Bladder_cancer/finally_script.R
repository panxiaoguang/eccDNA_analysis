
# Figure2A--data ----------------------------------------------------------

df <- openxlsx::read.xlsx("~/Project/WGS/allsamples.ecDNA_det.xlsx") %>%
  as_tibble()
### define amplicon type
amp_type <- df %>%
  mutate(
    type = amplicon_decomposition_class,
    type = if_else(`ecDNA+` == "Positive", "ecDNA", type),
    type = if_else(`BFB+` == "Positive", "BFB", type)
  )

get_sample_numbers <- function(x) {
  amp_type %>%
    filter(type == x) %>%
    pull(sample_name) %>%
    unique() %>%
    length()
}
five_class <- c(
  "ecDNA", "BFB", "Complex non-cyclic",
  "Linear amplification", "No amp/Invalid"
)
amp_type$type <- factor(amp_type$type, levels = five_class)
## define sample type

  
sample_tables<-tibble(class=five_class,numbers=sapply(five_class,get_sample_numbers))
sample_tables

# Figure2B---data ---------------------------------------------------------

sample_tables <- amp_type %>%
  arrange(sample_name, type) %>%
  mutate(type = as.character(type)) %>%
  group_by(sample_name) %>%
  tidyr::nest() %>%
  mutate(sampletype = purrr::map(data, function(x) {
    (x$type)[[1]]
  })) %>%
  select(sample_name, sampletype) %>%
  ungroup() %>%
  group_by(sampletype) %>%
  summarise(count = n())

Pvalue <- fisher.test(matrix(c(28, 45, 29, 3, 10, 5, 4, 4, 29, 23), nrow = 2))

Pvalue$p.value


# Figure2D,2E,2F -----data ------------------------------------------------

amp_class <- openxlsx::read.xlsx("下机信息（含barcode）/all_sum.xlsx") %>%
  as_tibble(0)

newClass <- amp_class %>%
  as_tibble() %>%
  filter(AmpliconType != "No amp/Invalid") %>%
  select(1:3, 7)
names(newClass) <- c("sample", "amp", "type", "length")
breakPoints <- read_tsv("graphs/all.detective.breakpoints.tsv")

breakPoints <- breakPoints %>%
  left_join(newClass, by = c("sample", "amp"))
wocao <- breakPoints %>%
  mutate(name = paste(sample, amp, sep = "_"), start = bp - 1) %>%
  select(chrom, start, bp, type)
wocao <- na.omit(wocao)

plotData <- breakPoints %>%
  na.omit() %>%
  group_by(sample, amp, type, length) %>%
  summarise(count = n()) %>%
  mutate(ratio = count / length * 1000000)

plotData <- newClass %>%
  left_join(plotData, by = c("sample", "amp", "type", "length")) %>%
  tidyr::replace_na(replace = list(count = 0, ratio = 0, logratio = 0))

amp_class %>%
  filter(AmpliconType != "No amp/Invalid") %>%
  rstatix::wilcox_test(AmplifiedIntervalSize ~ AmpliconType) %>%
  filter(group1 == "ecDNA" | group2 == "ecDNA")

amp_class %>%
  filter(AmpliconType != "No amp/Invalid") %>%
  rstatix::wilcox_test(AverageAmplifiedCopyCount ~ AmpliconType) %>%
  filter(group1 == "ecDNA" | group2 == "ecDNA")

plotData %>%
  rstatix::wilcox_test(ratio ~ type) %>%
  filter(group1 == "ecDNA" | group2 == "ecDNA")
## ecDNA_linear = 1.64 e-16 , ecDNA_BFB = 3.11e-1,ecDNA_complex = 4.4e-2
## ecDNA_linear = 3.89 e-22 , ecDNA_BFB = 1.2e-1,ecDNA_complex = 7.15e-10
## ecDNA_linear = 4.87 e-13 , ecDNA_BFB = 2.7e-2,ecDNA_complex = 3.9e-1


# figure2G-----plot -------------------------------------------------------

########################################## oncoprint####################################################
######## CCGA PART#####################################################################################
df <- openxlsx::read.xlsx("all_detective_gene_list.xlsx")
need_genes1 <- openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx", sheet = 1)
need_genes1 <- need_genes1 %>%
  as_tibble() %>%
  filter(stringr::str_detect(Role, "oncogene"))
df1 <- df %>%
  filter(gene %in% need_genes1$Gene)
mat <- df1 %>%
  as_tibble() %>%
  select(sample_name, feature, gene) %>%
  mutate(feature = stringr::str_remove(feature, "_\\d")) %>%
  tidyr::pivot_wider(
    names_from = "sample_name",
    values_from = "feature",
    values_fn = function(x) {
      paste(x, collapse = ";")
    }, values_fill = ""
  ) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.matrix()

needs <- df1 %>%
  group_by(gene) %>%
  summarise(ct = n()) %>%
  arrange(desc(ct)) %>%
  filter(ct > 2) %>%
  pull(gene)
needs <- needs[1:20]

mat2 <- mat[needs, ]
#############################################################################################
######################################## TCGA PART################################################
TCGA <- openxlsx::read.xlsx("TCGA_UBC.xlsx")
TCGAbed <- TCGA %>%
  as_tibble() %>%
  mutate(info = paste(sample_barcode,
    amplicon_index,
    amplicon_classification,
    sep = "&"
  )) %>%
  select(amplicon_intervals, info) %>%
  tidyr::separate_rows(amplicon_intervals, sep = ",") %>%
  tidyr::separate(amplicon_intervals, into = c("chrom", "start", "end"), sep = "[:-]") %>%
  mutate(chrom = paste0("chr", chrom))

tcga_anno <- read_tsv("TCGA_data.anno.bed", col_names = F)
mat3 <- tcga_anno %>%
  select(X4, X8) %>%
  tidyr::separate(X4, into = c("sample", "amp", "feature"), sep = "&") %>%
  mutate(feature = case_when(
    feature == "BFB" ~ "BFB",
    feature == "Circular" ~ "ecDNA",
    feature == "Heavily-rearranged" ~ "Complex non-cyclic",
    feature == "Linear" ~ "Linear amplification"
  )) %>%
  filter(X8 %in% (need_genes1$Gene)) %>%
  select(sample, feature, X8) %>%
  dplyr::rename(gene = X8) %>%
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = "feature",
    values_fn = function(x) {
      paste(x, collapse = ";")
    }, values_fill = ""
  ) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.matrix()

mat4 <- mat3[needs, ]
#############################################################################################
final_data <- cbind(mat2, mat4)
nima <- tibble(genes = rownames(final_data), num = as.character(apply(final_data, 1, function(x) {
  sum(x == "ecDNA")
})))

col <- c(
  BFB = "#436693", `Complex non-cyclic` = "#C58D65",
  ecDNA = "#BC4137", `Linear amplification` = "#A1B2C0",
  unknown = "#ACABAC"
)
############################################## anno part#############################################
anno_table <- openxlsx::read.xlsx("TCGA-CCGA.xlsx", sheet = 2)
anno_table <- anno_table %>%
  mutate(label = stringr::str_extract(Sample, "\\d+"))
anno_table <- anno_table %>%
  filter(label %in% (colnames(mat2))) %>%
  tidyr::replace_na(replace = list("N" = "No avalible", "M" = "No avalible")) %>%
  select(label, Survival, Gender, Age, N, M, `NMIBC/MIBC`, Grade)

anno_table2 <- openxlsx::read.xlsx("TCGA-CCGA.xlsx", sheet = 1)
colnames(anno_table2)[6] <- "NMIBC/MIBC"
anno_table2 <- anno_table2 %>%
  filter(sample_barcode %in% (colnames(mat4))) %>%
  tidyr::replace_na(replace = list("NMIBC/MIBC" = "No avalible", "N" = "No avalible", "M" = "No avalible")) %>%
  select(`sample_barcode`, Survival, Gender, Age, N, M, `NMIBC/MIBC`, Grade) %>%
  dplyr::rename(label = "sample_barcode") %>%
  mutate(Grade = stringr::str_remove(Grade, " "))

fin_anno <- bind_rows(anno_table, anno_table2)
fin_anno <- fin_anno %>%
  tibble::column_to_rownames(var = "label") %>%
  as.data.frame()

fin_anno <- fin_anno[colnames(final_data), ]
###############################################################################################
############################################ plot oncoprint#######################################
ht <- oncoPrint(final_data,
  alter_fun = function(x, y, w, h, v) {
    n <- sum(v)
    h <- h * 0.9
    if (n) {
      grid.rect(x, y - h * 0.5 + 1:n / n * h, w * 0.9, 1 / n * h,
        gp = gpar(fill = col[names(which(v))], lwd = 0.5, col = "white"), just = "top"
      )
    } else {
      grid.rect(x, y, w, h, gp = gpar(fill = "#E4E4E4", lwd = 0.5, col = "white"))
    }
  }, col = col,
  row_names_gp = gpar(fontface = "italic", fontsize = 8.5),
  column_split = c(rep("CCGA", 46), rep("TCGA", 68)),
  row_order = order(nima$num, decreasing = T),
  top_annotation =
    HeatmapAnnotation( # cbar = anno_oncoprint_barplot(),
      Age = fin_anno$Age,
      N = fin_anno$N,
      M = fin_anno$M,
      Survival = fin_anno$Survival,
      Gender = fin_anno$Gender,
      `NMIBC/MIBC` = fin_anno$`NMIBC/MIBC`,
      Grade = fin_anno$Grade,
      simple_anno_size = unit(0.3, "cm"),
      annotation_name_gp = gpar(fontsize = 8.5),
      col = list(
        Age = c("<=65" = "#B9DFFB", ">65" = "#68C84D"),
        N = c("N0" = "#F3F3F4", "N1" = "#ABDAE4", "N2" = "#4B95E9", "N3" = "#123294", "No avalible" = "#747070"),
        M = c("M0" = "#F3F3F4", "M1" = "#F09F37", "No avalible" = "#747070"),
        Survival = c("Alive" = "#F3F3F4", "Death" = "#010101"),
        Gender = c("FEMALE" = "#E93420", "MALE" = "#316DBB"),
        `NMIBC/MIBC` = c("MIBC" = "#AE2417", "NMIBC" = "#F3F3F4", "No avalible" = "#747070"),
        Grade = c("High" = "#4FADEB", "Low" = "#F3F3F4")
      )
    )
)
pdf("oncoprint5.pdf", width = 12.55, height = 4.37)
draw(ht)
dev.off()


# Figure3A-----plot -------------------------------------------------------

EPM <- openxlsx::read.xlsx("~/Project/Bladder/Bladder_circle_numbers.xlsx") %>%
  as_tibble()
EPM <- EPM %>% mutate(paired = stringr::str_extract(sample, "cBca_\\d+"))
EPM$paired <- factor(EPM$paired, levels = unique(stringr::str_sort(EPM$paired, numeric = T)))
NAT <- EPM %>% filter(group == "Normal")
TUM <- EPM %>% filter(group == "Tumour")
NAT <- NAT %>%
  arrange(paired) %>%
  mutate(label = seq(1, 80))
TUM <- TUM %>%
  arrange(paired) %>%
  mutate(label = seq(1, 80))
plotData <- bind_rows(NAT, TUM)
## one line to plot figures
ggplot(plotData, aes(x = label, y = EPM)) +
  geom_line(aes(group = paired),
    color = "#bdbdbd",
    size = 0.2
  ) +
  geom_point(aes(color = group)) +
  geom_smooth(aes(color = group, fill = group),
    linetype = "longdash",
    size = 0.7,
    method = "loess",
    alpha = 0.5
  ) +
  scale_color_manual(values = c(
    "Normal" = "#263272",
    "Tumour" = "#B83C3E"
  )) +
  scale_fill_manual(values = c(
    "Normal" = "#263272",
    "Tumour" = "#B83C3E"
  )) +
  theme_pubr() +
  xlab("Samples") +
  ylab("Number of eccDNAs per \nmillion mapped reads")


# Figure3B-----data -------------------------------------------------------

normals <- c("cBca_1N", "cBca_2N", "cBca_3N", "cBca_4N", "cBca_5N", "cBca_6N", "cBca_7N", "cBca_8N", "cBca_9N", "cBca_10N", "cBca_11N", "cBca_12N", "cBca_13N", "cBca_14N", "cBca_15N", "cBca_16N", "cBca_17N", "cBca_18N", "cBca_19T", "cBca_20T", "cBca_21N", "cBca_23N", "cBca_24N", "cBca_25N", "cBca_26N", "cBca_27N", "cBca_28N", "cBca_29N", "cBca_30N", "cBca_31N", "cBca_32N", "cBca_33N", "cBca_34N", "cBca_35N", "cBca_36N", "cBca_37N", "cBca_38N", "cBca_39T", "cBca_40N", "cBca_41N", "cBca_42N", "cBca_43N", "cBca_44N", "cBca_45N", "cBca_46N", "cBca_47N", "cBca_48N", "cBca_50N", "cBca_51N", "cBca_52N", "cBca_54N", "cBca_55N", "cBca_56N", "cBca_57N", "cBca_58N", "cBca_59N", "cBca_60N", "cBca_62N", "cBca_63N", "cBca_64N", "cBca_65N", "cBca_66N", "cBca_67N", "cBca_68N", "cBca_69N", "cBca_70N", "cBca_71N", "cBca_72N", "cBca_74N", "cBca_75N", "cBca_76N", "cBca_77N", "cBca_78N", "cBca_79N", "cBca_80N", "cBca_84N", "cBca_85N", "cBca_86N", "cBca_87N", "cBca_88N")

cases <- c("cBca_1T", "cBca_2T", "cBca_3T", "cBca_4T", "cBca_5T", "cBca_6T", "cBca_7T", "cBca_8T", "cBca_9T", "cBca_10T", "cBca_11T", "cBca_12T", "cBca_13T", "cBca_14T", "cBca_15T", "cBca_16T", "cBca_17T", "cBca_18T", "cBca_19N", "cBca_20N", "cBca_21T", "cBca_23T", "cBca_24T", "cBca_25T", "cBca_26T", "cBca_27T", "cBca_28T", "cBca_29T", "cBca_30T", "cBca_31T", "cBca_32T", "cBca_33T", "cBca_34T", "cBca_35T", "cBca_36T", "cBca_37T", "cBca_38T", "cBca_39N", "cBca_40T", "cBca_41T", "cBca_42T", "cBca_43T", "cBca_44T", "cBca_45T", "cBca_46T", "cBca_47T", "cBca_48T", "cBca_50T", "cBca_51T", "cBca_52T", "cBca_54T", "cBca_55T", "cBca_56T", "cBca_57T", "cBca_58T", "cBca_59T", "cBca_60T", "cBca_62T", "cBca_63T", "cBca_64T", "cBca_65T", "cBca_66T", "cBca_67T", "cBca_68T", "cBca_69T", "cBca_70T", "cBca_71T", "cBca_72T", "cBca_74T", "cBca_75T", "cBca_76T", "cBca_77T", "cBca_78T", "cBca_79T", "cBca_80T", "cBca_84T", "cBca_85T", "cBca_86T", "cBca_87T", "cBca_88T")
cal_total_eccs <- function(x) {
  read_tsv(stringr::str_glue("~/Project/Bladder/filter_eccs/{x}_circle_site.filter.tsv")) %>%
    nrow()
}

eccnumbers <- tibble(sample = c(cases, normals), eccs = sapply(c(cases, normals), cal_total_eccs))

cal_ingenes <- function(x) {
  read_tsv(stringr::str_glue("~/Project/Bladder/gene_anno/{x}.startAnno.bed"), col_names = F) %>%
    select(X4) %>%
    distinct(X4) %>%
    nrow()
}
genenumbers <- tibble(sample = c(cases, normals), eccs = sapply(c(cases, normals), cal_ingenes))
haha <- full_join(eccnumbers, genenumbers, by = "sample")
haha <- haha %>%
  mutate(geneP = eccs.y / eccs.x * 100)

calrepeat <- function(x) {
  df <- read_tsv(stringr::str_glue("~/Project/Bladder/repeatStat/repStats.{x}.tsv"))
  sum(df$ratio)
}
repnumbers <- tibble(sample = c(cases, normals), eccs = sapply(c(cases, normals), calrepeat))

all <- haha %>% full_join(repnumbers, by = "sample")


haha <- read_tsv("~/Project/Bladder/gene_vs_repeat.tsv")

# Figure3G-----plot -------------------------------------------------------

EPMs <- readRDS("../Bladder/gene_drops.RDS")
newEPM <- EPMs %>%
  tidyr::gather(samples, counts, -gene) %>%
  mutate(gp = if_else(stringr::str_ends(samples, "T"), "Tumor", "Normal"))
newEPM <- newEPM %>%
  tidyr::replace_na(replace = list(counts = 0)) %>%
  mutate(logc = log2(counts + 1)) %>%
  group_by(samples) %>%
  mutate(rank = dense_rank(desc(logc))) %>%
  ungroup()

ggplot(newEPM, aes(x = rank, y = logc)) +
  geom_line(aes(group = samples, color = gp)) +
  scale_color_manual(values = c("Tumor" = "#C8413B", "Normal" = "#364BBA")) +
  xlab("Genes ranks") +
  ylab("eccDNA density") +
  facet_wrap(. ~ gp) +
  theme_pubr() +
  theme(strip.background = element_blank())



# Figure3H-----plot -------------------------------------------------------

### first get abundance matrix
genes <- read_tsv("../qd-ECC4/S/ECC_report/FinallyData/bedFile/dbs/hg38.coding.bed", col_names = F)
names(genes) <- c("chr", "Start", "End", "gene")
genes <- genes %>%
  mutate(length = End - Start)
get_genes <- function(x, dbs = genes) {
  df <- read_tsv(stringr::str_glue("cgene_anno/{x}.startAnno.bed"), col_names = F)
  df <- df %>%
    select(1:4, 8) %>%
    mutate(ecc = paste(paste(X1, X2, sep = ":"), X3, sep = "-")) %>%
    group_by(X8) %>%
    distinct(ecc, .keep_all = T) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    rename(gene = X8) %>%
    left_join(dbs, by = "gene") %>%
    mutate(pct = count / length) %>%
    mutate(pct2 = pct / sum(pct) * (10^6)) %>%
    select(1, 8)
  names(df)[2] <- x
  df
}
hebing <- function(x, y) {
  full_join(x, y, by = "gene")
}
fin <- Reduce(hebing, lapply(c(cases, normals), get_genes))
fin <- fin %>%
  tidyr::gather(sample, value, -gene) %>%
  tidyr::replace_na(replace = list(value = 0)) %>%
  tidyr::pivot_wider(names_from = "sample", values_from = "value")
## calculate logfc
logfc <- fin %>%
  mutate(
    mean1 = rowMeans(across(cBca_1T:cBca_88T)),
    mean2 = rowMeans(across(cBca_1N:cBca_88N))
  ) %>%
  select(gene, mean1, mean2) %>%
  mutate(logfc = log2(mean1) - log2(mean2))
### different analysis
dffs <- fin %>%
  tidyr::gather(sample, TPM, -gene) %>%
  tidyr::replace_na(replace = list(TPM = 0)) %>%
  mutate(gp = if_else(sample %in% cases, "Case", "Normal")) %>%
  group_by(gene) %>%
  rstatix::pairwise_wilcox_test(TPM ~ gp, p.adjust.method = "BH")

dffs <- dffs %>% left_join(logfc, by = "gene")
### prepare data for heatmap
plotData <- fin %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.matrix()
nima <- dffs %>% filter(p.adj < 0.01, abs(logfc) > 0.5)
nima <- nima %>%
  arrange(logfc)
plt <- plotData[nima$gene, ]
plt[is.na(plt)] <- 0
plt <- log2(plt + 1)
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b"))
ht <- Heatmap(plt2,
  show_row_names = F,
  show_column_names = F,
  cluster_rows = F,
  cluster_columns = F,
  # col = col_fun,
  top_annotation = HeatmapAnnotation(
    group = c(rep("Case", 80), rep("Normal", 80)),
    col = list(group = c("Case" = "#e41a1c", "Normal" = "#377eb8"))
  ),
  name = "Log normalized \neccDNA count",
  row_split = c(rep("down regular", 326), rep("up regular", 1019))
)

# Figure3I------data ------------------------------------------------------
# proteinGenes<-read_tsv("~/Project/Bladder/dbs/hg38.coding.bed",col_names = F)%>%
#  pull(X4)
# newdb<-dbs%>%
#  tidyr::separate(X4,into = c("gene","type"),sep = ":")%>%
#  filter(gene%in%proteinGenes)%>%
#  tidyr::unite("X4",gene,type,sep = ":")
dbs <- read_tsv("~/Project/Bladder/dbs/hg38.genetic_elements_exceptCPG_fix.bed", col_names = F)
length_corr <- dbs %>%
  tidyr::separate(X4, into = c("gene", "type"), sep = ":") %>%
  mutate(length = X3 - X2) %>%
  group_by(type) %>%
  summarise(total_l = sum(length))
total_chrom_length <- read_tsv("~/Project/甲状腺癌/old/plasma/dbs/hg38.chromo.size", col_names = F)
total_chrom_length <- sum((total_chrom_length[1:24, ])$X2)
length_corr <- length_corr %>%
  mutate(pct = total_l / total_chrom_length)

cal_element <- function(x) {
  fs <- read_tsv(stringr::str_glue("~/Project/Bladder/elementAnno/{x}.startAnno.bed"), col_names = F)
  tongji <- fs %>%
    select(1, 2, 3, 8) %>%
    mutate(ecc = paste(paste(X1, X2, sep = ":"), X3, sep = "-")) %>%
    group_by(ecc) %>%
    distinct(X8, .keep_all = T) %>%
    tidyr::separate(X8, into = c("gene", "type"), sep = ":") %>%
    ungroup() %>%
    count(type)
  tongji$samples <- x
  tongji
}
fin <- do.call("rbind", lapply(samples, cal_element))
### pie plot data
pie_plot_data <- fin %>%
  group_by(type) %>%
  summarise(counts = sum(n)) %>%
  mutate(ratio = counts / sum(counts))

## correct by element length for barplot

barplot_data <- pie_plot_data %>%
  left_join(length_corr, by = "type") %>%
  mutate(enrichment = ratio / pct)

# Figure3J----plot --------------------------------------------------------


normals <- c("cBca_1N", "cBca_2N", "cBca_3N", "cBca_4N", "cBca_5N", "cBca_6N", "cBca_7N", "cBca_8N", "cBca_9N", "cBca_10N", "cBca_11N", "cBca_12N", "cBca_13N", "cBca_14N", "cBca_15N", "cBca_16N", "cBca_17N", "cBca_18N", "cBca_19T", "cBca_20T", "cBca_21N", "cBca_23N", "cBca_24N", "cBca_25N", "cBca_26N", "cBca_27N", "cBca_28N", "cBca_29N", "cBca_30N", "cBca_31N", "cBca_32N", "cBca_33N", "cBca_34N", "cBca_35N", "cBca_36N", "cBca_37N", "cBca_38N", "cBca_39T", "cBca_40N", "cBca_41N", "cBca_42N", "cBca_43N", "cBca_44N", "cBca_45N", "cBca_46N", "cBca_47N", "cBca_48N", "cBca_50N", "cBca_51N", "cBca_52N", "cBca_54N", "cBca_55N", "cBca_56N", "cBca_57N", "cBca_58N", "cBca_59N", "cBca_60N", "cBca_62N", "cBca_63N", "cBca_64N", "cBca_65N", "cBca_66N", "cBca_67N", "cBca_68N", "cBca_69N", "cBca_70N", "cBca_71N", "cBca_72N", "cBca_74N", "cBca_75N", "cBca_76N", "cBca_77N", "cBca_78N", "cBca_79N", "cBca_80N", "cBca_84N", "cBca_85N", "cBca_86N", "cBca_87N", "cBca_88N")

cases <- c("cBca_1T", "cBca_2T", "cBca_3T", "cBca_4T", "cBca_5T", "cBca_6T", "cBca_7T", "cBca_8T", "cBca_9T", "cBca_10T", "cBca_11T", "cBca_12T", "cBca_13T", "cBca_14T", "cBca_15T", "cBca_16T", "cBca_17T", "cBca_18T", "cBca_19N", "cBca_20N", "cBca_21T", "cBca_23T", "cBca_24T", "cBca_25T", "cBca_26T", "cBca_27T", "cBca_28T", "cBca_29T", "cBca_30T", "cBca_31T", "cBca_32T", "cBca_33T", "cBca_34T", "cBca_35T", "cBca_36T", "cBca_37T", "cBca_38T", "cBca_39N", "cBca_40T", "cBca_41T", "cBca_42T", "cBca_43T", "cBca_44T", "cBca_45T", "cBca_46T", "cBca_47T", "cBca_48T", "cBca_50T", "cBca_51T", "cBca_52T", "cBca_54T", "cBca_55T", "cBca_56T", "cBca_57T", "cBca_58T", "cBca_59T", "cBca_60T", "cBca_62T", "cBca_63T", "cBca_64T", "cBca_65T", "cBca_66T", "cBca_67T", "cBca_68T", "cBca_69T", "cBca_70T", "cBca_71T", "cBca_72T", "cBca_74T", "cBca_75T", "cBca_76T", "cBca_77T", "cBca_78T", "cBca_79T", "cBca_80T", "cBca_84T", "cBca_85T", "cBca_86T", "cBca_87T", "cBca_88T")


getLength <- function(x) {
  df <- read_tsv(stringr::str_glue("~/Project/Bladder/filter_eccs/{x}_circle_site.filter.tsv"))
  tibble(sample = x, length = df$length)
}

fin <- do.call("bind_rows", lapply(c(normals, cases), getLength))

fin <- fin %>%
  mutate(group = if_else(sample %in% cases, "Tumour", "Normal"))

ggplot(plotData, aes(x = length, y = ..count..)) +
  geom_density(aes(color = group)) +
  theme_prism(border = T) +
  xlab("The length distribution of eccDNA") +
  ylab("Count") +
  scale_color_manual(values = c("#3E68B2", "#A5303B"))


# Figure3K-----data -------------------------------------------------------

## catogray
haha2 <- fin %>%
  mutate(gps = cut(length,
    breaks = c(-1, 2000, 10000, Inf),
    labels = c("<2K", "2k~10K", ">10K")
  )) %>%
  group_by(sample, gps) %>%
  summarise(number = n()) %>%
  mutate(ratio = number / sum(number) * 100) %>%
  ungroup()

haha2 <- haha2 %>%
  mutate(group = if_else(sample %in% cases, "Tumour", "Normal"))


# Figure3L----plot --------------------------------------------------------


get_BAF <- function(x) {
  df1 <- read_tsv(stringr::str_glue("BAF/circleseq/AlleleCNV/{x}T_tumourBAF.txt")) %>%
    select(1, 4) %>%
    setNames(c("pos", "Circle-seq"))
  df2 <- read_tsv(stringr::str_glue("BAF/wgs/WGS_BAF/{x}T_tumourBAF.txt")) %>%
    select(1, 4) %>%
    setNames(c("pos", "WGS"))
  its <- df1 %>%
    inner_join(df2, by = "pos") %>%
    filter(WGS != 0) %>%
    filter(WGS != 1) %>%
    tidyr::gather(type, BAF, -pos)
  its
}

samples <- c(
  "1", "2", "3", "4", "5", "6", "7",
  "8", "9", "10", "11", "12", "13",
  "14", "15", "16", "17", "18", "19",
  "20", "21", "23", "24", "25", "26",
  "27", "28", "29", "30", "31", "32",
  "33", "34", "35", "36", "37", "38",
  "39", "40", "41", "42", "43", "44",
  "45", "46", "47", "48", "50",
  "51", "52", "54", "55", "56", "57",
  "58", "59", "60", "62", "63", "64",
  "65", "66", "67", "68", "69", "70",
  "71", "72", "74", "75", "76", "77",
  "78", "79", "80", "84", "85", "86",
  "87", "88"
)

plotData <- do.call("bind_rows", lapply(samples, get_BAF))

plotData <- do.call("bind_rows", lapply(c("11", "12"), get_BAF))

ggplot(plotData, aes(
  x = BAF,
  color = type,
  fill = type,
  y = ..count.. / sum(..count..) / 2
)) +
  geom_histogram(alpha = 0.7, position = "identity") +
  scale_color_manual(values = c("Circle-seq" = "#393A80", "WGS" = "#D1352B")) +
  scale_fill_manual(values = c("Circle-seq" = "#393A80", "WGS" = "#D1352B")) +
  xlab("B-allele frequency") +
  ylab("Fraction of heterozygous\nSNPs") +
  theme_pubr()


# Figure3M ------------data---------------------------------------------------


samples <- c("5N", "5T", "13N", "13T", "15N", "15T", "16N", "16T", "17N", "17T", "21N", "21T", "29N", "29T", "37N", "37T", "50N", "50T")

get_nsegment <- function(x) {
  test <- read_tsv(stringr::str_glue("~/Project/Bladder/longreads/eccDNAs/{x}.info.txt"))
  test <- test %>%
    filter(Nfullpass >= 2) %>%
    distinct(fragments, .keep_all = T)
  test$label <- seq(1, nrow(test))
  haha <- test %>%
    select(label, fragments) %>%
    tidyr::separate_rows(fragments, sep = "\\|") %>%
    mutate(orign = sapply(stringr::str_split(fragments, ":"), function(x) x[[1]])) %>%
    group_by(label) %>%
    summarise(count = n()) %>%
    group_by(count) %>%
    summarise(num = n())
  haha$sample <- x
  haha
}

plotData <- do.call("bind_rows", lapply(samples, get_nsegment))
plotData <- plotData %>%
  tidyr::pivot_wider(
    names_from = count,
    values_from = num, values_fill = 0
  )


# Figure4A-----plot -------------------------------------------------------

EPMs <- readRDS("~/Project/Bladder/gene_drops.RDS")
mat_ecc <- EPMs %>%
  select(1:81) %>%
  setNames(stringr::str_remove(stringr::str_remove(names(EPMs)[1:81], "cBca_"), "T")) %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.data.frame()

CNVs <- read_tsv("~/Project/WGS/CNV/annotate/all_CNV.tsv")
mat_cnv <- CNVs %>%
  setNames(stringr::str_remove(stringr::str_remove(names(CNVs), "Bca_"), "T")) %>%
  tibble::column_to_rownames(var = "Geneid") %>%
  as.data.frame()
##### process chrome location###

chrome_sourceOmics <- c(
  "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
  "12", "13", "14", "15", "16", "17", "18", "19", "20",
  "21", "22", "X", "Y"
)

chrome_targetOmics <- c(
  "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
  "12", "13", "14", "15", "16", "17", "18", "19", "20",
  "21", "22", "X", "Y"
)

####### Extract sub list#########
genelocate <- read_tsv("~/Project/Bladder/dbs/hg38.coding.bed", col_names = F)
genelocate <- genelocate %>%
  select(4, 1, 2, 3) %>%
  mutate(X1 = stringr::str_remove(X1, "chr")) %>%
  setNames(c("Symbol", "chrom", "start", "end")) %>%
  as.data.frame()

genelocate_sourceOmics <- genelocate[genelocate[, 2] %in%
  chrome_sourceOmics, ]
genelocate_targetOmics <- genelocate[genelocate[, 2] %in%
  chrome_targetOmics, ]

intG <- intersect(rownames(targetOmics), genelocate_targetOmics[, 1])

targetOmics <- targetOmics[intG, ]

source_gene <- rownames(sourceOmics)
source_gene_locate <- intersect(unique(genelocate_sourceOmics[, 1]), source_gene)
source_gene <- sourceOmics[source_gene_locate, ]
genelocate_sourceOmics <- genelocate_sourceOmics[genelocate_sourceOmics[, 1] %in% source_gene_locate, ]
genelocate_targetOmics <- genelocate_targetOmics[genelocate_targetOmics[, 1] %in% intG, ]
### Calculate the correlation between cna and other omics data######
corrArray <- calculateCorForTwoMatrices(source_gene, targetOmics, 0.01)
## functions#############
calculateChromLength <- function(chromLength, selectedChrom, genelocate) {
  chromLength <- chromLength[chromLength[, 1] %in% selectedChrom, , drop = FALSE]

  if (length(selectedChrom) == 1) {
    x <- 0
  } else {
    x <- c(0, chromLength[1:(nrow(chromLength) - 1), 2])
  }
  chromLength[, 3] <- cumsum(as.numeric(x))
  chromLength[, 4] <- cumsum(as.numeric(chromLength[, 2]))

  genelocate <- cbind(genelocate, 0, 0)

  colnames(genelocate)[5:6] <- c("finalstart", "finalend")

  for (i in c(1:nrow(genelocate))) {
    chr <- genelocate[i, 2]
    s <- genelocate[i, 3]
    e <- genelocate[i, 4]
    cs <- chromLength[chromLength[, 1] == chr, 3]
    genelocate[i, 5] <- s + cs
    genelocate[i, 6] <- e + cs
  }
  re <- list(chromLength = chromLength, genelocate = genelocate)
  return(re)
}
plotHeatMap <- function(corrArray, genelocate_sourceOmics,
                        chromLength_sourceOmics, genelocate_targetOmics, chromLength_targetOmics,
                        sourceOmicsName, targetOmicsName, dim = 1) {
  allChromlen_sourceOmics <-
    chromLength_sourceOmics[nrow(chromLength_sourceOmics), 4]
  allChromlen_targetOmics <-
    chromLength_targetOmics[nrow(chromLength_targetOmics), 4]

  if (dim == 1) {
    par(mar = c(4, 4, 4, 0))
  } else {
    par(mar = c(0, 4, 4, 0))
  }

  p <- which(corrArray != 0, arr.ind = TRUE)
  allcnagene <- rownames(corrArray)
  allovgene <- colnames(corrArray)

  la <- 1
  for (i in c(1:nrow(p))) {
    cnag <- allcnagene[p[i, 1]]
    ovg <- allovgene[p[i, 2]]
    cnagp <- genelocate_sourceOmics[genelocate_sourceOmics[, 1] == cnag, 5]
    ovgp <- genelocate_targetOmics[genelocate_targetOmics[, 1] == ovg, 5]

    if (length(cnagp) == 0 || length(ovgp) == 0) {
      next
    }

    cov <- corrArray[cnag, ovg]
    color <- ifelse(cov > 0, "#E63126", "#0932E3")

    if (la == 1) {
      if (dim == 1) {
        plot(cnagp, ovgp,
          main = paste(sourceOmicsName, "-", targetOmicsName,
            " correlation",
            sep = ""
          ), xlim = c(0, allChromlen_sourceOmics),
          ylim = c(0, allChromlen_targetOmics), xaxt = "n", yaxt = "n",
          frame.plot = FALSE, xlab = paste(sourceOmicsName,
            " chromosomal location",
            sep = ""
          ), ylab = paste(targetOmicsName,
            " chromosomal location",
            sep = ""
          ), pch = 20, col = color, cex = 0.2
        )
        axis(
          side = 1, at = (chromLength_sourceOmics[, 4] -
            chromLength_sourceOmics[, 2] / 2),
          labels = chromLength_sourceOmics[, 1]
        )
      } else {
        plot(cnagp, ovgp,
          main = paste(sourceOmicsName, "-",
            targetOmicsName, " correlation",
            sep = ""
          ),
          xlim = c(0, allChromlen_sourceOmics),
          ylim = c(0, allChromlen_targetOmics), xaxt = "n", yaxt = "n",
          frame.plot = FALSE, ylab = paste(targetOmicsName, " chromosomal
location", sep = ""),
          xlab = "", pch = 20, col = color, cex = 0.2
        )
      }
      axis(
        side = 2, at = (chromLength_targetOmics[, 4] -
          chromLength_targetOmics[, 2] / 2),
        labels = chromLength_targetOmics[, 1]
      )

      # abline(h=c(0,chromLength_targetOmics[,4]),v=c(0,chromLength_sourceOmics[,4]),
      #       col="gray",lty=3)
      la <- la + 1
    } else {
      for (u in seq_len(length(cnagp))) {
        for (v in seq_len(length(ovgp))) {
          points(cnagp[u], ovgp[v], pch = 20, col = color, cex = 0.2)
        }
      }
    }
  }
}

#####################
## Calculate the location of genes in the heatmap
chromLength <- read_tsv("../../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size", col_names = F)
chromLength <- chromLength %>%
  filter(X1 %in% c(paste0("chr", seq(1, 22)), "chrX", "chrY")) %>%
  mutate(X1 = stringr::str_remove(X1, "chr")) %>%
  setNames(c("V1", "V2"))
chromLength$V1 <- factor(chromLength$V1, levels = c(seq(1, 22), "X", "Y"))
chromLength <- chromLength %>%
  arrange(V1) %>%
  as.data.frame()

re <- calculateChromLength(chromLength, chrome_sourceOmics, genelocate_sourceOmics)
genelocate_sourceOmics <- re$genelocate
chromLength_sourceOmics <- re$chromLength
re <- calculateChromLength(chromLength, chrome_targetOmics, genelocate_targetOmics)
genelocate_targetOmics <- re$genelocate
chromLength_targetOmics <- re$chromLength
### plot
plotHeatMap(corrArray, genelocate_sourceOmics, chromLength_sourceOmics,
  genelocate_targetOmics, chromLength_targetOmics, "eccDNA",
  "mRNA",
  dim = 1
)


# Figure4B-----plot -------------------------------------------------------


chromSize <- read_tsv("../qd-ECC4/S/ECC_report/DataBase/hg38.chromo.size", col_names = F)
total_chrom_size <- chromSize %>%
  filter(X1 %in% c(paste0("chr", seq(1, 22)), "chrX", "chrY")) %>%
  summarise(total = sum(X2)) %>%
  pull(total)


get_stat <- function(x, total_chrom_size = 3088269832) {
  seed_region <- read_tsv(stringr::str_glue("~/Project/Bladder/seeds/{x}_AA_CNV_SEEDS.bed"), col_names = F)
  seed_length <- seed_region %>%
    mutate(length = X3 - X2) %>%
    summarise(total = sum(length)) %>%
    pull(total)
  inters <- read_tsv(stringr::str_glue("~/Project/Bladder/seeds/{x}.inters.bed"), col_names = F)
  eccs <- read_tsv(stringr::str_glue("~/Project/Bladder/bedFile/cBca_{x}T.ecc.bed"), col_names = F)
  if (nrow(inters) != nrow(eccs)) {
    cat(x, "may be have some problems!")
  }
  inters <- inters %>%
    filter(X1 %in% c(paste0("chr", seq(1, 22)), "chrX", "chrY")) %>%
    group_by(X5) %>%
    summarise(count = n()) %>%
    filter(X5 < 2) %>%
    arrange(X5) %>%
    mutate(length = c(total_chrom_size - seed_length, seed_length)) %>%
    mutate(
      count2 = count / length,
      count3 = count2 / sum(count2)
    ) %>%
    select(X5, count3) %>%
    setNames(c("type", "counts"))
  inters$sample <- x
  inters
}

samples <- c("1", "2", "5", "8", "9", "13", "14", "15", "16", "17", "21", "23", "24", "25", "28", "29", "30", "31", "32", "33", "34", "36", "37", "39", "40", "41", "42", "44", "46", "47", "48", "50", "51", "54", "56", "57", "60", "62", "63", "65", "67", "68", "69", "70", "71", "72", "74", "75", "76", "77", "79", "80", "84", "85", "86", "87", "88")

fin <- do.call("bind_rows", lapply(samples, function(x) get_stat(x)))

plotData <- fin %>%
  mutate(type = if_else(type == 0, "non-amp", "amp"))

# Figure4C-----plot -------------------------------------------------------


df <- read_tsv("PABPC1_geneRegion_readcoverage3.txt")
plotData <- df %>%
  mutate(coord = c(seq(100680816, 100727809, by = 50)[2:940], 100727809)) %>%
  tidyr::gather(sample, value, -coord) %>%
  mutate(group = if_else(sample %in% c("5T", "16T", "17T", "21T", "50T", "84T"), "PABPC1-amplified", "non-PABPC1-amplified"))

plotData$value2 <- log2(plotData$value + 1)

plotData3 <- plotData %>%
  group_by(coord, group) %>%
  summarise(value3 = list(mean_ci(value2))) %>%
  tidyr::unnest() %>%
  ungroup()

ggplot(plotData3, aes(x = coord, y = y, group = group)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#D3D2D3") +
  geom_line(aes(color = group)) +
  geom_vline(xintercept = 100685816, color = "grey") +
  geom_vline(xintercept = 100722809, color = "grey") +
  scale_color_manual(values = c("#34327F", "#CF3430")) +
  xlab("Genomic range") +
  ylab("Mean read coverage(log2)") +
  theme_pubr()


# Figure4E----plot --------------------------------------------------------

TPMs <- read_tsv("RNAseq/Protein_coding.tpms.txt")
## should filter
goodGenes <- TPMs %>%
  tidyr::gather(sample, TPM, -Geneid) %>%
  mutate(type = if_else(TPM > 5, "good", "bad")) %>%
  group_by(Geneid, type) %>%
  summarise(count = n()) %>%
  tidyr::pivot_wider(names_from = "type", values_from = count, values_fill = 0) %>%
  filter(good == 126) %>%
  pull(Geneid)
TPMs <- TPMs %>%
  filter(Geneid %in% goodGenes)
TPM2 <- TPMs %>%
  tidyr::gather(sample, TPM, -Geneid)
TPM3 <- TPM2 %>% filter(stringr::str_ends(sample, "T"))


## mean+sd
stst <- TPM3 %>%
  group_by(Geneid) %>%
  summarise(mean = mean(TPM), sd = sd(TPM))
## zscore
zscores <- TPM3 %>%
  left_join(stst, by = "Geneid") %>%
  mutate(zscore = (TPM - mean) / sd)

zscores <- zscores %>%
  dplyr::select(Geneid, sample, zscore) %>%
  dplyr::rename(TPM = zscore) %>%
  mutate(sample = stringr::str_extract(sample, "\\d+"))

eccgenes <- read_tsv("~/Project/Bladder/allgene20.bed", col_names = F)
eccgenes <- eccgenes %>%
  mutate(sample = stringr::str_extract(X4, "cBca_\\d+[NT]")) %>%
  dplyr::select(sample, X8) %>%
  distinct(sample, X8, .keep_all = T)
names(eccgenes)[2] <- "Geneid"

part1 <- eccgenes %>%
  filter(sample %in% c("cBca_19N", "cBca_20N", "cBca_39N"))
part2 <- eccgenes %>%
  filter(stringr::str_ends(sample, "T")) %>%
  filter(!(sample %in% c("cBca_19T", "cBca_20T", "cBca_39T")))
eccgenes <- bind_rows(part1, part2)

eccgenes <- eccgenes %>%
  mutate(sample = stringr::str_extract(sample, "\\d+")) %>%
  left_join(zscores, by = c("sample", "Geneid"))

ecgenes <- openxlsx::read.xlsx("all_detective_gene_list.xlsx") %>%
  as_tibble() %>%
  filter(stringr::str_starts(feature, "ecDNA")) %>%
  dplyr::select(sample_name, gene) %>%
  distinct(sample_name, gene, .keep_all = T) %>%
  dplyr::rename(sample = sample_name, Geneid = gene)
ecgenes$sample <- as.character(ecgenes$sample)
ecgenes <- ecgenes %>%
  left_join(zscores, by = c("sample", "Geneid")) %>%
  na.omit()
ecgenes <- ecgenes %>%
  mutate(ding = paste(sample, Geneid, sep = "-"))
eccgenes <- eccgenes %>%
  mutate(ding = paste(sample, Geneid, sep = "-"))
### use all focal-amp genes to filter
allGenes <- openxlsx::read.xlsx("all_detective_gene_list.xlsx") %>%
  as_tibble() %>%
  dplyr::select(sample_name, gene) %>%
  distinct(sample_name, gene, .keep_all = T) %>%
  dplyr::rename(sample = sample_name, Geneid = gene) %>%
  mutate(ding = paste(sample, Geneid, sep = "-"))
jiaoji <- intersect(allGenes$ding, eccgenes$ding)
eccgenes <- eccgenes %>%
  filter(!ding %in% jiaoji)
eccgenes$type <- "eccDNA"
ecgenes$type <- "ecDNA"

fin <- bind_rows(ecgenes, eccgenes)

fin <- fin %>%
  na.omit()

ggplot(fin, aes(x = TPM, color = type)) +
  geom_density() +
  scale_color_manual(values = c("#9E3735", "#48436D")) +
  theme_pubr() +
  xlab("RNA expression(Z-score)")

# Figure4F----plot --------------------------------------------------------

ecgenes <- ecgenes %>% filter(sample == "16")

plotData <- TPM3 %>%
  mutate(sample = stringr::str_extract(sample, "\\d+")) %>%
  filter(sample == "16") %>%
  mutate(rank = dense_rank(desc(TPM))) %>%
  mutate(ecgene = if_else(Geneid %in% c("PABPC1", "BRK1", "YWHAZ", "TADA3", "CTNNB1", "OXR1", "UBR5", "ANKRD46"), Geneid, "")) %>%
  mutate(col = if_else(ecgene == "", "no", "yes"))


ggplot(plotData, aes(x = rank, y = TPM)) +
  geom_line() +
  geom_point(aes(color = col)) +
  scale_y_log10() +
  scale_color_manual(values = list(yes = "red", no = "transparent")) +
  geom_text_repel(aes(label = ecgene), max.overlaps = 100000000000) +
  theme(legend.position = "none") +
  xlab("Genes ranked by expression level\n(TPM>5)") +
  ylab("TPM") +
  theme_pubr()



# Figure4G-----plot -------------------------------------------------------

####################### WGS##################################

df <- read_tsv("MYEOV/MYEOV.WGS.allele.txt")
df <- df %>%
  filter(Good_depth > 10)
wgs_alle <- df %>%
  tidyr::gather(alleleType, count, -`#CHR`, -POS, -Good_depth) %>%
  arrange(`#CHR`, POS, Good_depth, desc(count)) %>%
  group_by(`#CHR`, POS, Good_depth) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  mutate(tp = stringr::str_sub(alleleType, start = 7, end = 7)) %>%
  mutate(AF = count / Good_depth) %>%
  mutate(tp2 = rep(c("major", "minor"), times = 242))

data <- wgs_alle %>%
  select(`#CHR`, POS, AF, tp2) %>%
  tidyr::pivot_wider(names_from = "tp2", values_from = "AF") %>%
  arrange(POS) %>%
  as.data.frame()
gtrack <- GenomeAxisTrack(range = IRanges(
  start = 69253332,
  end = 70626244,
))

gr <- GRanges(
  seqnames = "chr8", strand = "*",
  ranges = IRanges(start = data$POS, width = 1),
  major = data[, 3], minor = data[, 4]
)

dtrack <- DataTrack(gr,
  name = "WGS AF",
  groups = c("major", "minor"),
  col = c("#c5362c", "#4475a7"),
  type = "p",
  yTicksAt = c(0, 0.5, 1),
  ylim = c(0, 1), grid = T, lty.grid = "dashed", lwd.grid = 0.5, h = 3
)

wgsCOV <- read_tsv("PABPC1/WGS.cov.bdg", col_names = F)

cr <- GRanges(
  seqnames = "chr8", strand = "*",
  ranges = IRanges(start = wgsCOV$X2, end = wgsCOV$X3),
  count = wgsCOV$X4
)
covtrack <- DataTrack(cr, name = "WGS cov", col = "grey", type = "histogram")

################################################################################
############################### RNAseq##########################################

df2 <- read_tsv("PABPC1/PABPC1.RNA.allele.txt")
df2 <- df2 %>%
  filter(Good_depth > 15) %>%
  mutate(
    Count_A = if_else(Count_A > 2, Count_A, 0),
    Count_C = if_else(Count_C > 2, Count_C, 0),
    Count_G = if_else(Count_G > 2, Count_G, 0),
    Count_T = if_else(Count_T > 2, Count_T, 0)
  ) %>%
  mutate(Good_depth = Count_A + Count_C + Count_G + Count_T)
RNA_alle <- df2 %>%
  tidyr::gather(alleleType, count, -`#CHR`, -POS, -Good_depth) %>%
  mutate(
    tp = stringr::str_sub(alleleType, start = 7, end = 7),
    AF = count / Good_depth
  ) %>%
  select(`#CHR`, POS, tp, AF) %>%
  dplyr::rename(RNAAF = AF)

data2 <- wgs_alle %>%
  inner_join(RNA_alle, by = c("#CHR", "POS", "tp")) %>%
  select(`#CHR`, POS, RNAAF, tp2) %>%
  tidyr::pivot_wider(names_from = "tp2", values_from = "RNAAF") %>%
  arrange(POS) %>%
  as.data.frame()


gr2 <- GRanges(
  seqnames = "chr8", strand = "*",
  ranges = IRanges(start = data2$POS, width = 1),
  major = data2[, 3], minor = data2[, 4]
)

dtrack2 <- DataTrack(gr2,
  name = "RNAcount",
  groups = c("major", "minor"),
  col = c("#c5362c", "#4475a7"),
  yTicksAt = c(0, 0.5, 1),
  ylim = c(0, 1), grid = T, lty.grid = "dashed", lwd.grid = 0.5, h = 3
)

rnaCOV <- read_tsv("PABPC1/RNA.cov.bdg", col_names = F)

cr2 <- GRanges(
  seqnames = "chr8", strand = "*",
  ranges = IRanges(start = rnaCOV$X2, end = rnaCOV$X3),
  count = rnaCOV$X4
)
covtrack2 <- DataTrack(cr2, name = "RNA cov", col = "grey", type = "histogram", ylim = c(0, 500))
##############################################################################
##################################### CIRCLE###################################
df3 <- read_tsv("PABPC1/PABPC1.CIRCLE.allele.txt")
df3 <- df3 %>%
  filter(Good_depth > 15) %>%
  mutate(
    Count_A = if_else(Count_A > 2, Count_A, 0),
    Count_C = if_else(Count_C > 2, Count_C, 0),
    Count_G = if_else(Count_G > 2, Count_G, 0),
    Count_T = if_else(Count_T > 2, Count_T, 0)
  ) %>%
  mutate(Good_depth = Count_A + Count_C + Count_G + Count_T)

circle_alle <- df3 %>%
  tidyr::gather(alleleType, count, -`#CHR`, -POS, -Good_depth) %>%
  mutate(
    tp = stringr::str_sub(alleleType, start = 7, end = 7),
    AF = count / Good_depth
  ) %>%
  select(`#CHR`, POS, tp, AF) %>%
  dplyr::rename(CIRCLEAF = AF)

data3 <- wgs_alle %>%
  inner_join(circle_alle, by = c("#CHR", "POS", "tp")) %>%
  select(`#CHR`, POS, CIRCLEAF, tp2) %>%
  tidyr::pivot_wider(names_from = "tp2", values_from = "CIRCLEAF") %>%
  arrange(POS) %>%
  as.data.frame()

gr3 <- GRanges(
  seqnames = "chr8", strand = "*",
  ranges = IRanges(start = data3$POS, width = 1),
  major = data3[, 3], minor = data3[, 4]
)

dtrack3 <- DataTrack(gr3,
  name = "CIRCLEcount", groups = c("major", "minor"), col = c("#c5362c", "#4475a7"), yTicksAt = c(0, 0.5, 1),
  ylim = c(0, 1), grid = T, lty.grid = "dashed", lwd.grid = 0.5, h = 3
)

circleCOV <- read_tsv("PABPC1/CIRCLE.cov.bdg", col_names = F)

cr3 <- GRanges(
  seqnames = "chr8", strand = "*",
  ranges = IRanges(start = circleCOV$X2, end = circleCOV$X3),
  count = circleCOV$X4
)
covtrack3 <- DataTrack(cr3, name = "CIRCLE cov", col = "grey", type = "histogram", ylim = c(0, 1000))

png("test_region.png", width = 6.24, height = 4.86, units = "in", res = 300)
plotTracks(list(covtrack, dtrack, covtrack3, dtrack3, covtrack2, dtrack2, gtrack),
  from = 99497918,
  to = 119476539,
  sizes = c(0.5, 1, 0.5, 1, 0.5, 1, 0.5), legend = F,
  background.title = "white", col.title = "black", col.axis = "black"
)
dev.off()


# Figure4H,4J,4K----- data-----------------------------------------------------

AA_derived_CNV <- openxlsx::read.xlsx("all_detective_gene_list.xlsx")
AA_derived_CNV$sample_name <- as.character(AA_derived_CNV$sample_name)

TPMs <- read_tsv("RNAseq/Protein_coding.tpms.txt")
allRNA_samples <- stringr::str_remove(stringr::str_remove(names(TPMs)[2:71], "Bca_"), "_T")

## should filter
goodGenes <- TPMs %>%
  tidyr::gather(sample, TPM, -Geneid) %>%
  mutate(type = if_else(TPM > 5, "good", "bad")) %>%
  group_by(Geneid, type) %>%
  summarise(count = n()) %>%
  tidyr::pivot_wider(names_from = "type", values_from = count, values_fill = 0) %>%
  filter(good == 126) %>%
  pull(Geneid)
TPMs <- TPMs %>%
  filter(Geneid %in% goodGenes)
TPM2 <- TPMs %>%
  tidyr::gather(sample, TPM, -Geneid)
TPM3 <- TPM2 %>% filter(stringr::str_ends(sample, "T"))
TPM_td <- TPM3 %>%
  dplyr::rename(sample_name = sample, gene = Geneid)
TPM_td$sample_name <- stringr::str_remove(stringr::str_remove(TPM_td$sample_name, "Bca_"), "_T")

higher_AA <- AA_derived_CNV %>%
  filter(sample_name %in% allRNA_samples)

gene_in_sample <- higher_AA %>%
  group_by(gene) %>%
  summarise(nohave = paste0(setdiff(allRNA_samples, sample_name), collapse = ","))

higher_AA <- higher_AA %>%
  left_join(TPM_td, by = c("sample_name", "gene"))
higher_AA <- na.omit(higher_AA)

ecc_amp <- higher_AA %>%
  filter(stringr::str_detect(feature, "ecDNA"))
non_ecc_amp <- higher_AA %>%
  filter(!stringr::str_detect(feature, "ecDNA"))

get_fold <- function(x) {
  noneed <- higher_AA %>%
    filter(gene == x) %>%
    pull(sample_name)
  tmp <- TPM_td %>%
    filter(gene == x) %>%
    filter(!(sample_name %in% noneed)) %>%
    summarise(noamp = mean(TPM))
  tmp$noamp
}

ecc_amp <- ecc_amp %>%
  mutate(noampTPM = purrr::map_dbl(gene, function(x) {
    get_fold(x)
  }))
non_ecc_amp <- non_ecc_amp %>%
  mutate(noampTPM = purrr::map_dbl(gene, function(x) {
    get_fold(x)
  }))
ecc_amp <- ecc_amp %>%
  mutate(FC = (TPM + 1) / (noampTPM + 1))
non_ecc_amp <- non_ecc_amp %>%
  mutate(FC = (TPM + 1) / (noampTPM + 1))

ecc_amp$type <- "ecDNA"
non_ecc_amp$type <- "others"
###### for all protein coding
plotData <- bind_rows(ecc_amp, non_ecc_amp)
###### for oncogenes
oncogenes <- openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx")
hao <- oncogenes %>%
  pull(Gene) %>%
  unique()
plotData2 <- plotData %>% filter(gene %in% hao)


# Figure4L ------data --------------------------------------------------


fin <- readRDS("~/Project/Bladder/gene_drops.RDS")
fin <- fin %>%
  dplyr::select(1:81) %>%
  tidyr::gather(sample, value, -gene) %>%
  dplyr::rename(Geneid = gene, EPM = value) %>%
  tidyr::replace_na(replace = list(EPM = 0))

TPMs <- read_tsv("RNAseq/Protein_coding.tpms.txt")
goodGenes <- TPMs %>%
  tidyr::gather(sample, TPM, -Geneid) %>%
  mutate(type = if_else(TPM > 5, "good", "bad")) %>%
  group_by(Geneid, type) %>%
  summarise(count = n()) %>%
  tidyr::pivot_wider(names_from = "type", values_from = count, values_fill = 0) %>%
  filter(good == 126) %>%
  pull(Geneid)

TPMs <- TPMs %>%
  filter(Geneid %in% goodGenes) %>%
  dplyr::select(1:71) %>%
  tidyr::gather(sample, TPM, -Geneid)

fin$sample <- stringr::str_remove(stringr::str_remove(fin$sample, "cBca_"), "T")
TPMs$sample <- stringr::str_remove(stringr::str_remove(TPMs$sample, "Bca_"), "_T")
fin <- inner_join(fin, TPMs, by = c("Geneid", "sample"))
rst <- fin %>%
  mutate(
    logT = log2(TPM + 1),
    logE = log2(EPM + 1)
  ) %>%
  group_by(Geneid) %>%
  rstatix::cor_test(vars = "logE", vars2 = "logT", method = "spearman")
# Figure4M-----plot -------------------------------------------------------

need_genes1 <- openxlsx::read.xlsx("Cancer_driver genes_Bladder cancer.xlsx", sheet = 1)
need_genes1 <- need_genes1 %>%
  as_tibble() %>%
  filter(stringr::str_detect(Role, "oncogene"))

tmp <- rst %>%
  mutate(rank = dense_rank(desc(cor))) %>%
  mutate(cosmic = if_else(Geneid %in% (need_genes1$Gene), "cosmic", "non-cosmic")) %>%
  mutate(tp = if_else(cosmic == "cosmic", "cosmic", tp))
tmp$tp <- factor(tmp$tp, levels = c("no", "yes", "cosmic"))

lizi <- tmp %>%
  filter(tp == "cosmic") %>%
  sample_frac(size = 0.02) %>%
  arrange(rank) %>%
  pull(Geneid)

tmp <- tmp %>%
  mutate(label = if_else(Geneid %in% c("KDM5A", lizi, "PRIM2"), Geneid, ""))

ggplot(tmp, aes(x = rank, y = cor)) +
  geom_point(size = 0.5) +
  geom_text_repel(aes(label = label), size = 3.4, max.overlaps = 10000000000) +
  facet_wrap(~tp) +
  theme_pubr()
ggsave("cor.generank.pdf", width = 12.24, height = 3.73)


# Figure6A,6B,6C-----data -------------------------------------------------

## exps
TPMs <- read_tsv("RNAseq/Protein_coding.tpms.txt")
TPM2 <- TPMs %>%
  tidyr::gather(sample, TPM, -Geneid)
TPM3 <- TPM2 %>% filter(stringr::str_ends(sample, "T"))
TPM3$sample <- stringr::str_remove(TPM3$sample, "_T")

## group
groupInfo <- read_tsv("class_sample.tsv")
names(groupInfo) <- c("sample", "type")
groupInfo <- groupInfo %>%
  mutate(type = if_else(type == "ecDNA", "ecDNA", "others"))
TPM3 <- TPM3 %>%
  left_join(groupInfo, by = "sample")

EPMS <- openxlsx::read.xlsx("~/Project/Bladder/Bladder_circle_numbers.xlsx") %>%
  as_tibble()
EPMs <- EPMS %>%
  filter(group == "Tumour") %>%
  dplyr::select(sample, EPM)
EPMs$sample <- stringr::str_remove(EPMs$sample, "c")
EPMs$sample <- stringr::str_remove(EPMs$sample, "T")
EPMs$sample <- stringr::str_remove(EPMs$sample, "N")
TPM3 <- TPM3 %>%
  left_join(EPMs, by = "sample")

TPM3 <- TPM3 %>%
  mutate(logT = log2(TPM + 1), logE = log2(EPM + 1))

goodGenes <- TPMs %>%
  tidyr::gather(sample, TPM, -Geneid) %>%
  filter(stringr::str_ends(sample, "T")) %>%
  mutate(type = if_else(TPM > 1, "good", "bad")) %>%
  group_by(Geneid, type) %>%
  summarise(count = n()) %>%
  tidyr::pivot_wider(names_from = "type", values_from = count, values_fill = 0) %>%
  filter(good == 70) %>%
  pull(Geneid)

diffGenes <- TPM3 %>%
  dplyr::filter(Geneid %in% goodGenes) %>%
  dplyr::group_by(Geneid) %>%
  rstatix::wilcox_test(logT ~ type)

diffCors2 <- TPM3 %>%
  dplyr::group_by(Geneid) %>%
  rstatix::cor_test(vars = "logT", vars2 = "logE", method = "spearman")

diffGene1 <- diffGenes %>%
  filter(p < 0.05)
diffCor1 <- diffCors %>%
  filter(p < 0.05)

jiaoji2 <- intersect(diffGene1$Geneid, diffCor1$Geneid)

zhuanhua <- bitr(jiaoji, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(zhuanhua2$ENTREZID, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, ont = "BP", readable = T)

eko <- enrichKEGG(zhuanhua2$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

# Figure6D----plot --------------------------------------------------------

genes <- c("LIG3", "LIG4", "POLM", "POLQ", "PRKDC", "BRCA1", "BRCA2", "MSH3")
exp_plot <- TPM3 %>%
  filter(Geneid %in% genes)

p1 <- ggplot(exp_plot %>% filter(Geneid == genes[1]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[1]} expression value")) +
  theme(legend.position = "none")

p2 <- ggplot(exp_plot %>% filter(Geneid == genes[2]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[2]} expression value")) +
  theme(legend.position = "none")

p3 <- ggplot(exp_plot %>% filter(Geneid == genes[3]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[3]} expression value")) +
  theme(legend.position = "none")

p4 <- ggplot(exp_plot %>% filter(Geneid == genes[4]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[4]} expression value")) +
  theme(legend.position = "none")

p5 <- ggplot(exp_plot %>% filter(Geneid == genes[5]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[5]} expression value")) +
  theme(legend.position = "none")

p6 <- ggplot(exp_plot %>% filter(Geneid == genes[6]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[6]} expression value")) +
  theme(legend.position = "none")

p7 <- ggplot(exp_plot %>% filter(Geneid == genes[7]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[7]} expression value")) +
  theme(legend.position = "none")

p8 <- ggplot(exp_plot %>% filter(Geneid == genes[8]), aes(x = type, y = logT)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  scale_fill_manual(values = c("others" = "#FBF3AA", "ecDNA" = "#F3B854")) +
  xlab("") +
  ylab(stringr::str_glue("{genes[8]} expression value")) +
  theme(legend.position = "none")

ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, nrow = 2, align = "v")

ggsave("somegenes.plot1.pdf", width = 6.85, height = 5.3)

# Figure6E-----plot -------------------------------------------------------


p1 <- ggplot(exp_plot %>% filter(Geneid == genes[1]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[1]} eccDNA value"))

p2 <- ggplot(exp_plot %>% filter(Geneid == genes[2]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[2]} eccDNA value"))

p3 <- ggplot(exp_plot %>% filter(Geneid == genes[3]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[3]} eccDNA value"))

p4 <- ggplot(exp_plot %>% filter(Geneid == genes[4]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[4]} eccDNA value"))

p5 <- ggplot(exp_plot %>% filter(Geneid == genes[5]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[5]} eccDNA value"))

p6 <- ggplot(exp_plot %>% filter(Geneid == genes[6]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[6]} eccDNA value"))

p7 <- ggplot(exp_plot %>% filter(Geneid == genes[7]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[7]} eccDNA value"))

p8 <- ggplot(exp_plot %>% filter(Geneid == genes[8]), aes(x = logT, y = logE)) +
  geom_point(color = "#F3B854", size = 2.1) +
  stat_smooth(geom = "line", method = "lm") +
  xlab("Expression value") +
  ylab(stringr::str_glue("{genes[8]} eccDNA value"))

ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, nrow = 2, align = "v")

# Figure7F -------------plot---------------------------------------------------


########################## diffGene_with_eccDNA(higher and lower)#####################################################
counts <- read_tsv("RNAseq/allSample.count.txt")
counts <- counts %>%
  tibble::column_to_rownames(var = "Geneid") %>%
  as.data.frame()
counts <- counts[, 1:70]
groupInfo <- read_tsv("~/Project/Bladder/EPM and clinical variables.tsv")
groupInfo <- groupInfo %>%
  select(1, 4)
colData <- tibble(Sample = colnames(counts))
colData <- colData %>% left_join(groupInfo, by = "Sample")
colData <- colData %>%
  mutate(group = if_else(`EPM Group` == "Low", "A", "B")) %>%
  select(Sample, group) %>%
  as.data.frame()
colData$group <- factor(colData$group, levels = c("A", "B"))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~group)
dds <- dds[rowSums(counts(dds)) > 5, ]
dds <- DESeq(dds)
res <- results(dds)
allresult <- as.data.frame(res) %>%
  tibble::rownames_to_column(var = "Geneid") %>%
  as_tibble()

allresult <- allresult %>%
  mutate(logp = -log10(padj))

allresult <- allresult %>%
  mutate(col = case_when(
    (log2FoldChange < (-1)) & (logp > 2) ~ "downregular",
    (log2FoldChange > 1) & (logp > 2) ~ "upregular",
    TRUE ~ "none"
  ))

geneNames <- bitr(allresult$Geneid, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
names(geneNames) <- c("Geneid", "id")
allresult <- allresult %>% inner_join(geneNames, by = "Geneid")
geneList <- allresult$log2FoldChange
names(geneList) <- allresult$id
geneList <- sort(geneList, decreasing = TRUE)

## download GSEA db
sub1 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
sub2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
sub3 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
sub4 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
C2_t2g <- rbind(sub1, sub2, sub3, sub4) %>%
  dplyr::select(gs_name, entrez_gene)
C2_t2n <- rbind(sub1, sub2, sub3, sub4) %>%
  dplyr::select(gs_name, gs_description)

em2 <- GSEA(geneList, TERM2GENE = C2_t2g, TERM2NAME = C2_t2n, minGSSize = 20, eps = 0)

rankvalues <- tibble(geneid = names(rev(geneList)), foldchange = rev(geneList), order = 1:18680)
needID <- c(6, 9, 11, 20, 22, 26, 30, 33, 34, 38, 43, 58)
getPlotData <- function(i) {
  tibble(geneid = em2@geneSets[[em2@result[needID[i], ]$ID]]) %>%
    left_join(rankvalues, by = "geneid") %>%
    mutate(y1 = i - 1, y2 = i - 1 + 0.9, col = if_else(foldchange > 0, "yes", "no")) %>%
    tidyr::drop_na()
}


haha <- do.call("rbind", lapply(1:12, getPlotData))
ggplot(haha, aes(x = order, y = y1)) +
  geom_segment(aes(xend = order, yend = y2, color = col)) +
  scale_color_manual(values = c("no" = "#5B799D", "yes" = "#CC726A")) +
  theme_void()

ggsave("GSEA_result.part3.png", width = 7.04, height = 2.37, dpi = 300)

rankvalues <- rankvalues %>%
  mutate(col = if_else(foldchange > 0, "yes", "no"))

ggplot(rankvalues, aes(x = order, y = foldchange)) +
  geom_area(aes(fill = col)) +
  scale_fill_manual(values = c("no" = "#5B799D", "yes" = "#CC726A")) +
  theme_pubr()
ggsave("GSEA_result.part2.pdf", width = 7.04, height = 4.14)


ggplot(rankvalues, aes(x = foldchange)) +
  geom_histogram(fill = "#5B799D", color = "black", bins = 100, size = 0.2) +
  scale_x_continuous(limits = c(-2, 2)) +
  theme(legend.position = "none") +
  theme_pubr() +
  ylab("Gene Density") +
  xlab("Correlation")
ggsave("GSEA_result.part1.pdf", width = 9.44, height = 4.14)
