source ("01.function.R")

# ref: reference allele, nef: non-reference allele
options(width = 120)
fdr_cutoff <- 0.05          # FDR cut-off 5%
cutoff.mat <- 0.70
cutoff.pat <- 0.30
depth <- 20
cv <- c("NIP", "KAS")
head <-
  c(
    "locus_name",
    "mat_Refnef_1",
    "pat_Refnef_1",
    "mat_Refnef_2",
    "pat_Refnef_2",
    "pat_Nefref_1",
    "mat_Nefref_1",
    "pat_Nefref_2",
    "mat_Nefref_2"
  )

  

######################
# load and prepare data
######################
# ref
anno <- read.table("../annotation/ALL_anno_class.txt", header = T)
seed <- read.table("./R/list/seed_cort_gene.ratio.txt", header = T)

#anno <- subset(anno, class == "coding")
anno <- anno[grep("LOC_", anno$locus_name), ]


# import count data
RefNef1 <- read.table("./reads_count2/NK_count.txt", header = T)
RefNef1 <- RefNef1[, c(-2:-6)]

RefNef2 <- read.table("../imprinting_emf2a/count/count.txt", header = T)
RefNef2 <- RefNef2[, c(-2:-6)]

NefRef1 <- read.table("./reads_count2/KN_count.txt", header = T)
NefRef1 <- NefRef1[, c(-2:-6)]

data1 <- cbind(RefNef1[, c(1, 2:5)], NefRef1[, c(2:5)])      # 3DAP
data2 <- cbind(RefNef1[, c(1, 6:9)], NefRef1[, c(6:9)])      # 5DAP
data3 <- cbind(RefNef2[, c(1, 2:5)], RefNef2[, c(6:9)])  # 7DAP
data4 <- cbind(RefNef1[, c(1, 14:17)], NefRef1[, c(14:17)])  # embryo
rm(RefNef1)
rm(NefRef1)
rm(RefNef2)

# Assign maternal and paternal depending on cross direction
mod_data1 <- modidyData.v3(data1)

mod_data1filt <- modidyData.SCfilt(data1)

rm(data1)
rm(data2)
rm(data3)
rm(data4)


######################
# create design matrix
#####################
cros = factor(c(rep("Refnef", 3), rep("Nefref", 3)))
mother = factor(rep(c("mother","father"), 3))
edgeR.design <- model.matrix( ~ cros + mother)

######################
# detect allele-specific reads
#####################
BH1 <- GLM_edgeR(mod_data1, "NK_3DAP", out_path)
BH1filt <- GLM_edgeR(mod_data1filt, "NK_3DAP.filt", out_path)

######################
# parental proportion by TMM normalized
#####################
MatRatePlot.filt(mod_data1, mod_data1filt, "NK_3DAP.filt", BH1, BH1filt, out_path, "NIP-KAS 3 DAP endosperm")
