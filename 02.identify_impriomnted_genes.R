source ("01.function.R")

fdr_cutoff <- 0.05
cutoff.mat <- 0.70
cutoff.pat <- 0.30
depth <- 20
cv <- c("NIP", "KAS")
head <- c("locus_name",
          "mat_Refnef_1","pat_Refnef_1","mat_Refnef_2","pat_Refnef_2",
          "pat_Nefref_1","mat_Nefref_1","pat_Nefref_2","mat_Nefref_2")

######################
# load and prepare data
######################
# ref
seed <- read.table("seed_cort_gene.ratio.txt", header = T)

# import read count data from featureCounts
RefNef <- read.table("NipKas_count.txt", header = T)[, c(-2:-6)]
NefRef <- read.table("KasNip_count.txt", header = T)[, c(-2:-6)]

data1 <- cbind(RefNef[, c(1, 2:5)], NefRef[, c(2:5)])      # 3DAP
data2 <- cbind(RefNef[, c(1, 6:9)], NefRef[, c(6:9)])      # 5DAP
data3 <- cbind(RefNef[, c(1, 10:13)], NefRef[, c(10:13)])  # 7DAP
data4 <- cbind(RefNef[, c(1, 14:17)], NefRef[, c(14:17)])  # embryo

# Assign maternal and paternal depending on cross direction
mod_data_maternal <- modidyData.mat(data1, depth, seed)
mod_data_paternal <- modidyData.pat(data1, depth)

# edgeR
cros = factor(c(rep("Refnef", 3), rep("Nefref", 3)))
mother = factor(rep(c("mother","father"), 3))
edgeR.design <- model.matrix( ~ cros + mother)

BH_maternal <- GLM_edgeR(mod_data_maternal)
BH_paternal <- GLM_edgeR(mod_data_paternal)

# filt by maternal rate
meg_lis <- findMEGs(mod_data1, BH1)
peg_lis <- findPEGs(mod_data1, BH1filt)
