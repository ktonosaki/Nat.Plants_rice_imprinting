library("data.table")
library("tidyverse")

file1 <- "NK_r1_mat.txt.gz" # maternal 
file2 <- "NK_r1_pat.txt.gz" # paternal
file3 <- "NK_r1_CpG_DMRs.bed"
FDR_cutoff <- 0.05
inforC <- 3
context <- "CpG"
diff_cutoff <- 0.2

# import data
data1 <- fread(file1, col.names = c("chr", "start", "end", "text", "metC.1", "unmetC.1", "Csite.1"))
data2 <- fread(file2, col.names = c("chr", "start", "end", "text", "metC.2", "unmetC.2", "Csite.2"))
data1$rate.1 <- data1$metC.1 / (data1$metC.1 + data1$unmetC.1)
data2$rate.2 <- data2$metC.2 / (data2$metC.2 + data2$unmetC.2)

data <- merge(subset(data1, text == context & Csite.1 >= inforC)[,-4],
              subset(data2, text == context & Csite.2 >= inforC)[,-4], by = c("chr", "start", "end"))

# head(data)

# mod 210709
#data$diff <- data %>% select(rate.1, rate.2) %>%
#                      summarize(rate.1 - rate.2)
data$diff <- data$rate.1 - data$rate.2
##

# fisher test
p.value <- data %>%
  select(metC.1, unmetC.1, metC.2, unmetC.2) %>%
  pmap_dbl(function(metC.1, unmetC.1, metC.2, unmetC.2){
    matrix <- matrix(c(metC.1, unmetC.1, metC.2, unmetC.2), 2)
    out <- fisher.test(matrix)
    return(signif(out$p.value, 3))})

# FDR
data$FDR <- p.adjust(p.value, method = "BH")
data$DMR[data$FDR < FDR_cutoff] <- TRUE
data$DMR[is.na(data$DMR)] <- FALSE
data$class[data$diff >= diff_cutoff] <- "hyper"
data$class[data$diff <= -diff_cutoff] <- "hypo"
data$DMR[is.na(data$class)] <- FALSE
DMRs <- subset(data, DMR == TRUE)

DMRs <- DMRs[, c("chr","start","end","diff","FDR","class")]

num_chr <- sort(unique(DMRs$chr))
if (mode(num_chr) == "numeric") {
  DMRs$chr <- formatC(as.numeric(DMRs$chr),width=2,flag=0)
}
#DMRs$chr <- formatC(as.numeric(DMRs$chr), width=2, flag = 0)
DMRs$diff <- round(DMRs$diff, 4)
options(scipen=2)
DMRs$FDR <- signif(DMRs$FDR, digits = 3)

# output data
fwrite(DMRs, file3, append = F, quote = F, sep = "\t", col.names = F, row.names = F)

print("")
print("end multiple Fisher exact test")
print("")
