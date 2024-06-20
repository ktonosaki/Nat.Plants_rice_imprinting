library("edgeR")

modidyData.mat <- function(data, depth, seed){
  names(data) <- head
  data1.ratio <- merge(data, seed[,c(1,4)], by = "locus_name", all =T)
  data1.ratio$ratio[is.na(data1.ratio$ratio)] <- 1
  data1.ratio <- na.omit(data1.ratio)

  row.names(data1.ratio) <- data1.ratio$locus_name
  data1.ratio$locus_name <- NULL
  data1.ratio$class <- NULL
  filter <- data.frame(row.names = row.names(data1.ratio))
  filter$Nefref_1 <- rowSums(data1.ratio[, grepl("Nefref_1", colnames(data1.ratio))]) >= depth
  filter$Refnef_1 <- rowSums(data1.ratio[, grepl("Refnef_1", colnames(data1.ratio))]) >= depth
  filter$Nefref_2 <- rowSums(data1.ratio[, grepl("Nefref_2", colnames(data1.ratio))]) >= depth
  filter$Refnef_2 <- rowSums(data1.ratio[, grepl("Refnef_2", colnames(data1.ratio))]) >= depth
  data1.ratio$filter <- apply(filter, 1, any)
  data.rm <- filter(data1.ratio, filter == TRUE)[, -ncol(data1.ratio)]

  # matrate apply
  data.rm2 <- data.rm
  data.rm2$mat_Refnef_1 <- data.rm$mat_Refnef_1 * data.rm$ratio
  data.rm2$mat_Refnef_2 <- data.rm$mat_Refnef_2 * data.rm$ratio
  data.rm2$mat_Nefref_1 <- data.rm$mat_Nefref_1 * data.rm$ratio
  data.rm2$mat_Nefref_2 <- data.rm$mat_Nefref_2 * data.rm$ratio
  mod_data1 <- data.rm2[, 1:8]

  return(mod_data1)
}

modidyData.pat <- function(data, depth){
  print(dim(data))
  names(data) <- head
  row.names(data) <- data$locus_name
  data$locus_name <- NULL
  data$class <- NULL
  filter <- data.frame(row.names = row.names(data))
  filter$Nefref_1 <- rowSums(data[,  grepl("Nefref_1", colnames(data))]) >= depth
  filter$Refnef_1 <- rowSums(data[, grepl("Refnef_1", colnames(data))]) >= depth
  filter$Nefref_2 <- rowSums(data[, grepl("Nefref_2", colnames(data))]) >= depth
  filter$Refnef_2 <- rowSums(data[, grepl("Refnef_2", colnames(data))]) >= depth
  data$filter <- apply(filter, 1, any)
  data.rm <- filter(data, filter == TRUE)[, -ncol(data)]
  print(dim(data.rm))
  return(data.rm)
}

GLM_edgeR <- function(mod_data){
  # GLM fit
  edgeR <- DGEList(counts=mod_data, group=mother)
  edgeR <- calcNormFactors(edgeR)
  edgeR <- estimateGLMCommonDisp(edgeR, edgeR.design)
  edgeR <- estimateGLMTrendedDisp(edgeR, edgeR.design)
  edgeR <- estimateGLMTagwiseDisp(edgeR, edgeR.design)
  edgeR.fit <- glmFit(edgeR, edgeR.design)
  edgeR.lrt <- glmLRT(edgeR.fit, coef= "mothermother")
  BH <- decideTestsDGE(edgeR.lrt, p=fdr_cutoff, adjust="BH")
  return(BH)
}


findMEGs <- function(Count.Data, BH){

  Count.Data.rate <- data.frame(Gene=row.names(Count.Data))
  Count.Data.rate$mat_RxN <- apply(Count.Data[, grepl("mat_Refnef", colnames(Count.Data))], 1, mean)
  Count.Data.rate$pat_RxN <- apply(Count.Data[, grepl("pat_Refnef", colnames(Count.Data))], 1, mean)
  Count.Data.rate$RxN_total <- with(Count.Data.rate, mat_RxN + pat_RxN)
  Count.Data.rate$mat_NxR <- apply(Count.Data[, grepl("mat_Nefref", colnames(Count.Data))], 1, mean)
  Count.Data.rate$pat_NxR <- apply(Count.Data[, grepl("pat_Nefref", colnames(Count.Data))], 1, mean)
  Count.Data.rate$NxR_total <- with(Count.Data.rate, mat_NxR + pat_NxR)

  Count.Data$pat_Refnef_1 <- Count.Data$pat_Refnef_1*2
  Count.Data$pat_Refnef_2 <- Count.Data$pat_Refnef_2*2
  Count.Data$pat_Nefref_1 <- Count.Data$pat_Nefref_1*2
  Count.Data$pat_Nefref_2 <- Count.Data$pat_Nefref_2*2

  Count.Data.rate$Refnef_1.rate <- with(Count.Data, mat_Refnef_1 / (mat_Refnef_1 + pat_Refnef_1))
  Count.Data.rate$Refnef_2.rate <- with(Count.Data, mat_Refnef_2  / (mat_Refnef_2 + pat_Refnef_2 ))
  Count.Data.rate$Nefref_1.rate <- with(Count.Data, mat_Nefref_1  / (pat_Nefref_1 + mat_Nefref_1 ))
  Count.Data.rate$Nefref_2.rate <- with(Count.Data, mat_Nefref_2  / (pat_Nefref_2 + mat_Nefref_2))
  Count.Data.rate$Refnef.rate <- with(Count.Data.rate, (Refnef_1.rate + Refnef_2.rate)/2)
  Count.Data.rate$Nefref.rate <- with(Count.Data.rate, (Nefref_1.rate + Nefref_2.rate)/2)

  # result table
  rate_list.filt <- cbind(Count.Data.rate, BH)
  rate_list.filt$locus_name <- row.names(rate_list.filt)
  rate_list_meg <- filter(rate_list.filt,
                               mothermother == 1,
                               Refnef_1.rate >= cutoff.mat,
                               Refnef_2.rate >= cutoff.mat,
                               Nefref_1.rate >= cutoff.mat,
                               Nefref_2.rate >= cutoff.mat)

  return(rate_list_meg)

}

findPEGs <- function(Count.Data, BH){

  Count.Data.rate <- data.frame(Gene=row.names(Count.Data))
  Count.Data.rate$mat_RxN <- apply(Count.Data[, grepl("mat_Refnef",
                                                      colnames(Count.Data))], 1, mean)
  Count.Data.rate$pat_RxN <- apply(Count.Data[, grepl("pat_Refnef",
                                                      colnames(Count.Data))], 1, mean)
  Count.Data.rate$RxN_total <- with(Count.Data.rate, mat_RxN + pat_RxN)

  Count.Data.rate$mat_NxR <- apply(Count.Data[, grepl("mat_Nefref",
                                                      colnames(Count.Data))], 1, mean)
  Count.Data.rate$pat_NxR <- apply(Count.Data[, grepl("pat_Nefref",
                                                      colnames(Count.Data))], 1, mean)
  Count.Data.rate$NxR_total <- with(Count.Data.rate, mat_NxR + pat_NxR)

  Count.Data$pat_Refnef_1 <- Count.Data$pat_Refnef_1*2
  Count.Data$pat_Refnef_2 <- Count.Data$pat_Refnef_2*2
  Count.Data$pat_Nefref_1 <- Count.Data$pat_Nefref_1*2
  Count.Data$pat_Nefref_2 <- Count.Data$pat_Nefref_2*2

  Count.Data.rate$Refnef_1.rate <- with(Count.Data, mat_Refnef_1 / (mat_Refnef_1 + pat_Refnef_1))
  Count.Data.rate$Refnef_2.rate <- with(Count.Data, mat_Refnef_2  / (mat_Refnef_2 + pat_Refnef_2 ))
  Count.Data.rate$Nefref_1.rate <- with(Count.Data, mat_Nefref_1  / (pat_Nefref_1 + mat_Nefref_1 ))
  Count.Data.rate$Nefref_2.rate <- with(Count.Data, mat_Nefref_2  / (pat_Nefref_2 + mat_Nefref_2))

  # result table
  rate_list <- cbind(Count.Data.rate, BH)
  rate_list$locus_name <- row.names(rate_list)

  rate_list_peg <- filter(rate_list, mothermother == -1,
                          Refnef_1.rate <= cutoff.pat,
                          Refnef_2.rate <= cutoff.pat,
                          Nefref_1.rate <= cutoff.pat,
                          Nefref_2.rate <= cutoff.pat)

  return(rate_list_peg)

}
