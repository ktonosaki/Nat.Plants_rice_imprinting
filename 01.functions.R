library("edgeR")
library("ggplot2")
library("ggsci")
library("ggplot2")
library("TCC")
library(Biobase)
library(dplyr)
library('plotrix')
library('reshape2')
library(svglite)


modidyData.v3 <- function(data){
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

modidyData.SCfilt <- function(data){
  names(data) <- head
  data1.ratio <- merge(data, seed[,c(1,4)], by = "locus_name", all =T)
  data1.ratio$ratio[is.na(data1.ratio$ratio)] <- 1
  data1.ratio <- na.omit(data1.ratio)
  #  data1.ratio <- subset(data1.ratio, ratio > 0 )

  print(dim(data1.ratio))
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
  print(dim(data.rm))

  data.rm2 <- data.rm

  # matrate apply
  data.rm2$mat_Refnef_1 <- data.rm$mat_Refnef_1 * data.rm$ratio
  data.rm2$mat_Refnef_2 <- data.rm$mat_Refnef_2 * data.rm$ratio
  data.rm2$mat_Nefref_1 <- data.rm$mat_Nefref_1 * data.rm$ratio
  data.rm2$mat_Nefref_2 <- data.rm$mat_Nefref_2 * data.rm$ratio
  mod_data1 <- data.rm2[, 1:8]
  return(mod_data1)
}


GLM_edgeR <- function(mod_data, cross, path){
  # GLM fit
  edgeR <- DGEList(counts=mod_data, group=mother)
  edgeR <- calcNormFactors(edgeR)
  edgeR <- estimateGLMCommonDisp(edgeR, edgeR.design)
  edgeR <- estimateGLMTrendedDisp(edgeR, edgeR.design)
  edgeR <- estimateGLMTagwiseDisp(edgeR, edgeR.design)
  edgeR.fit <- glmFit(edgeR, edgeR.design)
  edgeR.lrt <- glmLRT(edgeR.fit, coef= "mothermother")

  # BVC
  tiff(paste0(path, "/BCV_", cross, ".tiff"), pointsize=30, width=param_fig[1], height=param_fig[2])
  plotBCV(edgeR)
  dev.off()
  rm(edgeR)
  rm(edgeR.fit)

  # result table
  test_result <- as.data.frame(edgeR.lrt$table)
  BH <- decideTestsDGE(edgeR.lrt, p=fdr_cutoff, adjust="BH")
  BH_list <- cbind(test_result, BH)
  BH_list$QValue <- p.adjust(BH_list$PValue, method="BH")
  BH_list$locus_name <- row.names(BH_list)
  BH_list <- BH_list[, c(7, 1:4,6,5)]
  print(dim(BH_list))

  BH_list <- merge(BH_list, anno[, c(1,7)], by = "locus_name")
  print(dim(BH_list))

  outfile <- paste0(path, "/List_edgR_", cross, ".txt")
  write.table(BH_list, outfile, append=F, quote=F, row.names=F, sep ="\t")
  rm(outfile)

  # plot MA
  outfig <- paste0(path, "/MA_", cross, ".tiff")
  tiff(outfig, pointsize=30, width=param_fig[1], height=param_fig[2])
  tag <- rownames(BH_list)[BH_list$PValue < fdr_cutoff]
  plotSmear(edgeR.lrt, de.tags=tag, ylim = c(-15, 15))
  dev.off()

  rm(edgeR.lrt)
  rm(test_result)
  rm(tag)
  rm(BH_list)
  rm(outfig)

  return(BH)
}

MatRatePlot.filt <- function(data, data.filt, cross, BH.all, BH.filt, path, title){

Count.Data <- data
BH <- BH.all


#Count.Data <- mod_data1
#BH <- BH1

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
Count.Data.rate$Refnef.rate <- with(Count.Data.rate, (Refnef_1.rate + Refnef_2.rate)/2)
Count.Data.rate$Nefref.rate <- with(Count.Data.rate, (Nefref_1.rate + Nefref_2.rate)/2)

# result table
rate_list <- cbind(Count.Data.rate, BH)
rate_list$locus_name <- row.names(rate_list)
rate_list_meg <- filter(rate_list, mothermother == 1,
                        Refnef_1.rate >= cutoff.mat, Refnef_2.rate >= cutoff.mat, Nefref_1.rate >= cutoff.mat, Nefref_2.rate >= cutoff.mat)
rate_list_peg <- filter(rate_list, mothermother == -1,
                        Refnef_1.rate <= cutoff.pat, Refnef_2.rate <= cutoff.pat, Nefref_1.rate <= cutoff.pat, Nefref_2.rate <= cutoff.pat)
rate_list_Not <- filter(rate_list, !locus_name %in% c(rate_list_meg$locus_name, rate_list_peg$locus_name))
rate_list_Not$mothermother <- 0
rate_list <- rbind(rate_list_meg, rate_list_peg, rate_list_Not)

head(rate_list)
rm(Count.Data, Count.Data.rate)

## MEGs
Count.Data <- data.filt
BH <- BH.filt

#Count.Data <- mod_data1filt
#BH <- BH1filt

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

head(Count.Data.rate)

# Count.Data.rate$mean <- (Count.Data.rate$Refnef.rate + Count.Data.rate$Nefref.rate) / 2
# mat_rate <- Count.Data.rate[c("Gene", "mean")]
# names(mat_rate)[1] <- "locus_name"

# result table
rate_list.filt <- cbind(Count.Data.rate, BH)
rate_list.filt$locus_name <- row.names(rate_list.filt)
rate_list_meg.filt <- filter(rate_list.filt,
                             mothermother == 1,
                             Refnef_1.rate >= cutoff.mat,
                             Refnef_2.rate >= cutoff.mat,
                             Nefref_1.rate >= cutoff.mat,
                             Nefref_2.rate >= cutoff.mat)

rate_list2 <- merge(rate_list,
                    rate_list_meg.filt[,c(14,15)], by = "locus_name", all= T)
names(rate_list2)[15:16] <- c("mothermother","filt")

# NA replace 0
rate_list2$filt[is.na(rate_list2$filt)] <- 0

# replace by SC filter
rate_list2$mothermother[rate_list2$mothermother == 1 & rate_list2$filt == 0] <- 0
rate_list2$mothermother[rate_list2$mothermother == 0 & rate_list2$filt == 1] <- 1

outfile <- paste0(path, "/List_MatRate_", cross, ".txt")
write.table(rate_list2, outfile, append=F, quote=F, row.names=F, sep ="\t")


# plot
rate_list_meg <- subset(rate_list2, mothermother == 1)
rate_list_peg <- subset(rate_list2, mothermother == -1)
rate_list_Not <- subset(rate_list2, mothermother == 0)

rate_list_meg$class <- "MEG"
rate_list_peg$class <- "PEG"
rate_list_Not$class <- "Bi-allele"

ggplot(data=rate_list_Not) +
  geom_point(data=rate_list_Not,aes(x=Refnef.rate,
                                    y=Nefref.rate,
                                    color=class, alpha=class)) +
  geom_point(data=rate_list_meg,aes(x=Refnef.rate,
                                    y=Nefref.rate,
                                    color=class)) +
  geom_point(data=rate_list_peg,aes(x=Refnef.rate,
                                    y=Nefref.rate,
                                    color=class)) +
  theme_bw() +
  scale_color_manual(values=c("darkgrey","firebrick","steelblue")) +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1) ,
        axis.ticks=element_line(size = 1, colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(face="plain", color="black",size=14, family = "Helvetica"),
        axis.title=element_text(face="plain", color="black",size=14)) +
  labs(title=title) +
  xlab(paste0("Average maternal rate\n(", cv[1], " × ", cv[2], ")")) +
  ylab(paste0("Average maternal rate\n(", cv[2], " × ", cv[1], ")"))

ggsave(paste0(path, "/plot_", cross, ".eps"),
       dpi = 1200, width = 5, height = 5.1,
       device=cairo_ps)

}


MatRatePlot.New <- function(data, cross, BH, path){

  Count.Data <- data

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
  rate_list <- cbind(Count.Data.rate, BH)

  rate_list$locus_name <- row.names(rate_list)
  rate_list_meg <- filter(rate_list,
                          mothermother == 1,
                          Refnef_1.rate >= cutoff.mat, Refnef_2.rate >= cutoff.mat, Nefref_1.rate >= cutoff.mat, Nefref_2.rate >= cutoff.mat)
  rate_list_peg <- filter(rate_list,
                          mothermother == -1,
                          Refnef_1.rate <= cutoff.pat, Refnef_2.rate <= cutoff.pat, Nefref_1.rate <= cutoff.pat, Nefref_2.rate <= cutoff.pat)
  rate_list_Not <- filter(rate_list, !locus_name %in% c(rate_list_meg$locus_name, rate_list_peg$locus_name))
  rate_list_Not$mothermother <- 0
  rate_list <- rbind(rate_list_meg, rate_list_peg, rate_list_Not)

  outfile <- paste0(path, "/List_MatRate_", cross, ".txt")
  write.table(rate_list, outfile, append=F, quote=F, row.names=F, sep ="\t")


  rate_list_meg$class <- "MEG"
  rate_list_peg$class <- "PEG"
  rate_list_Not$class <- "Bi-allele"

  ggplot(data=rate_list_Not) +
    geom_point(data=rate_list_Not,aes(x=Refnef.rate,
                                      y=Nefref.rate,
                                      color=class, alpha=class)) +
    geom_point(data=rate_list_meg,aes(x=Refnef.rate,
                                      y=Nefref.rate,
                                      color=class)) +
    geom_point(data=rate_list_peg,aes(x=Refnef.rate,
                                      y=Nefref.rate,
                                      color=class)) +
    theme_bw() +
    scale_color_manual(values=c("darkgrey","firebrick","steelblue")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 1) ,
          axis.ticks=element_line(size = 1, colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color="black", size=12),
          axis.title=element_text(size=14,face="bold")) +
    xlab(paste0("Maternal rate in ", cv[1], " x ", cv[2])) +
    ylab(paste0("Maternal proportion in ", cv[2], " x ", cv[1]))

  ggsave(paste0(path, "/plot_", cross, ".tiff"),
         dpi = 300, width = 5.4, height = 4.8)
}

Imp.pie <- function(data, code, file){
  consistent <- length(grep("consistent", data$status))
  NK.acc <- length(grep("NK accession", data$status))
  NK.pos <- length(grep("NK possible", data$status))
  NQ.acc <- length(grep("NQ accession", data$status))
  NQ.pos <- length(grep("NQ possible", data$status))
  acc <- NK.acc + NQ.acc
  pos <- NK.pos + NQ.pos

  if(code == "MEG"){
    col <- list(cons = "#C01100",
                acc = "#ED7D31",
                pos = "#FFC000",
                Kas = "#A633D4",
                Q = "#DFB10F",
                NO = "white")
  } else if(code == "PEG") {
    col <- list(cons = "#0570C0",
                acc = "#5B9BD5",
                pos = "#9DC3E6",
                Kas = "#00ADFF",
                Q = "#01C165",
                NO = "white")
  }

  setEPS()
  postscript(file)
  par(lwd = 5)
  pie(c(consistent, NK.acc, NQ.acc, NK.pos, NQ.pos),
      radius=1,lty = "solid",　lwd = 2,
      col=as.character(col[c('NO', rep(c("Kas", "Q"), 2))]),
      border="white", labels='',init.angle = 90, clockwise = TRUE)
  par(new=T)
  pie(c(consistent, acc, pos),
      radius=.8,　
      col=as.character(col[c('cons',"acc", "pos")]),
      border="white", labels='',init.angle = 90, clockwise = TRUE)
  par(new=T)
  pie(c(consistent, acc, pos),
      radius=.4,　
      col="white",
      border="white", labels='',init.angle = 90, clockwise = TRUE)
  dev.off()

}





imp_status <- function(data_list, code){
  data_list$contami <- NULL
  data_list$ref <- NULL
  data_list$status[with(data_list, State_nk == code & State_nq == code)] <- paste0("consistent ",code)
  data_list$status[with(data_list, State_nk == code & State_nq == "no detect")] <- paste0("NK possible ",code)
  data_list$status[with(data_list, State_nk == code & State_nq == "biallele")] <- paste0("NK accession ",code)
  data_list$status[with(data_list, State_nk == "no detect" & State_nq == code)] <- paste0("NQ possible ",code)
  data_list$status[with(data_list, State_nk == "biallele" & State_nq == code)] <- paste0("NQ accession ",code)

  data_list$MatRate_nq[with(data_list, status == paste0("NK possible ",code))] <- -0.05
  data_list$MatRate_nk[with(data_list, status == paste0("NQ possible ",code))] <- -0.05

  return(data_list)
}




plot_imp <- function(data_non, data_meg, data_peg, stage){

  hoge0 <- data_non
  hoge1 <- data_meg
  hoge2 <- data_peg
  cross <- stage

  hoge0 <- subset(hoge0, MatRate_nk != "-" & MatRate_nq != "-")
  hoge0$MatRate_nk <- as.numeric(hoge0$MatRate_nk)
  hoge0$MatRate_nq <- as.numeric(hoge0$MatRate_nq)
  #  hoge1 <- subset(hoge1, MatRate_nk != "-" & MatRate_nq != "-")
  hoge1$MatRate_nk <- as.numeric(hoge1$MatRate_nk)
  hoge1$MatRate_nq <- as.numeric(hoge1$MatRate_nq)
  #  hoge2 <- subset(hoge2, MatRate_nk != "-" & MatRate_nq != "-")
  hoge2$MatRate_nk <- as.numeric(hoge2$MatRate_nk)
  hoge2$MatRate_nq <- as.numeric(hoge2$MatRate_nq)

  ggplot(data=hoge0) +
    geom_hline(yintercept=0,linetype="dashed", colour="black") +
    geom_vline(xintercept=0,linetype="dashed", colour="black") +
    geom_point(data=hoge0, aes(x=MatRate_nk, y=MatRate_nq, color=status, alpha = 0.2)) +
    geom_point(data=hoge1, aes(x=MatRate_nk, y=MatRate_nq, color=status, alpha = 0.3)) +
    geom_point(data=hoge2, aes(x=MatRate_nk, y=MatRate_nq, color=status, alpha = 0.3)) +
    #    xlim(0,1) + ylim(0, 1) +
    theme_bw() +
    scale_color_manual(values=c("darkgrey",
                                "firebrick", "#0168b3",
                                "darkorchid","#01A6FF", "darkorchid", "#01A6FF",
                                "goldenrod", "#20B85D", "goldenrod", "#20B85D")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", size = 1) ,
          axis.ticks=element_line(size = 1, colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color="black", size=20),
          axis.title=element_text(size=14,face="bold")) +
    xlab(paste0("maternal proportion of Nipponbare and Kasalath combination")) +
    ylab(paste0("maternal proportion of Nipponbare and 93-11 combination"))
  ggsave(paste0("identify_gene/compare/plot_", cross, ".tiff"),
         dpi = 300, width = 6.4, height = 5.8)

  ggplot(data=hoge0) +
    geom_point(data=hoge0, aes(x=MatRate_nk, y=MatRate_nq, color=status)) +
    geom_point(data=hoge1, aes(x=MatRate_nk, y=MatRate_nq, color=status)) +
    geom_point(data=hoge2, aes(x=MatRate_nk, y=MatRate_nq, color=status)) +
    #    xlim(0,1) + ylim(0, 1) +
    theme_bw() +
    scale_color_manual(values=c("darkgrey",
                                "firebrick", "#0168b3",
                                "darkorchid","#01A6FF", "darkorchid", "#01A6FF",
                                "goldenrod", "#20B85D", "goldenrod", "#20B85D")) +
    theme(legend.position = 'right',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks=element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold")) +
    xlab(paste0("maternal proportion of Nipponbare and Kasalath combination")) +
    ylab(paste0("maternal proportion of Nipponbare and 93-11 combination"))
  ggsave(paste0("identify_gene/compare/plot_regend_", cross, ".tiff"),
         dpi = 300, width = 5.4, height = 4.8)
}
