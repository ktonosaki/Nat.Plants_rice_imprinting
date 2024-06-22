library(DoubletFinder)
library(Seurat)
library(hdf5r)
library(tidyverse)


# create counts matrix
data_path = "/after_biobender"

samp <- c("Chro-5dEN-r1-1", "Chro-5dEN-r1-2", "Chro-5dEN-r2-1", "Chro-5dEN-r2-2")

# 01.import
f_list <- list(NULL)
num <- seq(1,length(samp))
for (i in num){
  file <- paste0(data_path, "/", samp[i], ".output_filtered_seurat.h5")
  data <- Read10X_h5(filename = file, use.names = TRUE)
  data <- CreateSeuratObject(counts = data, min.cells = 5)
  data$mitoPercent <- PercentageFeatureSet(data, pattern = "^gene:")
  data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA)
  data$sample <- samp[i]
  f_list[[i]] <- data
}

# data check
for (i in num){
  VlnPlot(f_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "mitoPercent", "log10GenesPerUMI"), ncol = 4)
  ggplot2::ggsave(paste0(fig_path, "/vlnplot/", samp[i], "_pre.pdf"), device = "pdf",
                dpi = 1200, width = 10, height = 8)
}

# 02.QC
f_list.filt <- list(NULL)
for (i in num){
  lim <- quantile(f_list[[i]]$nFeature_RNA)
  data <- subset(f_list[[i]], subset = nFeature_RNA > 100 &
                   nFeature_RNA < 10000 &
                   mitoPercent < 5 &
                   log10GenesPerUMI > 0.8)
  f_list.filt[[i]] <- data
}

# 03.pre-process standard workflow
Prep.doubletF <- function(pbmc.seurat.filtered){
  pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
  ElbowPlot(pbmc.seurat.filtered)
  pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:50)
  pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:50)
  return(pbmc.seurat.filtered)
}

f_list.doub_pre <- list(NULL)
for (i in num){
  data <- Prep.doubletF(f_list.filt[[i]])
  f_list.doub_pre[[i]] <- data
}

# 04.doblet finder
f_list.doub <- list(NULL)
for (i in num){
    pbmc.seurat.filtered <- f_list.doub_pre[[i]]
    ## pK Identification (no ground-truth)
    sweep.res.list_pbmc <- paramSweep_v3(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)
    sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
    bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

    pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
      filter(BCmetric == max(BCmetric))
    pK <- as.numeric(as.character(pK[2]))

    ## Homotypic Doublet Proportion Estimate
    annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)                  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    # run doubletFinder
    pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered,
                                             PCs = 1:20,
                                             pN = 0.25,
                                             pK = pK,
                                             nExp = nExp_poi.adj,
                                             reuse.pANN = FALSE, sct = FALSE)
    pbmc.seurat.filtered@meta.data$doublet <- pbmc.seurat.filtered@meta.data[,10]
    f_list.doub[[i]] <- pbmc.seurat.filtered

    df <- data.frame(cellname=names(pbmc.seurat.filtered@active.ident),
                     doublet=pbmc.seurat.filtered@meta.data$doublet)
    write.table(df,
                paste0("out_doubletFinder/", samp[i], "_doublet.txt"), quote = F, sep = "\t")
}
