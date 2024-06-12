library(Seurat)
library(hdf5r)
library(tidyverse)
library(ggplot2)

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
  data$sample <- samp[i]
  f_list[[i]] <- data
}

# 02.QC
f_list.filt <- list(NULL)
for (i in num){
  lim <- quantile(f_list[[i]]$nFeature_RNA)[2]
  lim <- 1000
  data <- subset(f_list[[i]], subset = nFeature_RNA > 100 &
                   nFeature_RNA < 10000 )
  f_list.filt[[i]] <- data
}

## 03.remove doublet
f_list.rm <- NULL
for (i in num){
  pbmc <- f_list.filt[[i]]
  doublet_file <- read.table(paste0("/home/epigenome/snRNA/R_Seurat/doublet/", samp[i], "_doublet.txt"), header = T, stringsAsFactors = F)
  doublet <- data.frame(cellname = colnames(pbmc), doublet = "Singlet")

  id1 <- doublet_file$doublet[match(doublet$cellname, doublet_file$cellname)] == "Doublet"
  id1[is.na(id1)] <- FALSE
  doublet$doublet[id1] <- "Doublet"
  sum(doublet$cellname == colnames(pbmc))
  pbmc[["doublet_filter"]] <- doublet$doublet
  pbmc <- subset(pbmc, subset = doublet_filter == "Singlet")
  f_list.rm[[i]] <- pbmc
}


## 04.merge object
pbmc.comb1 <- merge(f_list.rm[[1]], y = f_list.rm[[2]],
                    add.cell.ids = c("tech1", "tech2"),
                    project = "5d_r1",
                    merge.data = TRUE)
pbmc.comb2 <- merge(f_list.rm[[3]], y = f_list.rm[[4]],
                    add.cell.ids = c("tech1", "tech2"),
                    project = "5d_r2",
                    merge.data = TRUE)

rep1.nc <- pbmc.comb1$nCount_RNA
rep2.nc <- pbmc.comb2$nCount_RNA
d1 <- density(rep1.nc)
d2 <- density(rep2.nc)
dat <- data.frame(x = c(d1$x, d2$x),
                  y = c(d1$y, d2$y),
                  repn = rep(c("5d-1", "5d-2"),
                             times = c(length(d1$x), length(d2$x))))
s1d <- as.numeric(quantile(rep1.nc, 0.1))
ggplot(dat, aes(x=x, y=y, color=repn)) +
  geom_line() + theme_classic() + xlab("Total mumber of UMIs") + ylab("Density") +
  scale_x_continuous(expand = c(0,0), limits = c(min(dat$x), max(dat$x)*1.1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(dat$y)*1.1)) +
  annotate("text", x = s1d + 0.2*max(dat$x),
           y = max(dat$y)*1.05,
           label= paste0("Thresh = ", round(s1d, 2)),
           color = "black") + xlim(0, 5000)
ggsave("./Figure/hist_count.pdf", device = "pdf",
       dpi = 600, width = 5, height = 4)

# 05.SCT normalize
ifnb.list <- list(five.r1.1 = f_list.rm[[1]], five.r1.2 = f_list.rm[[2]],
                  five.r2.1 = f_list.rm[[3]], five.r2.2 = f_list.rm[[4]])

ifnb.list.norm <- lapply(X = ifnb.list, FUN = function(x) {
  x <- SCTransform(object = x, vst.flavor = "v2",
                   assay = 'RNA',
                   new.assay.name = 'SCT',
                   vars.to.regress = c('nFeature_RNA', 'nCount_RNA'),
                   return.only.var.genes = F) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
})


# 06. INTERGRATION
#### Perform INTERGRATION using pearson residuals
features <- SelectIntegrationFeatures(object.list = ifnb.list.norm,
                                      nfeatures = 2000)
ifnb.list.norm <- PrepSCTIntegration(object.list = ifnb.list.norm,
                                     anchor.features = features,
                                     verbose = T)
anchors <- FindIntegrationAnchors(object.list = ifnb.list.norm,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  reduction="rpca",
                                  verbose = T)
all_genes <- lapply(ifnb.list.norm, row.names) %>% Reduce(intersect, .)
combined.sct <- IntegrateData(anchorset = anchors,
                              normalization.method = "SCT",
                              features.to.integrate = all_genes,
                              verbose = T,
                              new.assay.name = "Integrated")

## add samples info
samples <- c("Chro-5dEN-r1", "Chro-5dEN-r2")
combined.sct@meta.data$sample2 <- NULL
for (fa in samples){
combined.sct@meta.data$sample2[grepl(fa, combined.sct@meta.data$sample)] <- fa
}

# 08. integrated analysis
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:50, verbose = TRUE)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:50)
combined.sct <- FindClusters(combined.sct, resolution = 0.4)


DimPlot(combined.sct, reduction = "umap", label = T)
ggsave(paste("./Figure/UMAP/DimPlot_umap_all.5dap.pdf"), device = "pdf",
       dpi = 600, width = 5, height = 4)

write.table(combined.sct$seurat_clusters,
            "pbmc.combined_cluster.txt", quote = F, sep = "\t")
