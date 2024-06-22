#!/usr/bin/env Rscript

rm(list = ls())
args <- commandArgs(TRUE)


in_file <- "NK_K27me3_pat.bdg"
out <- "NK_K27me3_pat_z-score"

# import
header <- c("chr", "start","end","rpgc")
file <- read.csv(file = in_file, sep = "\t", header = F, col.names = header)

# mod file
rpgc_mean <- mean(file$rpgc)
rpgc_variance <- sd(file$rpgc)
file$z <- (file$rpgc - rpgc_mean) / rpgc_variance
file <- file[, c("chr", "start","end","z")]
# file$z[file$z < 0] <- 0

if (mode(file$chr) == "numeric") {
  file$chr <- formatC(as.numeric(file$chr),width=2,flag=0)
}

write.table(x = file, paste0(out, ".bdg"),
            append=F, quote=F, sep = "\t", row.names = F, col.names = F)
