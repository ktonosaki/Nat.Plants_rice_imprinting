library("data.table")
library("tidyverse")
library("R.utils")
library("furrr")

core <- 16
file <- "NK_met_fies.txt"
out <- "out_dir"
coverage <- 5
FDR_cutoff <- 0.01
Pt <- "Pt"

# control chromosome
dir.create(out)

# paramaters
options(scipen=100)
result.table <- NULL


# import data
print("")
print(paste0("reads data ", file))
data <- fread(file)
print("")
names(data) <- c("chr", "pos", "strand","met", "unmet", "context")
#data$cov <- data %>% select(met,unmet) %>% summarize(met+unmet)
data$cov <- data %>% select(met,unmet) %>% reframe(met+unmet)
data <- filter(data, cov >= coverage)

# conversion rate
Pt.P <- data %>%
  select(chr, met, cov) %>%
  filter(chr == Pt) %>%
  reframe(mean(met/cov)) %>%
  as.numeric()
print("")
print(paste0("Non conversion rate is ", Pt.P))
print("")
write.table(Pt.P, paste0(out, "/NonConversionRate_on_", Pt, ".txt"), append = F, quote = F, sep = "\t", col.names = F, row.names = F)

# binomial test by conversion rate in Pt
data.chr <- filter(data, chr != Pt)
num_chr <- sort(unique(data.chr$chr))


# run test by each chromosome
plan("multisession", workers = core)
for (i in 1:length(num_chr)){
  print("")
  print(paste0("Binominal test for Methylated cytosine of chr ", num_chr[i]))
  print("")
  data.chr.ex <- filter(data.chr, chr == num_chr[i])
  p.value <- future_map2(.x = data.chr.ex$met,
                         .y = data.chr.ex$cov,
                         .f = function(x,y){binom.test(x, y, Pt.P)$p.value},
                         .progress = TRUE) %>%
    as.vector(mode = "numeric")
  print("")
  print(paste0("Binominal test Ended for chr ", num_chr[i]))
  print("")
  data.chr.ex$FDR <- p.adjust(p.value, method = "BH")
  data.chr.ex$exact[data.chr.ex$FDR < FDR_cutoff] <- TRUE
  data.chr.ex$exact[is.na(data.chr.ex$exact)] <- FALSE
  data.chr.ex$sig.met <- data.chr.ex$met
  data.chr.ex$sig.met[data.chr.ex$exact == FALSE] <- 0
  data.chr.ex$end <- data.chr.ex$pos + 1
  data.chr.ex <- data.chr.ex %>% select(context, chr, pos, end, strand, sig.met, unmet, cov)

  if (mode(num_chr[i]) == "numeric") {
    data.chr.ex$chr <- formatC(as.numeric(data.chr.ex$chr),width=2,flag=0)
  }

  print(paste0("Save temp files for chr ", num_chr[i]))
  print("")
  fwrite(data.chr.ex, paste0(out, "/binom.data_", num_chr[i], ".txt.gz"), append = F, quote = F, sep = "\t", col.names = F, row.names = F)
  print("")
  rm(data.chr.ex)
}
