# C. Ahern
# 09/06/2023
# Using ALDEx2 on my qiime2 data

library(tidySummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(BiocManager)
library(Biostrings)
library(mia)
library(ALDEx2)
library(qiime2R)
library(readxl)
BiocManager::install("ALDEx2")
setwd("/Users/colleenahern/Documents/BASF")

# Import qiime2 results and metadata info
pswd <- qza_to_phyloseq(features="Analysis2/16S-table-noplant-rarefied-10000_filtered.qza",
                        tree="Analysis2/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "Analysis2/16S-rep-seqs-taxonomy.qza", 
                        metadata= "Analysis2/metadata2.txt")

df_env <- read.delim("/Users/colleenahern/Documents/BASF/Analysis2/metadata2.txt") %>%
  arrange(sample.id) %>% column_to_rownames('sample.id')

# Subset conditions of interest for ALDEx2 analysis
dfenvsub1 <- df_env[df_env$Day == 90,]
dfenvsub2 <- dfenvsub1[dfenvsub1$Substrate == "Oligomer" | dfenvsub1$Substrate == "Copolymer" | dfenvsub1$Substrate == "HDPE" | dfenvsub1$Substrate == "Blank",] # results in 76

dfenvsub2$Category <- 0
dfenvsub2$Category[c(1:9)] <- "No Biodegradation"
dfenvsub2$Category[c(10:12)] <- "Biodegradation"

df <- as.data.frame(pswd@otu_table) # %>% column_to_rownames('argtype')

dfsub2 = as.data.frame(df[,rownames(dfenvsub2)])

condit_mm = data.frame("Substrate" = as_factor(dfenvsub2$Substrate),
                       "Day" = as_factor(dfenvsub2$Day),
                       "Replicate" = dfenvsub2$Replicate,
                       "Group" = dfenvsub2$Group,
                       "Category" = as_factor(dfenvsub2$Category))

rownames(condit_mm) = colnames(dfsub2) 
mm = model.matrix(~Category, condit_mm)

r_df = as.matrix(round(dfsub2)) # Rounding data for clr

drx = aldex.clr(dfsub2, dfenvsub2$Category)    # this one worked for subsetting for ttest

z_tt <- aldex.ttest(drx, paired.test = FALSE, verbose = TRUE)
z_effect <- aldex.effect(drx, CI = TRUE, verbose = TRUE)
aldex_out <- data.frame(z_tt, z_effect)
par(mfrow = c(1, 2))

aldex.plot(aldex_out,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

# aldex_out %>%
#   rownames_to_column(var = "Genus") %>%
#   # here we choose the wilcoxon output rather than t-test output
#   filter(wi.eBH <= 0.05)  %>%
#   dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
#   knitr::kable()

aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  #filter(abs(effect) >= 1 & abs(diff.btw >2))  %>%
  filter(abs(effect) >= 1)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()


aldex_out4 <- aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  #filter(abs(effect) >= 1 & abs(diff.btw >2))  %>%
  filter(abs(effect) >= 1)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) 

# subset OTUs from aldex_out2
df_ed <- dfsub2[aldex_out4$Genus,]

drxed = aldex.clr(df_ed, dfenvsub2$Category)  

z_tted <- aldex.ttest(drxed, paired.test = FALSE, verbose = TRUE)
z_effected <- aldex.effect(drxed, CI = TRUE, verbose = TRUE)
aldex_outed <- data.frame(z_tted, z_effected)
par(mfrow = c(1, 2))

aldex.plot(aldex_outed,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_outed,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

aldex_outed %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()

aldex_out3 <- aldex_outed %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap)

taxonomyca <- pswd@tax_table
taxonomyca <- as.data.frame(taxonomyca@.Data)
taxonomyca <- rownames_to_column(taxonomyca, "ID")

library("dplyr")
aldex_out3 <- aldex_out3 %>%
  dplyr::rename("ID" = "Genus")

aldex_outed3 <- merge(taxonomyca, aldex_out3, "ID")

taxonomyca$otu <- paste0("otu", 1:nrow(otu_countsca))  # MAKE SURE YOUR OTU NUMBERS MATCH UP WITH THE ORIGINAL COMPLICATED NAMES
taxonomyca <- taxonomyca %>%
  select(otu, everything()) %>%
  rename_all(tolower)
