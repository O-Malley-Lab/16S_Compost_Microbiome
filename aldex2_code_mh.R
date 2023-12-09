# C. Ahern 
# Using ALDEx2 on my qiime2s data

library(tidySummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(BiocManager)
library(Biostrings)
library(mia)
library(ALDEx2)
library(qiime2R)
library(readxl)

# Load data
otu_tab = read_excel("/Users/colleenahern/Downloads/checking_daa/otu_table.xlsx")
sampleinfo = read_excel("/Users/colleenahern/Downloads/checking_daa/metadata.xlsx") %>%
  column_to_rownames("sample.id")
df_env <- read_excel("/Users/colleenahern/Downloads/checking_daa/metadata.xlsx") %>%
  arrange(sample.id) %>% column_to_rownames('sample.id')

# Grouping sampleinfo
sampleinfo = sampleinfo %>%
  mutate(newgroup = ifelse(Day==0,"ctrl",ifelse(grepl("Cellulose", Substrate), "Cellulose", 
                                                ifelse(Substrate=="Oligomer","Oligomer","Other")))) %>% 
  mutate(newgroup2 = ifelse(Day==0, "t0", Substrate))

# Prepare phyloseq object
tax_tab = otu_tab %>% dplyr::select(ID, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames("ID") %>% as.matrix() %>% tax_table()
counts_table = otu_tab[c(1,9:41)] %>% column_to_rownames("ID")
sample_tab = sample_data(sampleinfo)

# Build phyloseq object
physeq_obj = phyloseq(otu_table(counts_table, taxa_are_rows = T), tax_tab, sample_tab)

# Subset Only the desired conditions I want from the sample info: 
#I want to compare Oligomer Day 90 vs. [Copolymer Day 90 + HDPE Day 90 + Blank Day 90]
dfenvsub1 <- df_env[df_env$Day == 90,]
dfenvsub2 <- dfenvsub1[dfenvsub1$Substrate == "Oligomer" | dfenvsub1$Substrate == "Copolymer" | dfenvsub1$Substrate == "HDPE" | dfenvsub1$Substrate == "Blank",] # results in 76

# Create a pseudo category to group together [Copolymer Day 90 + HDPE Day 90 + Blank Day 90]
dfenvsub2$Category <- 0
dfenvsub2$Category[c(1:9)] <- "No Biodegradation"
dfenvsub2$Category[c(10:12)] <- "Biodegradation"

# Load OTU table
df <- as.data.frame(physeq_obj@otu_table) 

# Subset only the desired conditions from the OTU table
dfsub2 = as.data.frame(df[,rownames(dfenvsub2)])

# Perform ALDEx2 analysis on this OTU table
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

# If you look at the aldex_out data frame, you can see that none of the Wilcoxin p values ("wi.eBH") are < 0.05. 
# The creators of Aldex2 recommend looking at the effect size ("effect") as well. I have multiple taxa that have an effect size >= 1,
# which is a cutoff they recommend
aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()    

# So if I subset out the taxa that have an effect size >=1 even if the p value is >0.05 and I re-run the analysis, many of the p values 
# of these taxa are now <+= 0.05. But I don't know if this subsetting would be considered valid. 
aldex_out_ed <- as.data.frame(aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  filter(abs(effect) >= 1)  %>%   # subsets taxa with effect size >= 1
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap))

# Subset OTU table data for the taxa with effect size >=1 from aldex_out
df_ed <- dfsub2[aldex_out_ed$Genus,]

# Re-run ALDEx2 analysis
drx_ed = aldex.clr(df_ed, dfenvsub2$Category)  

z_tt_ed <- aldex.ttest(drx_ed, paired.test = FALSE, verbose = TRUE)
z_effect_ed <- aldex.effect(drx_ed, CI = TRUE, verbose = TRUE)
aldex_out_ed <- data.frame(z_tt_ed, z_effected)
par(mfrow = c(1, 2))

aldex.plot(aldex_out_ed,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out_ed,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

# Now many of the wi.eBH values are >= 0.05. But is this okay?
aldex_out_ed2 <- as.data.frame(aldex_out_ed %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap))
