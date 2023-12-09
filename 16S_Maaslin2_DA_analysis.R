# C. Ahern
# 09/10/2023
# Using Maaslin2 for 16S differential abundance analysis 

BiocManager::install(version = "3.10")
BiocManager::install("Maaslin2")
install.packages("remotes")
remotes::install_github("david-barnett/microViz")
library(tidyverse)
library(readxl)
library(Maaslin2)
library(phyloseq)
library(microViz) # https://david-barnett.github.io/microViz
# library(ALDEx2)

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

# Maaslin2 analysis
mas_1 <- Maaslin2(
  input_data = dfsub2,
  input_metadata = dfenvsub2,
  output = "/Users/colleenahern/Downloads/checking_daa/colres", ### TOCHANGE
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("Category"),
  #reference = c("newgroup,ctrl"), 
  correction = "BH",
  standardize = FALSE,
  cores = 2)
mas_res_df <- mas_1$results

mas_res_df

# ####################### transformed data
# ### transform and aggregate 
# physeq_filt = physeq_taxfix %>% 
#   # keep only taxa belonging to genera that have over 100 counts in at least 10% of samples
#   tax_filter(min_prevalence = 0.1, undetected = 100, tax_level = "Genus") %>%
#   # aggregate counts at genus-level & transform with robust CLR transformation
#   tax_transform(trans = "rclr", rank = "Genus")
# ### maaslin with new data
# fit_data = Maaslin2(
#   input_data = data.frame(otu_table(physeq_filt)),
#   input_metadata = data.frame(sample_data(physeq_filt)),
#   output = "/home/mherold/Work/Various/Ahern_diffOTus/maaslin_genus_out",
#   normalization  = "NONE",
#   fixed_effects  = c("newgroup"),
#   reference      = c("newgroup,ctrl"))  
# fdr_mas2 <- fit_data$results %>%
#   dplyr::filter(qval < 0.05)  %>%
#   filter(value == "Oligomer") %>% 
#   arrange(desc(coef))
# 
# physeq_filt = physeq_taxfix %>% 
#   # keep only taxa belonging to genera that have over 100 counts in at least 10% of samples
#   tax_filter(min_prevalence = 0.1, undetected = 100, tax_level = "Genus") %>%
#   # aggregate counts at genus-level & transform with robust CLR transformation
#   tax_transform(trans = "rclr", rank = "Genus")
# ### maaslin with new data
# fit_data = Maaslin2(
#   input_data = data.frame(otu_table(physeq_filt)),
#   input_metadata = data.frame(sample_data(physeq_filt)),
#   output = "/home/mherold/Work/Various/Ahern_diffOTus/maaslin_genusSubsT0_out",
#   normalization  = "NONE",
#   fixed_effects  = c("newgroup2"),
#   reference      = c("newgroup2,t0"))  
# fdr_mas2_3 <- fit_data$results %>%
#   dplyr::filter(qval < 0.05)  %>%
#   # filter(value == "Oligomer" | value == "Copolymer") %>% 
#   filter(value == "Oligomer") %>% 
#   arrange(desc(coef))
# 
# 
# fit_data = Maaslin2(
#   input_data = data.frame(otu_table(physeq_filt)),
#   input_metadata = data.frame(sample_data(physeq_filt)),
#   output = "/home/mherold/Work/Various/Ahern_diffOTus/maaslin_genusSubsDay_out",
#   normalization  = "NONE",
#   fixed_effects  = c("Substrate", "Day"),
#   reference      = c("Substrate,Blank"))  
# fdr_mas2_5 <- fit_data$results %>%
#   dplyr::filter(qval < 0.05)  %>%
#   filter(value == "Oligomer") %>% 
#   arrange(desc(coef))
# 
# ########### try aldex
# # Generate Monte Carlo samples of the Dirichlet distribution foreach sample.↪
# # Convert each instance using the centered log-ratio transform.
# # This is the input for all further analyses.
# # set.seed(123)
# # x <- aldex.clr(data.frame(otu_table(physeq_taxfix)), conds=sampleinfo$newgroup)
# # x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)
# # x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# 
# ################ remove t=0 days
# physeq_filt_noctrl = physeq_taxfix %>% 
#   ps_filter(newgroup!="ctrl") %>%
#   # keep only taxa belonging to genera that have over 100 counts in at least 20% of samples
#   tax_filter(min_prevalence = 0.1, undetected = 100, tax_level = "Genus") %>%
#   # aggregate counts at genus-level & transform with robust CLR transformation
#   tax_transform(trans = "rclr", rank = "Genus")
# fit_data = Maaslin2(
#   input_data = data.frame(otu_table(physeq_filt_noctrl)),
#   input_metadata = data.frame(sample_data(physeq_filt_noctrl)),
#   output = "/home/mherold/Work/Various/Ahern_diffOTus/maaslin_genusNOCTRL_out",
#   min_abundance = 0.0,
#   min_prevalence = 0.0,
#   normalization  = "NONE",
#   fixed_effects  = c("newgroup"),
#   reference = c("newgroup,other")) 
# fit_data$results %>% head()
# fdr_mas3 <- fit_data$results %>%
#   dplyr::filter(qval < 0.05)  %>%
#   filter(value == "Oligomer")
# 
# ######### Individual abundances viz
# physeq_taxfix %>%
#   comp_barplot("Genus", n_taxa = 15, merge_other = FALSE, label = NULL) +
#   facet_wrap(vars(newgroup), scales = "free") + # scales = "free" is IMPORTANT!
#   coord_flip() +
#   theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))
# 
# #### selected genera
# featurelist = 
#   fdr_mas2_3 %>% # CHANGE HERE
#   filter(abs(coef) > 2) %>%
#   mutate(G_name = gsub("\\.\\.","\\: ", feature)) %>% 
#   mutate(G_name = gsub("\\.Genus"," Genus", G_name)) %>%
#   mutate(G_name = gsub("\\.","\\-", G_name)) %>%
#   pull(G_name)
# 
# 
# # featurelist = c("G: d__Bacteria Genus", featurelist) # add individual
# 
# plot_data <- physeq_taxfix %>%
#   tax_fix() %>%
#   tax_transform("compositional", rank = "Genus") %>%
#   # tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>%
#   ps_get() %>%
#   ps_otu2samdat( featurelist) %>% # adds abundance as sample data!
#   samdat_tbl()
# 
# plot_data  = plot_data %>%
#   # filter(Day==90) %>%
#   pivot_longer(featurelist, names_to="genus", values_to="abundance")
# 
# plot_data %>% 
#   ggplot(aes(x = newgroup2, y = abundance)) +
#   geom_boxplot(width = 0.5, colour = "grey35", outlier.shape = NA) +
#   geom_jitter(aes(colour = Substrate), width = 0.2, alpha = 0.5) +
#   # scale_y_continuous(
#   #   breaks = log2(1 / 2^(0:13)),
#   #   labels = function(x) paste0(100 * round(2^x, digits = 5), "%"),
#   #   limits = c(log2(0.00005), log2(0.25))
#   # ) +
#   theme_bw() +
#   # facet_grid(Day~ genus, scales="free")
#   facet_wrap(~ genus, scales="free")

