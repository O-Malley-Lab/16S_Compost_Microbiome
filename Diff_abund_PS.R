# Differential abundance - Pat Schloss
# C. Ahern
# 08/23/2023

library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)

set.seed(19760620)

# metadata <- read_excel("raw_data/schubert.metadata.xlsx", na="NA") %>%
#   select(sample_id, disease_stat) %>%
#   drop_na(disease_stat)

metadata <- read_excel("/Users/colleenahern/Documents/BASF/minimalR-raw_data-master/schubert.metadata.xlsx", na="NA") %>%
  select(sample_id, disease_stat) %>%
  drop_na(disease_stat)

metadataca <- read.delim("/Users/colleenahern/Documents/BASF/qiime2/metadata2.txt") %>%
  rename("sample_id" = "sample.id")


otu_counts <- read_delim("/Users/colleenahern/Downloads/minimalR-raw_data-0.3/schubert.subsample.shared.txt") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample_id = Group) %>%
  pivot_longer(-sample_id, names_to="otu", values_to = "count")


otu_counts <- read_tsv("/Users/colleenahern/Documents/BASF/minimalR-raw_data-master/schubert.subsample.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample_id = Group) %>%
  pivot_longer(-sample_id, names_to="otu", values_to = "count")

setwd("/Users/colleenahern/Documents/BASF/qiime2/qiime2downloads")

library(qiime2R)
psme <- qza_to_phyloseq(features="Analysis2/16S-table-noplant-rarefied-10000_filtered.qza",
                        tree="Analysis2/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "Analysis2/16S-rep-seqs-taxonomy.qza", 
                        metadata= "Analysis2/metadata2.txt")

otu_countsca <- as.matrix(psme@otu_table)
otu_countsca <- as.data.frame(otu_countsca@.Data)
colnames(otu_countsca) <- gsub("\\-16S.*", "", colnames(psme@otu_table))
otu_countsca$otu <- paste0("otu", 1:nrow(otu_countsca))
otu_countsca <- otu_countsca %>%
  select(otu, everything()) %>%
  gather(sample_id, count, 2:34) %>%
  select(sample_id, everything())

nseqs_per_sample <- otu_counts %>%
  group_by(sample_id) %>%
  summarize(N = sum(count), .groups = "drop") %>%
  count(N) %>%
  pull(N)

nseqs_per_sampleca <- otu_countsca %>%
  group_by(sample_id) %>%
  summarize(N = sum(count), .groups = "drop") %>%
  count(N) %>%
  pull(N)

nseqs_per_sampleca <- min(nseqs_per_sampleca)
  
stopifnot(length(nseqs_per_sample) == 1)
stopifnot(length(nseqs_per_sampleca) == 1)

lod <- 100* 1/nseqs_per_sample
lodca <- 100* 1/nseqs_per_sampleca

# taxonomy <- read_tsv("raw_data/schubert.cons.taxonomy") %>%
taxonomy <- read_tsv("/Users/colleenahern/Documents/BASF/minimalR-raw_data-master/schubert.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")

taxonomyca <- psme@tax_table
taxonomyca <- as.data.frame(taxonomyca@.Data)
taxonomyca$otu <- paste0("otu", 1:nrow(otu_countsca))  # MAKE SURE YOUR OTU NUMBERS MATCH UP WITH THE ORIGINAL COMPLICATED NAMES
taxonomyca <- taxonomyca %>%
  select(otu, everything()) %>%
  rename_all(tolower)
  
otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(
    c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
    names_to="level",
    values_to="taxon") #%>%
  mutate(disease_stat = factor(disease_stat,
                               levels=c("Case",
                                        "DiarrhealControl",
                                        "NonDiarrhealControl")))

metadataca$sample_id <- gsub("\\-16S.*", "", metadataca$sample_id)

otu_rel_abundca <- inner_join(metadataca, otu_countsca, by="sample_id") %>%
  inner_join(., taxonomyca, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(
    c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
    names_to="level",
    values_to="taxon") %>%
  mutate(Group = factor(Group,
                               levels=c("Blank Day 0",
                                        "Cellulose Day 0",
                                        "Copolymer Day 0",
                                        "HDPE Day 0",
                                        "Oligomer Day 0",
                                        "Big Cellulose Day 90",
                                        "Blank Day 90",
                                        "Cellulose Day 90",
                                        "Copolymer Day 90",
                                        "HDPE Day 90",
                                        "Oligomer Day 90")))


taxon_rel_abund <- otu_rel_abund %>%
  filter(level=="genus") %>%
  group_by(disease_stat, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified<br>*\\1*"),
         taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"),
         taxon = str_replace_all(taxon,
                                 "_", " "))

taxon_rel_abundca <- otu_rel_abundca %>%
  filter(level=="genus") %>%
  group_by(Group, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified<br>*\\1*"),
         taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"),
         taxon = str_replace_all(taxon,
                                 "_", " "))

library(ggtext)

taxon_pool <- taxon_rel_abund %>%
  group_by(disease_stat, taxon) %>%
  summarize(median=median(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(median) < 1,
            median = median(median),
            .groups="drop")

taxon_poolca <- taxon_rel_abundca %>%
  group_by(Group, taxon) %>%
  summarize(median=median(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(median) < 1,
            median = median(median),
            .groups="drop")

inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, disease_stat, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            median = min(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE)) %>%
  mutate(rel_abund = if_else(rel_abund==0,
                            2/3 * lod,
                            rel_abund)) %>%
  ggplot(aes(y=taxon, x=rel_abund, color=disease_stat)) +
  geom_vline(xintercept = lod, size = 0.2) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.6)) +
  # scale_x_log10() +
  coord_trans(x = "log10") +
  scale_x_continuous(limits = c(NA, 100),
                     breaks = c(0.1, 1, 10, 100),
                     labels = c(0.1, 1, 10, 100)) +
  scale_color_manual(name=NULL,
                     breaks=c("NonDiarrhealControl",
                              "DiarrhealControl",
                              "Case"),
                     labels=c("Healthy",
                              "Diarrhea,<br>*C. difficile* negative",
                              "Diarrhea,<br>*C. difficile* positive"),
                     values=c("gray", "blue", "red")) +
  labs(y=NULL,
       x="Relative  Abundance (%)") +
  theme_classic() +
  theme(axis.text.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown(),
        # legend.position = c(0.8, 0.6),
        legend.background = element_rect(color="black", fill = NA),
        legend.margin = margin(t=-5, r=3, b=3)
  )


inner_join(taxon_rel_abundca, taxon_poolca, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, Group, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            median = min(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE)) %>%
  mutate(rel_abund = if_else(rel_abund==0,
                             2/3 * lodca,
                             rel_abund)) %>%
  ggplot(aes(y=taxon, x=rel_abund, color=Group)) +
  geom_vline(xintercept = lodca, size = 0.2) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.6)) +
  # scale_x_log10() +
  coord_trans(x = "log10") +
  scale_x_continuous(limits = c(NA, 100),
                     breaks = c(0.1, 1, 10, 100),
                     labels = c(0.1, 1, 10, 100)) +
  scale_color_manual(name=NULL,
                     breaks=c("Blank Day 0",
                              "Cellulose Day 0",
                              "Copolymer Day 0",
                              "HDPE Day 0",
                              "Oligomer Day 0",
                              "Big Cellulose Day 90",
                              "Blank Day 90",
                              "Cellulose Day 90",
                              "Copolymer Day 90",
                              "HDPE Day 90",
                              "Oligomer Day 90"),
                     labels=c("Blank Day 0",
                              "Cellulose Day 0",
                              "Copolymer Day 0",
                              "HDPE Day 0",
                              "Oligomer Day 0",
                              "Big Cellulose Day 90",
                              "Blank Day 90",
                              "Cellulose Day 90",
                              "Copolymer Day 90",
                              "HDPE Day 90",
                              "Oligomer Day 90"),
                     values=c("coral1","goldenrod3", "gold","darkgreen", "green", "turquoise4", "turquoise", "blue","deepskyblue", "deeppink2", "deeppink4")) +
  labs(y=NULL,
       x="Relative  Abundance (%)") +
  theme_classic() +
  theme(axis.text.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown(),
        # legend.position = c(0.8, 0.6),
        legend.background = element_rect(color="black", fill = NA),
        legend.margin = margin(t=-5, r=3, b=3)
  )


# my data but taking out all day 0 data

metadataca90 <- metadataca[!grepl("Day 0", metadataca$Group),]

otu_countsca <- as.matrix(psme@otu_table)
otu_countsca <- as.data.frame(otu_countsca@.Data)
colnames(otu_countsca) <- gsub("\\-16S.*", "", colnames(psme@otu_table))
otu_countsca$otu <- paste0("otu", 1:nrow(otu_countsca))
otu_countsca <- otu_countsca %>%
  select(otu, everything()) %>%
  gather(sample_id, count, 2:34) %>%
  select(sample_id, everything())

otu_countsca90 <- otu_countsca[!grepl("t0", otu_countsca$sample_id),]

nseqs_per_sampleca90 <- otu_countsca90 %>%
  group_by(sample_id) %>%
  summarize(N = sum(count), .groups = "drop") %>%
  count(N) %>%
  pull(N)

nseqs_per_sampleca90 <- min(nseqs_per_sampleca90)

stopifnot(length(nseqs_per_sampleca90) == 1)

lodca90 <- 100* 1/nseqs_per_sampleca90

taxonomyca90 <- taxonomyca

otu_rel_abundca90 <- inner_join(metadataca90, otu_countsca90, by="sample_id") %>%
  inner_join(., taxonomyca90, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(
    c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
    names_to="level",
    values_to="taxon") %>%
  mutate(Group = factor(Group,
                        levels=c("Big Cellulose Day 90",
                                 "Blank Day 90",
                                 "Cellulose Day 90",
                                 "Copolymer Day 90",
                                 "HDPE Day 90",
                                 "Oligomer Day 90")))


taxon_rel_abundca90 <- otu_rel_abundca90 %>%
  filter(level=="genus") %>%
  group_by(Group, sample_id, taxon) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "(.*)_unclassified", "Unclassified<br>*\\1*"),
         taxon = str_replace(taxon,
                             "^([^<]*)$", "*\\1*"),
         taxon = str_replace_all(taxon,
                                 "_", " "))

library(ggtext)

taxon_poolca90 <- taxon_rel_abundca90 %>%
  group_by(Group, taxon) %>%
  summarize(median=median(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(median) < 2,
            median = median(median),
            .groups="drop")

inner_join(taxon_rel_abundca90, taxon_poolca90, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, Group, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            median = min(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE)) %>%
  mutate(rel_abund = if_else(rel_abund==0,
                             2/3 * lodca90,
                             rel_abund)) %>%
  ggplot(aes(y=taxon, x=rel_abund, color=Group)) +
  geom_vline(xintercept = lodca, size = 0.2) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.6)) +
  # scale_x_log10() +
  coord_trans(x = "log10") +
  scale_x_continuous(limits = c(NA, 100),
                     breaks = c(0.1, 1, 10, 100),
                     labels = c(0.1, 1, 10, 100)) +
  scale_color_manual(name=NULL,
                     breaks=c( "Big Cellulose Day 90",
                              "Blank Day 90",
                              "Cellulose Day 90",
                              "Copolymer Day 90",
                              "HDPE Day 90",
                              "Oligomer Day 90"),
                     labels=c("Big Cellulose Day 90",
                              "Blank Day 90",
                              "Cellulose Day 90",
                              "Copolymer Day 90",
                              "HDPE Day 90",
                              "Oligomer Day 90"),
                     values=c("coral1","goldenrod3", "green", "turquoise", "blue","deeppink2")) +
  labs(y=NULL,
       x="Relative  Abundance (%)") +
  theme_classic() +
  theme(axis.text.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown(),
        # legend.position = c(0.8, 0.6),
        legend.background = element_rect(color="black", fill = NA),
        legend.margin = margin(t=-5, r=3, b=3)
  )

# ggsave("schubert_genus.tiff", width=5, height=6).
