library(tidyverse)
library(data.table)
library(UpSetR)
library(wesanderson)
library(ggsignif)

select <- dplyr::select
rename <- dplyr::rename

# TPM Expression analysis of SD-98 genes in CORTECON dataset -------------------

## 1. Reading and pre-processing input files -----------------------------------

# Germinal zones human developing neocortex
fietz <- 
  fread("tpms/allSamples_Fietz_TPM.tsv") %>% 
  select(!c(median, mean, sd)) %>% 
  pivot_longer(-gene, names_to = "run", values_to = "counts") %>%
  mutate(rep = gsub('[[:digit:]]+', '', run)) %>%
  group_by(gene, rep) %>%
  summarize(median_tpm = median(counts))
  
# Cortical neural progenitor cells 
florio <- 
  fread("tpms/allSamples_Florio_TPM.tsv") %>% 
  select(!c(median, mean, sd)) %>%
  pivot_longer(-gene, names_to = "run", values_to = "counts") %>%
  mutate(rep = gsub('_rep[[:digit:]]+', '', run)) %>%
  group_by(gene, rep) %>%
  summarize(median_tpm = median(counts))

# SH-SY5Y neuroblastoma cells 
pezzini <- 
  fread("tpms/allSamples_Pezzini_TPM.tsv") %>%
  select(!c(median, mean, sd)) %>%
  pivot_longer(-gene, names_to = "run", values_to = "counts") %>%
  group_by(gene) %>%
  summarize(median_tpm = median(counts)) %>%
  rename("SH_SY5Y" = median_tpm)

# BrainSpan
brainspan <- 
  fread("tpms/allSamples_brainspan_TPM.tsv") %>%
  select(!c(median, mean, sd)) %>%
  pivot_longer(-gene, names_to = "run", values_to = "counts")

metadata_brainspan <- 
  read.table(file = "metadata/Brainspan_metadata.txt", sep="\t", fill=TRUE, 
             strip.white=TRUE, quote="", header=TRUE, na.strings=c("","NA"))

# only developmental samples
brainspan2 <- 
  left_join(brainspan, metadata_brainspan, by = c("run" = "run")) %>%
  filter(!(developmental_stage %in% c("Adulthood", "Adolescence", 
                                      "Early childhood", "Late childhood",
                                      "Early infancy", "Late infancy"))) %>%
  drop_na() %>%
  group_by(gene, developmental_stage) %>%
  summarize(median_tpm = median(counts))

# including adult samples
brainspan3 <-
  left_join(brainspan, metadata_brainspan, by = c("run" = "run")) %>%
  drop_na() %>%
  group_by(gene, developmental_stage) %>%
  summarize(median_tpm = median(counts))

# CORTECON
cortecon <- 
  fread("tpms/allSamples_vandeLeemput_TPM.tsv") %>%
  select(!c(median, mean, sd)) %>%
  pivot_longer(-gene, names_to = "run", values_to = "counts") %>%
  mutate(run = str_remove(run, "_1"))

metadata_cortecon <- 
  read.table(file = 'metadata/CORTECON_metadata.tsv', 
             sep = '\t', header = TRUE) %>%
  select(Run, RUN_ID, Stage) %>% 
  rename("run" = Run, "dataset" = RUN_ID, "condition" = Stage)

cortecon2 <- 
  left_join(cortecon, metadata_cortecon, by = c("run" = "run")) %>%
  filter(dataset == "Cortecon") %>%
  group_by(gene, condition) %>%
  summarize(median_tpm = median(counts))

fietz_all <- fietz %>% pivot_wider(values_from = median_tpm, names_from = rep)
florio_all <- florio %>% pivot_wider(values_from = median_tpm, names_from = rep)
brainspan2_all <- brainspan2 %>% pivot_wider(values_from = median_tpm, names_from = developmental_stage)
cortecon2_all <- cortecon2 %>% pivot_wider(values_from = median_tpm, names_from = condition)

list_df_all <- list(fietz_all, florio_all, pezzini, brainspan2_all, cortecon2_all)
all <- purrr::reduce(list_df_all, dplyr::left_join, by = c('gene')) # alternatively, we can use plyr::join_all
all <- all %>% mutate(number_expressed = rowSums(across(where(is.numeric), ~ . >= 1)))

bt <- fread("ID2Biotype.tsv", col.names = c("gene", "biotype"))
all_annot <- all %>% left_join(bt, by="gene")
all_filt <- all_annot %>% filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype))

dim(all) # 62216
dim(all_filt) # 23395
all_filt %>% filter(number_expressed>=1) %>% dim() # 18918

## 2. SD-98 genes --------------------------------------------------------------

# Selecting SD98 genes only

genes_sd98 <- read.table(file = 'metadata/sd98_genes.tsv', sep = '\t', header = TRUE)

fietz_sd98 <- 
  left_join(fietz, genes_sd98, by = c("gene" = "GeneID")) %>% 
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = rep)
fietz_sd98 <- fietz_sd98 %>% select("gene", "GeneName", "Biotype", "parCNClass", "VZ", "ISVZ", "OSVZ", "CP") # reordering

florio_sd98 <- left_join(florio, genes_sd98, by = c("gene" = "GeneID")) %>% 
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = rep)
florio_sd98 <- florio_sd98 %>% select("gene", "GeneName", "Biotype", "parCNClass", "aRG", "bRG", "N") # reordering

pezzini_sd98 <- left_join(pezzini, genes_sd98, by = c("gene" = "GeneID")) %>%
  drop_na()
pezzini_sd98 <- pezzini_sd98 %>% select("gene", "GeneName", "Biotype", "parCNClass", "SH_SY5Y") # reordering

brainspan2_sd98 <- left_join(brainspan2, genes_sd98, by = c("gene" = "GeneID")) %>%
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = developmental_stage)
brainspan2_sd98 <- brainspan2_sd98 %>% select("gene", "GeneName", "Biotype", "parCNClass", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal")

cortecon2_sd98 <- left_join(cortecon2, genes_sd98, by = c("gene" = "GeneID")) %>%
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = condition)
cortecon2_sd98 <- cortecon2_sd98 %>% select("gene", "GeneName", "Biotype", "parCNClass", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")

list_df <- list(fietz_sd98, florio_sd98, pezzini_sd98, brainspan2_sd98, cortecon2_sd98)
all_sd98 <- purrr::reduce(list_df, dplyr::left_join, by = c('gene', 'GeneName', 'Biotype', "parCNClass")) # alternatively, we can use plyr::join_all

fwrite(all_sd98, "results/sd98_genes_tpms.tsv", sep = "\t", quote = F, na = "-")

brainspan3_sd98 <- left_join(brainspan3, genes_sd98, by = c("gene" = "GeneID")) %>%
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = developmental_stage)
brainspan3_sd98 <- brainspan3_sd98 %>% select("gene", "GeneName", "Biotype", "parCNClass", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Early infancy", "Late infancy", "Early childhood", "Late childhood", "Adolescence", "Adulthood")

brainspan3_sd98 %>%
  mutate(across(where(is.numeric), ~ if_else(. < 1, "-", as.character(.)))) %>%
  fwrite("results/sd98_genes_tpms.full_brainspan.tsv", sep = "\t", quote = F, na = "-")

## 3. Bar plot SD98 genes ------------------------------------------------------

# All SD98 genes

# Note: 209 SD98 genes were present in the TPM files
genes_sd98[!genes_sd98$GeneID %in% fietz_sd98$gene,] %>% dim()

all_sd98_any <- all_sd98
all_sd98_any$Any <- apply(all_sd98[, 5:21], 1, max)

# All
all_sd98_any %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>%
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 900) +
  theme_minimal()
ggsave("plots/expressed.all.pdf")

# CN~2
all_sd98_any %>%
  filter(parCNClass != "Polymorphic") %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>% 
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 900) +
  theme_minimal()
ggsave("plots/expressed.all_cn2.pdf")

# Protein coding only

# All
all_sd98_any %>%
  filter(Biotype == "protein_coding") %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>%
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 900) +
  theme_minimal()
ggsave("plots/expressed.protein_coding.pdf")

# CN~2
all_sd98_any %>%
  filter(Biotype == "protein_coding") %>%
  filter(parCNClass != "Polymorphic") %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>% 
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 900) +
  theme_minimal()
ggsave("plots/expressed.protein_coding_cn2.pdf")

## 4. Human duplicated gene families -------------------------------------------

genes_hd <- read.table(file = 'metadata/human_duplicated_genes.tsv', sep = '\t', header = TRUE)

fietz_hd <- 
  left_join(fietz, genes_hd, by = c("gene" = "GeneID")) %>% 
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = rep)
fietz_hd <- fietz_hd %>% select("gene", "GeneName", "Biotype", "parCNClass", "VZ", "ISVZ", "OSVZ", "CP") # reordering

florio_hd <- left_join(florio, genes_hd, by = c("gene" = "GeneID")) %>% 
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = rep)
florio_hd <- florio_hd %>% select("gene", "GeneName", "Biotype", "parCNClass", "aRG", "bRG", "N") # reordering

pezzini_hd <- left_join(pezzini, genes_hd, by = c("gene" = "GeneID")) %>%
  drop_na()
pezzini_hd <- pezzini_hd %>% select("gene", "GeneName", "Biotype", "parCNClass", "SH_SY5Y") # reordering

brainspan2_hd <- left_join(brainspan2, genes_hd, by = c("gene" = "GeneID")) %>%
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = developmental_stage)
brainspan2_hd <- brainspan2_hd %>% select("gene", "GeneName", "Biotype", "parCNClass", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal")

cortecon2_hd <- left_join(cortecon2, genes_hd, by = c("gene" = "GeneID")) %>%
  drop_na() %>%
  pivot_wider(values_from = median_tpm, names_from = condition)
cortecon2_hd <- cortecon2_hd %>% select("gene", "GeneName", "Biotype", "parCNClass", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")

list_df <- list(fietz_hd, florio_hd, pezzini_hd, brainspan2_hd, cortecon2_hd)
all_hd <- purrr::reduce(list_df, dplyr::left_join, by = c('gene', 'GeneName', 'Biotype', "parCNClass")) # alternatively, we can use plyr::join_all

fwrite(all_hd, "results/human_duplicated_genes_tpms.tsv", sep = "\t", quote = F, na = "-")

# Brainspan for all stages for Table S1

## 5. Bar plot human duplicated genes ------------------------------------------ 

# All human duplicated genes

all_hd_any <- all_hd
all_hd_any$Any <- apply(all_hd[, 5:21], 1, max)

# All
all_hd_any %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>%
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 500) +
  theme_minimal()
ggsave("plots/human_duplicated.all.pdf")

# CN~2
all_hd_any %>%
  filter(parCNClass != "Polymorphic") %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>% 
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 500) +
  theme_minimal()
ggsave("plots/human_duplicated.all_cn_constrained.pdf")

# Protein coding only

# All
all_hd_any %>%
  filter(Biotype == "protein_coding") %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>%
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 500) +
  theme_minimal()
ggsave("plots/human_duplicated.protein_coding.pdf")

# CN~2
all_hd_any %>%
  filter(Biotype == "protein_coding") %>%
  filter(parCNClass != "Polymorphic") %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  group_by(name) %>% 
  mutate(name = factor(name, levels = rev(c("Any", "VZ", "ISVZ", "OSVZ", "CP", "aRG", "bRG", "N", "SH_SY5Y", "Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal", "Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")))) %>%
  summarize(expressed_genes = sum(value >= 1)) %>%
  ggplot(aes(x = expressed_genes, y = name)) +
  geom_bar(stat="identity") +
  xlim(0, 500) +
  theme_minimal()
ggsave("plots/human_duplicated.protein_coding_cn_constrained.pdf")

## 6. Violin plot --------------------------------------------------------------

# CORTECON all genes
all_hd %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  filter(name %in% c("Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")) %>%
  mutate(name = factor(name, levels = c("Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation"))) %>%
  filter(value > 1) %>%
  ggplot(aes(x = interaction(parCNClass, name), y = log2(value))) + 
  geom_violin(aes(fill = name), color = NA) +
  geom_boxplot(width = 0.25) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.1) +
  scale_fill_manual(values = wes_palette("AsteroidCity1", 5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("plots/volin_expressed_cortecon_all.pdf", width = 8, height = 5)

# CORTECON protein coding only
all_hd %>%
  filter(Biotype == "protein_coding") %>%
  pivot_longer(!c("gene", "GeneName", "Biotype", "parCNClass")) %>%
  filter(name %in% c("Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation")) %>%
  mutate(name = factor(name, levels = c("Pluripotency", "Neural Differentation", "Cortical Specification", "Deep Layer Formation", "Upper Layer Formation"))) %>%
  filter(value > 1) %>%
  ggplot(aes(x = interaction(parCNClass, name), y = log2(value))) + 
  geom_violin(aes(fill = name), color = NA) +
  geom_boxplot(width = 0.25) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.1) +
  scale_color_manual(values = wes_palette("AsteroidCity1", 5)) +
  scale_fill_manual(values = wes_palette("AsteroidCity1", 5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
ggsave("plots/volin_expressed_cortecon.pdf", width = 8, height = 5)

## 7. Upset plot ---------------------------------------------------------------

all_hd_binary <- 
  all_hd %>%
  select(!c(GeneName, Biotype, parCNClass)) %>%
  column_to_rownames("gene") %>%
  mutate(across(1:17, function(x) case_when(x >= 1 ~ 1, x < 1 ~ 0, is.na(x) ~ 0)))

pdf("plots/upset_plot.pdf", width = 10, height = 7)
upset(all_hd_binary, 
      order.by = "freq",
      nsets = 20, 
      nintersects = 40,
      keep.order = T, 
      sets = colnames(all_hd_binary))
dev.off()

all_hd_binary_filt <- 
  all_hd_binary %>%
  select(!Pluripotency)

pdf("plots/upset_plot_filt.pdf", width = 10, height = 7)
upset(all_hd_binary_filt, 
      order.by = "freq",
      nsets = 20, 
      nintersects = 40,
      keep.order = T, 
      sets = colnames(all_hd_binary_filt))
dev.off()

## 8. Expressed genes CN=2 -----------------------------------------------------

all_hd$expressed <- rowSums(all_hd[,5:21] >= 1)

# all
all_hd %>% dim() # 1584
all_hd %>% filter(expressed != 0) %>% dim() # 889
889/1584*100 # 56.12374

all_hd %>% filter(parCNClass != "Polymorphic") %>% dim() # 1166
all_hd %>% filter(parCNClass != "Polymorphic")%>%
  filter(expressed != 0) %>% dim() # 741
741/1166*100 # 63.5506

all_hd %>% filter(parCNClass == "Fixed") %>% dim() # 315
all_hd %>% filter(parCNClass == "Fixed")%>%
  filter(expressed != 0) %>% dim() # 248
248/315*100 # 78.73016

# protein coding
all_hd %>% filter(Biotype == "protein_coding") %>% dim() # 625
all_hd %>% filter(Biotype == "protein_coding") %>%
  filter(expressed != 0) %>% dim() # 436
436/625*100 # 69.76

all_hd %>% filter(Biotype == "protein_coding") %>% 
  filter(parCNClass != "Polymorphic") %>% dim() # 472
all_hd %>% filter(Biotype == "protein_coding") %>% 
  filter(parCNClass != "Polymorphic") %>%
  filter(expressed != 0) %>% dim() # 380
380/472*100 # 80.508476

all_hd %>% filter(Biotype == "protein_coding") %>% 
  filter(parCNClass == "Fixed") %>% dim() # 161
all_hd %>% filter(Biotype == "protein_coding") %>% 
  filter(parCNClass == "Fixed") %>%
  filter(expressed != 0) %>% dim() # 151
151/161*100 #  93.78882

# -> enrichment of brain genes among CN fixed SD98

## 9. Expression of interesting genes ------------------------------------------

# FRMPD2: CHM13_G0006295
# FRMPD2B: CHM13_G0006264

all_sd98 %>% 
  select(!c("GeneName", "Biotype", "parCNClass")) %>%
  filter(gene == "CHM13_G0006295" | gene == "CHM13_G0006264") %>%
  pivot_longer(-gene, names_to = "stage", values_to = "TPM") %>%
  mutate(stage = factor(stage, levels = colnames(all_sd98))) %>%
  ggplot(aes(x = stage, y = TPM, color = gene)) +
    geom_point(size = 2) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/FRMPD2.pdf", width = 5, height = 5)

all_sd98 %>% 
  select(!c("GeneName", "Biotype", "parCNClass")) %>%
  filter(gene == "CHM13_G0006295" | gene == "CHM13_G0006264") %>%
  pivot_longer(-gene, names_to = "stage", values_to = "TPM") %>%
  mutate(stage = factor(stage, levels = colnames(all_sd98))) %>%
  ggplot(aes(x = stage, y = TPM, color = gene)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## 10. Lymphoblastoid cell lines -----------------------------------------------

# reading Pickrell data
bt <- fread("ID2Biotype.tsv", col.names = c("gene", "biotype"))

pickrell <- 
  fread("tpms/allSamples_Pickrell_TPM.tsv") %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>%
  mutate(sample = str_remove(sample, "_.*")) %>%
  # consolidating multiple runs per sample
  group_by(gene, sample) %>%
  summarize(median_tpm = median(counts)) %>%
  # aggregating multiple samples per gene
  group_by(gene) %>%
  summarize(LCL = median(median_tpm)) #%>%
  # adding biotype information
  #left_join(bt, by="gene")

# selecting human-duplicated genes only
genes_hd <- read.table(file = 'metadata/human_duplicated_genes.tsv', sep = '\t', header = TRUE)

pickrell_hd <- 
  left_join(pickrell, genes_hd, by = c("gene" = "GeneID")) %>%
  drop_na()

genes_sd98 <- read.table(file = 'metadata/sd98_genes.tsv', sep = '\t', header = TRUE)

pickrell_sd98 <- 
  left_join(pickrell, genes_sd98, by = c("gene" = "GeneID")) %>%
  drop_na()

# combining with cortecon for plot
list_cortecon_pickrell <- list(cortecon2_hd, pickrell_hd)
cortecon_pickrell <- purrr::reduce(list_cortecon_pickrell, dplyr::left_join, by = c('gene', 'GeneName', 'Biotype', "parCNClass")) %>% drop_na()

# Protein coding CORTECON+Pickrell w/significance
my_comparisons <- list(
  c("Fixed", "Nearly-Fixed"),
  c("Fixed", "Polymorphic"),
  c("Nearly-Fixed", "Polymorphic")
)

df_for_plot <- cortecon_pickrell %>%
  pivot_longer(
    cols = -c(gene, GeneName, Biotype, parCNClass, ENSEMBL),
    names_to = "name",
    values_to = "val"
  ) %>%
  mutate(name = factor(name, levels = c("Pluripotency", 
                                        "Neural Differentation", 
                                        "Cortical Specification", 
                                        "Deep Layer Formation", 
                                        "Upper Layer Formation",
                                        "LCL"))) %>%
  filter(val > 1) %>%
  mutate(log_val = log2(val))

ggplot(df_for_plot, aes(x = parCNClass, y = log_val)) +
  geom_violin(aes(fill = name), color = NA) +
  geom_boxplot(width = 0.25, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.5) +
  facet_wrap(~name, ncol = 6) +
  scale_fill_manual(values = c(wes_palette("AsteroidCity1", 5), wes_palette("AsteroidCity3", 4)[4])) +
  theme_classic() +
  geom_signif(
    comparisons = my_comparisons,
    test = "wilcox.test",          # or "t.test"
    map_signif_level = TRUE,       # convert p to asterisks
    step_increase = 0.07,          # vertical spacing between brackets
    textsize = 3
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/volin_expressed_cortecon_plusLCL.pdf", width = 12, height = 5)

## 11. Brain expression enrichment ---------------------------------------------

### 11.1 SD98 genes enrichment -------------------------------------------------

N <- dim(all_filt)[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- all_filt %>% filter(number_expressed>=1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- dim(genes_sd98)[1] # size of drawn subset (total SD98 protein coding genes and unprocessed pseudogenes)
x <- all_sd98_any %>% filter(Any>=1) %>% dim() %>% .[1] # successes in the subset (total SD98 protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

### 11.2 Human duplicated genes enrichment -------------------------------------

N <- dim(all_filt)[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- all_filt %>% filter(number_expressed>=1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- dim(genes_hd)[1] # size of drawn subset (total HD protein coding genes and unprocessed pseudogenes)
x <- all_hd_any %>% filter(Any>=1) %>% dim() %>% .[1] # successes in the subset (total HD protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

### 11.3 Human duplicated CN constraint genes enrichment -----------------------

N <- dim(all_filt)[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- all_filt %>% filter(number_expressed>=1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- genes_hd %>% filter(parCNClass != "Polymorphic") %>% dim() %>% .[1] # size of drawn subset (total HD protein coding genes and unprocessed pseudogenes)
x <- all_hd_any %>% filter(parCNClass != "Polymorphic") %>% filter(Any>=1) %>% dim() %>% .[1] # successes in the subset (total HD protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

## 12. LCL enrichment ----------------------------------------------------------

bt <- fread("ID2Biotype.tsv", col.names = c("gene", "biotype"))

pickrell2 <- 
  fread("tpms/allSamples_Pickrell_TPM.tsv") %>% 
  dplyr::select(!c(median, mean, sd)) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>%
  mutate(sample = str_remove(sample, "_.*")) %>%
  # consolidating multiple runs per sample
  group_by(gene, sample) %>%
  summarize(median_tpm = median(counts)) %>%
  # aggregating multiple samples per gene
  group_by(gene) %>%
  summarize(LCL = median(median_tpm)) %>%
  # adding biotype information
  left_join(bt, by="gene")

### 12.1. SD98 genes enrichment ------------------------------------------------

N <- bt %>% filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>% dim() %>% .[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- pickrell2 %>% filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>% filter(LCL >= 1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- dim(genes_sd98)[1] # size of drawn subset (total SD98 protein coding genes and unprocessed pseudogenes)
x <- pickrell_sd98 %>% filter(LCL>=1) %>% dim() %>% .[1] # successes in the subset (total SD98 protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

### 12.2. Human duplicated genes enrichment ------------------------------------

N <- bt %>% filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>% dim() %>% .[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- pickrell2 %>% filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>% filter(LCL >= 1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- dim(genes_hd)[1] # size of drawn subset (total SD98 protein coding genes and unprocessed pseudogenes)
x <- pickrell_hd %>% filter(LCL>=1) %>% dim() %>% .[1] # successes in the subset (total SD98 protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

### 12.3. Human duplicated CN-constrained genes enrichment ---------------------

N <- bt %>% filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>% dim() %>% .[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- pickrell2 %>% filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>% filter(LCL >= 1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- genes_hd %>% filter(parCNClass != "Polymorphic") %>% dim() %>% .[1] # size of drawn subset (total HD protein coding genes and unprocessed pseudogenes)
x <- pickrell_hd %>% filter(parCNClass != "Polymorphic") %>% filter(LCL>=1) %>% dim() %>% .[1] # successes in the subset (total HD protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

## 13. SH-SY5Y enrichment  -----------------------------------------------------

### 13.1. SD98 genes enrichment ------------------------------------------------

N <- dim(all_filt)[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- all_filt %>% filter(SH_SY5Y>=1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- dim(genes_sd98)[1] # size of drawn subset (total SD98 protein coding genes and unprocessed pseudogenes)
x <- all_sd98_any %>% filter(SH_SY5Y>=1) %>% dim() %>% .[1] # successes in the subset (total SD98 protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

### 13.2. Human_duplicated genes enrichment ------------------------------------

N <- dim(all_filt)[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- all_filt %>% filter(SH_SY5Y>=1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- dim(genes_hd)[1] # size of drawn subset (total SD98 protein coding genes and unprocessed pseudogenes)
x <- all_hd_any %>% filter(SH_SY5Y>=1) %>% dim() %>% .[1] # successes in the subset (total SD98 protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value

### 13.3. Human_duplicated CN-constrained genes enrichment ---------------------

N <- dim(all_filt)[1] # total population size (total protein coding genes and unprocessed pseudogenes)
K <- all_filt %>% filter(SH_SY5Y>=1) %>% dim() %>% .[1] # total successes in population (total protein coding genes and unprocessed pseudogenes expressed)
n <- genes_hd %>% filter(parCNClass != "Polymorphic") %>% dim() %>% .[1] # size of drawn subset (total HD protein coding genes and unprocessed pseudogenes)
x <- all_hd_any %>% filter(parCNClass != "Polymorphic") %>% filter(SH_SY5Y>=1) %>% dim() %>% .[1] # successes in the subset (total SD98 protein coding genes and unprocessed pseudogenes expressed)

p_value <- phyper(
  q = x - 1,   # number of observed successes minus 1
  m = K,       # total successes in population
  n = N - K,   # total failures in population
  k = n,       # size of drawn subset
  lower.tail = TRUE, log.p = FALSE
)

p_value


## 14. Tissue expression comparison --------------------------------------------

# brain datasets
brain_counts <- all_filt %>% 
  dplyr::select(!number_expressed) %>%
  filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>%
  mutate(human_duplicated = ifelse(gene %in% genes_hd$GeneID, "yes", "no")) %>%
  pivot_longer(!c(gene, human_duplicated, biotype), names_to = "dataset", values_to = "tpm") %>%
  group_by(dataset) %>%
  summarize(all_expressed = sum(tpm >= 1), hd_expressed = sum(tpm >= 1 & human_duplicated == "yes")) %>%
  mutate(tissue = "Brain")

# lcl datasets
lcl_counts <- 
  fread("tpms/allSamples_Pickrell_TPM.tsv") %>% 
  dplyr::select(!c(median, mean, sd)) %>%
  left_join(bt, by="gene") %>%
  filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>%
  pivot_longer(-c(gene, biotype), names_to = "dataset", values_to = "counts") %>%
  mutate(dataset = str_remove(dataset, "_.*")) %>%
  # consolidating multiple runs per sample
  group_by(gene, dataset) %>%
  summarize(tpm = median(counts)) %>%
  mutate(human_duplicated = ifelse(gene %in% genes_hd$GeneID, "yes", "no")) %>%
  group_by(dataset) %>%
  summarize(all_expressed = sum(tpm >= 1), hd_expressed = sum(tpm >= 1 & human_duplicated == "yes")) %>%
  mutate(tissue = "LCL")

# other tissues
other_metadata <- fread("metadata/Blake_metadata.tsv", header = FALSE, col.names = c("dataset", "tissue_id", "tissue"))
other_counts <-
  fread("tpms/allSamples_blake_TPM.tsv") %>%
  dplyr::select(!c(median, mean, sd)) %>%
  left_join(bt, by="gene") %>%
  filter(biotype == "protein_coding" | grepl("unprocessed_pseudogene", biotype)) %>%
  pivot_longer(-c(gene, biotype), names_to = "dataset", values_to = "tpm") %>%
  left_join(other_metadata, by = "dataset") %>%
  mutate(human_duplicated = ifelse(gene %in% genes_hd$GeneID, "yes", "no")) %>%
  group_by(dataset, tissue) %>%
  summarize(all_expressed = sum(tpm >= 1), hd_expressed = sum(tpm >= 1 & human_duplicated == "yes")) %>%
  dplyr::select(dataset, all_expressed, hd_expressed, tissue)

combined_counts <- rbind(brain_counts, lcl_counts, other_counts)

# plot
library(ggrepel)
combined_counts %>%
  mutate(dataset = ifelse(tissue == "Brain", dataset, NA)) %>%
  pivot_longer(!c(dataset, tissue), names_to = "gene_group", values_to = "expressed") %>%
  ggplot(aes(y = expressed, x = tissue)) + 
    geom_boxplot(outlier.shape = NA, width = 0.8) + 
    #geom_label_repel(aes(label = dataset), max.overlaps = 1000, force = 20, nudge_y = 1) +
    geom_jitter(width = 0.25) + 
    #scale_x_discrete(labels = c("Brain", "Heart", "Kidney", "LCL", "Liver", "Lung")) + 
    facet_wrap(~gene_group, 
               scale = "free_y", 
               labeller = as_labeller(c("all_expressed" = "Protein coding + unprocessed pseudogenes", "hd_expressed" = "Human-duplicated genes"))
               ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 14, face = "bold")
    ) + 
  xlab("Tissues") +
  ylab("No. genes expressed (TPM>1)")
ggsave("boxplot_expression.pdf", width = 10)

# Brainspan Geschwind WGCNA genes ----------------------------------------------

# Geschwind filtered genes that were present at an RPKM of 1 in 80% of the 
# samples from at least one neocortical region at one major temporal epoch, 
# resulting in 22,084 coding and non-coding transcripts.

# We filtered Brainspain genes focusing on Geshwind samples using only prenatal
# samples from the frontal cortex (DFC + VFC + MFC + OFC)

brainspan <-
  fread("tpms/allSamples_brainspan_TPM.tsv") %>%
  select(!c(median, mean, sd)) %>%
  pivot_longer(-gene, names_to = "run", values_to = "counts")

metadata_brainspan <- 
  read.table(file = "metadata/Brainspan_metadata.txt", sep="\t", fill=TRUE, 
             strip.white=TRUE, quote="", header=TRUE, na.strings=c("","NA"))

brainspan_wgcna <-
  left_join(brainspan, metadata_brainspan, by = c("run" = "run")) %>%
  filter(geschwind == "yes") %>%
  filter(body_site %in% c("DFC", "VFC", "MFC", "OFC")) %>%
  filter(developmental_stage %in% c("Early prenatal", "Early mid-prenatal", "Late mid-prenatal")) %>%
  dplyr::select(gene, run, counts, developmental_stage, body_site) 

genes_wgcna <- brainspan_wgcna %>%
  group_by(gene, developmental_stage, body_site) %>%
  summarize(percent_over_1 = sum(counts > 1)/n()) %>%
  filter(percent_over_1 >= 0.8) %>%
  pull(gene) %>% unique()

write(genes_wgcna, "../03_WGCNA/BrainSpan/brainspan_wgcna_genes.txt", ncolumns = 1)

# FRMP2 ------------------------------------------------------------------------

library(janitor)

tpms <- read.table(file = 'TPMs/allSamples_brainspan_TPM.tsv', sep = '\t', header = TRUE)
tpms <- t(tpms)
tpms2 <- row_to_names(tpms, row_number = 1)
tpms2 <- as.data.frame(tpms2)

frmpd2_exp <- tpms2['CHM13_G0006295'] %>% 
  rownames_to_column('run') %>% 
  mutate(run = factor(run, levels = order)) %>%
  mutate(CHM13_G0006295 = as.numeric(CHM13_G0006295))
frmpd2_exp <- left_join(frmpd2_exp, traits, by = "run") %>%
  filter(developmental_stage != "Adulthood" & 
           developmental_stage != "Adolescence" &
           developmental_stage != "Early childhood" &
           developmental_stage != "Late childhood") %>% 
  mutate(run = factor(run, levels = order)) 

ggplot(frmpd2_exp, aes(x = run, y = CHM13_G0006295, color = developmental_stage)) +
  geom_point() + facet_wrap(~ site) +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())


fam72_exp <- tpms2['CHM13_G0004511'] %>% 
  rownames_to_column('run') %>% 
  mutate(run = factor(run, levels = order)) %>%
  mutate(CHM13_G0006295 = as.numeric(CHM13_G0004511))
fam72_exp <- left_join(fam72_exp, traits, by = "run") %>%
  filter(developmental_stage != "Adulthood" & 
           developmental_stage != "Adolescence" &
           developmental_stage != "Early childhood" &
           developmental_stage != "Late childhood") %>% 
  mutate(run = factor(run, levels = order)) 

ggplot(fam72_exp, aes(x = run, y = CHM13_G0006295, color = developmental_stage)) +
  geom_point() + facet_wrap(~ site) +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

frmpd2_exp <- norm.counts['CHM13_G0006295'] %>% 
  rownames_to_column('run') %>% 
  mutate(run = factor(run, levels = order)) 
frmpd2_exp <- left_join(frmpd2_exp, traits, by = "run") %>%
  filter(developmental_stage != "Adulthood" & 
           developmental_stage != "Adolescence" &
           developmental_stage != "Early childhood" &
           developmental_stage != "Late childhood") %>% 
  mutate(run = factor(run, levels = order)) 

ggplot(frmpd2_exp, aes(x = run, y = CHM13_G0006295, color = developmental_stage)) +
  geom_point() + facet_wrap(~ site) +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())
