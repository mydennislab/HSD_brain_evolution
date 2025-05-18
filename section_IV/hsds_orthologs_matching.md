# Analysis of expression patterns of HSD genes in model organisms (mouse and zebrafish)

# STEP 0: load libraries and color palettes
```{r}
library(Seurat)
library(SeuratDisk)
library(stringr)
library(tidyverse)
library(GenomicFeatures)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(preprocessCore)
library(ggbiplot)
```

# STEP 1: gather gene count matrices per developmental timepoint in each model organisms.

## Mouse:
```{r}
# Source of data: PRJNA637987 (LaManno et al 2021 - Nature) [data focused on mouse brain development]

# Convert from loom to Seurat:
mouse.loom <- Connect(filename="mouse_data/dev_all.loom", mode="r")
#mouse.sc <- as.Seurat(mouse.loom) # this is giving an error so will extract manually
mouse.counts <- mouse.loom[["matrix"]][,] # full gene x cell count table
dim(mouse.counts) # 292495 cells x 31053 genes
colnames(mouse.counts) <- make.unique(mouse.loom[["row_attrs"]][["Gene"]][],sep = "-") # gene names as columns
rownames(mouse.counts) <- make.unique(mouse.loom[["col_attrs"]][["CellID"]][],sep="-") # cell ID as rows
mouse.sc <- CreateSeuratObject(counts= t(mouse.counts), project = "mouse")
mouse.sc$Age <- mouse.loom[["col_attrs"]][["Age"]][]
mouse.sc$Cluster <- mouse.loom[["col_attrs"]][["ClusterName"]][]
mouse.sc$Sample <- mouse.loom[["col_attrs"]][["SampleID"]][]
# names(x = mouse.loom[['col_attrs']]) # this trick helps to know which other attributes can be found in the loom file (same as the meta.data in Seurat)
mouse.loom$close_all()

# Aggregate expression by timepoint:
# Timepoints: e7.0, e8.0, e8.5, e9.0, e10.0, e11.0, e12.0, e12.5, e13.0, e13.5, e14.0, e14.5, e15.0, e15.5, e16.0, e16.25, e16.5, e17.0, e17.5, e18.0
mouse.sc$Age <- factor(mouse.sc$Age, levels = c("e7.0", "e8.0", "e8.5", "e9.0", "e10.0", "e11.0", "e12.0", "e12.5", "e13.0", "e13.5", "e14.0", "e14.5", "e15.0", "e15.5", "e16.0", "e16.25", "e16.5", "e17.0", "e17.5", "e18.0")) # order timepoints

mouse.sc.pseudo <- AggregateExpression(mouse.sc, assays="RNA", group.by=c("Age", "Sample"), return.seurat = F)
mouse.sc.pseudo.df <- as.data.frame(mouse.sc.pseudo$RNA)
write.csv(mouse.sc.pseudo.df, file="mouse_data/mouse.sc.pseudo.df.csv")

mouse.metadata <- data.frame(FullID= colnames(mouse.sc.pseudo.df),
                             Age= str_extract(colnames(mouse.sc.pseudo.df), "^[^_]+")
                             )
write.csv(mouse.metadata, file="mouse_data/mouse.metadata.csv", row.names = F)
```

## Zebrafish:
```{r}
# Source of data: GSE158142 (Raj et al. 2020 - Neuron) [data focused on zebrafish brain development]

# Load each Seurat object per timepoint:
# Timepoints: 12 hpf (6s), 14 hpf (10s), 16 hpf (14s), 18 hpf (18s), 20 hpf, 24 hpf, 36 hpf, 2 dpf, 3 dpf, 5 dpf, 8 dpf, 15 dpf.
z12hpf <- readRDS("zfish_data/GSE158142_zf6s_cc_filt.cluster.rds")
z14hpf <- readRDS("zfish_data/GSE158142_zf10s_cc_filt.cluster.rds")
z16hpf <- readRDS("zfish_data/GSE158142_zf14s_cc_filt.cluster.rds")
z18hpf <- readRDS("zfish_data/GSE158142_zf18s_cc_filt.cluster.rds")
z20hpf <- readRDS("zfish_data/GSE158142_zf20hpf_cc_filt.cluster.rds")
z24hpf <- readRDS("zfish_data/GSE158142_zf24hpf_cc_filt.cluster.rds")
z36hpf <- readRDS("zfish_data/GSE158142_zf36hpf_cc_filt.cluster.rds")
z2dpf <- readRDS("zfish_data/GSE158142_zf2dpf_cc_filt.cluster.rds")
z3dpf <- readRDS("zfish_data/GSE158142_zf3dpf_cc_filt.cluster.rds")
z5dpf <- readRDS("zfish_data/GSE158142_zf5dpf_cc_filt.cluster.rds")
z8dpf <- readRDS("zfish_data/GSE158142_zf8dpf_cc_filt.cluster4.rds")
z15dpf <- readRDS("zfish_data/GSE158142_zf15dpf_PCAALL.rds")

z12hpf <- UpdateSeuratObject(z12hpf)
z14hpf <- UpdateSeuratObject(z14hpf)
z16hpf <- UpdateSeuratObject(z16hpf)
z18hpf <- UpdateSeuratObject(z18hpf)
z20hpf <- UpdateSeuratObject(z20hpf)
z24hpf <- UpdateSeuratObject(z24hpf)
z36hpf <- UpdateSeuratObject(z36hpf)
z2dpf <- UpdateSeuratObject(z2dpf)
z3dpf <- UpdateSeuratObject(z3dpf)
z5dpf <- UpdateSeuratObject(z5dpf)
z8dpf <- UpdateSeuratObject(z8dpf)
z15dpf <- UpdateSeuratObject(z15dpf)

# Merge into one big seurat object:
zfish.sc <- merge(z12hpf, y=c(z14hpf,z16hpf,z18hpf,z20hpf,z24hpf,z36hpf,z2dpf,z3dpf,z5dpf,z8dpf,z15dpf))

# Add timepoint variable to metadata:
split_names <- strsplit(zfish.sc$orig.ident, "-")
zfish.sc$Age <- sapply(split_names, function(x) x[[2]])
zfish.sc$Age <- factor(zfish.sc$Age, levels=c("12h", "14h", "16h", "18h", "20h", "24h", "36h", "2d", "3d", "5d", "8d", "15d"))

# Aggregate expression by timepoint:
zfish.sc.pseudo <- AggregateExpression(zfish.sc, assays="RNA", group.by=c("Age", "orig.ident"), return.seurat = F)
zfish.sc.pseudo.df <- as.data.frame(zfish.sc.pseudo$RNA)
write.csv(zfish.sc.pseudo.df, file="zfish_data/zfish.sc.pseudo.df.csv")

zfish.metadata <- data.frame(FullID= colnames(zfish.sc.pseudo.df),
                             Age= str_extract(colnames(zfish.sc.pseudo.df), "^[^_]+")
                             )
write.csv(zfish.metadata, file="zfish_data/zfish.metadata.csv", row.names = F)

```

## Add human orthologs to matrices:
```{r}
# Load human ortholog lists:
human_mouse_ortho <- read.csv("files/human_mouse_orthologs_all.csv")
human_zfish_ortho <- read.csv("files/human_zfish_orthologs_all.csv")

# Merge 
mouse.sc.pseudo.df.merged <- merge(mouse.sc.pseudo.df, human_mouse_ortho, by.x = "row.names", by.y = "Mouse_GeneName", all.x = TRUE)
mouse.sc.pseudo.df.merged <- mouse.sc.pseudo.df.merged[complete.cases(mouse.sc.pseudo.df.merged$Human_GeneName), ] # 19,949 gene with 1:1 ortholog kept
rownames(mouse.sc.pseudo.df.merged) <- make.unique(mouse.sc.pseudo.df.merged$Row.names, sep="-")
mouse.sc.pseudo.df.merged <- mouse.sc.pseudo.df.merged[,-1]

zfish.sc.pseudo.df.merged <- merge(zfish.sc.pseudo.df, human_zfish_ortho, by.x = "row.names", by.y = "Zebrafish_GeneName", all.x = TRUE)
zfish.sc.pseudo.df.merged <- zfish.sc.pseudo.df.merged[complete.cases(zfish.sc.pseudo.df.merged$Human_GeneName), ] # 16,910 gene with 1:1 ortholog kept
rownames(zfish.sc.pseudo.df.merged) <- make.unique(zfish.sc.pseudo.df.merged$Row.names, sep="-")
zfish.sc.pseudo.df.merged <- zfish.sc.pseudo.df.merged[,-1]
```

# STEP 2: Match developmental age between model organisms and human using the BrainSpan data.

## Load functions needed:
```{r}
paper.pca <- function(g1, pc.1=1, pc.2=3, my.title){
  g1 <- g1 + geom_point(size = 3) + 
    coord_fixed() + 
    xlab(paste0("PC", pc.1, ": ", round(percentVar[pc.1] * 100), "% variance")) + 
    ylab(paste0("PC", pc.2, ": ", round(percentVar[pc.2] * 100), "% variance")) + 
    ggtitle(my.title) + NULL
  return(g1)
}

paper.ggbiplot2 <- function (pcobj,
                              choices = 1:2, scale = 1, pc.biplot = TRUE, 
                              obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                              ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                              alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                              varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE) {
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else {
    stop("Expected a object of class prcomp")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  df.u <- df.u * nobs.factor
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = muted("red"))
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar, angle = angle, hjust = hjust), 
                       color = "darkred", size = varname.size)
  }
  return(df.u)
}

```

## Obtain RPKM values from mouse/zfish.
```{r}
## Mouse:
mouse.metadata$Age <- factor(mouse.metadata$Age, levels=c("e7.0", "e8.0", "e8.5", "e9.0", "e10.0", "e11.0", "e12.0", "e12.5", "e13.0", "e13.5", "e14.0", "e14.5", "e15.0", "e15.5", "e16.0", "e16.25", "e16.5", "e17.0", "e17.5", "e18.0"))
mouse_timepoints <- levels(mouse.metadata$Age)
mouse.samples <- mouse.metadata$FullID
#mouse.samples <- gsub("-(?=[^-]*$)", ".", mouse.samples, perl = TRUE) # change the - for a . to match colnames of dataframe

# Get CPM values:
mouse.cpm <- cpm(mouse.sc.pseudo.df.merged[,colnames(mouse.sc.pseudo.df.merged) %in% mouse.samples], log = FALSE) %>% as.data.frame()
mouse.cpm$Human_GeneName <- mouse.sc.pseudo.df.merged$Human_GeneName

## Zebrafish:
zfish.metadata$Age <- factor(zfish.metadata$Age, levels=c("12h", "14h", "16h", "18h", "20h", "24h", "36h", "2d", "3d", "5d", "8d", "15d"))
zfish_timepoints <- levels(zfish.metadata$Age)
zfish.samples <- zfish.metadata$FullID

# Get CPM values:
zfish.cpm <- cpm(zfish.sc.pseudo.df.merged[,colnames(zfish.sc.pseudo.df.merged) %in% zfish.samples], log = FALSE) %>% as.data.frame()
zfish.cpm$Human_GeneName <- zfish.sc.pseudo.df.merged$Human_GeneName

```

## Process BrainSpan data.
```{r}
# Load and process the BrainSpan data:
brain_span <- read.csv("brainspan_RNAseq/expression_matrix.csv", header = F)
brain_span_meta_row <- read.csv("brainspan_RNAseq/rows_metadata.csv")
brain_span_meta_col <- read.csv("brainspan_RNAseq/columns_metadata.csv")
rownames(brain_span) <- make.names(brain_span_meta_row$gene_symbol, unique = T)
brain_span <- brain_span[,-1]
full_id <- gsub(" ", "", paste0(brain_span_meta_col$structure_acronym, "_",
                                brain_span_meta_col$age, "_", brain_span_meta_col$donor_id))
colnames(brain_span) <-  brain_span_meta_col$FullID <- full_id

bs_metadata <- brain_span_meta_col
bs_metadata$tmpage <- as.numeric(ifelse(bs_metadata$age %in% 
                                          grep("pcw", bs_metadata$age, value = T),
                                        gsub("pcw", "", bs_metadata$age), 10000))
bs_metadata$AgeGroup <- "Other"
bs_metadata[bs_metadata$tmpage < 13, "AgeGroup"] <- "Early Fetal"
bs_metadata[bs_metadata$tmpage >=13 & bs_metadata$tmpage <19, "AgeGroup"] <- "Early Mid Fetal"
bs_metadata[bs_metadata$tmpage >=19 & bs_metadata$tmpage <24, "AgeGroup"] <- "Late Mid Fetal"
bs_metadata[bs_metadata$tmpage >=24 & bs_metadata$tmpage <38, "AgeGroup"] <- "Late Fetal"
bs_metadata[bs_metadata$age == "4 mos", "AgeGroup"] <- "Neonatal and early infancy"
bs_metadata[bs_metadata$age %in% c("6 mos","10 mos", "12 mos"), "AgeGroup"] <- "Late infancy"
bs_metadata$tmpage <- as.numeric(ifelse(bs_metadata$age %in% 
                                          grep("yrs", bs_metadata$age, value = T),
                                        gsub("yrs", "", bs_metadata$age), 10000))
bs_metadata[bs_metadata$tmpage < 6 & bs_metadata$tmpage >=1 , "AgeGroup"] <- "Early Childhood"
bs_metadata[bs_metadata$tmpage < 12 & bs_metadata$tmpage >=6 , "AgeGroup"] <- "Middle and Late Childhood"
bs_metadata[bs_metadata$tmpage < 20 & bs_metadata$tmpage >=12 , "AgeGroup"] <- "Adolescence"
bs_metadata[bs_metadata$tmpage < 40 & bs_metadata$tmpage >=20 , "AgeGroup"] <- "Young adulthood"
bs_metadata[bs_metadata$tmpage < 1000 & bs_metadata$tmpage >=40, "AgeGroup"] <- "Adulthood"
agelevels=c("Early Fetal", "Early Mid Fetal", "Late Mid Fetal", "Late Fetal", "Neonatal and early infancy", "Late infancy", "Early Childhood", "Middle and Late Childhood", "Adolescence", "Young adulthood", "Adulthood")
bs_metadata$AgeGroup <- factor(bs_metadata$AgeGroup, levels=agelevels)

bs_metadata$Struc_Group <- "Other"
bs_metadata[bs_metadata$structure_acronym %in% c("DFC","VFC","MFC","OFC","M1C-S1C","M1C","S1C"), ]$Struc_Group <- "PFC-MSC"
bs_metadata[bs_metadata$structure_acronym %in% c("V1C","IPC"), ]$Struc_Group <- "IPC-V1C"
bs_metadata[bs_metadata$structure_acronym %in% c("HIP","STR"), ]$Struc_Group <- "HIP-STR"
bs_metadata[bs_metadata$structure_acronym %in% c("MD","CBC"), ]$Struc_Group <- "MD-CBC"

# Filter low expressing genes from BrainSpan (RPKM<3):
brain_span_sub2 <- brain_span[rowSums(brain_span) > 0,]
brain_span_sub3 <- brain_span_sub2[rowSums(brain_span_sub2 > 3) >= 5,]
```

## Merge human and model organisms data.
```{r}
## Collapse human gene names with more than 1 ortholog by obtaining their median value:

mouse.dup <- mouse.cpm$Human_GeneName[duplicated(mouse.cpm$Human_GeneName)] # 3210 human genes are found multiple times
pt1 <- mouse.cpm %>% 
  group_by(Human_GeneName) %>% 
  summarise_at(vars(all_of(mouse.samples)), median)
pt2 <- mouse.cpm %>%
  tibble::rownames_to_column() %>% 
  group_by(Human_GeneName) %>% 
  summarise_at(vars("rowname"), toString)
mouse.collapsed <- merge(pt1, pt2, by="Human_GeneName")
rownames(mouse.collapsed) <- mouse.collapsed$Human_GeneName # make rows the HUMAN genes
mouse.collapsed <- subset(mouse.collapsed, select=mouse.samples)
dim(mouse.collapsed) # 16,739 genes

zfish.dup <- zfish.cpm$Human_GeneName[duplicated(zfish.cpm$Human_GeneName)] # 4,310 human genes are found multiple times
pt1 <- zfish.cpm %>% 
  group_by(Human_GeneName) %>% 
  summarise_at(vars(all_of(zfish.samples)), median)
pt2 <- zfish.cpm %>%
  tibble::rownames_to_column() %>% 
  group_by(Human_GeneName) %>% 
  summarise_at(vars("rowname"), toString)
zfish.collapsed <- merge(pt1, pt2, by="Human_GeneName")
rownames(zfish.collapsed) <- zfish.collapsed$Human_GeneName # make rows the HUMAN genes
zfish.collapsed <- subset(zfish.collapsed, select=zfish.samples)
dim(zfish.collapsed) # 12,599 genes

## Merge all datasets (human/mouse/zebrafish):

# Metadata:
metadata.mouse.zfish <- merge(mouse.metadata, zfish.metadata, all=T)
colnames(metadata.mouse.zfish) <- c("FullID", "Age")  # match colname for Sample in BrainSpan data
metadata.complete <- merge(bs_metadata, metadata.mouse.zfish, by="FullID", all=T)
metadata.complete$Experiment <- ifelse(metadata.complete$FullID %in% mouse.samples,
                                       "Mouse",
                                       ifelse(metadata.complete$FullID %in% zfish.samples,
                                              "Zebrafish",
                                              "BrainSpan"))

metadata.complete$AgeGroup <- as.factor(ifelse(is.na(metadata.complete$AgeGroup),
                                       as.character(metadata.complete$Age),
                                       as.character(metadata.complete$AgeGroup)))

alllevels <- c(as.character(agelevels), levels(mouse.metadata$Age), levels(zfish.metadata$Age))
metadata.complete$AgeGroup <- factor(metadata.complete$AgeGroup, levels=alllevels)

# CPM values:
genes_in_common <- Reduce(intersect, list(rownames(brain_span_sub3), rownames(mouse.collapsed), rownames(zfish.collapsed))) # 9,238 common genes
brain_span.common <- brain_span_sub3[genes_in_common,]
mouse.collapsed.common <- mouse.collapsed[genes_in_common,]
zfish.collapsed.common <- zfish.collapsed[genes_in_common,]

merged_temp <- merge(brain_span.common, mouse.collapsed.common[, mouse.samples], by = "row.names", all = TRUE)
rownames(merged_temp) <- merged_temp$Row.names
merged_temp <- merged_temp[,-1]
complete.cpm <- merge(merged_temp, 
                      zfish.collapsed.common[,zfish.samples],
                      by="row.names", all=T)
rownames(complete.cpm) <- complete.cpm$Row.names
complete.cpm <- complete.cpm[,-1]

## Normalize data:
complete.cpm.norm <- normalize.quantiles(as.matrix(complete.cpm), copy = TRUE)
rownames(complete.cpm.norm) <- rownames(complete.cpm)
colnames(complete.cpm.norm) <- colnames(complete.cpm)

## Keep genes with variance > 0:
genes.to.keep <- t(complete.cpm.norm[as.numeric(which(apply(complete.cpm.norm, 1, var) > 0 )), ]) # 9,238 genes, all kept

## Obtain log2 values:
complete.cpm.norm.log <- log2(complete.cpm.norm+1) # add a 1 to avoid Inf conflicts

## Perform exploratory PCAs per organism:

# BrainSpan:
BS.data <- t(complete.cpm.norm.log[, colnames(complete.cpm.norm.log) %in% bs_metadata$FullID])

pca.BS <- prcomp(BS.data, center = T, scale. = T)
percentVar.BS <- pca.BS$sdev^2/sum(pca.BS$sdev^2)
BS.tmp.df <- data.frame(PC1 = pca.BS$x[, 1], PC2 = pca.BS$x[, 2])
BS.tmp.df2 <- merge(BS.tmp.df, bs_metadata, all=TRUE, by.x="row.names", by.y="FullID")
p.PCA.BS <- BS.tmp.df2 %>%
  ggplot(aes(x=PC1, y=PC2, color=AgeGroup)) + geom_point(size = 3) + 
    coord_fixed() + theme_bw() +
    xlab(paste0("PC1 ",round(percentVar.BS[1] * 100), "%")) + 
    ylab(paste0("PC2 ",round(percentVar.BS[2] * 100), "%")) + 
    ggtitle("BrainSpan (common genes)") + scale_color_viridis_d()
ggsave("PCA_BrainSpan.pdf", p.PCA.BS, height = 8, width= 8)

# Mouse:
mouse.data <- t(complete.cpm.norm.log[, colnames(complete.cpm.norm.log) %in% mouse.samples])

pca.mouse <- prcomp(mouse.data, center = T, scale. = T)
percentVar.mouse <- pca.mouse$sdev^2/sum(pca.mouse$sdev^2)
mouse.tmp.df <- data.frame(PC1 = pca.mouse$x[, 1], PC2 = pca.mouse$x[, 2])
mouse.tmp.df2 <- merge(mouse.tmp.df, mouse.metadata, all=TRUE, by.x="row.names", by.y="Sample")
p.PCA.mouse <- mouse.tmp.df2 %>%
  ggplot(aes(x=PC1, y=PC2, color=Age)) + geom_point(size = 3) + 
    coord_fixed() + theme_bw() +
    xlab(paste0("PC1 ", round(percentVar.BS[1] * 100), "%")) + 
    ylab(paste0("PC2 ", round(percentVar.BS[2] * 100), "%")) + 
    ggtitle("Mouse (common genes)") + scale_color_viridis_d()
ggsave("PCA_mouse.pdf", p.PCA.mouse, height = 8, width= 8)

# Zebrafish:
zfish.data <- t(complete.cpm.norm.log[, colnames(complete.cpm.norm.log) %in% zfish.samples])

pca.zfish <- prcomp(zfish.data, center = T, scale. = T)
percentVar.zfish <- pca.zfish$sdev^2/sum(pca.zfish$sdev^2)
zfish.tmp.df <- data.frame(PC1 = pca.zfish$x[, 1], PC2 = pca.zfish$x[, 2])
zfish.tmp.df2 <- merge(zfish.tmp.df, zfish.metadata, all=TRUE, by.x="row.names", by.y="Sample")
p.PCA.zfish <- zfish.tmp.df2 %>%
  ggplot(aes(x=PC1, y=PC2, color=Age)) + geom_point(size = 3) + 
    coord_fixed() + theme_bw() +
    xlab(paste0("PC1 ", round(percentVar.BS[1] * 100), "%")) + 
    ylab(paste0("PC2 ", round(percentVar.BS[2] * 100), "%")) + 
    ggtitle("Zebrafish (common genes)") + scale_color_viridis_d()
ggsave("PCA_zebrafish.pdf", p.PCA.zfish, height = 8, width= 8)

## Predict PC coordinates for mouse and zebrafish data in human PC space:
models.data <- rbind(mouse.data, zfish.data)
predicted <- scale(models.data, center= pca.BS$center, scale = pca.BS$scale) %*% pca.BS$rotation
pca.combined <- pca.BS
pca.combined$x <- rbind(pca.combined$x, predicted)
label_order <- metadata.complete$AgeGroup

## Extract common PC values from human PC space:
pc.values <- paper.ggbiplot2(pca.combined, var.axes=FALSE) 
colnames(pc.values) <- c("PC1","PC2")

## Add PC values to metadata:
pc.values_reordered <- pc.values[match(metadata.complete$FullID, rownames(pc.values)), ]
metadata.complete.final <- bind_cols(metadata.complete, pc.values_reordered)
```

## Make final plots.
```{r}

## Figure 1: PCA with embedded samples:

agelevels= c("Early Fetal", "Early Mid Fetal", "Late Mid Fetal", "Late Fetal",
             "Neonatal and early infancy", "Late infancy", "Early Childhood",
             "Middle and Late Childhood", "Adolescence", "Young adulthood",
             "Adulthood")
agelevelsNames= c("Early Fetal (<12 pcw)", "Early Mid Fetal (13-17 pcw)", "Late Mid Fetal (19-21 pcw)", "Late Fetal (24-37 pcw)",
                  "Neonatal and early infancy (4 months)", "Late infancy (10 months)", "Early Childhood (1-4 years)",
                  "Middle and Late Childhood (8-11 years)", "Adolescence (13-19 years)", "Young adulthood (21-37 years)",
                  "Adulthood (>40 years)")

metadata.complete.final$AgeGroup_sum <- ifelse(metadata.complete.final$AgeGroup %in% agelevels,
                                               as.character(metadata.complete.final$AgeGroup), 
                                               as.character(metadata.complete.final$Experiment)) # make a summarize value
metadata.complete.final$AgeGroup_sum <- factor(metadata.complete.final$AgeGroup_sum, levels = c("Early Fetal", "Early Mid Fetal", "Late Mid Fetal", "Late Fetal", "Neonatal and early infancy", "Late infancy", "Early Childhood", "Middle and Late Childhood", "Adolescence", "Young adulthood","Adulthood", "Mouse", "Zebrafish"))

my.pal <- c(scales::viridis_pal()(length(agelevels)), "#D2042D", "#3F00FF")

p1 <- ggplot(data = metadata.complete.final, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = AgeGroup_sum, shape = Experiment), size = 3) +
  scale_color_manual(values = my.pal) +
  scale_shape_manual(values = c("BrainSpan"=16,"Mouse" = 17, "Zebrafish" = 18)) +
  geom_label_repel(data = subset(metadata.complete.final, Experiment != "BrainSpan"),
                   aes(label = Age, fill = Experiment),
                   show.legend = FALSE,
                   nudge_x = 0.2, nudge_y = 0.5, 
                   color = "black", 
                   size = 3, 
                   box.padding = 0.25, max.overlaps = 20) +
  scale_fill_manual(values = c("Mouse" = "red", "Zebrafish" = "blue")) +
  xlab(paste0("PC1 ", round(percentVar.BS[1] * 100), "%")) + 
  ylab(paste0("PC2 ", round(percentVar.BS[2] * 100), "%")) +
  theme_bw()

ggsave("PCA_combined_all.pdf",
       p1,
       height=8, width=8)

## Figure 2: PCA with BrainSpan timepoints as blocks and zebrafish samples predicted:

mouse.palette <- colorRampPalette(c("#FFEBE6", "#6E0E0A"))(n = length(mouse_timepoints))
zfish.palette <- colorRampPalette(c("#E0F3FF", "#08233F"))(n = length(zfish_timepoints))

BSlimits <- data.frame(AgeGroup=NA, xmin=NA, xmax=NA)
for (i in agelevels){
  metad2 <- metadata.complete.final[metadata.complete.final$Experiment == "BrainSpan",]
    if (i == "Early Fetal"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Early Mid Fetal"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Late Mid Fetal"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Late Fetal"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Neonatal and early infancy"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Late infancy"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Early Childhood"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Middle and Late Childhood"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Adolescence"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Young adulthood"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else if(i == "Adulthood"){vect <- metad2$PC1[metad2$AgeGroup == i]}
    else {vect <- as.numeric(metadata.complete.final[metadata.complete.final$AgeGroup ==i, "PC1"])}
    my.mean <- mean(vect)
    my.sd <- sd(vect)
    me <- (qt(0.975, (length(vect)-1)) * my.sd)/sqrt(length(vect))
    temp <- c(i, my.mean-me, my.mean+me)
    BSlimits <- rbind(BSlimits, temp)
  }

BSlimits <- BSlimits[!is.na(BSlimits$AgeGroup),] # remove NA row
BSlimits$xmin <- as.numeric(BSlimits$xmin)
BSlimits$xmax <- as.numeric(BSlimits$xmax)
BSlimits$AgeGroup <- factor(BSlimits$AgeGroup, levels=agelevels)

metadata.complete.final$manual_jitter <- runif(nrow(metadata.complete.final), min=0, max=1)

data_to_plot <- metadata.complete.final[metadata.complete.final$Experiment != "BrainSpan",]

p2 <- ggplot(data = data_to_plot, aes(x = PC1, y = manual_jitter, color = AgeGroup, shape = Experiment)) +
  # geom rect first so it's underneath
  geom_rect(data = BSlimits, inherit.aes = FALSE,  
            aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = AgeGroup), alpha = 0.5) +
  scale_fill_manual(name = "Human Developmental Stage (95% CI)",
                    values = my.pal,
                    labels = agelevelsNames) +
  # points with color mapped to AgeGroup and shape mapped to Experiment
  geom_point(size = 3) +
  scale_color_manual(name = "Model organisms", values = c(mouse.palette, zfish.palette), drop = TRUE) +
  scale_shape_manual(name = "Experiment", values = c("Mouse" = 17, "Zebrafish" = 16)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical") +
  xlab("PC1")


ggsave("PCA_combined_all_rectangular.pdf",
       p2,
       height=8, width=10)

```

# STEP 3: Extract expression data for HSD genes.
```{r}
# Load list of HSDs orthologs:
mouse_hsd_orthologs <- read.csv("files/hsds_mouse_orthologs.csv", head=F) # 503 orthologous genes
zfish_hsd_orthologs <- read.csv("files/hsds_zfish_orthologs.csv", head=F) # 393 orthologous genes

# Subset expression values:
genes_shared <- unique(rownames(mouse.cpm)[rownames(mouse.cpm) %in% mouse_hsd_orthologs$V1])
mouse.cpm.hsds <- mouse.cpm[genes_shared,] # 271 HSD ortholog genes with expression values in mouse
write.csv(as.data.frame(mouse.cpm.hsds), file="mouse.cpm.hsds.csv")

genes_shared <- unique(rownames(zfish.cpm)[rownames(zfish.cpm) %in% zfish_hsd_orthologs$V1])
zfish.cpm.hsds <- zfish.cpm[genes_shared,] # 221 HSD ortholog genes with expression values in zebrafish
write.csv(as.data.frame(zfish.cpm.hsds), file="zfish.cpm.hsds.csv")
```

# STEP 4: Make a heatmap with expression of HSD genes and their orthologs across brain development.

## Assign model organisms timepoints to their closest human developmental stages.
```{r}
models.PC.medians <- read.csv("files/models.PC1.means.csv", head=T, stringsAsFactors = F)

human_stage <- NULL
for(i in 1:nrow(models.PC.medians)){
  value <- models.PC.medians$Average.PC1[i]
  tmp <-  BSlimits$AgeGroup[which.min(pmin(pmax(BSlimits$xmin, value), pmin(BSlimits$xmax, value)))]
  human_stage[i] <- tmp
}

```


## Plot the summarized expression in model organisms by repeating same ortholog across gene family
```{r}
# Load processed dataset:
df.hsds <- read.csv("files/TableS11_new.csv", head=T, stringsAsFactors = F)

# Scale each sub-matrix by gene across development:
human_timepoints <- c("Early.prenatal", "Early.mid.prenatal",
                      "Late.mid.prenatal", "Late.prenatal")
mouse_timepoints <- c("Mouse_EarlyMidFetal", "Mouse_LateFetal")
zfish_timepoints <- c("Zebrafish_EarlyFetal", "Zebrafish_EarlyMidFetal", "Zebrafish_LateFetal")

metadata_cols <- df.hsds[,1:8]
scaled_human <- t(apply(df.hsds[,human_timepoints],1,scale))
scaled_mouse <- t(apply(df.hsds[,mouse_timepoints],1,scale))
scaled_zfish <- t(apply(df.hsds[,zfish_timepoints],1,scale))

df.hsds.scaled <- cbind(metadata_cols,scaled_human,scaled_mouse,scaled_zfish)
colnames(df.hsds.scaled) <- c(colnames(metadata_cols), human_timepoints, mouse_timepoints, zfish_timepoints)

# Make heatmap:
df.hsds$To.plot_fixed <- as.factor(df.hsds$To.plot_fixed)
gene_families_to_plot <- unique(df.hsds$Family.ID[df.hsds$To.plot_fixed=="Yes"])

data.to.plot <- df.hsds.scaled[df.hsds$Family.ID%in%gene_families_to_plot,c(9:17)]
data.to.plot$gene <- df.hsds.scaled$Gene.Name[df.hsds$Family.ID%in%gene_families_to_plot]
data.to.plot.long <- pivot_longer(data.to.plot, cols = -gene, names_to = "Timepoint", values_to = "Expression") %>% as.data.frame()
data.to.plot.long$Timepoint <- factor(data.to.plot.long$Timepoint)
data.to.plot.long$gene <- factor(data.to.plot.long$gene, levels = unique(data.to.plot.long$gene))
data.to.plot.long$gene <- factor(data.to.plot.long$gene, levels = rev(levels(data.to.plot.long$gene))) # reverse order for plotting

p <- ggplot(data.to.plot.long, aes(x = Timepoint, y = gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", na.value = "gray") +
  labs(x = "Timepoint", y = "Gene") +
  theme_classic() + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("HSDs_ortho_heatmap_summarized_models_repeated_per_family.pdf", p, height=30, width=8)

```
