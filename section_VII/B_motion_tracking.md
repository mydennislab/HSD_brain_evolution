# Behavior analysis

```{r}
library(ggplot2)
library(dunn.test)
library(tidyverse)
library(patchwork)
library(MASS)

## Load and pre-process data:

df <- read.csv("Motion_tracking.csv", head=T, stringsAsFactors = T)
str(df)
df$N_events <- as.numeric(df$N_events)
df$N_HSE <- as.numeric(df$N_HSE)
df$PTZ <- as.factor(df$PTZ)
df$Plate <- as.factor(df$Plate)

## Decide threshold for dead fish based on low movement:

# threshold: 1 sd away from mean in TotalMovement:
threshold <- mean(df$TotalMovement[df$PTZ=="0"]) - (1*sd(df$TotalMovement[df$PTZ=="0"]))

# Plot TotalMovement with threshold:
ggplot(df, aes(x= TotalMovement)) + geom_histogram(bins = 100) + 
  facet_wrap(~PTZ) + geom_vline(xintercept = threshold, linetype="dashed", color="red") + theme_minimal()

# Filter database based on threshold and get summary data:
df_filt <- df[df$TotalMovement >= threshold & df$Genotype!="",]
df_filt <- droplevels(df_filt)

## Test if there are differences between controls:

df_ctrl <- df_filt[df_filt$Group=="control",]
df_ctrl <- droplevels(df_ctrl)
summary(df_ctrl)

ctrl1 <- ggplot(df_ctrl, aes(x= Genotype, y=TotalMovement, color= Genotype)) +
  geom_boxplot() + geom_jitter(width=0.1) + theme_minimal() +
  facet_wrap(~PTZ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Total movement (mm)") + theme(legend.position = "none")

ctrl2 <- ggplot(df_ctrl, aes(x= Genotype, y=MaxSpeed, color= Genotype)) +
  geom_boxplot() + geom_jitter(width=0.1) + theme_minimal() +
  facet_wrap(~PTZ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Max speed (mm/s)") + theme(legend.position = "none")

ctrl3 <- ggplot(df_ctrl, aes(x= Genotype, y=N_HSE_permin, color= Genotype)) +
  geom_boxplot() + geom_jitter(width=0.1) + theme_minimal() +
  facet_wrap(~PTZ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("HSE per min") + theme(legend.position = "none")

ctrl1 / ctrl2 / ctrl3

dunn.test(x= df_ctrl$TotalMovement[df_ctrl$PTZ=="0"],
          g= df_ctrl$Genotype[df_ctrl$PTZ=="0"], method = "bh")
# No differences: ControlMO vs eGFP: 0.0735, ControlMO vs scramble: 0.2546, scramble vs eGFP: 0.0519

dunn.test(x= df_ctrl$TotalMovement[df_ctrl$PTZ=="2.5"],
          g= df_ctrl$Genotype[df_ctrl$PTZ=="2.5"], method = "bh")
# Diff vs controlMO and the others, so will merge scrambled and eGFP and not controlMO (for arhgap11MO control)

dunn.test(x= df_ctrl$MaxSpeed[df_ctrl$PTZ=="0"],
          g= df_ctrl$Genotype[df_ctrl$PTZ=="0"], method = "bh")
# No differences: ControlMO vs eGFP: 0.4444, ControlMO vs scramble: 0.3334, scramble vs eGFP: 0.2710

dunn.test(x= df_ctrl$MaxSpeed[df_ctrl$PTZ=="2.5"],
          g= df_ctrl$Genotype[df_ctrl$PTZ=="2.5"], method = "bh")
# Diff vs controlMO and the others, so will merge scrambled and eGFP and not controlMO (for arhgap11MO control)

dunn.test(x= df_ctrl$N_HSE_permin[df_ctrl$PTZ=="0"],
          g= df_ctrl$Genotype[df_ctrl$PTZ=="0"], method = "bh")
# No differences: ControlMO vs eGFP: 0.3708, ControlMO vs scramble: 0.4269, scramble vs eGFP: 0.1425

dunn.test(x= df_ctrl$N_HSE_permin[df_ctrl$PTZ=="2.5"],
          g= df_ctrl$Genotype[df_ctrl$PTZ=="2.5"], method = "bh")
# Diff vs controlMO and the others, so will merge scrambled and eGFP and not controlMO (for arhgap11MO control)


## Re-order Genotype variable for plots:
df_filt$Genotype2 <- factor(df_filt$Genotype2, levels = c("control", "controlMO",
                                                    "srgap2KO", "SRGAP2A", "SRGAP2C",
                                                    "arhgap11MO", "ARHGAP11A", "ARHGAP11B",
                                                    "gpr89KO", "GPR89B",
                                                    "npy4rKO", "NPY4R",
                                                    "frmpd2KO", "FRMPD2A", "FRMPD2B",
                                                    "fam72KO", "FAM72B",
                                                    "ncf1KO", "NCF1A", "NCF1C",
                                                    "hydinKO",
                                                    "pdzk1KO", "PDZK1", "PDZK1P1",
                                                    "ptpn20KO", "PTPN20", "PTPN20CP"))

## Summary of data:

summary(df_filt$Genotype2[df$PTZ=="0"])
# control: 699
# controlMO: 47
# srgap2KO: 39, SRGAP2A: 39, SRGAP2C: 39
# arhgap11MO: 64, ARHGAP11A: 60, ARHGAP11B: 68
# gpr89KO: 156, GPR89B: 146
# npy4rKO: 126, NPY4R: 74
# frmpd2KO: 54, FRMPD2A: 48, FRMPD2B: 65
# fam72KO: 40, FAM72B: 27
# ncf1KO: 46, NCF1A: 59, NCF1C: 50
# hydinKO: 92
# pdzk1KO: 45, PDZK1: 76, PDZK1P1: 54
# ptpn20KO: 33, PTPN20: 36, PTPN20CP: 44

summary(df_filt$Genotype2[df$PTZ=="2.5"])
# control: 693
# controlMO: 43
# srgap2KO: 51, SRGAP2A: 50, SRGAP2C: 79
# arhgap11MO: 64, ARHGAP11A: 40, ARHGAP11B: 63
# gpr89KO: 161, GPR89B: 143
# npy4rKO: 100, NPY4R: 63
# frmpd2KO: 84, FRMPD2A: 60, FRMPD2B: 71
# fam72KO: 46, FAM72B: 16
# ncf1KO: 52, NCF1A: 27, NCF1C: 59
# hydinKO: 62
# pdzk1KO: 44, PDZK1: 54, PDZK1P1: 44
# ptpn20KO: 34, PTPN20: 41, PTPN20CP: 46

## General plots and tests:

# Total movement:

fig1 <- ggplot(df_filt, aes(x= Genotype2, y= TotalMovement, color= Genotype2)) +
  geom_boxplot() + geom_jitter(width=0.1) + theme_minimal() +
  facet_wrap(~PTZ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Total movement (mm)") + theme(legend.position = "none")

TotMov0 <- as.data.frame(dunn.test(x= df_filt$TotalMovement[df_filt$PTZ=="0"],
          g= df_filt$Genotype2[df_filt$PTZ=="0"], method = "bh"))

aggregate(df_filt$TotalMovement[df_filt$PTZ=="0"], 
          list(df_filt$Genotype2[df_filt$PTZ=="0"]), FUN=median) 

TotMov2.5 <- as.data.frame(dunn.test(x= df_filt$TotalMovement[df_filt$PTZ=="2.5"],
                                   g= df_filt$Genotype2[df_filt$PTZ=="2.5"], method = "bh"))

aggregate(df_filt$TotalMovement[df_filt$PTZ=="2.5"], 
          list(df_filt$Genotype2[df_filt$PTZ=="2.5"]), FUN=median)

# Max speed:

fig2 <- ggplot(df_filt, aes(x= Genotype2, y= MaxSpeed, color= Genotype2)) +
  geom_boxplot() + geom_jitter(width=0.1) + theme_minimal() +
  facet_wrap(~PTZ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("Max speed (mm/s)") + theme(legend.position = "none")

MaxSpeed0 <- as.data.frame(dunn.test(x= df_filt$MaxSpeed[df_filt$PTZ=="0"],
                                   g= df_filt$Genotype2[df_filt$PTZ=="0"], method = "bh"))

aggregate(df_filt$MaxSpeed[df_filt$PTZ=="0"], 
          list(df_filt$Genotype2[df_filt$PTZ=="0"]), FUN=median) 

MaxSpeed2.5 <- as.data.frame(dunn.test(x= df_filt$MaxSpeed[df_filt$PTZ=="2.5"],
                                     g= df_filt$Genotype2[df_filt$PTZ=="2.5"], method = "bh"))

aggregate(df_filt$MaxSpeed[df_filt$PTZ=="2.5"], 
          list(df_filt$Genotype2[df_filt$PTZ=="2.5"]), FUN=median)

# HSE per min:

fig3 <- ggplot(df_filt, aes(x= Genotype2, y= N_HSE_permin, color= Genotype2)) +
  geom_boxplot() + geom_jitter(width=0.1) + theme_minimal() +
  facet_wrap(~PTZ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("HSE per min") + theme(legend.position = "none")

HSE0 <- as.data.frame(dunn.test(x= df_filt$N_HSE_permin[df_filt$PTZ=="0"],
                                     g= df_filt$Genotype2[df_filt$PTZ=="0"], method = "bh"))

aggregate(df_filt$N_HSE_permin[df_filt$PTZ=="0"], 
          list(df_filt$Genotype2[df_filt$PTZ=="0"]), FUN=median) 

HSE2.5 <- as.data.frame(dunn.test(x= df_filt$N_HSE_permin[df_filt$PTZ=="2.5"],
                                       g= df_filt$Genotype2[df_filt$PTZ=="2.5"], method = "bh"))

aggregate(df_filt$N_HSE_permin[df_filt$PTZ=="2.5"], 
          list(df_filt$Genotype2[df_filt$PTZ=="2.5"]), FUN=median)

# Merged figures:

fig1 / fig2 / fig3

```
