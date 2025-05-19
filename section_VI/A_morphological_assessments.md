# Morphometric data analysis

```{r}
# Libraries:
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(mvnormtest)
library(dunn.test)
library(car)
library(ggbiplot)
library(patchwork)
library(DescTools)
library(Hmisc)
library(corrplot)
library(rstatix)
library(ggpubr)

`%notin%` <- Negate(`%in%`)

## Load and pre-process data:
df <- read.csv("Morphometric_data.csv", head=T, stringsAsFactors = T)
df$Plate <- as.factor(df$Plate)

## Filter images with general issues (dead, decomposing):
df <- df[df$GeneralIssues=="No",]
df <- df[df$Truncated=="No",]
df <- df[df$Genotype != "uninjected",]

## Add head area minus the eye area:
df$HeadArea_noEyes <- (df$AreaBeforeXMin_contourDV_yolkDV_um2) - (df$eye1DV_regionpropsArea_um2 + df$eye2DV_regionpropsArea_um2)

## Quon Lab flagged SRGAP2, GPR89, FAM72, and FRMPD2 as showing distinguishable differences to 

## Extract genotypes/features of interest:
df_filt <- df[df$Genotype %in% c("scrambled", "eGFP", "srgap2KO", "SRGAP2C", "fam72KO", "FAM72B", "frmpd2KO", "FRMPD2B", "gpr89KO", "GPR89B"),] %>%
  dplyr::select(DateVAST, Plate, Well, Group, Genotype, Genotype2, Age,
                contourDV_regionpropsLengthOfCentralLine_um,
                AreaBeforeXMin_contourDV_yolkDV_um2,
                MeanOfArea_eye1DV_eye2DV_um2,
                HeadArea_noEyes,
                HeadTrunkAngle_angle.1) %>%
  droplevels()

colnames(df_filt) <- c("DateVAST", "Plate", "Well", "Group", "Genotype", "Genotype2", "Age",
                           "Length", "HeadArea", "MeanEyeArea", "HeadArea_noEyes","HeadTrunkAngle")


## Separate the data into 3 and 5 dpf for the selected variables of interest:
df_3 <- df_filt[df_filt$Age=="3",] %>%
  droplevels()
# N= 1012 larvae at 3 dpf

df_5 <- df_filt[df_filt$Age=="5",] %>%
  droplevels()
# N= 813 larvae at 5 dpf

## Filtering outliers potentially due to errors in shape-definition by FishInspecto (+/- 1.5*IQR):

# Find outlier values for each variable:
outliers.length <- df_3$Length[df_3$Length > quantile(df_3$Length,na.rm = T)[4] + 2*IQR(df_3$Length,na.rm = T) | df_3$Length < quantile(df_3$Length,na.rm = T)[2] - 2*IQR(df_3$Length,na.rm = T)]
outliers.headarea <- df_3$HeadArea[df_3$HeadArea > quantile(df_3$HeadArea,na.rm = T)[4] + 2*IQR(df_3$HeadArea,na.rm = T) | df_3$HeadArea < quantile(df_3$HeadArea,na.rm = T)[2] - 2*IQR(df_3$HeadArea,na.rm = T)]
outliers.eyearea <- df_3$MeanEyeArea[df_3$MeanEyeArea > quantile(df_3$MeanEyeArea,na.rm = T)[4] + 2*IQR(df_3$MeanEyeArea,na.rm = T) | df_3$MeanEyeArea < quantile(df_3$MeanEyeArea,na.rm = T)[2] - 2*IQR(df_3$MeanEyeArea,na.rm = T)]
outliers.headangle <- df_3$HeadTrunkAngle[df_3$HeadTrunkAngle > quantile(df_3$HeadTrunkAngle,na.rm = T)[4] + 2*IQR(df_3$HeadTrunkAngle,na.rm = T) | df_3$HeadTrunkAngle < quantile(df_3$HeadTrunkAngle,na.rm = T)[2] - 2*IQR(df_3$HeadTrunkAngle,na.rm = T)]

outliers.length <- outliers.length[!is.na(outliers.length)] # n= 30
outliers.headarea <- outliers.headarea[!is.na(outliers.headarea)] # n= 14
outliers.eyearea <- outliers.eyearea[!is.na(outliers.eyearea)] # n= 35
outliers.headangle <- outliers.headangle[!is.na(outliers.headangle)] # n= 30

df_3_filt <- df_3[df_3$Length %notin% outliers.length & 
                    df_3$HeadArea %notin% outliers.headarea &
                    df_3$MeanEyeArea %notin% outliers.eyearea &
                    df_3$HeadTrunkAngle %notin% outliers.headangle,]
# N= 936
mean(table(df_3_filt$Genotype)) # 93.4 larvae per group
sd(table(df_3_filt$Genotype)) # 70.56786 larvae per group

outliers.length <- df_5$Length[df_5$Length > quantile(df_5$Length,na.rm = T)[4] + 2*IQR(df_5$Length,na.rm = T) | df_5$Length < quantile(df_5$Length,na.rm = T)[2] - 2*IQR(df_5$Length,na.rm = T)]
outliers.headarea <- df_5$HeadArea[df_5$HeadArea > quantile(df_5$HeadArea,na.rm = T)[4] + 2*IQR(df_5$HeadArea,na.rm = T) | df_5$HeadArea < quantile(df_5$HeadArea,na.rm = T)[2] - 2*IQR(df_5$HeadArea,na.rm = T)]
outliers.eyearea <- df_5$MeanEyeArea[df_5$MeanEyeArea > quantile(df_5$MeanEyeArea,na.rm = T)[4] + 2*IQR(df_5$MeanEyeArea,na.rm = T) | df_5$MeanEyeArea < quantile(df_5$MeanEyeArea,na.rm = T)[2] - 2*IQR(df_5$MeanEyeArea,na.rm = T)]
outliers.headangle <- df_5$HeadTrunkAngle[df_5$HeadTrunkAngle > quantile(df_5$HeadTrunkAngle,na.rm = T)[4] + 2*IQR(df_5$HeadTrunkAngle,na.rm = T) | df_5$HeadTrunkAngle < quantile(df_5$HeadTrunkAngle,na.rm = T)[2] - 2*IQR(df_5$HeadTrunkAngle,na.rm = T)]

outliers.length <- outliers.length[!is.na(outliers.length)] # n= 24
outliers.headarea <- outliers.headarea[!is.na(outliers.headarea)] # n= 5
outliers.eyearea <- outliers.eyearea[!is.na(outliers.eyearea)] # n= 13
outliers.headangle <- outliers.headangle[!is.na(outliers.headangle)] # n= 6

df_5_filt <- df_5[df_5$Length %notin% outliers.length & 
                    df_5$HeadArea %notin% outliers.headarea &
                    df_5$MeanEyeArea %notin% outliers.eyearea &
                    df_5$HeadTrunkAngle %notin% outliers.headangle,]
# N= 770
mean(table(df_5_filt$Genotype)) # 77.3 larvae per group
sd(table(df_5_filt$Genotype)) # 69.48389 larvae per group

## General evaluation of variables at each timepoint:

# Normality:

# 3dpf:
shapiro.test(df_3_filt$Length) # p-value = 3.289e-06
hist(df_3_filt$Length)
ggqqplot(df_3_filt$Length)

shapiro.test(df_3_filt$HeadArea) # p-value = 0.01071
hist(df_3_filt$HeadArea)
ggqqplot(df_3_filt$HeadArea)

shapiro.test(df_3_filt$MeanEyeArea) # p-value = 9.46e-05
hist(df_3_filt$MeanEyeArea)
ggqqplot(df_3_filt$MeanEyeArea)

shapiro.test(df_3_filt$HeadTrunkAngle) # p-value = 1.022e-05
hist(df_3_filt$HeadTrunkAngle)
ggqqplot(df_3_filt$HeadTrunkAngle)

shapiro.test(df_3_filt$HeadArea_noEyes) # p-value = 0.003659
hist(df_3_filt$HeadArea_noEyes)
ggqqplot(df_3_filt$HeadArea_noEyes)

# 5dpf:
shapiro.test(df_5_filt$Length) # p-value = 5.657e-14
hist(df_5_filt$Length)
ggqqplot(df_5_filt$Length)

shapiro.test(df_5_filt$HeadArea) # p-value = 0.04936
hist(df_5_filt$HeadArea)
ggqqplot(df_5_filt$HeadArea)

shapiro.test(df_5_filt$MeanEyeArea) # p-value = 6.343e-12
hist(df_5_filt$MeanEyeArea)
ggqqplot(df_5_filt$MeanEyeArea)

shapiro.test(df_5_filt$HeadTrunkAngle) # p-value = 0.138
hist(df_5_filt$HeadTrunkAngle)
ggqqplot(df_5_filt$HeadTrunkAngle)

shapiro.test(df_5_filt$HeadArea_noEyes) # p-value = 0.1504
hist(df_5_filt$HeadArea_noEyes)
ggqqplot(df_5_filt$HeadArea_noEyes)

## Statistical tests:

## 3 dpf:

# Are controls different?
wilcox.test(df_3_filt$Length[df_3_filt$Genotype=="eGFP"], df_3_filt$Length[df_3_filt$Genotype=="scrambled"]) # W = 25636, p-value = 0.1871
wilcox.test(df_3_filt$HeadArea[df_3_filt$Genotype=="eGFP"], df_3_filt$HeadArea[df_3_filt$Genotype=="scrambled"]) # W = 23216, p-value = 0.6644
wilcox.test(df_3_filt$MeanEyeArea[df_3_filt$Genotype=="eGFP"], df_3_filt$MeanEyeArea[df_3_filt$Genotype=="scrambled"]) # W = 24858, p-value = 0.4644
wilcox.test(df_3_filt$HeadTrunkAngle[df_3_filt$Genotype=="eGFP"], df_3_filt$HeadTrunkAngle[df_3_filt$Genotype=="scrambled"]) # W = 23916, p-value = 0.133
wilcox.test(df_3_filt$HeadArea_noEyes[df_3_filt$Genotype=="eGFP"], df_3_filt$HeadArea_noEyes[df_3_filt$Genotype=="scrambled"]) # W = 23134, p-value = 0.6197

# Extract the data for each genotypic group (because we need to account for their proper controls)

# SRGAP2:
df_3_srgap2 <- df_3_filt[df_3_filt$DateVAST %in% names(which(summary(df_3_filt[df_3_filt$Group=="SRGAP2",]$DateVAST) > 0)),]
df_3_srgap2 <- df_3_srgap2[df_3_srgap2$Group == "control" | df_3_srgap2$Group == "SRGAP2",]
df_3_srgap2 <- droplevels(df_3_srgap2)
df_3_srgap2$Genotype2 <- relevel(df_3_srgap2$Genotype2, ref = "control")
table(df_3_srgap2$Genotype2)
#control  SRGAP2C srgap2KO 
#53       44       43 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_3_srgap2_rank <- df_3_srgap2
df_3_srgap2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_3_srgap2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

srgap2_3_length <- lm(Length ~ Genotype2 + DateVAST + Plate, data = df_3_srgap2_rank)
summary(srgap2_3_length) # F-statistic: 7.393 on 5 and 134 DF,  p-value: 3.762e-06,	Adjusted R-squared:  0.187 
mean(df_3_srgap2$Length[df_3_srgap2$Genotype2=="control"])
mean(df_3_srgap2$Length[df_3_srgap2$Genotype2=="srgap2KO"])
mean(df_3_srgap2$Length[df_3_srgap2$Genotype2=="SRGAP2C"])
sd(df_3_srgap2$Length[df_3_srgap2$Genotype2=="control"])
sd(df_3_srgap2$Length[df_3_srgap2$Genotype2=="srgap2KO"])
sd(df_3_srgap2$Length[df_3_srgap2$Genotype2=="SRGAP2C"])

srgap2_3_headarea <- lm(HeadArea ~ Genotype2 + DateVAST + Plate, data = df_3_srgap2_rank)
summary(srgap2_3_headarea) # F-statistic: 19.34 on 5 and 134 DF,  p-value: 1.846e-14, Adjusted R-squared:  0.3975 
mean(df_3_srgap2$HeadArea[df_3_srgap2$Genotype2=="control"])
mean(df_3_srgap2$HeadArea[df_3_srgap2$Genotype2=="srgap2KO"])
mean(df_3_srgap2$HeadArea[df_3_srgap2$Genotype2=="SRGAP2C"])
sd(df_3_srgap2$HeadArea[df_3_srgap2$Genotype2=="control"])
sd(df_3_srgap2$HeadArea[df_3_srgap2$Genotype2=="srgap2KO"])
sd(df_3_srgap2$HeadArea[df_3_srgap2$Genotype2=="SRGAP2C"])

srgap2_3_eyearea <- lm(MeanEyeArea ~ Genotype2 + DateVAST + Plate, data = df_3_srgap2_rank)
summary(srgap2_3_eyearea) # F-statistic: 14.96 on 5 and 134 DF,  p-value: 1.178e-11, Adjusted R-squared:  0.3343 
mean(df_3_srgap2$MeanEyeArea[df_3_srgap2$Genotype2=="control"])
mean(df_3_srgap2$MeanEyeArea[df_3_srgap2$Genotype2=="srgap2KO"])
mean(df_3_srgap2$MeanEyeArea[df_3_srgap2$Genotype2=="SRGAP2C"])
sd(df_3_srgap2$MeanEyeArea[df_3_srgap2$Genotype2=="control"])
sd(df_3_srgap2$MeanEyeArea[df_3_srgap2$Genotype2=="srgap2KO"])
sd(df_3_srgap2$MeanEyeArea[df_3_srgap2$Genotype2=="SRGAP2C"])

srgap2_3_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2 + DateVAST + Plate, data = df_3_srgap2_rank)
summary(srgap2_3_headtrunkangle) # F-statistic: 5.852 on 5 and 134 DF,  p-value: 6.367e-05, Adjusted R-squared:  0.1486
mean(df_3_srgap2$HeadTrunkAngle[df_3_srgap2$Genotype2=="control"], na.rm = T)
mean(df_3_srgap2$HeadTrunkAngle[df_3_srgap2$Genotype2=="srgap2KO"], na.rm = T)
mean(df_3_srgap2$HeadTrunkAngle[df_3_srgap2$Genotype2=="SRGAP2C"], na.rm = T)
sd(df_3_srgap2$HeadTrunkAngle[df_3_srgap2$Genotype2=="control"], na.rm = T)
sd(df_3_srgap2$HeadTrunkAngle[df_3_srgap2$Genotype2=="srgap2KO"], na.rm = T)
sd(df_3_srgap2$HeadTrunkAngle[df_3_srgap2$Genotype2=="SRGAP2C"], na.rm = T)

srgap2_3_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2 + DateVAST + Plate, data = df_3_srgap2_rank)
summary(srgap2_3_headareanoeyes) # F-statistic: 18.92 on 5 and 134 DF,  p-value: 3.326e-14, Adjusted R-squared:  0.392
mean(df_3_srgap2$HeadArea_noEyes[df_3_srgap2$Genotype2=="control"], na.rm = T)
mean(df_3_srgap2$HeadArea_noEyes[df_3_srgap2$Genotype2=="srgap2KO"], na.rm = T)
mean(df_3_srgap2$HeadArea_noEyes[df_3_srgap2$Genotype2=="SRGAP2C"], na.rm = T)
sd(df_3_srgap2$HeadArea_noEyes[df_3_srgap2$Genotype2=="control"], na.rm = T)
sd(df_3_srgap2$HeadArea_noEyes[df_3_srgap2$Genotype2=="srgap2KO"], na.rm = T)
sd(df_3_srgap2$HeadArea_noEyes[df_3_srgap2$Genotype2=="SRGAP2C"], na.rm = T)

# GPR89:
df_3_gpr89 <- df_3_filt[df_3_filt$DateVAST %in% names(which(summary(df_3_filt[df_3_filt$Group=="GPR89",]$DateVAST) > 0)),]
df_3_gpr89 <- df_3_gpr89[df_3_gpr89$Group == "control" | df_3_gpr89$Group == "GPR89",]
df_3_gpr89 <- droplevels(df_3_gpr89)
df_3_gpr89$Genotype2 <- relevel(df_3_gpr89$Genotype2, ref = "control")
table(df_3_gpr89$Genotype2)
#control  GPR89B gpr89KO 
#175       62       79 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_3_gpr89_rank <- df_3_gpr89
df_3_gpr89_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_3_gpr89_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

gpr89_3_length <- lm(Length ~ Genotype2 + DateVAST + Plate, data = df_3_gpr89_rank)
summary(gpr89_3_length) # F-statistic: 18.75 on 7 and 308 DF,  p-value: < 2.2e-16,	Adjusted R-squared:  0.2829 
mean(df_3_gpr89$Length[df_3_gpr89$Genotype2=="control"])
mean(df_3_gpr89$Length[df_3_gpr89$Genotype2=="gpr89KO"])
mean(df_3_gpr89$Length[df_3_gpr89$Genotype2=="GPR89B"])
sd(df_3_gpr89$Length[df_3_gpr89$Genotype2=="control"])
sd(df_3_gpr89$Length[df_3_gpr89$Genotype2=="gpr89KO"])
sd(df_3_gpr89$Length[df_3_gpr89$Genotype2=="GPR89B"])

gpr89_3_headarea <- lm(HeadArea ~ Genotype2 + DateVAST + Plate, data = df_3_gpr89_rank)
summary(gpr89_3_headarea) # F-statistic: 32.31 on 7 and 308 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.4103 
mean(df_3_gpr89$HeadArea[df_3_gpr89$Genotype2=="control"], na.rm = T)
mean(df_3_gpr89$HeadArea[df_3_gpr89$Genotype2=="gpr89KO"], na.rm = T)
mean(df_3_srgap2$HeadArea[df_3_gpr89$Genotype2=="GPR89B"], na.rm = T)
sd(df_3_gpr89$HeadArea[df_3_gpr89$Genotype2=="control"], na.rm = T)
sd(df_3_gpr89$HeadArea[df_3_gpr89$Genotype2=="gpr89KO"], na.rm = T)
sd(df_3_gpr89$HeadArea[df_3_gpr89$Genotype2=="GPR89B"], na.rm = T)

gpr89_3_eyearea <- lm(MeanEyeArea ~ Genotype2 + DateVAST + Plate, data = df_3_gpr89_rank)
summary(gpr89_3_eyearea) # F-statistic: 14.6 on 7 and 308 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.2321 
mean(df_3_gpr89$MeanEyeArea[df_3_gpr89$Genotype2=="control"])
mean(df_3_gpr89$MeanEyeArea[df_3_gpr89$Genotype2=="gpr89KO"])
mean(df_3_gpr89$MeanEyeArea[df_3_gpr89$Genotype2=="GPR89B"])
sd(df_3_gpr89$MeanEyeArea[df_3_gpr89$Genotype2=="control"])
sd(df_3_gpr89$MeanEyeArea[df_3_gpr89$Genotype2=="gpr89KO"])
sd(df_3_gpr89$MeanEyeArea[df_3_gpr89$Genotype2=="GPR89B"])

gpr89_3_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2 + DateVAST + Plate, data = df_3_gpr89_rank)
summary(gpr89_3_headtrunkangle) # F-statistic:  14.6 on 7 and 308 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.2321
mean(df_3_gpr89$HeadTrunkAngle[df_3_gpr89$Genotype2=="control"], na.rm = T)
mean(df_3_gpr89$HeadTrunkAngle[df_3_gpr89$Genotype2=="gpr89KO"], na.rm = T)
mean(df_3_gpr89$HeadTrunkAngle[df_3_gpr89$Genotype2=="GPR89B"], na.rm = T)
sd(df_3_gpr89$HeadTrunkAngle[df_3_gpr89$Genotype2=="control"], na.rm = T)
sd(df_3_gpr89$HeadTrunkAngle[df_3_gpr89$Genotype2=="gpr89KO"], na.rm = T)
sd(df_3_gpr89$HeadTrunkAngle[df_3_gpr89$Genotype2=="GPR89B"], na.rm = T)

gpr89_3_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2 + DateVAST + Plate, data = df_3_gpr89_rank)
summary(gpr89_3_headareanoeyes) # F-statistic: 32.18 on 7 and 308 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.4093
mean(df_3_gpr89$HeadArea_noEyes[df_3_gpr89$Genotype2=="control"], na.rm = T)
mean(df_3_gpr89$HeadArea_noEyes[df_3_gpr89$Genotype2=="gpr89KO"], na.rm = T)
mean(df_3_gpr89$HeadArea_noEyes[df_3_gpr89$Genotype2=="GPR89B"], na.rm = T)
sd(df_3_gpr89$HeadArea_noEyes[df_3_gpr89$Genotype2=="control"], na.rm = T)
sd(df_3_gpr89$HeadArea_noEyes[df_3_gpr89$Genotype2=="gpr89KO"], na.rm = T)
sd(df_3_gpr89$HeadArea_noEyes[df_3_gpr89$Genotype2=="GPR89B"], na.rm = T)

# FAM72:
df_3_fam72 <- df_3_filt[df_3_filt$DateVAST %in% names(which(summary(df_3_filt[df_3_filt$Group=="FAM72",]$DateVAST) > 0)),]
df_3_fam72 <- df_3_fam72[df_3_fam72$Group == "control" | df_3_fam72$Group == "FAM72",]
df_3_fam72 <- droplevels(df_3_fam72)
df_3_fam72$Genotype2 <- relevel(df_3_fam72$Genotype2, ref = "control")
table(df_3_fam72$Genotype2)
#control  FAM72B fam72KO 
#95       27       66 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_3_fam72_rank <- df_3_fam72
df_3_fam72_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_3_fam72_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

fam72_3_length <- lm(Length ~ Genotype2 + DateVAST + Plate, data = df_3_fam72_rank)
summary(fam72_3_length) # F-statistic: 1.57 on 7 and 181 DF,  p-value: < 2.2e-16,	Adjusted R-squared:  0.4337 
mean(df_3_fam72$Length[df_3_fam72$Genotype2=="control"])
mean(df_3_fam72$Length[df_3_fam72$Genotype2=="fam72KO"])
mean(df_3_fam72$Length[df_3_fam72$Genotype2=="FAM72B"])
sd(df_3_fam72$Length[df_3_fam72$Genotype2=="control"])
sd(df_3_fam72$Length[df_3_fam72$Genotype2=="fam72KO"])
sd(df_3_fam72$Length[df_3_fam72$Genotype2=="FAM72B"])

fam72_3_headarea <- lm(HeadArea ~ Genotype2 + DateVAST + Plate, data = df_3_fam72_rank)
summary(fam72_3_headarea) # F-statistic: 18.93 on 7 and 181 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.4003 
mean(df_3_fam72$HeadArea[df_3_fam72$Genotype2=="control"], na.rm = T)
mean(df_3_fam72$HeadArea[df_3_fam72$Genotype2=="fam72KO"], na.rm = T)
mean(df_3_fam72$HeadArea[df_3_fam72$Genotype2=="FAM72B"], na.rm = T)
sd(df_3_fam72$HeadArea[df_3_fam72$Genotype2=="control"], na.rm = T)
sd(df_3_fam72$HeadArea[df_3_fam72$Genotype2=="fam72KO"], na.rm = T)
sd(df_3_fam72$HeadArea[df_3_fam72$Genotype2=="FAM72B"], na.rm = T)

fam72_3_eyearea <- lm(MeanEyeArea ~ Genotype2 + DateVAST + Plate, data = df_3_fam72_rank)
summary(fam72_3_eyearea) # F-statistic: 16.35 on 7 and 181 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.3637 
mean(df_3_fam72$MeanEyeArea[df_3_fam72$Genotype2=="control"])
mean(df_3_fam72$MeanEyeArea[df_3_fam72$Genotype2=="fam72KO"])
mean(df_3_fam72$MeanEyeArea[df_3_fam72$Genotype2=="FAM72B"], na.rm=T)
sd(df_3_fam72$MeanEyeArea[df_3_fam72$Genotype2=="control"])
sd(df_3_fam72$MeanEyeArea[df_3_fam72$Genotype2=="fam72KO"])
sd(df_3_fam72$MeanEyeArea[df_3_fam72$Genotype2=="FAM72B"], na.rm = T)

fam72_3_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2 + DateVAST + Plate, data = df_3_fam72_rank)
summary(fam72_3_headtrunkangle) # F-statistic: 6.192 on 7 and 181 DF,  p-value: 1.651e-06, Adjusted R-squared:  0.162
mean(df_3_fam72$HeadTrunkAngle[df_3_fam72$Genotype2=="control"], na.rm = T)
mean(df_3_fam72$HeadTrunkAngle[df_3_fam72$Genotype2=="fam72KO"], na.rm = T)
mean(df_3_fam72$HeadTrunkAngle[df_3_fam72$Genotype2=="FAM72B"], na.rm = T)
sd(df_3_fam72$HeadTrunkAngle[df_3_fam72$Genotype2=="control"], na.rm = T)
sd(df_3_fam72$HeadTrunkAngle[df_3_fam72$Genotype2=="fam72KO"], na.rm = T)
sd(df_3_fam72$HeadTrunkAngle[df_3_fam72$Genotype2=="FAM72B"], na.rm = T)

fam72_3_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2 + DateVAST + Plate, data = df_3_fam72_rank)
summary(fam72_3_headareanoeyes) # F-statistic: 16.51 on 7 and 181 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.3661
mean(df_3_fam72$HeadArea_noEyes[df_3_fam72$Genotype2=="control"], na.rm = T)
mean(df_3_fam72$HeadArea_noEyes[df_3_fam72$Genotype2=="fam72KO"], na.rm = T)
mean(df_3_fam72$HeadArea_noEyes[df_3_fam72$Genotype2=="FAM72B"], na.rm = T)
sd(df_3_fam72$HeadArea_noEyes[df_3_fam72$Genotype2=="control"], na.rm = T)
sd(df_3_fam72$HeadArea_noEyes[df_3_fam72$Genotype2=="fam72KO"], na.rm = T)
sd(df_3_fam72$HeadArea_noEyes[df_3_fam72$Genotype2=="FAM72B"], na.rm = T)

# FRMPD2:
df_3_frmpd2 <- df_3_filt[df_3_filt$DateVAST %in% names(which(summary(df_3_filt[df_3_filt$Group=="FRMPD2",]$DateVAST) > 0)),]
df_3_frmpd2 <- df_3_frmpd2[df_3_frmpd2$Group == "control" | df_3_frmpd2$Group == "FRMPD2",]
df_3_frmpd2 <- droplevels(df_3_frmpd2)
df_3_frmpd2$Genotype2 <- relevel(df_3_frmpd2$Genotype2, ref = "control")
table(df_3_frmpd2$Genotype2)
#control  FRMPD2B frmpd2KO 
#143       80       93 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_3_frmpd2_rank <- df_3_frmpd2
df_3_frmpd2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_3_frmpd2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

fam72_3_length <- lm(Length ~ Genotype2 + DateVAST + Plate, data = df_3_frmpd2_rank)
summary(fam72_3_length) # F-statistic: 8.663 on 7 and 308 DF,  p-value: 1.036e-09,	Adjusted R-squared:  0.1455 
mean(df_3_frmpd2$Length[df_3_frmpd2$Genotype2=="control"])
mean(df_3_frmpd2$Length[df_3_frmpd2$Genotype2=="frmpd2KO"])
mean(df_3_frmpd2$Length[df_3_frmpd2$Genotype2=="FRMPD2B"])
sd(df_3_frmpd2$Length[df_3_frmpd2$Genotype2=="control"])
sd(df_3_frmpd2$Length[df_3_frmpd2$Genotype2=="frmpd2KO"])
sd(df_3_frmpd2$Length[df_3_frmpd2$Genotype2=="FRMPD2B"])

frmpd2_3_headarea <- lm(HeadArea ~ Genotype2 + DateVAST + Plate, data = df_3_frmpd2_rank)
summary(frmpd2_3_headarea) # F-statistic: 19.27 on 7 and 308 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.2888 
mean(df_3_frmpd2$HeadArea[df_3_frmpd2$Genotype2=="control"], na.rm = T)
mean(df_3_frmpd2$HeadArea[df_3_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
mean(df_3_frmpd2$HeadArea[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)
sd(df_3_frmpd2$HeadArea[df_3_frmpd2$Genotype2=="control"], na.rm = T)
sd(df_3_frmpd2$HeadArea[df_3_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
sd(df_3_frmpd2$HeadArea[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)

frmpd2_3_eyearea <- lm(MeanEyeArea ~ Genotype2 + DateVAST + Plate, data = df_3_frmpd2_rank)
summary(frmpd2_3_eyearea) # F-statistic: 8.066 on 7 and 308 DF,  p-value: 5.201e-09, Adjusted R-squared:  0.1357 
mean(df_3_frmpd2$MeanEyeArea[df_3_frmpd2$Genotype2=="control"])
mean(df_3_frmpd2$MeanEyeArea[df_3_frmpd2$Genotype2=="frmpd2KO"])
mean(df_3_frmpd2$MeanEyeArea[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm=T)
sd(df_3_frmpd2$MeanEyeArea[df_3_frmpd2$Genotype2=="control"])
sd(df_3_frmpd2$MeanEyeArea[df_3_frmpd2$Genotype2=="frmpd2KO"])
sd(df_3_frmpd2$MeanEyeArea[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)

frmpd2_3_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2 + DateVAST + Plate, data = df_3_frmpd2_rank)
summary(frmpd2_3_headtrunkangle) # F-statistic: 7.213 on 7 and 308 DF,  p-value: 5.261e-08, Adjusted R-squared:  0.1213
mean(df_3_frmpd2$HeadTrunkAngle[df_3_frmpd2$Genotype2=="control"], na.rm = T)
mean(df_3_frmpd2$HeadTrunkAngle[df_3_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
mean(df_3_frmpd2$HeadTrunkAngle[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)
sd(df_3_frmpd2$HeadTrunkAngle[df_3_frmpd2$Genotype2=="control"], na.rm = T)
sd(df_3_frmpd2$HeadTrunkAngle[df_3_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
sd(df_3_frmpd2$HeadTrunkAngle[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)

frmpd2_3_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2 + DateVAST + Plate, data = df_3_frmpd2_rank)
summary(frmpd2_3_headareanoeyes) # F-statistic: 19.44 on 7 and 308 DF,  p-value: < 2.2e-16, Adjusted R-squared:  0.2907
mean(df_3_frmpd2$HeadArea_noEyes[df_3_frmpd2$Genotype2=="control"], na.rm = T)
mean(df_3_frmpd2$HeadArea_noEyes[df_3_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
mean(df_3_frmpd2$HeadArea_noEyes[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)
sd(df_3_frmpd2$HeadArea_noEyes[df_3_frmpd2$Genotype2=="control"], na.rm = T)
sd(df_3_frmpd2$HeadArea_noEyes[df_3_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
sd(df_3_frmpd2$HeadArea_noEyes[df_3_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)

## 5 dpf:

# Are controls different?
wilcox.test(df_5_filt$Length[df_5_filt$Genotype=="eGFP"], df_5_filt$Length[df_5_filt$Genotype=="scrambled"]) # W = 21050, p-value = 0.6263
wilcox.test(df_5_filt$HeadArea[df_5_filt$Genotype=="eGFP"], df_5_filt$HeadArea[df_5_filt$Genotype=="scrambled"]) # W = 22054, p-value = 0.1811
wilcox.test(df_5_filt$MeanEyeArea[df_5_filt$Genotype=="eGFP"], df_5_filt$MeanEyeArea[df_5_filt$Genotype=="scrambled"]) # W = 21568, p-value = 0.3547
wilcox.test(df_5_filt$HeadTrunkAngle[df_5_filt$Genotype=="eGFP"], df_5_filt$HeadTrunkAngle[df_5_filt$Genotype=="scrambled"]) # W = 17950, p-value = 0.1079
wilcox.test(df_5_filt$HeadArea_noEyes[df_5_filt$Genotype=="eGFP"], df_5_filt$HeadArea_noEyes[df_5_filt$Genotype=="scrambled"]) # W = 22092, p-value = 0.1707

# Extract the data for each genotypic group (because we need to account for their proper controls)

# SRGAP2:
df_5_srgap2 <- df_5_filt[df_5_filt$DateVAST %in% names(which(summary(df_5_filt[df_5_filt$Group=="SRGAP2",]$DateVAST) > 0)),]
df_5_srgap2 <- df_5_srgap2[df_5_srgap2$Group == "control" | df_5_srgap2$Group == "SRGAP2",]
df_5_srgap2 <- droplevels(df_5_srgap2)
df_5_srgap2$Genotype2 <- relevel(df_5_srgap2$Genotype2, ref = "control")
table(df_5_srgap2$Genotype2)
#control  SRGAP2C srgap2KO 
#81       42       33 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_5_srgap2_rank <- df_5_srgap2
df_5_srgap2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_5_srgap2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

srgap2_5_length <- lm(Length ~ Genotype2 + Plate, data = df_5_srgap2_rank)
summary(srgap2_5_length) # F-statistic: 10.72 on 4 and 151 DF,  p-value: 1.131e-07,	Adjusted R-squared:  0.2005 
mean(df_5_srgap2$Length[df_5_srgap2$Genotype2=="control"])
mean(df_5_srgap2$Length[df_5_srgap2$Genotype2=="srgap2KO"])
mean(df_5_srgap2$Length[df_5_srgap2$Genotype2=="SRGAP2C"])
sd(df_5_srgap2$Length[df_5_srgap2$Genotype2=="control"])
sd(df_5_srgap2$Length[df_5_srgap2$Genotype2=="srgap2KO"])
sd(df_5_srgap2$Length[df_5_srgap2$Genotype2=="SRGAP2C"])

srgap2_5_headarea <- lm(HeadArea ~ Genotype2 + Plate, data = df_5_srgap2_rank)
summary(srgap2_5_headarea) # F-statistic: 3.909 on 4 and 149 DF,  p-value: 0.00477, Adjusted R-squared:  0.07067 
mean(df_5_srgap2$HeadArea[df_5_srgap2$Genotype2=="control"])
mean(df_5_srgap2$HeadArea[df_5_srgap2$Genotype2=="srgap2KO"], na.rm=T)
mean(df_5_srgap2$HeadArea[df_5_srgap2$Genotype2=="SRGAP2C"])
sd(df_5_srgap2$HeadArea[df_5_srgap2$Genotype2=="control"])
sd(df_5_srgap2$HeadArea[df_5_srgap2$Genotype2=="srgap2KO"], na.rm=T)
sd(df_5_srgap2$HeadArea[df_5_srgap2$Genotype2=="SRGAP2C"])

srgap2_5_eyearea <- lm(MeanEyeArea ~ Genotype2 + Plate, data = df_5_srgap2_rank)
summary(srgap2_5_eyearea) # F-statistic: 2.957 on 4 and 149 DF,  p-value: 0.02187, Adjusted R-squared:  0.04868 
mean(df_5_srgap2$MeanEyeArea[df_5_srgap2$Genotype2=="control"])
mean(df_5_srgap2$MeanEyeArea[df_5_srgap2$Genotype2=="srgap2KO"])
mean(df_5_srgap2$MeanEyeArea[df_5_srgap2$Genotype2=="SRGAP2C"])
sd(df_5_srgap2$MeanEyeArea[df_5_srgap2$Genotype2=="control"])
sd(df_5_srgap2$MeanEyeArea[df_5_srgap2$Genotype2=="srgap2KO"])
sd(df_5_srgap2$MeanEyeArea[df_5_srgap2$Genotype2=="SRGAP2C"])

srgap2_5_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2 + Plate, data = df_5_srgap2_rank)
summary(srgap2_5_headtrunkangle) # F-statistic: 1.545 on 4 and 149 DF,  p-value: 0.1922, Adjusted R-squared:  0.01404
mean(df_5_srgap2$HeadTrunkAngle[df_5_srgap2$Genotype2=="control"], na.rm = T)
mean(df_5_srgap2$HeadTrunkAngle[df_5_srgap2$Genotype2=="srgap2KO"], na.rm = T)
mean(df_5_srgap2$HeadTrunkAngle[df_5_srgap2$Genotype2=="SRGAP2C"], na.rm = T)
sd(df_5_srgap2$HeadTrunkAngle[df_5_srgap2$Genotype2=="control"], na.rm = T)
sd(df_5_srgap2$HeadTrunkAngle[df_5_srgap2$Genotype2=="srgap2KO"], na.rm = T)
sd(df_5_srgap2$HeadTrunkAngle[df_5_srgap2$Genotype2=="SRGAP2C"], na.rm = T)

srgap2_5_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2 + Plate, data = df_5_srgap2_rank)
summary(srgap2_5_headareanoeyes) # F-statistic: 3.935 on 4 and 149 DF,  p-value: 0.004577, Adjusted R-squared:  0.07125
mean(df_5_srgap2$HeadArea_noEyes[df_5_srgap2$Genotype2=="control"], na.rm = T)
mean(df_5_srgap2$HeadArea_noEyes[df_5_srgap2$Genotype2=="srgap2KO"], na.rm = T)
mean(df_5_srgap2$HeadArea_noEyes[df_5_srgap2$Genotype2=="SRGAP2C"], na.rm = T)
sd(df_5_srgap2$HeadArea_noEyes[df_5_srgap2$Genotype2=="control"], na.rm = T)
sd(df_5_srgap2$HeadArea_noEyes[df_5_srgap2$Genotype2=="srgap2KO"], na.rm = T)
sd(df_5_srgap2$HeadArea_noEyes[df_5_srgap2$Genotype2=="SRGAP2C"], na.rm = T)

# GPR89:
df_5_gpr89 <- df_5_filt[df_5_filt$DateVAST %in% names(which(summary(df_5_filt[df_5_filt$Group=="GPR89",]$DateVAST) > 0)),]
df_5_gpr89 <- df_5_gpr89[df_5_gpr89$Group == "control" | df_5_gpr89$Group == "GPR89",]
df_5_gpr89 <- droplevels(df_5_gpr89)
df_5_gpr89$Genotype2 <- relevel(df_5_gpr89$Genotype2, ref = "control")
table(df_5_gpr89$Genotype2)
#control  GPR89B gpr89KO 
#37       60       65 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_5_gpr89_rank <- df_5_gpr89
df_5_gpr89_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_5_gpr89_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

gpr89_5_length <- lm(Length ~ Genotype2 + Plate, data = df_5_gpr89_rank)
summary(gpr89_5_length) # F-statistic: 1.594 on 4 and 157 DF,  p-value: 0.1784,	Adjusted R-squared:  0.01455 
mean(df_5_gpr89$Length[df_5_gpr89$Genotype2=="control"])
mean(df_5_gpr89$Length[df_5_gpr89$Genotype2=="gpr89KO"])
mean(df_5_gpr89$Length[df_5_gpr89$Genotype2=="GPR89B"])
sd(df_5_gpr89$Length[df_5_gpr89$Genotype2=="control"])
sd(df_5_gpr89$Length[df_5_gpr89$Genotype2=="gpr89KO"])
sd(df_5_gpr89$Length[df_5_gpr89$Genotype2=="GPR89B"])

gpr89_5_headarea <- lm(HeadArea ~ Genotype2 + Plate, data = df_5_gpr89_rank)
summary(gpr89_5_headarea) # F-statistic: 2.448 on 4 and 157 DF,  p-value: 0.04852, Adjusted R-squared:  0.03474 
mean(df_5_gpr89$HeadArea[df_5_gpr89$Genotype2=="control"], na.rm = T)
mean(df_5_gpr89$HeadArea[df_5_gpr89$Genotype2=="gpr89KO"], na.rm = T)
mean(df_5_srgap2$HeadArea[df_5_gpr89$Genotype2=="GPR89B"], na.rm = T)
sd(df_5_gpr89$HeadArea[df_5_gpr89$Genotype2=="control"], na.rm = T)
sd(df_5_gpr89$HeadArea[df_5_gpr89$Genotype2=="gpr89KO"], na.rm = T)
sd(df_5_gpr89$HeadArea[df_5_gpr89$Genotype2=="GPR89B"], na.rm = T)

gpr89_5_eyearea <- lm(MeanEyeArea ~ Genotype2 + Plate, data = df_5_gpr89_rank)
summary(gpr89_5_eyearea) # F-statistic:  4.29 on 4 and 157 DF,  p-value: 0.002541, Adjusted R-squared:  0.07556 
mean(df_5_gpr89$MeanEyeArea[df_5_gpr89$Genotype2=="control"])
mean(df_5_gpr89$MeanEyeArea[df_5_gpr89$Genotype2=="gpr89KO"])
mean(df_5_gpr89$MeanEyeArea[df_5_gpr89$Genotype2=="GPR89B"])
sd(df_5_gpr89$MeanEyeArea[df_5_gpr89$Genotype2=="control"])
sd(df_5_gpr89$MeanEyeArea[df_5_gpr89$Genotype2=="gpr89KO"])
sd(df_5_gpr89$MeanEyeArea[df_5_gpr89$Genotype2=="GPR89B"])

gpr89_5_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2 + Plate, data = df_5_gpr89_rank)
summary(gpr89_5_headtrunkangle) # F-statistic: 1.852 on 4 and 157 DF,  p-value: 0.1215, Adjusted R-squared:  0.02074
mean(df_5_gpr89$HeadTrunkAngle[df_5_gpr89$Genotype2=="control"], na.rm = T)
mean(df_5_gpr89$HeadTrunkAngle[df_5_gpr89$Genotype2=="gpr89KO"], na.rm = T)
mean(df_5_gpr89$HeadTrunkAngle[df_5_gpr89$Genotype2=="GPR89B"], na.rm = T)
sd(df_5_gpr89$HeadTrunkAngle[df_5_gpr89$Genotype2=="control"], na.rm = T)
sd(df_5_gpr89$HeadTrunkAngle[df_5_gpr89$Genotype2=="gpr89KO"], na.rm = T)
sd(df_5_gpr89$HeadTrunkAngle[df_5_gpr89$Genotype2=="GPR89B"], na.rm = T)

gpr89_5_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2 + Plate, data = df_5_gpr89_rank)
summary(gpr89_5_headareanoeyes) # F-statistic: 2.461 on 4 and 157 DF,  p-value: 0.0476, Adjusted R-squared:  0.03502
mean(df_5_gpr89$HeadArea_noEyes[df_5_gpr89$Genotype2=="control"], na.rm = T)
mean(df_5_gpr89$HeadArea_noEyes[df_5_gpr89$Genotype2=="gpr89KO"], na.rm = T)
mean(df_5_gpr89$HeadArea_noEyes[df_5_gpr89$Genotype2=="GPR89B"], na.rm = T)
sd(df_5_gpr89$HeadArea_noEyes[df_5_gpr89$Genotype2=="control"], na.rm = T)
sd(df_5_gpr89$HeadArea_noEyes[df_5_gpr89$Genotype2=="gpr89KO"], na.rm = T)
sd(df_5_gpr89$HeadArea_noEyes[df_5_gpr89$Genotype2=="GPR89B"], na.rm = T)

# FAM72:
df_5_fam72 <- df_5_filt[df_5_filt$DateVAST %in% names(which(summary(df_5_filt[df_5_filt$Group=="FAM72",]$DateVAST) > 0)),]
df_5_fam72 <- df_5_fam72[df_5_fam72$Group == "control" | df_5_fam72$Group == "FAM72",]
df_5_fam72 <- droplevels(df_5_fam72)
df_5_fam72$Genotype2 <- relevel(df_5_fam72$Genotype2, ref = "control")
table(df_5_fam72$Genotype2)
#control  FAM72B fam72KO 
#33       10       38 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_5_fam72_rank <- df_5_fam72
df_5_fam72_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_5_fam72_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

fam72_5_length <- lm(Length ~ Genotype2, data = df_5_fam72_rank)
summary(fam72_5_length) # F-statistic: 6.083 on 2 and 78 DF,  p-value: 0.003507,	Adjusted R-squared:  0.1128 
mean(df_5_fam72$Length[df_5_fam72$Genotype2=="control"])
mean(df_5_fam72$Length[df_5_fam72$Genotype2=="fam72KO"])
mean(df_5_fam72$Length[df_5_fam72$Genotype2=="FAM72B"])
sd(df_5_fam72$Length[df_5_fam72$Genotype2=="control"])
sd(df_5_fam72$Length[df_5_fam72$Genotype2=="fam72KO"])
sd(df_5_fam72$Length[df_5_fam72$Genotype2=="FAM72B"])

fam72_5_headarea <- lm(HeadArea ~ Genotype2, data = df_5_fam72_rank)
summary(fam72_5_headarea) # F-statistic: 1.604 on 2 and 78 DF,  p-value: 0.2077, Adjusted R-squared:  0.01487 
mean(df_5_fam72$HeadArea[df_5_fam72$Genotype2=="control"], na.rm = T)
mean(df_5_fam72$HeadArea[df_5_fam72$Genotype2=="fam72KO"], na.rm = T)
mean(df_5_fam72$HeadArea[df_5_fam72$Genotype2=="FAM72B"], na.rm = T)
sd(df_5_fam72$HeadArea[df_5_fam72$Genotype2=="control"], na.rm = T)
sd(df_5_fam72$HeadArea[df_5_fam72$Genotype2=="fam72KO"], na.rm = T)
sd(df_5_fam72$HeadArea[df_5_fam72$Genotype2=="FAM72B"], na.rm = T)

fam72_5_eyearea <- lm(MeanEyeArea ~ Genotype2, data = df_5_fam72_rank)
summary(fam72_5_eyearea) # F-statistic: 4.601 on 2 and 78 DF,  p-value: 0.01292, Adjusted R-squared:  0.08258 
mean(df_5_fam72$MeanEyeArea[df_5_fam72$Genotype2=="control"])
mean(df_5_fam72$MeanEyeArea[df_5_fam72$Genotype2=="fam72KO"])
mean(df_5_fam72$MeanEyeArea[df_5_fam72$Genotype2=="FAM72B"], na.rm=T)
sd(df_5_fam72$MeanEyeArea[df_5_fam72$Genotype2=="control"])
sd(df_5_fam72$MeanEyeArea[df_5_fam72$Genotype2=="fam72KO"])
sd(df_5_fam72$MeanEyeArea[df_5_fam72$Genotype2=="FAM72B"], na.rm = T)

fam72_5_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2, data = df_5_fam72_rank)
summary(fam72_5_headtrunkangle) # F-statistic: 6.455 on 2 and 78 DF,  p-value: 0.002547, Adjusted R-squared:  0.12
mean(df_5_fam72$HeadTrunkAngle[df_5_fam72$Genotype2=="control"], na.rm = T)
mean(df_5_fam72$HeadTrunkAngle[df_5_fam72$Genotype2=="fam72KO"], na.rm = T)
mean(df_5_fam72$HeadTrunkAngle[df_5_fam72$Genotype2=="FAM72B"], na.rm = T)
sd(df_5_fam72$HeadTrunkAngle[df_5_fam72$Genotype2=="control"], na.rm = T)
sd(df_5_fam72$HeadTrunkAngle[df_5_fam72$Genotype2=="fam72KO"], na.rm = T)
sd(df_5_fam72$HeadTrunkAngle[df_5_fam72$Genotype2=="FAM72B"], na.rm = T)

fam72_5_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2, data = df_5_fam72_rank)
summary(fam72_5_headareanoeyes) # F-statistic: 1.414 on 2 and 78 DF,  p-value: 0.2494, Adjusted R-squared:  0.01024
mean(df_5_fam72$HeadArea_noEyes[df_5_fam72$Genotype2=="control"], na.rm = T)
mean(df_5_fam72$HeadArea_noEyes[df_5_fam72$Genotype2=="fam72KO"], na.rm = T)
mean(df_5_fam72$HeadArea_noEyes[df_5_fam72$Genotype2=="FAM72B"], na.rm = T)
sd(df_5_fam72$HeadArea_noEyes[df_5_fam72$Genotype2=="control"], na.rm = T)
sd(df_5_fam72$HeadArea_noEyes[df_5_fam72$Genotype2=="fam72KO"], na.rm = T)
sd(df_5_fam72$HeadArea_noEyes[df_5_fam72$Genotype2=="FAM72B"], na.rm = T)

# FRMPD2:
df_5_frmpd2 <- df_5_filt[df_5_filt$DateVAST %in% names(which(summary(df_5_filt[df_5_filt$Group=="FRMPD2",]$DateVAST) > 0)),]
df_5_frmpd2 <- df_5_frmpd2[df_5_frmpd2$Group == "control" | df_5_frmpd2$Group == "FRMPD2",]
df_5_frmpd2 <- droplevels(df_5_frmpd2)
df_5_frmpd2$Genotype2 <- relevel(df_5_frmpd2$Genotype2, ref = "control")
table(df_5_frmpd2$Genotype2)
#control  FRMPD2B frmpd2KO 
#32       61       55 

# Because of non-normality, will use an ANCOVA with rank-transformation:
df_5_frmpd2_rank <- df_5_frmpd2
df_5_frmpd2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")] <- lapply(df_5_frmpd2_rank[,c("Length", "HeadArea", "MeanEyeArea", "HeadTrunkAngle", "HeadArea_noEyes")], rank)

fam72_5_length <- lm(Length ~ Genotype2 + Plate, data = df_5_frmpd2_rank)
summary(fam72_5_length) # F-statistic: 6.099 on 4 and 143 DF,  p-value: 0.0001467,	Adjusted R-squared:  0.1218 
mean(df_5_frmpd2$Length[df_5_frmpd2$Genotype2=="control"])
mean(df_5_frmpd2$Length[df_5_frmpd2$Genotype2=="frmpd2KO"])
mean(df_5_frmpd2$Length[df_5_frmpd2$Genotype2=="FRMPD2B"])
sd(df_5_frmpd2$Length[df_5_frmpd2$Genotype2=="control"])
sd(df_5_frmpd2$Length[df_5_frmpd2$Genotype2=="frmpd2KO"])
sd(df_5_frmpd2$Length[df_5_frmpd2$Genotype2=="FRMPD2B"])

frmpd2_5_headarea <- lm(HeadArea ~ Genotype2 + Plate, data = df_5_frmpd2_rank)
summary(frmpd2_5_headarea) # F-statistic: 10.12 on 4 and 143 DF,  p-value: 3.034e-07, Adjusted R-squared:  0.1989 
mean(df_5_frmpd2$HeadArea[df_5_frmpd2$Genotype2=="control"], na.rm = T)
mean(df_5_frmpd2$HeadArea[df_5_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
mean(df_5_frmpd2$HeadArea[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)
sd(df_5_frmpd2$HeadArea[df_5_frmpd2$Genotype2=="control"], na.rm = T)
sd(df_5_frmpd2$HeadArea[df_5_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
sd(df_5_frmpd2$HeadArea[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)

frmpd2_5_eyearea <- lm(MeanEyeArea ~ Genotype2 + Plate, data = df_5_frmpd2_rank)
summary(frmpd2_5_eyearea) # F-statistic: 9.045 on 4 and 143 DF,  p-value: 1.53e-06, Adjusted R-squared:  0.1796 
mean(df_5_frmpd2$MeanEyeArea[df_5_frmpd2$Genotype2=="control"])
mean(df_5_frmpd2$MeanEyeArea[df_5_frmpd2$Genotype2=="frmpd2KO"])
mean(df_5_frmpd2$MeanEyeArea[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm=T)
sd(df_5_frmpd2$MeanEyeArea[df_5_frmpd2$Genotype2=="control"])
sd(df_5_frmpd2$MeanEyeArea[df_5_frmpd2$Genotype2=="frmpd2KO"])
sd(df_5_frmpd2$MeanEyeArea[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)

frmpd2_5_headtrunkangle <- lm(HeadTrunkAngle ~ Genotype2 + Plate, data = df_5_frmpd2_rank)
summary(frmpd2_5_headtrunkangle) # F-statistic: 2.844 on 4 and 143 DF,  p-value: 0.02631, Adjusted R-squared:  0.04778
mean(df_5_frmpd2$HeadTrunkAngle[df_5_frmpd2$Genotype2=="control"], na.rm = T)
mean(df_5_frmpd2$HeadTrunkAngle[df_5_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
mean(df_5_frmpd2$HeadTrunkAngle[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)
sd(df_5_frmpd2$HeadTrunkAngle[df_5_frmpd2$Genotype2=="control"], na.rm = T)
sd(df_5_frmpd2$HeadTrunkAngle[df_5_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
sd(df_5_frmpd2$HeadTrunkAngle[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)

frmpd2_5_headareanoeyes <- lm(HeadArea_noEyes ~ Genotype2 + Plate, data = df_5_frmpd2_rank)
summary(frmpd2_5_headareanoeyes) # F-statistic: 9.884 on 4 and 143 DF,  p-value: 4.336e-07, Adjusted R-squared:  0.1947
mean(df_5_frmpd2$HeadArea_noEyes[df_5_frmpd2$Genotype2=="control"], na.rm = T)
mean(df_5_frmpd2$HeadArea_noEyes[df_5_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
mean(df_5_frmpd2$HeadArea_noEyes[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)
sd(df_5_frmpd2$HeadArea_noEyes[df_5_frmpd2$Genotype2=="control"], na.rm = T)
sd(df_5_frmpd2$HeadArea_noEyes[df_5_frmpd2$Genotype2=="frmpd2KO"], na.rm = T)
sd(df_5_frmpd2$HeadArea_noEyes[df_5_frmpd2$Genotype2=="FRMPD2B"], na.rm = T)


### Plots:

# Detailed plots for GPR89 and FRMPD2:

# Load pre-processed summary data and perform p-value adjustment:
sum_df <- read.csv("Final_results_summary.csv", head=T)
sum_df$p.value.adj <- p.adjust(sum_df$p.value, method="BH")

# Head area plots for the individual genes FRMPD2 and GPR89:

df_3_gpr89 <- df_3_filt[df_3_filt$DateVAST %in% names(which(summary(df_3_filt[df_3_filt$Group=="GPR89",]$DateVAST) > 0)),]
df_3_gpr89 <- df_3_gpr89[df_3_gpr89$Group %in% c("control", "GPR89") ,]
df_3_gpr89 <- droplevels(df_3_gpr89)
df_3_gpr89$Genotype2 <- factor(df_3_gpr89$Genotype2, levels = c("control", "GPR89B", "gpr89KO"))
table(df_3_gpr89$Genotype2)
#control   GPR89B  gpr89KO 
#175       62       79      

p_3_gpr89 <- ggplot(data= df_3_gpr89, aes(x= Genotype2, y= HeadArea_noEyes, col= Genotype2)) + 
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.1) + geom_jitter(width=0.05, alpha=0.5) + 
  theme_bw() + xlab("") + ylab("Head area (um2)") + theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = c("black", "darkblue", "blue"))

ggsave("HeadArea_individual_plots_3dpf_GPR89.pdf",
       p_3_gpr89, height = 3, width = 4)

df_3_frmpd2 <- df_3_filt[df_3_filt$DateVAST %in% names(which(summary(df_3_filt[df_3_filt$Group=="FRMPD2",]$DateVAST) > 0)),]
df_3_frmpd2 <- df_3_frmpd2[df_3_frmpd2$Group %in% c("control", "FRMPD2") ,]
df_3_frmpd2 <- droplevels(df_3_frmpd2)
df_3_frmpd2$Genotype2 <- factor(df_3_frmpd2$Genotype2, levels = c("control", "FRMPD2B", "frmpd2KO"))
table(df_3_frmpd2$Genotype2)
#control   FRMPD2B  frmpd2KO 
#143       80       93      

p_3_frmpd2 <- ggplot(data= df_3_frmpd2, aes(x= Genotype2, y= HeadArea_noEyes, col= Genotype2)) + 
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.1) + geom_jitter(width=0.05, alpha=0.5) + 
  theme_bw() + xlab("") + ylab("Head area (um2)") + theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = c("black", "darkblue", "blue"))

ggsave("HeadArea_individual_plots_3dpf_FRMPD2.pdf",
       p_3_frmpd2, height = 3, width = 4)

df_ind_5 <- df_5_filt[df_5_filt$DateVAST %in% names(which(summary(df_5_filt[df_5_filt$Group=="GPR89" | df_5_filt$Group=="FRMPD2",]$DateVAST) > 0)),]
df_ind_5 <- df_ind_5[df_ind_5$Group %in% c("control", "FRMPD2", "GPR89") ,]
df_ind_5 <- droplevels(df_ind_5)
df_ind_5$Genotype2 <- factor(df_ind_5$Genotype2, levels = c("control", "GPR89B", "gpr89KO", "FRMPD2B", "frmpd2KO"))
table(df_ind_5$Genotype2)
#control   GPR89B  gpr89KO  FRMPD2B frmpd2KO 
#69       60       65       62       55 

p_ind5 <- ggplot(data= df_ind_5, aes(x= Genotype2, y= HeadArea_noEyes, col= Genotype2)) + 
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.1) + geom_jitter(width=0.05, alpha=0.5) + 
  theme_bw() + xlab("") + ylab("Head area (um2)") + theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = c("black", "darkred", "red", "darkblue", "blue"))

# Heatmap with all results:

heatmap_3dpf <- ggplot(sum_df[sum_df$Age=="3",], aes(x = Genotype, y = Feature, fill = delta)) +
  geom_tile(aes(fill = delta), color = "white") +
  geom_text(data = subset(sum_df[sum_df$Age=="3",], p.value.adj < 0.1), aes(label = "*"), color = "black", size= 10) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x= "", y="", fill = "Delta") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heatmap_5dpf <- ggplot(sum_df[sum_df$Age=="5",], aes(x = Genotype, y = Feature, fill = delta)) +
  geom_tile(aes(fill = delta), color = "white") +
  geom_text(data = subset(sum_df[sum_df$Age=="5",], p.value.adj < 0.1), aes(label = "*"), color = "black", size= 10) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x= "", y="", fill = "Delta") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Heatmap_combined.pdf",
       heatmap_3dpf + heatmap_5dpf, height = 6, width = 12)



# Plots for all:
genotypes <- c("control", "GPR89B", "gpr89KO",
               "FRMPD2B", "frmpd2KO",
               "FAM72B", "fam72KO",
               "SRGAP2C", "srgap2KO")
df_3_filt2 <- df_3[df_3$Genotype2 %in% genotypes,]
df_3_filt2 <- droplevels(df_3_filt2)
df_3_filt2$Genotype2 <- factor(df_3_filt2$Genotype2, 
                         levels = c("control", "GPR89B", "gpr89KO",
                                    "FRMPD2B", "frmpd2KO",
                                    "FAM72B", "fam72KO",
                                    "SRGAP2C", "srgap2KO"))
table(df_3_filt2$Genotype2)
p_3_head <- ggplot(data= df_3_filt2, 
                   aes(x= Genotype2, y= HeadArea_noEyes, col= Genotype2)) + 
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.1) + 
  geom_jitter(width=0.05, alpha=0.5) + 
  theme_bw() + xlab("") + ylab("Head area (um2)") + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_3_length <- ggplot(data= df_3_filt2, 
                   aes(x= Genotype2, y= Length, col= Genotype2)) + 
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.1) + 
  geom_jitter(width=0.05, alpha=0.5) + 
  theme_bw() + xlab("") + ylab("Length (um)") + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

df_5_filt2 <- df_5[df_5$Genotype2 %in% genotypes,]
df_5_filt2 <- droplevels(df_5_filt2)
df_5_filt2$Genotype2 <- factor(df_5_filt2$Genotype2, 
                               levels = c("control", "GPR89B", "gpr89KO",
                                          "FRMPD2B", "frmpd2KO",
                                          "FAM72B", "fam72KO",
                                          "SRGAP2C", "srgap2KO"))
table(df_5_filt2$Genotype2)
p_5_head <- ggplot(data= df_5_filt2, 
                   aes(x= Genotype2, y= HeadArea_noEyes, col= Genotype2)) + 
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.1) + 
  geom_jitter(width=0.05, alpha=0.5) + 
  theme_bw() + xlab("") + ylab("Head area (um2)") + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_5_length <- ggplot(data= df_5_filt2, 
                     aes(x= Genotype2, y= Length, col= Genotype2)) + 
  geom_boxplot(outlier.shape = NA, fill="grey", alpha=0.1) + 
  geom_jitter(width=0.05, alpha=0.5) + 
  theme_bw() + xlab("") + ylab("Length (um)") + 
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

```
