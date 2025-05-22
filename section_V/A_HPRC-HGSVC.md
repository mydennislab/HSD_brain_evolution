# HGSVC and HGSVC variants in pHSDs

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc
```

## 1. Data download

### 1.1 HPRC data download

Fasta assemblies were downloaded in `/share/dennislab/projects/HPRC/assemblies/`.

List of available assemblies:
```bash
/share/dennislab/programs/agc-2.1_x64-linux/agc listset /share/dennislab/projects/HPRC/assemblies/HPRC-yr1.agc | grep -v "CHM13Y" > hprc_yr1_samples.txt
# downloading missing sample
/programs/agc-2.1_x64-linux/agc getset HPRC-yr1.agc NA18906.2 > NA18906.2.fa
```

### 1.2 HGSVC data download

We downloaded HGSVC fully-phased genome assemblies from: https://zenodo.org/record/5607680#.Y-rd0-zMIwQ

The files were located in `/share/dennislab/datasets/data/HGSVC`

## 2. Variant calling

We mapped HPRC and HGSVC assemblies to T2T-CHM13v1.0, filter alignments per locus of interest (`pHSD_regions.tsv`), and selected the longest alignment per locus. These steps were automated in a Snakefile:
```bash
conda activate assembly

/share/dennislab/users/dcsoto/Miniconda3/envs/smk/bin/snakemake \
--snakefile asm-preprocessing.smk \
--config reference=00_reference/chm13.draft_v1.0.plusY.fasta \
regions=../coords/pHSDs_regions.tsv samples=samples.tsv \
-p -j 60s
```

Then, we performed variant calling per sample per locus, using  the longest alignment per region:
```bash
mkdir -p 03_variants_pHSDs

cat ../coords/pHSDs_regions.tsv | xargs -n2 -P10 bash -c '
for sample in $(cut -f1 samples.tsv); do
    echo "Procesing sample "$sample" in locus "$0
    REFERENCE="00_reference/chm13.draft_v1.0.plusY.fasta"

    # calling variants per sample
    htsbox pileup -q 0 -evcf $REFERENCE 02_alignments_pHSDs/${sample}.1.$0.paf.filt.bam 02_alignments_pHSDs/${sample}.2.$0.paf.filt.bam > 03_variants_pHSDs/${sample}.pair.$0.vcf

    # obtaining diploid calls
    dipcall-aux.js vcfpair -s ${sample} -a 03_variants_pHSDs/${sample}.pair.$0.vcf \
    | bcftools norm -Ov -m-any \
    | bcftools norm -Ov -d exact \
    | bcftools norm -Ov -m-any --fasta-ref $REFERENCE --check-ref w \
    | bcftools sort -m 4G \
    | bgzip -c > 03_variants_pHSDs/${sample}.dip.$0.vcf.gz
    tabix 03_variants_pHSDs/${sample}.dip.$0.vcf.gz
done'
```

Then, we identified fully callable samples per regions (samples where a locus of interest was fully covered by one contig per haplotype):
```bash
# callable regions per sample per locus
cat ../coords/pHSDs_regions.tsv | xargs -n2 -P10 bash -c '
for sample in $(cut -f1 samples.tsv); do
    echo "Procesing sample "$sample" in locus "$0
    REFERENCE="00_reference/chm13.draft_v1.0.plusY.fasta"

    # obtaining callable regions per sample
    cat 02_alignments_pHSDs/${sample}.1.$0.paf.filt.var | grep ^R | cut -f2- > 02_alignments_pHSDs/${sample}.1.$0.paf.filt.bed
    cat 02_alignments_pHSDs/${sample}.2.$0.paf.filt.var | grep ^R | cut -f2- > 02_alignments_pHSDs/${sample}.2.$0.paf.filt.bed
    bedtools intersect -a 02_alignments_pHSDs/${sample}.1.$0.paf.filt.bed -b 02_alignments_pHSDs/${sample}.2.$0.paf.filt.bed \
    > 03_variants_pHSDs/${sample}.dip.$0.bed
done'

# merging callable regions
cat ../coords/pHSDs_regions.bed | xargs -n4 -P10 bash -c '
for sample in $(cut -f1 samples.tsv); do

    bedtools intersect -f 1 -wao -a ../coords/pHSDs_regions.bed -b 03_variants_pHSDs/${sample}.dip.$3.bed \
    | grep "$3" | awk -v chr=$0 -v start=$1 -v end=$2 -v locus=$3 -v var=$sample "{{if(\$5!=\".\"){{print chr\"\t\"start\"\t\"end\"\t\"locus\"\t\"var}} }}" 
    
done' | bedtools sort | bedtools groupby -g 1,2,3,4 -c 5 -o collapse > pHSDs_regions.fully_covered.bed
```

> Note: FAM72C_SRGAP2D was not covered by any individual.

We merged VCF files per locus from fully callable samples only:
```bash
mkdir -p 04_final_VCFs_pHSDs

awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' pHSDs_regions.fully_covered.bed | xargs -n3 -P1 bash -c '
    echo "Processing locus "$1
    vcfs=$(for i in $(echo $2 | sed "s/,/ /g"); do printf "03_variants_pHSDs/${i}.dip.$1.vcf.gz "; done)

    bcftools merge -r $0 --missing-to-ref $vcfs \
    | bcftools norm -m+ \
    | bcftools norm -m- \
    | bcftools norm -f 00_reference/chm13.draft_v1.0.plusY.fasta \
    | bcftools view -f 'PASS,.' \
    | awk "{{if(\$0 ~ \"^#\" || \$5 != \"*\") {{print}} }}" \
    > 04_final_VCFs_pHSDs/$1.vcf
    bgzip -f 04_final_VCFs_pHSDs/$1.vcf; tabix -f 04_final_VCFs_pHSDs/$1.vcf.gz
'
```

## 3. HPRC+HGSVC analysis

### 3.1 Removing duplicated samples

We moved forward with HPRC assemblies plus the nine HGSVC-only assemblies.

Specifically, we removed the following samples from files:
- HG00733_HGSVC
- HG02818_HGSVC
- HG03486_HGSVC
- NA19240_HGSVC
- NA24385_HGSVC

```bash
mkdir -p 05_merged_VCFs_pHSDs

awk '{print $1":"$2"-"$3"\t"$4}' ../coords/pHSDs_regions.bed | grep -v "FAM72C_SRGAP2D" | xargs -n2 -P1 bash -c '
    bcftools view --force-samples -s ^HG00733_HGSVC,HG02818_HGSVC,HG03486_HGSVC,NA19240_HGSVC,NA24385_HGSVC -r $0 04_final_VCFs_pHSDs/$1.vcf.gz \
    | bcftools view --min-ac 1 \
    | bgzip -c > 05_merged_VCFs_pHSDs/$1.vcf.gz 
    tabix -f 05_merged_VCFs_pHSDs/$1.vcf.gz  
'
```

### 3.2 Variant counts

```bash
mkdir -p 06_variant_counts
```

Counting number of variants in SDs
```bash
# overall variants in SD regions
awk '{print $1":"$2"-"$3"\t"$4}' ../coords/pHSDs_regions.SD.bed \
| grep -v "FAM72C_SRGAP2D" \
| xargs -n2 -P1 bash -c '
    printf $1"\t"

    bcftools view -r $0 05_merged_VCFs_pHSDs/$1.vcf.gz \
    | bcftools stats \
    | grep "^SN" \
    | grep "number of samples\|number of records\|number of SNPs\|number of indels" \
    | cut -f4 \
    | tr "\n" "\t"
    printf "\n"
' > 06_variant_counts/hprc_hgsvc.SD.total_counts.tsv
```

Variant density per individual:
```bash
# heterozygous variants in SD per sample
awk '{print $1":"$2"-"$3"\t"$4}' ../coords/pHSDs_regions.SD.bed \
| grep -v "FAM72C_SRGAP2D" \
| xargs -n2 -P1 bash -c '

for smp in `bcftools query -l 05_merged_VCFs_pHSDs/$1.vcf.gz`
do
    printf "$1\t"
    printf "$smp\t"

    bcftools view -s ${smp} -c1 -r $0 05_merged_VCFs_pHSDs/$1.vcf.gz \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | grep -v "1|1\|1/1\|0/0" \
    | bcftools stats | grep "^SN" | grep "number of SNPs" | cut -d" " -f3- | cut -f2
done' > 06_variant_counts/hprc_hgsvc.SD.het_counts.tsv

# heterozygous variants in matching SD per sample
awk '{print $1":"$2"-"$3"\t"$5"\t"$6}' ../coords/pHSDs_txs.SD.fix.bed \
| grep -v "FAM72C_SRGAP2D" \
| xargs -n3 -P1 bash -c '

for smp in `bcftools query -l 05_merged_VCFs_pHSDs/$2.vcf.gz`
do
    printf "$1\t"
    printf "$smp\t"

    bcftools view -s ${smp} -c1 -r $0 05_merged_VCFs_pHSDs/$2.vcf.gz \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | grep -v "1|1\|1/1\|0/0" \
    | bcftools stats | grep "^SN" | grep "number of SNPs" | cut -d" " -f3- | cut -f2
done' > 06_variant_counts/hprc_hgsvc.SD_matching.het_counts.tsv
```

### 3.3 Matching individuals

Then, we compared derived paralogs to their ancestral counterparts using population genetics test statistic. To do so, we subsetted the VCF files to contain only matching individuals between paralogs, meaning individuals that have contigs fully covering all the members of a gene family, and filtered biallelic SNPs only.  

Obtaining matching individuals per paralog family:
```python
import pandas as pd

def get_family(x):
    if "ARHGAP11" in x: return "ARHGAP11"
    if "SRGAP2" in x: return "SRGAP2"
    if "GPR89" in x: return "GPR89"
    if "GPRIN2" in x: "GPRIN2"
    if "PTPN20" in x: "PTPN20"
    if "ROCK1" in x: return "ROCK1"
    if "CD8B" in x: return "CD8B"
    if "CFC1" in x: return "CFC1"
    if "HYDIN" in x: return "HYDIN"
    if "DUSP22" in x: return "DUSP22"
    else: return "NA"

def common_list_of_lists(lst):
    temp = set(lst[0]).intersection(*lst)
    return list(temp)

to_remove = ["HG00733_HGSVC", "HG02818_HGSVC", "HG03486_HGSVC", "NA19240_HGSVC", "NA24385_HGSVC"]

df = pd.read_csv('pHSDs_regions.fully_covered.bed', names=["chr","start","end","locus", "samples"], sep="\t")
df['family'] = [get_family(x) for x in df['locus']]
df['samples'] = [x.split(',') for x in df['samples']]
df['samples'] = [list(set(x) - set(to_remove)) for x in df['samples']]

df2 = df.groupby('family')['samples'].apply(list).reset_index().rename(columns={"samples": "family_samples"})
df3 = pd.merge(df, df2, how='left', on='family')
df3["common_samples"] = [list(set.intersection(*map(set, x))) for x in df3['family_samples']]
df3['common_samples'] = [','.join(x) for x in df3['common_samples']]

df3.drop(columns=['samples', 'family', 'family_samples']).to_csv('pHSDs_regions.fully_covered.matching.bed', sep="\t", header=False, index=False)
```

Per region, we selected only variant sites from matching individuals and in SD space:
```bash
mkdir -p 07_matching_VCF

# all matching individuals
awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' pHSDs_regions.fully_covered.matching.bed | xargs -n3 -P1 bash -c '
    echo "Processing locus "$1
    vcfs=$(for i in $(echo $2 | sed "s/,/ /g"); do printf "03_variants_pHSDs/${i}.dip.$1.vcf.gz "; done)

    bcftools merge -r $0 --missing-to-ref $vcfs \
    | bcftools norm -m+ \
    | bcftools norm -m- \
    | bcftools norm -f 00_reference/chm13.draft_v1.0.plusY.fasta \
    | bcftools view -f 'PASS,.' \
    | awk "{{if(\$0 ~ \"^#\" || \$5 != \"*\") {{print}} }}" \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | bgzip -c \
    > 07_matching_VCF/$1.SNPs.biallelic.vcf.gz
    tabix -f 07_matching_VCF/$1.SNPs.biallelic.vcf.gz
'

# intersecting matching individuals with SD space
bedtools intersect -a pHSDs_regions.fully_covered.matching.bed -b /share/dennislab/users/dcsoto/ms_hsd/sd-98/coordinates/chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed > pHSDs_regions.fully_covered.matching.SD.bed

awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' pHSDs_regions.fully_covered.matching.SD.bed | xargs -n3 -P1 bash -c '
    echo "Processing locus "$1
    vcfs=$(for i in $(echo $2 | sed "s/,/ /g"); do printf "03_variants_pHSDs/${i}.dip.$1.vcf.gz "; done)

    bcftools merge -r $0 --missing-to-ref $vcfs \
    | bcftools norm -m+ \
    | bcftools norm -m- \
    | bcftools norm -f 00_reference/chm13.draft_v1.0.plusY.fasta \
    | bcftools view -f 'PASS,.' \
    | awk "{{if(\$0 ~ \"^#\" || \$5 != \"*\") {{print}} }}" \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | bgzip -c \
    > 07_matching_VCF/$1.SD.SNPs.biallelic.vcf.gz
    tabix -f 07_matching_VCF/$1.SD.SNPs.biallelic.vcf.gz
'
```

### 3.4 Nucleotide diversity

```bash
mkdir -p 08_PopGenome_Pi
```

Pi per locus:
```r
# module load R/4.3.1
library(PopGenome)
library(data.table)
library(tidyverse)

regions <- fread("pHSDs_regions.fully_covered.matching.bed", col.names = c("chr", "start", "end", "locus", "samples"))
regions <- regions[,1:4]

results <- data.frame()
windows <- data.frame()

for(i in 1:nrow(regions)) {
    chr = regions[i,1][[1]]
    start = regions[i,2][[1]]
    end = regions[i,3][[1]]
    locus = regions[i,4][[1]]

    # overall pi
    GENOME.class <- readVCF(paste("07_matching_VCF/",locus,".SNPs.biallelic.vcf.gz", sep = ""), 
                            tid=chr, from=start, to=end, numcols=10000000)
    GENOME.class <- F_ST.stats(GENOME.class)

    df1 <- data.frame(
        locus = locus, 
        pi =  as.vector(GENOME.class@Pi/GENOME.class@n.sites)) 
    results <- rbind(results, df1)

    # sliding window
    GENOME.class.slide <- sliding.window.transform(GENOME.class, width=15000, jump=1000, type=2, whole.data=TRUE)
    GENOME.class.slide <- F_ST.stats(GENOME.class.slide)
    
    df2 <- data.frame(
        locus = locus,
        regions = GENOME.class.slide@region.names, 
        pi = as.vector(GENOME.class.slide@Pi/GENOME.class.slide@n.sites)
    )    
    windows <- rbind(windows, df2)

    #GENOME.class.slide <- diversity.stats(GENOME.class.slide, pi = TRUE)
    #GENOME.class.slide@nuc.diversity.within
}

fwrite(results, "08_PopGenome_Pi/hprc_hgsvc.overall_pi.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(windows, "08_PopGenome_Pi/hprc_hgsvc.windows_pi.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
```

### 3.5 Tajima's D

```bash
mkdir -p 09_PopGenome_Tajima
```

We calculated overall Tajima's D across the duplicated portion of each gene and compared to their corresponding paralog using PopGenome, and adding +-5kbp to each gene.

Tajima's D per locus:
```r
# module load R/4.3.1
library(PopGenome)
library(data.table)
library(tidyverse)

regions <- fread("pHSDs_regions.fully_covered.matching.bed", col.names = c("chr", "start", "end", "locus", "samples"))
regions <- regions[,1:4]

results <- data.frame()
windows <- data.frame()

for(i in 1:nrow(regions)) {
    chr = regions[i,1][[1]]
    start = regions[i,2][[1]]
    end = regions[i,3][[1]]
    locus = regions[i,4][[1]]

    # overall pi
    GENOME.class <- readVCF(paste("07_matching_VCF/",locus,".SNPs.biallelic.vcf.gz", sep = ""), 
                            tid=chr, from=start, to=end, numcols=10000000)
    GENOME.class <- neutrality.stats(GENOME.class) 

    df1 <- data.frame(
        locus = locus, 
        tajimas =  as.vector(GENOME.class@Tajima.D)) 
    results <- rbind(results, df1)

    # sliding window
    GENOME.class.slide <- sliding.window.transform(GENOME.class, width=15000, jump=1000, type=2, whole.data=TRUE)
    GENOME.class.slide <- neutrality.stats(GENOME.class.slide)
    
    df2 <- data.frame(
        locus = locus,
        regions = GENOME.class.slide@region.names, 
        tajimas = as.vector(GENOME.class.slide@Tajima.D)
    )    
    windows <- rbind(windows, df2)
}

fwrite(results, "09_PopGenome_Tajima/hprc_hgsvc.tajimasd.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(windows, "09_PopGenome_Tajima/hprc_hgsvc.windows_tajimasd.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
```

## 4. Extended Pi and Tajima's D

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc_extended
```

We repeated previous analysis, this time calling variants across selected paralogs loci and their surroundings. Coordinates were saved in file `extended.tsv`.
* GPR89A: chr1:144449050-145149169 | chr1:144297593-144800048
* GPR89B: chr1:146179939-146881662 | chr1:146376921-146880413
* FRMPD2A: chr10:48648324-49349022 | chr10:48776474-49290971
* FRMPD2B: chr10:47472911-48174201 | chr10:47562310-48062365
* CD8B: chr2:86719878-86911142
* CD8B2: chr2:106815307-107015852 
* SRGAP2C: chr1:121075201-121457259

Mapping HPRC and HGSVC assemblies to T2T-CHM13v1.0, filter alignments per locus of interest, and selected the longest alignment per locus:
```bash
conda activate assembly

/share/dennislab/users/dcsoto/Miniconda3/envs/smk/bin/snakemake \
--snakefile asm-preprocessing.smk \
--config reference=00_reference/chm13.draft_v1.0.plusY.fasta \
regions=extended.tsv samples=samples.tsv \
-p -j 60
```

Variant calling per sample per locus, using  the longest alignment per region:
```bash
mkdir -p 03_variants_pHSDs

cat extended.tsv | xargs -n2 -P10 bash -c '
for sample in $(cut -f1 samples.tsv); do
    echo "Procesing sample "$sample" in locus "$0
    REFERENCE="00_reference/chm13.draft_v1.0.plusY.fasta"

    # calling variants per sample
    htsbox pileup -q 0 -evcf $REFERENCE 02_alignments_pHSDs/${sample}.1.$0.paf.filt.bam 02_alignments_pHSDs/${sample}.2.$0.paf.filt.bam > 03_variants_pHSDs/${sample}.pair.$0.vcf

    # obtaining diploid calls
    dipcall-aux.js vcfpair -s ${sample} -a 03_variants_pHSDs/${sample}.pair.$0.vcf \
    | bcftools norm -Ov -m-any \
    | bcftools norm -Ov -d exact \
    | bcftools norm -Ov -m-any --fasta-ref $REFERENCE --check-ref w \
    | bcftools sort -m 4G \
    | bgzip -c > 03_variants_pHSDs/${sample}.dip.$0.vcf.gz
    tabix 03_variants_pHSDs/${sample}.dip.$0.vcf.gz
done'
```

Then, we identified fully callable samples per regions (samples where a locus of interest was fully covered by one contig per haplotype):
```bash
# callable regions per sample per locus
cat extended.tsv | xargs -n2 -P10 bash -c '
for sample in $(cut -f1 samples.tsv); do
    echo "Procesing sample "$sample" in locus "$0
    REFERENCE="00_reference/chm13.draft_v1.0.plusY.fasta"

    # obtaining callable regions per sample
    cat 02_alignments_pHSDs/${sample}.1.$0.paf.filt.var | grep ^R | cut -f2- > 02_alignments_pHSDs/${sample}.1.$0.paf.filt.bed
    cat 02_alignments_pHSDs/${sample}.2.$0.paf.filt.var | grep ^R | cut -f2- > 02_alignments_pHSDs/${sample}.2.$0.paf.filt.bed
    bedtools intersect -a 02_alignments_pHSDs/${sample}.1.$0.paf.filt.bed -b 02_alignments_pHSDs/${sample}.2.$0.paf.filt.bed \
    > 03_variants_pHSDs/${sample}.dip.$0.bed
done'

# merging callable regions
cat extended.bed | xargs -n4 -P10 bash -c '
for sample in $(cut -f1 samples.tsv); do

    bedtools intersect -f 1 -wao -a extended.bed -b 03_variants_pHSDs/${sample}.dip.$3.bed \
    | grep "$3" | awk -v chr=$0 -v start=$1 -v end=$2 -v locus=$3 -v var=$sample "{{if(\$5!=\".\"){{print chr\"\t\"start\"\t\"end\"\t\"locus\"\t\"var}} }}" 
      
done' | bedtools sort | bedtools groupby -g 1,2,3,4 -c 5 -o collapse > pHSDs_regions.fully_covered.bed
```

We merged VCF files per locus from fully callable samples only:
```bash
mkdir -p 04_final_VCFs_pHSDs

awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' pHSDs_regions.fully_covered.bed | xargs -n3 -P1 bash -c '
    echo "Processing locus "$1
    vcfs=$(for i in $(echo $2 | sed "s/,/ /g"); do printf "03_variants_pHSDs/${i}.dip.$1.vcf.gz "; done)

    bcftools merge -r $0 --missing-to-ref $vcfs \
    | bcftools norm -m+ \
    | bcftools norm -m- \
    | bcftools norm -f 00_reference/chm13.draft_v1.0.plusY.fasta \
    | bcftools view -f 'PASS,.' \
    | awk "{{if(\$0 ~ \"^#\" || \$5 != \"*\") {{print}} }}" \
    > 04_final_VCFs_pHSDs/$1.vcf
    bgzip -f 04_final_VCFs_pHSDs/$1.vcf; tabix -f 04_final_VCFs_pHSDs/$1.vcf.gz
'
```

Removing redundant samples and selecting biallelic SNPs only:
```bash
mkdir -p 05_merged_VCFs_pHSDs

awk '{print $1":"$2"-"$3"\t"$4}' extended.bed | xargs -n2 -P1 bash -c '
    bcftools view --force-samples -s ^HG00733_HGSVC,HG02818_HGSVC,HG03486_HGSVC,NA19240_HGSVC,NA24385_HGSVC -r $0 04_final_VCFs_pHSDs/$1.vcf.gz \
    | bcftools view --min-ac 1 \
    | bgzip -c > 05_merged_VCFs_pHSDs/$1.vcf.gz 
    tabix -f 05_merged_VCFs_pHSDs/$1.vcf.gz  

    cat 05_merged_VCFs_pHSDs/$1.vcf.gz \
    | bcftools norm -m+ \
    | bcftools norm -m- \
    | bcftools norm -f 00_reference/chm13.draft_v1.0.plusY.fasta \
    | bcftools view -f 'PASS,.' \
    | awk "{{if(\$0 ~ \"^#\" || \$5 != \"*\") {{print}} }}" \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | bgzip -c \
    > 05_merged_VCFs_pHSDs/$1.SNPs.biallelic.vcf.gz
    tabix -f 05_merged_VCFs_pHSDs/$1.SNPs.biallelic.vcf.gz
'
```

Pi per locus:
```r
# module load R/4.3.1
library(PopGenome)
library(data.table)
library(tidyverse)

regions <- fread("pHSDs_regions.fully_covered.bed", col.names = c("chr", "start", "end", "locus", "samples"), header = F, sep = "\t")
regions <- regions[,1:4]

results <- data.frame()
windows <- data.frame()

for(i in 1:nrow(regions)) {
    chr = regions[i,1][[1]]
    start = regions[i,2][[1]]
    end = regions[i,3][[1]]
    locus = regions[i,4][[1]]

    # overall pi
    GENOME.class <- readVCF(paste("05_merged_VCFs_pHSDs/",locus,".SNPs.biallelic.vcf.gz", sep = ""), 
                            tid=chr, from=start, to=end, numcols=10000000)
    GENOME.class <- F_ST.stats(GENOME.class)

    df1 <- data.frame(
        locus = locus, 
        pi =  as.vector(GENOME.class@Pi/GENOME.class@n.sites)) 
    results <- rbind(results, df1)
s
    # sliding window
    GENOME.class.slide <- sliding.window.transform(GENOME.class, width=20000, jump=1000, type=2, whole.data=TRUE)
    GENOME.class.slide <- F_ST.stats(GENOME.class.slide)
    
    df2 <- data.frame(
        locus = locus,
        regions = GENOME.class.slide@region.names, 
        pi = as.vector(GENOME.class.slide@Pi/GENOME.class.slide@n.sites)
    )    
    windows <- rbind(windows, df2)

    #GENOME.class.slide <- diversity.stats(GENOME.class.slide, pi = TRUE)
    #GENOME.class.slide@nuc.diversity.within
}

fwrite(results, "extended.hprc_hgsvc.overall_pi.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(windows, "extended.hprc_hgsvc.windows_pi.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
```

Tajima's D per locus:
```r
# module load R/4.3.1
library(PopGenome)
library(data.table)
library(tidyverse)

regions <- fread("pHSDs_regions.fully_covered.bed", col.names = c("chr", "start", "end", "locus", "samples"), header = F, sep = "\t")
regions <- regions[,1:4]

results <- data.frame()
windows <- data.frame()

for(i in 1:nrow(regions)) {
    chr = regions[i,1][[1]]
    start = regions[i,2][[1]]
    end = regions[i,3][[1]]
    locus = regions[i,4][[1]]

    # overall pi
    GENOME.class <- readVCF(paste("05_merged_VCFs_pHSDs/",locus,".SNPs.biallelic.vcf.gz", sep = ""), 
                            tid=chr, from=start, to=end, numcols=10000000)
    GENOME.class <- neutrality.stats(GENOME.class)

    df1 <- data.frame(
        locus = locus, 
        pi =  as.vector(GENOME.class@Tajima.D)) 
    results <- rbind(results, df1)

    # sliding window
    GENOME.class.slide <- sliding.window.transform(GENOME.class, width=20000, jump=1000, type=2, whole.data=TRUE)
    GENOME.class.slide <- neutrality.stats(GENOME.class.slide)
    
    df2 <- data.frame(
        locus = locus,
        regions = GENOME.class.slide@region.names, 
        pi = as.vector(GENOME.class.slide@Tajima.D)
    )    
    windows <- rbind(windows, df2)

    #GENOME.class.slide <- diversity.stats(GENOME.class.slide, pi = TRUE)
    #GENOME.class.slide@nuc.diversity.within
}

fwrite(results, "extended.hprc_hgsvc.overall_tajimad.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(windows, "extended.hprc_hgsvc.windows_tajimad.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
```

Genes in extended regions:
```bash
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.bed -b <(grep "GPR89A" extended.bed) | grep "protein_coding"
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.bed -b <(grep "GPR89B" extended.bed) | grep "protein_coding"
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.bed -b <(grep "FRMPD2A" extended.bed) | grep "protein_coding"
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.bed -b <(grep "FRMPD2B" extended.bed) | grep "protein_coding"
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.bed -b <(grep -w "CD8B" extended.bed) | grep "protein_coding"
bedtools intersect -a /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.bed -b <(grep "CD8B2" extended.bed) | grep "protein_coding"
```
