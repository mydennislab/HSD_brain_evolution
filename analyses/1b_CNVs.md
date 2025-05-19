# Copy-number estimation

## 1. 1KGP data

### 1.1 Metadata

**1KGP**

ENA projects associated with the 1KGP data are the following:
- Tab delimited file list for the 2504 panel or 2504 ENA Study at https://www.ebi.ac.uk/ena/data/view/PRJEB31736
- Tab delimited file list for the 698 related samples or 698 ENA Study at https://www.ebi.ac.uk/ena/browser/view/PRJEB36890

The FTP site contains complete metadata about both dataset (2,504 and 698 samples).

**HGDP**

High coverage sequencing of the Human Genome Diversity Project (HGDP) Cell Line Panel samples on the Illumina X10.

- 35x coverage of 929 samples described in Bergstrom et al. 2020
- Raw read alignments are available from the European Nucleotide Archive under study accession no. [PRJEB6463](https://www.ebi.ac.uk/ena/browser/view/PRJEB6463)
- Metadata available here http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/

**GIAB**

All GIAB metadata is here: https://github.com/genome-in-a-bottle

Data indexes are located here: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/copy-number/metadata

for pop in ACB ASW ESN GWD LWK MSL YRI; do
    grep -w ${pop} 1kgp_2504.all_pops.tsv > 1kgp_2504.AFR_${pop}.tsv
done
```

### 1.2 Data download and QC

```bash
cd /share/dennislab/databases/data/quickmer-t2t/1kg_high
```

We linked 1KG 30x coverage files:
```bash
cd /share/dennislab/databases/data/quickmer-t2t/1kg_high/data
for file in $(ls /share/datasets/dennislab/1kg/*cram); do ln -s $file; done
ls *cram | xargs -n1 -P1 bash -c 'samtools index -@ 20 $0'
```

We re-downloaded corrupted samples in data_fix folder:
```bash
cd data_fix
grep -Ff corrupted.txt urls.txt | sed 's/^/wget /g' > download.sh
bash download.sh
ls *cram | xargs -n1 -P1 bash -c 'samtools index -@ 20 $0'
```

We also included CHM13 Illumina reads as a control:
```bash
cat /share/dennislab/databases/data/CHM13/illumina/SRR3986881_1.fastq.gz  /share/dennislab/databases/data/CHM13/illumina/SRR3986881_2.fastq.gz > CHM13.fastq.gz
```

```bash
cd /share/dennislab/databases/data/quickmer-t2t/1kg_high/data
ls *cram | xargs -n1 -P5 bash -c 'md5sum $0' > ../data_1kgp_md5.txt

cd /share/dennislab/databases/data/quickmer-t2t/1kg_high/data_missing
ls *cram | xargs -n1 -P5 bash -c 'md5sum $0' > ../data_missing_1kgp_md5.txt
```

## 2. Population CN

### 2.1 QuicKmer2

#### 2.1.1 1KGP

```bash
cd /share/dennislab/databases/data/quickmer-t2t/1kg_high
```

> "The -c option specifies a bedfile of regions that are not expected to be copy-number variable across the analyzed samples. Sex chromosomes, unplaced chromosome sequences, known segmental duplications and known copy-number variants should be excluded."

Preparing reference:
```bash
mkdir -p reference
cd reference
ln -s /share/dennislab/projects/t2t/assembly/v1.0/chm13.draft_v1.0.fasta
cd ..
```

Preparing control regions:
```bash
mkdir -p control
cd control

# Downloading SegDups
wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/sedefSegDups/chm13.draft_v1.0_plus38Y.SDs.bed.bb
bigBedToBed chm13.draft_v1.0_plus38Y.SDs.bed.bb chm13.draft_v1.0_plus38Y.SDs.bed
cut -f1-3 chm13.draft_v1.0_plus38Y.SDs.bed | grep -v "^chrY" | bedtools sort -faidx chm13.draft_v1.0.fasta.fai > chm13.draft_v1.0.SDs.bed

# Generating control file
bedtools complement -i chm13.draft_v1.0.SDs.bed -g chm13.draft_v1.0.fasta.fai | grep -v "chrM\|chrX" > chm13.draft_v1.0.control.bed

cd ..
```

Running QuickMer2:
```bash
conda activate variants # python3 with pandas and argparse, and samtools + bedToBigBed
export PATH=/share/dennislab/databases/data/quickmer-t2t/QuicK-mer2/:$PATH

# Runnning data
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile quickmer2.smk \
--config reference=reference/chm13.draft_v1.0.fasta control=control/chm13.draft_v1.0.control.bed \
input=data outdir=results_w500 distance=1 window=500 \
-p -j 20

# Runnning data_missing
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile quickmer2.smk \
--config reference=reference/chm13.draft_v1.0.fasta control=control/chm13.draft_v1.0.control.bed \
input=data_missing outdir=results_w500 distance=1 window=500 \
-p -j 24

# Runnning data_capture
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile quickmer2.smk \
--config reference=reference/chm13.draft_v1.0.fasta control=control/chm13.draft_v1.0.control.bed \
input=data_capture outdir=results_w500 distance=1 window=500 \
-p -j 24
```

#### 2.1.2 archaics

```bash
conda activate variants # python3 with pandas and argparse, and samtools + bedToBigBed
export PATH=/share/dennislab/databases/data/quickmer-t2t/QuicK-mer2/:$PATH

# Runnning archaics
cd /share/dennislab/databases/data/quickmer-t2t/archaics

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake  \
--snakefile quickmer2.smk \
--config reference=reference/chm13.draft_v1.0.fasta control=reference/chm13.draft_v1.0.control.bed \
input=data outdir=results_w500 distance=1 window=500 \
-p -j 20
```

#### 2.1.3 HPRC

Running QuicKmer2 in HPRC in T2T-CHM13v1.0:
```bash
cd /share/dennislab/databases/data/quickmer-t2t/hprc_yr1

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake  \
--snakefile quickmer2.smk \
--config reference=reference/chm13.draft_v1.0.plusY.fasta control=control/chm13.draft_v1.0.control.bed \
input=data outdir=results_w500 distance=1 window=500 \
-p -j 20

/share/dennislab/users/dcsoto/Miniconda3/envs/smk/bin/snakemake \
--snakefile quickmer2.smk \
--config reference=reference/chm13.draft_v1.0.plusY.fasta control=control/chm13.draft_v1.0.control.bed \
input=data_fastq outdir=results_w500 distance=1 window=500 \
-p -j 20
```

Running QuicKmer2 in HPRC in T2T-CHM13v2.0:
```bash
cd /share/dennislab/databases/data/quickmer-t2t/hprc_yr1_v2.0

conda activate quickmer2 # python3 with pandas and argparse, and samtools + bedToBigBed
export PATH=/share/dennislab/databases/data/quickmer-t2t/QuicK-mer2/:$PATH

/share/dennislab/users/dcsoto/Miniconda3/envs/smk/bin/snakemake \
--snakefile quickmer2.smk \
--config reference=reference/chm13v2.0.fa control=reference/chm13v2.0.control.bed \
input=data outdir=results_w500 distance=1 window=500 \
-p -j 20
```

### 2.2 CN genotyping

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/copy-number
```

**SDs-98**

First, we calculated copy number of HSD genes in the 1KGP unrelated 2,504 samples.

Linking 1KGP 2,504 unrelated individuals:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/copy-number/CN_1KGP/
for smp in $(cat samples_1kgp.txt); do ln -s /share/dennislab/databases/data/quickmer-t2t/1kg_high/results_w500/CN/$smp.CN.bed . ; done
```

Copy-number genotyping of HSDs (SD-98) genes in 1KGP unrelated individuals:
```bash
conda activate python3

split -l 100 -d <(cat coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.genes.protein_coding.bed coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.genes.unprocessed_pseudogene.bed) part_

for num in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20; do
    sbatch --partition=production --job-name="cn_${num}" --nodes=1 --ntasks=16 --mem=1G --time=1-12:00:00 --wrap="/share/dennislab/users/dcsoto/Miniconda3/envs/python3/bin/python3 genotype_cn_parallel.py --path CN_1KGP --genes part_${num} --output part_${num}.cn.tsv -t 16"
done

for num in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20; do
    sbatch --partition=production --job-name="count_${num}" --nodes=1 --ntasks=1 --mem=1G --time=2-00:00:00 --wrap="/share/dennislab/users/dcsoto/Miniconda3/envs/python3/bin/python3 count_windows.py --input CN_1KGP/HG00096.CN.bed --genes part_${num} --output windows_${num}.tsv"
done
```

Merging loci:
```r
# module load R/4.0.1
files <- list.files(pattern = ".cn.tsv$")
list_df <- lapply(files, read.table, header=TRUE)
joined_df <- do.call("rbind", list_df)
write.table(joined_df, "quickmer_t2t.1kgp_2504.sd98.tsv", sep="\t", quote=FALSE, row.names = FALSE)

files <- list.files(pattern = "^windows_")
list_df <- lapply(files, read.table, header=TRUE)
joined_df <- do.call("rbind", list_df)
write.table(joined_df, "quickmer_t2t.1kgp_2504.sd98_windows.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

CN genotyping of archaics:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/03_CNVs/CN_archaics/
ln -s /share/dennislab/databases/data/quickmer-t2t/archaics/results_w500/CN/*.CN.bed .

cat coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.genes.protein_coding.bed coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.genes.unprocessed_pseudogene.bed > coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.genes.protein_coding_pseudogenes.bed

/share/dennislab/users/dcsoto/Miniconda3/envs/python3/bin/python3 genotype_cn_parallel.py --path CN_archaics --genes coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.genes.protein_coding_pseudogenes.bed --output archaics.cn.tsv -t 30
```

Obtaining SD-98 with CN=2 in >90% of 1KGP individuals:
```r
library(tidyverse)
library(data.table)

df <- fread("quickmer_t2t.1kgp_2504.sd98.tsv")
df2 <- df %>% pivot_longer(starts_with(c("NA","HG")), names_to="sample", values_to="copy_number")
df3 <- df2 %>% group_by(chrom, chromStart, chromEnd, geneName) %>% summarize(CN_2=sum(copy_number >= 1.5 & copy_number < 2.5))
df4 <- df3 %>% filter(CN_2>=2258)

fwrite(df4, file = "SD-98_CN_2.bed", sep = "\t", col.names=FALSE)
```

Annotating exons:
```bash
bedtools intersect -F 1 -wao -a <(cut -f1-4 quickmer_t2t.1kgp_2504.sd98.tsv) -b /share/dennislab/projects/hsd_genes/T2T/annotation/CHM13.combined.v4.exonsOnly.bed | awk '{if($5!="."){print}}' | sed 's/;/\t/g' | cut -f1-4,11 | awk '{if($4==$5){print}}' | cut -f1-4 | sort | uniq > quickmer_t2t.1kgp_2504.sd98.isec_exon.tsv
```

Annotating HSD:
```bash
bedtools intersect -F 1 -wao -a <(cut -f1-4 quickmer_t2t.1kgp_2504.sd98.tsv) -b /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/HSD/hsds.bed | awk '{if($5!="."){print}}' | cut -f1-3 | sort | uniq >quickmer_t2t.1kgp_2504.sd98.isec_hsds.tsv
```

**pHSDs**

We also calculated the copy number of pHSDs in the 1KGP unrelated dataset plus additional unrelated individuals from the 1KGP extended dataset that will be included in pop gen analyses. The combined dataset will be used for estimation of CN=2 individuals using KDE.

Linking intended capture regions:
```bash
ln -s /share/dennislab/projects/hsd_genes/T2T/regions/intended_capture.t2t.SDs.bed .
```

Copy number genotyping of pHSDs in 1KGP unrelated individuals:
```bash
conda activate python3
python3 genotype_cn_parallel.py --path CN_1KGP --genes intended_capture.t2t.SDs.bed --output quickmer_t2t.1kgp_2504.phsds.tsv -t 30
```

Linking unrelated capture individuals from 1KGP (including extended dataset):
```bash
cd CN_1KGP_extra/
for smp in $(cat samples_extra.txt); do ln -s /share/dennislab/databases/data/quickmer-t2t/1kg_high/results_w500/CN/$smp.CN.bed . ; done
cd ..
```

Copy number genotyping of pHSDs capture individuals from the extended dataset:
```bash
conda activate python3
python3 genotype_cn_parallel.py --path CN_1KGP_extra --genes intended_capture.t2t.SDs.bed --output quickmer_t2t.1kgp_extra.phsds.tsv -t 30
```

Merging samples:
```R
library(tidyverse)
df1 <- read.table("quickmer_t2t.1kgp_2504.phsds.tsv", header=TRUE)
df2 <- read.table("quickmer_t2t.1kgp_extra.phsds.tsv", header=TRUE) %>% select(!c(chrom, chromStart, chromEnd, geneName))
df <- cbind(df1,df2)

write.table(df, "quickmer_t2t.1kgp_4kde.phsds.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

We also genotyped the unique section of GYPA for other project:
```bash
conda activate python3
python3 genotype_cn_parallel.py --path CN_1KGP --genes GYPA.unique.bed --output GYPA.unique.1kgp_2504.tsv -t 120
```

## 3. SSC CN

We genotyped copy number of pHSDs and SD-98 in trios from the Simons Simplex Collection.

### 3.1 SSC QuicK-mer2

We obtained copy-number estimates from QuickMer2 using T2T-CHM13v1.1 from the Turner Lab. The files were downloaded in `/share/dennislab/databases/data/quickmer-t2t/ssc_v1.1/`. We observed that the chromosomes had a different naming convention in these files, so we had to rename them.

```bash
for file in /share/dennislab/databases/data/quickmer-t2t/ssc_v1.1/*bed.gz; do
    filename=$(basename -s .qm2.CN.1k.bed.gz $file)
    echo $filename
    gunzip -c /share/dennislab/databases/data/quickmer-t2t/ssc_v1.1/$filename.qm2.CN.1k.bed.gz | sed 's/^.*chromosome_/chr/g' | grep -v "mitochondrion" > CN/$filename.qm2.CN.1k.bed
done
```

### 3.2 Coordinate liftover

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/ssc/CNV_t2tv1.1
```

We had to lift the coordinates of SD-98 and pHSDs from v1.0 to v1.1 using the LiftOver chain provided in (GitHuB)[https://github.com/marbl/CHM13]. 

Lifting over pHSDs coordinates from v1.0 to v1.1:
```bash
conda activate variants

# we added the original SD98 region in T2T-CHM13-v1.0 coordinates to the name
liftOver \
<(cat /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed | cut -f1-4 | awk '{print $0";"$1":"$2"-"$3}') \
v1.0_to_v1.1_rdna.chain \
protein_coding_pseudogenes.t2t.v1.1.bed \
protein_coding_pseudogenes.t2t.v1.1.unmapped
```

### 3.3 Copy number genotyping

**SD-98 CN genotyping**

First, we genotyped SD-98, in batches:
```bash
split -l 25 -d protein_coding_pseudogenes.t2t.v1.1.bed part_

conda activate python3
for num in 0{1..9} {10..73}; do
    python3 genotype_cn_parallel.py --path /share/dennislab/databases/data/quickmer-t2t/ssc_v1.1/CN --genes part_${num} --output part_${num}.cn.tsv -t 30
done
```

Merging loci:
```r
# module load R/4.3.1
files <- list.files(pattern = ".cn.tsv$")
list_df <- lapply(files, read.table, header=TRUE)
joined_df <- do.call("rbind", list_df)
write.table(joined_df, "SD98.t2t-chm13v1.1.SSC.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

```bash
sed -i 's/_SSC_T2Tref//g' SD98.t2t-chm13v1.1.SSC.tsv
```

**pHSDs CN genotyping**

Then, we genotyped CN in pHSDs:
```bash
conda activate python3

python3 genotype_cn_parallel.py --path /share/dennislab/databases/data/quickmer-t2t/ssc_v1.1/CN --genes intended_capture.t2t.SDs.v1.1.bed --output intended_capture.t2t.SDs.v1.1.SSC_CN.tsv -t 120

sed -i 's/_SSC_T2Tref//g' intended_capture.t2t.SDs.v1.1.SSC_CN.tsv
```

We also genotyped a unique region in 1q21.1 to distinguish individuals with the 1q21.1 microdeletion/microduplication syndrome:
- hg19: chr1:146578257-147312088
- T2T-CHM13 v1.0 and 1.1: chr1:146890147-147702540

```bash
python3 genotype_cn_parallel.py --path /share/dennislab/databases/data/quickmer-t2t/ssc_v1.1/CN --genes 1q21_unique.bed --output 1q21_unique.SSC_CN.tsv -t 120
sed -i 's/_SSC_T2Tref//g' 1q21_unique.SSC_CN.tsv
```

We also genotyped other genes of interest:
```bash
python3 genotype_cn_parallel.py --path /share/dennislab/databases/data/quickmer-t2t/ssc_v1.1/CN --genes YTHDF2_YTHDC1.bed --output YTHDF2_YTHDC1.SSC_CN.tsv -t 120
```

## 4. UK BioBank CNVs

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/uk-biobank
```

### 4.1 Downloading data and metadata

> From UK BioBank CNV study: The file  `cnv.assoc.release_v2.20190513.tar.gz` contains two types of summary statistics files: one for rare variant burden tests, the other for common variant associations. Files are named as: `ukb24983_v2.[PHENO_ID].[GWAS_TYPE].PHENO1.glm.[REGRESSION_TYPE]` where `PHENO_ID` is the GBE phenotype ID (as desribed in in the index file above); `GWAS_TYPE` is cnv or cnv-burden, for CNV GWAS or gene-level burden tests; and `REGRESSION_TYPE` is linear for quantitative traits and logistic for binary outcomes.

```bash
wget https://biobankengine.stanford.edu/static/cnvdata/pheno_index_v2.20190513.txt
wget https://biobankengine.stanford.edu/static/cnvdata/cnv.ukb.reference_v2.20190626.txt
wget https://biobankengine.stanford.edu/static/cnvdata/cnv.assoc.release_v2.20190513.tar.gz
```

> Importantly, UK BioBank CNV dataset is in Hg19 coordinates.

### 4.2 Gene burden of rare CNVs

We first look for an excess of rare CNVs in a phenotype cohort in pHSD genes.

We pre-processed original gene burden test results:
```bash
for file in $(find gene_burden_test -name ukb24983_v2.*linear); do
    filename=$(basename $file)
    echo $filename;
    awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $file | tail -n +2 > gene_burden_linear/$filename 
done

for file in $(find gene_burden_test -name ukb24983_v2.*hybrid); do
    filename=$(basename $file)
    echo $filename;
    awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' $file | tail -n +2 > gene_burden_logistic/$filename 
done
```

We then searched for pHSDs gene names with a significant burden of rare CNVs:
```bash
for file in gene_burden_linear/*; do
    filename=$(basename $file)
    awk 'FNR==NR {a[$1]; next} FNR> 1 && $1 in a' phsds.txt <(cut -f4,10,13 $file) | \
    awk -v NAME=$filename '{if($3<=0.05){print NAME"\t"$0}}'
done > phsds.gene_burden_linear_p005.tsv

for file in gene_burden_logistic/*; do
    filename=$(basename $file)
    awk 'FNR==NR {a[$1]; next} FNR> 1 && $1 in a' phsds.txt <(cut -f4,13,14 $file) | \
    awk -v NAME=$filename '{if($3<=0.05){print NAME"\t"$0}}'
done > phsds.gene_burden_logistic_p005.tsv
```

### 4.3 Common CNV associations

The UK Biobank CNVs are in Hg19. Therefore, we had to liftover the SD-98 genes and intended capture regions coordinates from T2T-CHM13 to Hg19.

Downloading chains
```bash
# chain from v1.0 to hg38
wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/hg38Lastz/t2t-chm13-v1.0.hg38.over.chain.gz
gunzip t2t-chm13-v1.0.hg38.over.chain.gz

# chain from hg38 to hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
gunzip hg38ToHg19.over.chain.gz
```

Pre-processing original gwas files and selected significant associations (p-value < 0.05).
```bash
# linear
# ID OBS_CT BETA SE T_STAT P
for file in $(find gwas -name ukb24983_v2.*linear); do
    filename=$(basename $file)
    echo $filename;
    cut -f3,8- $file | tail -n +2 | awk 'BEGIN{OFS="\t"}{sub("-","\t",$1); print}' | sed 's/:/\t/g;s/_/\t/g' | awk '{if($9<=0.05){print}}' > gwas_linear/$filename 
done

# logistic
# ID OBS_CT OR SE Z_STAT P
for file in $(find gwas -name ukb24983_v2.*hybrid); do
    filename=$(basename $file)
    echo $filename;
    cut -f3,9- $file | tail -n +2 | awk 'BEGIN{OFS="\t"}{sub("-","\t",$1); print}'  | sed 's/:/\t/g;s/_/\t/g' | awk '{if($9<=0.05){print}}' > gwas_logistic/$filename 
done
```

Finding common CNV associations in pHSDs:
```bash
conda activate variants

# liftover
liftOver /share/dennislab/projects/hsd_genes/T2T/regions/intended_capture.t2t.fix.bed t2t-chm13-v1.0.hg38.over.chain intended_capture.t2t.fix.hg38.bed intended_capture.t2t.fix.bed_unmapped
liftOver intended_capture.t2t.fix.hg38.bed hg38ToHg19.over.chain intended_capture.t2t.fix.hg19.bed intended_capture.t2t.fix.hg38.bed_unmapped

sed -i 's/^chr//g' intended_capture.t2t.fix.hg19.bed

# Linear
for file in gwas_linear/*; do
    filename=$(basename $file)
    bedtools intersect -F 0.5 -wao -a $file -b intended_capture.t2t.fix.hg19.bed | awk '{if($10!="."){print}}' | awk -v NAME=$filename '{print NAME"\t"$0}'
done > intended_capture.gwas_linear.tsv

# Logistic
for file in gwas_logistic/*; do
    filename=$(basename $file)
    bedtools intersect -F 0.5 -wao -a $file -b intended_capture.t2t.fix.hg19.bed | awk '{if($10!="."){print}}' | awk -v NAME=$filename '{print NAME"\t"$0}'
done > intended_capture.gwas_logistic.tsv
```

Finding common CNV and brain measurement associations in SD98 genes:
```bash
# we added the original SD98 region in T2T-CHM13-v1.0 coordinates to the name
liftOver <(cat /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.SD-98.protein_coding.bed /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.SD-98.unprocessed_pseudogenes.bed | cut -f1-4 | awk '{print $0";"$1":"$2"-"$3}') t2t-chm13-v1.0.hg38.over.chain protein_coding_pseudogenes.hg38.bed protein_coding_pseudogenes.t2t.bed_unmapped

liftOver protein_coding_pseudogenes.hg38.bed hg38ToHg19.over.chain protein_coding_pseudogenes.hg19.bed protein_coding_pseudogenes.hg38.bed_unmapped
sed -i 's/^chr//g' protein_coding_pseudogenes.hg19.bed

# obtaining genes failed to liftover
cat <(grep -v "^#" protein_coding_pseudogenes.t2t.bed_unmapped | cut -f4) <(grep -v "^#" protein_coding_pseudogenes.hg38.bed_unmapped | cut -f4) | sort | uniq > protein_coding_pseudogenes.unmapped.txt

# Logistic
for file in gwas_logistic/*; do
    filename=$(basename $file)
    bedtools intersect -F 0.5 -wao -a $file -b protein_coding_pseudogenes.hg19.bed | awk '{if($10!="."){print}}' | awk -v NAME=$filename '{print NAME"\t"$0}'
done > sd98_genes.gwas_logistic.tsv

# Linear
for file in gwas_linear/*; do
    filename=$(basename $file)
    bedtools intersect -F 0.5 -wao -a $file -b protein_coding_pseudogenes.hg19.bed | awk '{if($10!="."){print}}' | awk -v NAME=$filename '{print NAME"\t"$0}'
done > sd98_genes.gwas_linear.tsv

# selecting brain measurements only
cat sd98_genes.gwas_linear.tsv | sed 's/ukb24983_v2.//g' | sed 's/.cnv.PHENO1.glm.linear//g' \
> sd98_genes.gwas_linear.brain_measurement.tsv

# separate duplications
awk '{if($5=="+"){print}}' sd98_genes.gwas_linear.brain_measurement.tsv \
| cut -f1-4,7,14 \
| awk 'BEGIN{OFS="\t"}{print $1,$2":"$3"-"$4";"$5,$6}' \
> sd98_genes.gwas_linear.brain_measurement.dup.tsv

# separate deletions
awk '{if($5=="-"){print}}' sd98_genes.gwas_linear.brain_measurement.tsv \
| cut -f1-4,7,14 \
| awk 'BEGIN{OFS="\t"}{print $1,$2":"$3"-"$4";"$5,$6}' \
> sd98_genes.gwas_linear.brain_measurement.del.tsv
```
