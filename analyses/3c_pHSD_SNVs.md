# Combined cHiFi + HPRC/HGSVC SNVs

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/combined
```

## 1. Combining datasets 

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/combined/01_combined
```

- Final cHiFi dataset all: 
  - `/share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/04_unrelated/cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz`
- Final HPRC/HGSVC dataset per pHSD:
  - `/share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc/05_merged_VCFs_pHSDs/*.vcf.gz`

We coupled HPRC/HGSVC and cHiFi SNVs for downstream analyses:
```bash
# separating chifi VCF
awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed | xargs -n2 -P1 bash -c '
   bcftools view -r $0 /share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/04_unrelated/cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz \
   | bgzip -c > $1.chifi.vcf.gz
   tabix -f $1.chifi.vcf.gz
'

# merging HPRC/HGSVC with cHiFi
awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed | grep -v 'FAM72C_SRGAP2D' | xargs -n2 -P1 bash -c '
   bcftools merge $1.chifi.vcf.gz ../../hprc_hgsvc/05_merged_VCFs_pHSDs/$1.vcf.gz \
   | bcftools norm -m+ \
   | bcftools norm -m- \
   | bcftools norm -f ../../chifi/00_reference/chm13.draft_v1.0.plusY.fasta \
   | bcftools view -f 'PASS,.' \
   | awk "{{if(\$0 ~ \"^#\" || \$5 != \"*\") {{print}} }}" \
   | bgzip -c \
   | bcftools +fill-tags - -- -t AN,AC,AF \
   | bgzip -c \
   > $1.combined.vcf.gz
   tabix $1.combined.vcf.gz
'

# copying FAM72C/SRGAP2D
cp FAM72C_SRGAP2D.chifi.vcf.gz FAM72C_SRGAP2D.combined.vcf.gz
tabix -f FAM72C_SRGAP2D.combined.vcf.gz
```

## 2. Variant counts

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/combined/02_variant_counts
```

Counting number of variants in SDs
```bash
# overall variants in SD regions
awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed \
| xargs -n2 -P1 bash -c '
    printf $1"\t"

    bcftools view -r $0 ../01_combined/$1.combined.vcf.gz \
    | bcftools stats \
    | grep "^SN" \
    | grep "number of samples\|number of records\|number of SNPs\|number of indels" \
    | cut -f4 \
    | tr "\n" "\t"
    printf "\n"
' > combined.SD.total_counts.tsv
```

Variant density per individual:
```bash
# heterozygous variants in SD per sample
awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed \
| xargs -n2 -P1 bash -c '

for smp in `bcftools query -l ../01_combined/$1.combined.vcf.gz`
do
    printf "$1\t"
    printf "$smp\t"

    bcftools view -s ${smp} -c1 -r $0 ../01_combined/$1.combined.vcf.gz \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | grep -v "1|1\|1/1\|0/0" \
    | bcftools stats | grep "^SN" | grep "number of SNPs" | cut -d" " -f3- | cut -f2
done' > combined.SD.het_counts.tsv
```

## 3. Variant effect predictor

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/combined/03_vep
```

For this analysis, we did not consider the following:
- HYDIN2: not canonical transcript selected

We also consider canonical transcripts and matching SD regions. These coordinates were manually curated and saved in `/share/dennislab/users/dcsoto/ms_hsd/phsds/coords/pHSDs_txs.SD.fix.bed`

Running VEP:
```bash
conda activate vep

awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed | xargs -n2 -P1 bash -c '
  vep --force_overwrite \
  -i ../01_combined/$1.combined.vcf.gz \
  -o $1.vep.vcf --vcf \
  --gff /share/dennislab/users/dcsoto/ms_hsd/phsds/coords/CHM13.combined.v4.locus.srt.biotype.pHSDs.sort.gff3.gz \
  --fasta /share/dennislab/projects/hsd_genes/T2T/varcalling/reference/chm13.draft_v1.0.fasta
'

for file in *vcf; do bgzip -f $file; tabix -f $file.gz; done
```

Obtaining variant tables:
```bash
awk '{print $1":"$2"-"$3"\t"$4"\t"$5"\t"$6}' ../../coords/pHSDs_txs.SD.fix.bed | xargs -n4 -P1 bash -c '
  coords=$0; tx=$1; paralog=$2; locus=$3;

  bcftools view -r ${coords} ${locus}.vep.vcf.gz \
  | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\n" \
  | sed "s/|/\t/g" | cut -f1-6,8,9,13 | awk -v paralog=${paralog} "{print paralog\"\t\"\$0}" \
  | grep ${tx} \
  | grep "HIGH\|MODERATE\|LOW"
' > high_moderate_low_vep.txt
```

## 4. pN/pS

We identified signatures of natural selection using the ratio pN/pS (relative abundance of non-synonymous and synonymous polymorphisms), which measures the direct effect of natural selection removing slightly deleterious non-synonymous variants.

#### 3.3.1 Synonymous and nonsynonymous sites

Using Ka/Ks calculation from the from seqinr package, Alfie previously obtained the number of synonymous and non-synonymous sites per pHSD:
- L0: frequency of non-synonymous sites.
- L2: frequency of 2-fold synonymous sites. 
- L4: frequency of 4-fold synonymous sites.

> Note: A position of a codon is said to be a n-fold degenerate site if only n of four possible nucleotides (A, C, G, T) at this position specify the same amino acid.

Let Ks be the number of (synonymous) substitutions per synonymous site and Ka the number of (nonsynonymous) substitutions per nonsynonymous site. Following the convention to count one third of a 2-fold degenerate sites as synonymous and 2/3 as nonsynonymous, we obtain:
- Ks = (L2A2 + L4K4) / (L2/3 + L4)
- Ka = (L2B2 + LoKo) / (2L2/3 + Lo) 

So, the number of synonymous sites and nonsynonymous sites are:
- synonymous sites: (L2/3 + L4)
- nonsynonymous sites: (2*L2/3 + L0) 

```r
library(data.table)
library(tidyverse)

dat_kaks <- readRDS("hsd_kaks.rds")
genes <- names(dat_kaks)
results <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(results) <- c('gene', 'paralog', 'syn_sites', 'nonsyn_sites')

for(i in 1:length(genes)) {
  gene = genes[i]
  
  paralogs = c("human_A", "human_B")
  if(gene == "SRGAP2"){ paralogs = c("human_A", "human_B", "human_C")}
  
  for(h in 1:length(paralogs)){
    paralog <- paralogs[h]

    L0 <- as.matrix(dat_kaks[[gene]]$l0)['chimp', paralog]
    L2 <- as.matrix(dat_kaks[[gene]]$l2)['chimp', paralog]
    L4 <- as.matrix(dat_kaks[[gene]]$l4)['chimp', paralog]
    
    syn_sites <- (L2/3 + L4)
    nonsyn_sites <- (2*L2/3 + L0) 
    
    new_results <- data.frame('gene' = gene, 'paralog' = paralog, 'syn_sites' = syn_sites, 'nonsyn_sites' = nonsyn_sites)
    results <- rbind(results, new_results)
  }
}

fwrite(results, "pHSDs.sites.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = TRUE)
```

#### 3.3.2 Synonymous and nonsynonymous substitutions

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/combined/04_pnps
```

We annotated variant consequence for biallelic SNPs for MKT test.

Obtaining biallelic SNPs:
```bash
awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed | xargs -n2 -P1 bash -c '
  bcftools view -r $0 ../01_combined/$1.combined.vcf.gz \
  | bcftools norm -m- \
  | bcftools view --exclude-types indels \
  | bcftools norm -m+ \
  | bcftools view --max-alleles 2 \
  | bgzip -c \
  > $1.combined.SNPs.biallelic.vcf.gz
  tabix -f $1.combined.SNPs.biallelic.vcf.gz
'
```

Running VEP on biallelic SNPs:
```bash
conda activate vep

awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed | xargs -n2 -P1 bash -c '
  vep --force_overwrite \
  -i $1.combined.SNPs.biallelic.vcf.gz \
  -o $1.SNPs.biallelic.vep.vcf --vcf \
  --gff /share/dennislab/users/dcsoto/ms_hsd/phsds/coords/CHM13.combined.v4.locus.srt.biotype.pHSDs.sort.gff3.gz \
  --fasta /share/dennislab/projects/hsd_genes/T2T/varcalling/reference/chm13.draft_v1.0.fasta
'

for file in *vcf; do bgzip -f $file; tabix -f $file.gz; done
```

Obtaining variant tables:
```bash
awk '{print $1":"$2"-"$3"\t"$4"\t"$5"\t"$6}' ../../coords/pHSDs_txs.SD.fix.bed | xargs -n4 -P1 bash -c '
  coords=$0; tx=$1; paralog=$2; locus=$3;

  bcftools view -r ${coords} ${locus}.SNPs.biallelic.vep.vcf.gz \
  | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\n" \
  | sed "s/|/\t/g" | cut -f1-6,8,9,13 | awk -v paralog=${paralog} "{print paralog\"\t\"\$0}" \
  | grep ${tx} \
  | grep "HIGH\|MODERATE\|LOW"
' > high_moderate_low_vep_biallelic.txt
```
