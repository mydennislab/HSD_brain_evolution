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

## 4. CD8B balancing selection analysis

### 4.1 Haplotype networks

Colors:
* AFR <- "#EDA93A"
* EUR <- "#3F53A4"
* EAS <- "#417F3C"
* SAS <- "#73297C"
* AMR <- "#E13126"

Mapping PanTro6 to reference:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc/99_CD8B/pantro6

conda activate variants

minimap2 -t 32 -a --eqx --cs -x asm5 --secondary=no -s 25000 -K 8G /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc/00_reference/chm13.draft_v1.0.plusY.fasta panTro6.fa | samtools view -@ 4 -Sb | samtools sort -@ 4 > pantro6.CD8B.bam
samtools index pantro6.CD8B.bam
```

```bash
# CD8B network 1 
# chr2:86824967-86844966
mkdir network1
ls ../02_alignments_pHSDs/*CD8B.paf.filt.bam | xargs -n1 -P20 bash -c '
    sample=$(basename -s .CD8B.paf.filt.bam $0)
    echo $sample
    samtools consensus -m simple -r chr2:86824967-86844966 -f fasta $0 | sed "s/>chr2/>"$sample"/g" > network1/${sample}_CD8B_network1.fa
'

samtools consensus -m simple -r chr2:86824967-86844966 -f fasta pantro6/pantro6.CD8B.bam | sed "s/>chr2/>pantro6/g" > network1/pantro6_CD8B_network1.fa
cat network1/*_CD8B_network1.fa > network1/all.CD8B_network1.fa

# CD8B network 2
# chr2:86827380-86833378
mkdir network2
ls ../02_alignments_pHSDs/*CD8B.paf.filt.bam | xargs -n1 -P20 bash -c '
    sample=$(basename -s .CD8B.paf.filt.bam $0)
    echo $sample
    samtools consensus -m simple -r chr2:86827380-86833378 -f fasta $0 | sed "s/>chr2/>"$sample"/g" > network2/${sample}_CD8B_network2.fa
'

samtools consensus -m simple -r chr2:86827380-86833378 -f fasta pantro6/pantro6.CD8B.bam | sed "s/>chr2/>pantro6/g" > network2/pantro6_CD8B_network2.fa
cat network2/*_CD8B_network2.fa > network2/all.CD8B_network2.fa

# CD8B network 3
# chr2:86865880-86871878
mkdir network3
ls ../02_alignments_pHSDs/*CD8B.paf.filt.bam | xargs -n1 -P20 bash -c '
    sample=$(basename -s .CD8B.paf.filt.bam $0)
    echo $sample
    samtools consensus -m simple -r chr2:86865880-86871878 -f fasta $0 | sed "s/>chr2/>"$sample"/g" > network3/${sample}_CD8B_network3.fa
'

samtools consensus -m simple -r chr2:86865880-86871878 -f fasta pantro6/pantro6.CD8B.bam | sed "s/>chr2/>pantro6/g" > network3/pantro6_CD8B_network3.fa
cat network3/*_CD8B_network3.fa > network3/all.CD8B_network3.fa

# CD8B network 4
# chr2:86845128-86851128
mkdir network4
ls ../02_alignments_pHSDs/*CD8B.paf.filt.bam | xargs -n1 -P20 bash -c '
    sample=$(basename -s .CD8B.paf.filt.bam $0)
    echo $sample
    samtools consensus -m simple -r chr2:86845128-86851128 -f fasta $0 | sed "s/>chr2/>"$sample"/g" > network4/${sample}_CD8B_network4.fa
'

samtools consensus -m simple -r chr2:86845128-86851128 -f fasta pantro6/pantro6.CD8B.bam | sed "s/>chr2/>pantro6/g" > network4/pantro6_CD8B_network4.fa
cat network4/*_CD8B_network4.fa > network4/all.CD8B_network4.fa
```

Merging all sequences in CD8B positive Tajima's D window:

```bash
cat *_CD8B_window.fa > all.CD8B_window.fa
```

CD8B and CD8B2 coordinates:
- CD8B: chr2:86840067-86863799
- CD8B2: chr2:106948802-106972551
- CD8B positive Taj D window: chr2:86825000-86850000

- CD8B: 	
CD8B chr2 86822879 86847878

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc/99_CD8B

module load samtools/1.15.1

ls ../02_alignments_pHSDs/*CD8B.paf.filt.bam | xargs -n1 -P20 bash -c '
    sample=$(basename -s .CD8B.paf.filt.bam $0)
    echo $sample
    samtools consensus -m simple -r chr2:86840067-86863799 -f fasta $0 | sed "s/>chr2/>"$sample"/g" > ${sample}_CD8B.fa
'

ls ../02_alignments_pHSDs/*CD8B2.paf.filt.bam | xargs -n1 -P20 bash -c '
    sample=$(basename -s .CD8B2.paf.filt.bam $0)
    echo $sample
    samtools consensus -m simple -r chr2:106948802-106972551 -f fasta $0 | sed "s/>chr2/>"$sample"/g" > ${sample}_CD8B2.fa
'


ls ../02_alignments_pHSDs/*CD8B.paf.filt.bam | xargs -n1 -P20 bash -c '
    sample=$(basename -s .CD8B.paf.filt.bam $0)
    echo $sample
    samtools consensus -m simple -r chr2:86825000-86850000 -f fasta $0 | sed "s/>chr2/>"$sample"/g" > ${sample}_CD8B_window.fa
'
```

We generated a multi-fasta file from individuals containing a contig fully covering CD8B and CD8B2:

```bash
grep -w "CD8B" ../../hprc_hgsvc/pHSDs_regions.fully_covered.SD.tsv > CD8B_fully_covered.tsv
grep -w "CD8B2" ../../hprc_hgsvc/pHSDs_regions.fully_covered.SD.tsv > CD8B2_fully_covered.tsv

awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' CD8B_fully_covered.tsv | xargs -n3 -P1 bash -c '
    echo "Locus: "$1

    fastas=$(for i in $(echo $2 | sed "s/,/ /g"); do printf $i".1_${1}.fa "$i".2_${1}.fa "; done)
    cat $fastas > all.$1.fa
'

awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' CD8B2_fully_covered.tsv | xargs -n3 -P1 bash -c '
    echo "Locus: "$1

    fastas=$(for i in $(echo $2 | sed "s/,/ /g"); do printf $i".1_${1}.fa "$i".2_${1}.fa "; done)
    cat $fastas > all.$1.fa
'
```

```bash
samtools consensus -m simple -r chr2:86825000-86850000 -f fasta pantro6.CD8B.bam | sed "s/>chr2/>pantro6/g" > ../pantro6_CD8B_window.fa
```

### 4.2 CD8B Tajima's D 

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc_extended
```

Separating populations:
```bash
conda activate variants

for pop in AFR EUR EAS SAS AMR; do
    cat 05_merged_VCFs_pHSDs/CD8B.SNPs.biallelic.vcf.gz \
    | bcftools view -S hprc_${pop}.tsv --force-samples \
    | bcftools +fill-tags - -- -t AN,AC,AF \
    | bcftools view -c1 \
    | bcftools view -i 'F_MISSING<0.1' \
    | bgzip -c \
    > 05_merged_VCFs_pHSDs/CD8B.SNPs.biallelic.${pop}.vcf.gz 
    tabix -f 05_merged_VCFs_pHSDs/CD8B.SNPs.biallelic.${pop}.vcf.gz 
done

for pop in AFR EUR EAS SAS AMR; do
    cat 05_merged_VCFs_pHSDs/CD8B2.SNPs.biallelic.vcf.gz \
    | bcftools view -S hprc_${pop}.tsv --force-samples \
    | bcftools +fill-tags - -- -t AN,AC,AF \
    | bcftools view -c1 \
    | bcftools view -i 'F_MISSING<0.1' \
    | bgzip -c \
    > 05_merged_VCFs_pHSDs/CD8B2.SNPs.biallelic.${pop}.vcf.gz 
    tabix -f 05_merged_VCFs_pHSDs/CD8B2.SNPs.biallelic.${pop}.vcf.gz 
done

# mid frequency variants CD8B
bcftools view -i 'AF>0.4 && AF<0.6' 05_merged_VCFs_pHSDs/CD8B.SNPs.biallelic.AFR.vcf.gz -Oz -o CD8B.SNPs.biallelic.AFR.maf_0.4.vcf.gz
bcftools view -i 'AF>0.3 && AF<0.7' 05_merged_VCFs_pHSDs/CD8B.SNPs.biallelic.AFR.vcf.gz -Oz -o CD8B.SNPs.biallelic.AFR.maf_0.3.vcf.gz 

bcftools view -i 'AF>0.4 && AF<0.6' 05_merged_VCFs_pHSDs/CD8B.SNPs.biallelic.AMR.vcf.gz -Oz -o CD8B.SNPs.biallelic.AMR.maf_0.4.vcf.gz
bcftools view -i 'AF>0.3 && AF<0.7' 05_merged_VCFs_pHSDs/CD8B.SNPs.biallelic.AMR.vcf.gz -Oz -o CD8B.SNPs.biallelic.AMR.maf_0.3.vcf.gz 

# mid frequency variants CD8B2
bcftools view -i 'AF>0.3 && AF<0.7' 05_merged_VCFs_pHSDs/CD8B2.SNPs.biallelic.AFR.vcf.gz -Oz -o CD8B2.SNPs.biallelic.AFR.maf_0.3.vcf.gz 
bcftools view -i 'AF>0.3 && AF<0.7' 05_merged_VCFs_pHSDs/CD8B2.SNPs.biallelic.AMR.vcf.gz -Oz -o CD8B2.SNPs.biallelic.AMR.maf_0.3.vcf.gz 
```

Tajima's D calculations per population and overall:
```r
# module load R/4.3.1
library(PopGenome)
library(data.table)
library(tidyverse)

regions <- fread("pHSDs_regions.fully_covered.bed", col.names = c("chr", "start", "end", "locus", "samples"), header = F, sep = "\t")
regions <- regions[,1:4]
regions <- regions[regions$locus %in% c("CD8B","CD8B2"),]

# Per population

results <- data.frame(
  locus = character(),
  d = numeric(),
  population = character(),
  stringsAsFactors = FALSE
)
windows <- data.frame()

for(pop in c("AFR", "EUR", "SAS", "EAS", "AMR")){
    for(i in 1:nrow(regions)) {
        chr = regions[i,1][[1]]
        start = regions[i,2][[1]]
        end = regions[i,3][[1]]
        locus = regions[i,4][[1]]

        # overall pi
        GENOME.class <- readVCF(paste("05_merged_VCFs_pHSDs/",locus,".SNPs.biallelic.",pop,".vcf.gz", sep = ""), 
                                tid=chr, from=start, to=end, numcols=10000000)
        GENOME.class <- neutrality.stats(GENOME.class)

        df1 <- data.frame(
            locus = locus, 
            d =  as.vector(GENOME.class@Tajima.D),
            population = pop) 
        print(df1)
        results <- rbind(results, df1)

        # sliding window
        GENOME.class.slide <- sliding.window.transform(GENOME.class, width=6000, jump=500, type=2, whole.data=TRUE)
        GENOME.class.slide <- neutrality.stats(GENOME.class.slide)
        
        df2 <- data.frame(
            locus = locus,
            regions = GENOME.class.slide@region.names, 
            d = as.vector(GENOME.class.slide@Tajima.D),
            population = pop
        )    
        windows <- rbind(windows, df2)
    }    
}

fwrite(results, "extended.hprc_hgsvc.overall_tajimad.CD8B_per_pop.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(windows, "extended.hprc_hgsvc.windows_tajimad.CD8B_per_pop.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)

# All populations
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
    GENOME.class.slide <- sliding.window.transform(GENOME.class, width=6000, jump=500, type=2, whole.data=TRUE)
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

fwrite(results, "extended.hprc_hgsvc.overall_tajimad.CD8B_all_pops.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(windows, "extended.hprc_hgsvc.windows_tajimad.CD8B_all_pops.tsv", quote = FALSE,  sep = "\t", row.names = FALSE, col.names = FALSE)
```
