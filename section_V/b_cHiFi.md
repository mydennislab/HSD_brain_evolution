# cHiFi variants

We performed capture sequencing of ~200 individuals from the 1KGP project, HGDP, and PGP/GIAB.

## 0. Mean MAPQ across pHSDs

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/mean_mapq
```

```bash
files=$(ls /share/dennislab/projects/hsd_genes/T2T/varcalling/results/alignments/*bam | cut -f2 | tr '\n' ' ')
samtools merge finalBamFile.bam $files
samtools index -@ 10 finalBamFile.bam
```

Calculating median MAPQ per 1kbp window and generating browser track:
```bash
bedtools makewindows -b /share/dennislab/projects/hsd_genes/T2T/regions/intended_capture.t2t.fix.bed -w 1000 -s 1000 -i src > intended_capture.t2t.w1000.bed

cat intended_capture.t2t.w1000.bed | awk '{print $1":"$2+1"-"$3"\t"$4}' | xargs -n2 -P60 bash -c '
  meanMAPQ=$(samtools view finalBamFile.bam $0 | cut -f5 | awk "{{ sum+=\$0 }} END {{ if( NR>0 ) {{ print sum/NR }} else {{ print 0 }} }}"  )
  printf $0"\t"$1"\t"$meanMAPQ"\n" 
' | sed 's/:/\t/g;s/-/\t/g' | awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5}' > intended_capture.t2t.w1000.meanMAPQ.bed

bedtools sort -i intended_capture.t2t.w1000.meanMAPQ.bed > intended_capture.t2t.w1000.meanMAPQ.srt.bed

cut -f1-3,5 intended_capture.t2t.w1000.meanMAPQ.srt.bed | awk 'BEGIN{OFS="\t"}{
  if($4>=0 && $4<10){$5=0; $6="+"; $7=$2; $8=$3; $9="215,48,39"}
  else if($4>=10 && $4<20){$5=0; $6="+"; $7=$2; $8=$3; $9="252,141,89"}
  else if($4>=20 && $4<30){$5=0; $6="+"; $7=$2; $8=$3; $9="254,224,139"}
  else if($4>=30 && $4<40){$5=0; $6="+"; $7=$2; $8=$3; $9="217,239,139"}
  else if($4>=40 && $4<50){$5=0; $6="+"; $7=$2; $8=$3; $9="145,207,96"}
  else{$5=0; $6="+"; $7=$2; $8=$3; $9="26,152,80"}
  ;print $0}' > intended_capture.t2t.w1000.meanMAPQ.srt.color.bed
```

Getting mean mapq by locus:
```bash
bedtools groupby -i intended_capture.t2t.w1000.meanMAPQ.srt.color.bed -g 1 -c 4 -o median
bedtools groupby -i intended_capture.t2t.w1000.meanMAPQ.srt.bed -g 4 -c 5 -o median
```

## 1. Variant calling in T2T-CHM13

```bash
cd /share/dennislab/projects/hsd_genes/T2T/varcalling
```

We first called variants on each sample and then perform joint genotyping with GATK.

### 1.0 Pre-processing

Raw reads were previously pre-processed by demultiplexing runs and removing PCR duplicates.

### 1.1 Preparing merged runs

There were a total of 192 individuals successfully sequenced once and 8 individuals successfully sequenced twice. Runs of these individuals were written in `runs_to_merge.tsv`. Then, runs were merged.

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/T2T/merged_runs/
cat runs_to_merge.tsv | xargs -n 3 -P 1 bash -c 'samtools merge "$0".bam "$1" "$2"'
```

### 1.2 SNV calling per sample

Pre-processed BAM files from single or merged runs (for individuals sequenced more than once) were saved in `input_files.tsv`. Note: BAM files must be indexed.

Preparing reference:
```bash
mkdir -p reference
cd reference
ln -s /share/dennislab/projects/t2t/assembly/v1.0/chm13.draft_v1.0.fasta
cd ..
```

Running SNV/SV calling wiht different MAPQs:
```bash
conda activate chifi-sv
module load gatk/4.2.0.0

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake \
--snakefile var-calling-chifi.smk \
--config reference=reference/chm13.draft_v1.0.fasta files=input_files.tsv roi=../regions/intended_capture.t2t.fix.bed \
-p -j 40
```

We added two flags to GATK HaplotypeCaller to try different MAPQ thresholds based on [this](https://gatk.broadinstitute.org/hc/en-us/community/posts/6282475621403-Find-variants-with-HaplotypeCaller-on-MAPQ0-regions) link:
> * --mapping-quality-threshold
> * --mapping-quality-threshold-for-genotyping

## 2. cHiFi joint genotypes

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/
```

### 2.1 Joint genotyping

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/01_genotyping
```

> --standard-min-confidence-threshold-for-calling: Only variant sites with QUAL equal or greater than this threshold will be called. Note that when HaplotypeCaller is used in GVCF mode (using either -ERC GVCF or -ERC BP_RESOLUTION) the call threshold is automatically set to zero. Call confidence thresholding will then be performed in the subsequent GenotypeGVCFs command. Note that the default was changed from 10.0 to 30.0 in version 4.1.0.0.

```bash
module load gatk/4.2.0.0

for mapq in 1 2 5 10 20; do
  (
  echo $mapq
  variants=$(for file in $(ls /share/dennislab/projects/hsd_genes/T2T/varcalling/results/snvs/*_$mapq.gvcf); do echo "--variant "$file; done | tr '\n' ' ')
  gatk CombineGVCFs --pedigree ../pedigree.ped -R ../00_reference/chm13.draft_v1.0.plusY.fasta -O cohort_200.MAPQ_${mapq}.vcf.gz $variants
  ) &
done

for mapq in 1 2 5 10 20; do
  (
  for conf in 0 10 20 30; do
    gatk GenotypeGVCFs --pedigree ../pedigree.ped -R ../00_reference/chm13.draft_v1.0.plusY.fasta -V cohort_200.MAPQ_${mapq}.vcf.gz -O cohort_200.MAPQ_${mapq}.genotypes_${conf}.vcf.gz --standard-min-confidence-threshold-for-calling ${conf} 
  done
  ) &
done

for file in *genotypes*.vcf.gz; do printf $file"\t"; bcftools stats $file | grep "^SN" | grep "records" | cut -f4; done
```

### 2.2 Variant hard filtering

For hard filtering, we followed this instructions from GATK website: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering

> As with VQSR, hard-filter SNPs and indels separately. As of this writing, SelectVariants subsets SNP-only records, indel-only records or mixed-type, i.e. SNP and indel alternate alleles in the same record, separately. Therefore, when subsetting to SNP-only or indel-only records, mixed-type records are excluded. To avoid the loss of mixed-type variants, break up the multiallelic records into biallelic records before proceeding with the following subsetting.
> The * character was introduced in VCF 4.2 to indicate that the allele is missing due to an upstream deletion. More information on spanning deletion allele: https://gatk.broadinstitute.org/hc/en-us/articles/360035531912-Spanning-or-overlapping-deletions-allele-

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/02_filtering
```

Splitting multiallelic records and removing spanning deletions:
```bash
conda activate variants

for mapq in 1 2 5 10 20; do
  for conf in 0 10 20 30; do
    bcftools norm -m- ../01_genotyping/cohort_200.MAPQ_${mapq}.genotypes_${conf}.vcf.gz \
    | bcftools view -e'ALT="*"' \
    | bcftools norm -f ../00_reference/chm13.draft_v1.0.plusY.fasta \
    > cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.vcf
  done
done
```

Separating SNPs and indels:
```bash
module load gatk/4.2.0.0

for mapq in 1 2 5 10 20; do
  for conf in 0 10 20 30; do
    gatk SelectVariants \
        -V cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.vcf \
        -select-type SNP \
        -O cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.vcf 

    gatk SelectVariants \
        -V cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.vcf \
        -select-type INDEL \
        -O cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.indels.vcf
  done
done
```

Obtaining calls at different GQ and DP hard filters for SNPs and indels:
```bash
for mapq in 1 2 5 10 20; do
  for conf in 0 10 20 30; do
    for gq in 0 20 50 70; do
      for dp in 0 4 8 12 16; do
        echo "$mapq $conf $gq $dp"
      done
    done
  done
done | xargs -P 20 -n 4 bash -c '
  mapq=$1
  conf=$2
  gq=$3
  dp=$4

  echo Processing: mapq $mapq, conf $conf, gq $gq, dp $dp

  vcftools --minGQ $gq --minDP $dp --vcf cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.vcf --recode --recode-INFO-all --stdout \
  | bcftools +fill-tags - -- -t AN,AC \
  | bcftools view -i "AC>0" \
  | bgzip -c > cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.filt_${gq}_${dp}.vcf.gz

  tabix -f cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.filt_${gq}_${dp}.vcf.gz

' _

for mapq in 1 2 5 10 20; do
  for conf in 0 10 20 30; do
    for gq in 0 20 50 70; do
      for dp in 0 4 8 12 16; do
        echo "$mapq $conf $gq $dp"
      done
    done
  done
done | xargs -P 20 -n 4 bash -c '
  mapq=$1
  conf=$2
  gq=$3
  dp=$4

  echo Processing: mapq $mapq, conf $conf, gq $gq, dp $dp

  vcftools --minGQ $gq --minDP $dp --vcf cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.indels.vcf --recode --recode-INFO-all --stdout \
  | bcftools +fill-tags - -- -t AN,AC \
  | bcftools view -i "AC>0" \
  | bgzip -c > cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.indels.filt_${gq}_${dp}.vcf.gz

  tabix -f cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.indels.filt_${gq}_${dp}.vcf.gz

' _
```

### 2.3 Variant benchmarking

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/03_benchmarking
```

#### 2.3.1 Total variants

```bash
for gq in 0 20 50 70; do
  for dp in 0 4 8 12 16; do
    for mapq in 1 2 5 10 20; do
      for conf in 0 10 20 30; do
        gunzip -c ../02_filtering/cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.filt_${gq}_${dp}.vcf.gz \
        | bcftools norm -m+ \
        | bcftools view --max-alleles 2 \
        | bcftools view --exclude-types indels \
        | bcftools stats | grep "^SN" | grep "SNPs" | cut -f4 \
        | sed "s/^/$gq $dp $mapq $conf /g"
      done
    done
  done
done > total_snps.tsv
```

#### 2.3.2 Inbreeding coefficient

> InbreedingCoeff: Excess heterozygotes defined by an inbreeding coefficient < -0.3

```bash
for gq in 0 20 50 70; do
  for dp in 0 4 8 12 16; do
    for mapq in 1 2 5 10 20; do
      for conf in 0 10 20 30; do
        gunzip -c ../02_filtering/cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.filt_${gq}_${dp}.vcf.gz \
        | bcftools norm -m+ \
        | bcftools view --max-alleles 2 \
        | bcftools view --exclude-types indels \
        | bcftools query -f '%INFO/InbreedingCoeff\n' \
        | sed "s/^/$gq $dp $mapq $conf /g"
      done
    done
  done
done > inbreeding_coefficients.tsv
```

#### 2.3.3 Mendelian concordance

```bash
# preparing T2T-reference
conda activate variants
rtg format -o chm13.draft_v1.0.noalt_SDF /share/dennislab/projects/hsd_genes/T2T/varcalling/reference/chm13.draft_v1.0.fasta

for gq in 0 20 50 70; do
  for dp in 0 4 8 12 16; do
    for mapq in 1 2 5 10 20; do
      for conf in 0 10 20 30; do
        gunzip -c ../02_filtering/cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.filt_${gq}_${dp}.vcf.gz \
        | bcftools view -S trio_samples.txt -Ov \
        | bcftools norm -m+ \
        | bcftools view --max-alleles 2 \
        | bcftools view --exclude-types indels \
        | bcftools view -i 'F_MISSING=0' \
        | rtg mendelian -i - -t chm13.draft_v1.0.noalt_SDF --pedigree=../pedigree.ped \
        | grep "^Concordance" | sed "s/^/$gq $dp $mapq $conf /g"
      done
    done
  done
done | tr -s " " | cut -d" " -f1,2,3,4,6,12 | sed 's/://g' | sed "s/(//g" | sed "s/)//g" > mendelian_concordance.tsv
```

#### 2.3.4 HG002

We compared variants from cHiFi joint genotyping to HPRC variants for individual HG002.

```bash
conda activate vcfeval

# all regions
awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' /share/dennislab/projects/hsd_genes/T2T/regions/intended_capture.t2t.fix.bed | sed 's/;/_/g' | xargs -n2 -P22 bash -c '
echo "Locus: "$1
  for gq in 0 20 50 70; do
    for dp in 0 4 8 12 16; do
      for mapq in 1 2 5 10 20; do
        for conf in 0 10 20 30; do
          rtg vcfeval \
          --sample="HG002_HPRC","NA24385" \
          --region=$0 \
          --baseline /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc/07_merged_VCFs_pHSDs/$1.bisnps.all.vcf.gz \
          --calls ../02_filtering/cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.filt_${gq}_${dp}.biallelic.vcf.gz \
          -t chm13.draft_v1.0.noalt_SDF \
          -o vcfeval_all/$1_MAPQ${mapq}_CONF${conf}_DP${dp}_GQ${gq}
        done
      done
    done
  done   
'

# SD only
awk '{print $1":"$2"-"$3"\t"$4"\t"$5}' /share/dennislab/projects/hsd_genes/T2T/regions/intended_capture.t2t.fix.SDs.bed | sed 's/;/_/g' | xargs -n2 -P22 bash -c '
echo "Locus: "$1
  for gq in 0 20 50 70; do
    for dp in 0 4 8 12 16; do
      for mapq in 1 2 5 10 20; do
        for conf in 0 10 20 30; do
          rtg vcfeval \
          --sample="HG002_HPRC","NA24385" \
          --region=$0 \
          --baseline /share/dennislab/users/dcsoto/ms_hsd/phsds/hprc_hgsvc/07_merged_VCFs_pHSDs/$1.bisnps.all.vcf.gz \
          --calls ../02_filtering/cohort_200.MAPQ_${mapq}.genotypes_${conf}.norm.snps.filt_${gq}_${dp}.biallelic.vcf.gz \
          -t t2t-chm13v1.0.noalt_SDF \
          -o vcfeval_SD/$1_MAPQ${mapq}_CONF${conf}_DP${dp}_GQ${gq}
        done
      done
    done
  done   
'

find vcfeval_all -name summary.txt | xargs -i bash -c "echo {}; cat {}" | grep -v '^-\|^Th' | grep '^vcfeval\|None\|available' | awk '{c="\n"} NR%2 {c=" "} {printf("%s%s", $0, c)}' | tr -s ' ' | sed 's/\/summary.txt//g' > summary_vcfeval_all.tsv
cat text_to_replace.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0/$1/g" summary_vcfeval_all.tsv'

find vcfeval_SD -name summary.txt | xargs -i bash -c "echo {}; cat {}" | grep -v '^-\|^Th' | grep '^vcfeval\|None\|available' | awk '{c="\n"} NR%2 {c=" "} {printf("%s%s", $0, c)}' | tr -s ' ' | sed 's/\/summary.txt//g' > summary_vcfeval_SD.tsv
cat text_to_replace.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0/$1/g" summary_vcfeval_SD.tsv'
```

## 3. cHiFi unrelated dataset

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/04_unrelated
```

### 3.1 Selecting unrelated samples

Merging SNPs and indels from selected thresholds for unrelated samples only:
```bash
bcftools concat -a ../02_filtering/cohort_200.MAPQ_20.genotypes_30.norm.snps.filt_20_8.vcf.gz ../02_filtering/cohort_200.MAPQ_20.genotypes_30.norm.indels.filt_20_8.vcf.gz \
| bcftools view -S unrelated.txt -c1 \
| bcftools norm -m+ \
| bcftools norm -m- \
| bcftools norm -f ../00_reference/chm13.draft_v1.0.plusY.fasta \
| bcftools view -f 'PASS,.' \
| awk "{{if(\$0 ~ \"^#\" || \$5 != \"*\") {{print}} }}" \
| bgzip -c \
> cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz
tabix -f cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz
```

### 3.2 Variant counts

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/phsds/chifi/05_variant_counts
```

Counting number of variants in SDs
```bash
# overall variants in SD regions
awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed \
| xargs -n2 -P1 bash -c '
    printf $1"\t"

    bcftools view -r $0 ../04_unrelated/cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz \
    | bcftools stats \
    | grep "^SN" \
    | grep "number of samples\|number of records\|number of SNPs\|number of indels" \
    | cut -f4 \
    | tr "\n" "\t"
    printf "\n"
' > chifi.SD.total_counts.tsv
```

Variant density per individual:
```bash
# heterozygous variants in SD per sample
awk '{print $1":"$2"-"$3"\t"$4}' ../../coords/pHSDs_regions.SD.bed \
| xargs -n2 -P1 bash -c '

for smp in `bcftools query -l ../04_unrelated/cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz`
do
    printf "$1\t"
    printf "$smp\t"

    bcftools view -s ${smp} -c1 -r $0 ../04_unrelated/cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | grep -v "1|1\|1/1\|0/0" \
    | bcftools stats | grep "^SN" | grep "number of SNPs" | cut -d" " -f3- | cut -f2
done' > chifi.SD.het_counts.tsv

# heterozygous variants in matching SD per sample
awk '{print $1":"$2"-"$3"\t"$5"\t"$6}' ../../coords/pHSDs_txs.SD.fix.bed \
| xargs -n3 -P1 bash -c '

for smp in `bcftools query -l ../04_unrelated/cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz`
do
    printf "$1\t"
    printf "$smp\t"

    bcftools view -s ${smp} -c1 -r $0 ../04_unrelated/cohort_200.MAPQ_20.genotypes_30.norm.all.filt_20_8.unrelated.vcf.gz \
    | bcftools view --exclude-types indels \
    | bcftools norm -m+ \
    | bcftools view --max-alleles 2 \
    | grep -v "1|1\|1/1\|0/0" \
    | bcftools stats | grep "^SN" | grep "number of SNPs" | cut -d" " -f3- | cut -f2
done' > chifi.SD_matching.het_counts.tsv
```
