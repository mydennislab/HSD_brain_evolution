# Analysis of WGS SNVs

## 1. SRS accessibility

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/coordinates

# combined mask stats
awk '{sum+=$3-$2}END{print sum}' combined_mask.bed # 2671196397
bedtools sort -i combined_mask.bed | bedtools merge | awk '{sum+=$3-$2}END{print sum}' # 2671196397 -> no overlapping features
grep -v "^chrM\|chrX\|chrY" combined_mask.bed | awk '{sum+=$3-$2}END{print sum}' # 2518828751

# accessible regions
bedtools intersect -a chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed -b combined_mask.bed > chm13.draft_v1.0_plus38Y.SDs.merged.auto.accessible.bed
bedtools intersect -a chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed -b combined_mask.bed > chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.accessible.bed
bedtools intersect -a chm13.draft_v1.0_plus38Y.non-SDs_no-censat.merged.auto.bed -b combined_mask.bed > chm13.draft_v1.0_plus38Y.non-SDs_no-censat.merged.auto.accessible.bed

awk '{sum+=$3-$2}END{print sum}' chm13.draft_v1.0_plus38Y.SDs.merged.auto.accessible.bed # 68,258,070
awk '{sum+=$3-$2}END{print sum}' chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.accessible.bed # 10,623,540
awk '{sum+=$3-$2}END{print sum}' chm13.draft_v1.0_plus38Y.non-SDs_no-censat.merged.auto.accessible.bed # 2,439,661,620
```

SD-98 genes in accessible regions:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes

bedtools intersect -wao -a CHM13.combined.v4.genes.SD-98.nr.protein_coding.body.bed -b ../coordinates/combined_mask.bed | awk '{print $4"\t"$3-$2"\t"$8}' | bedtools groupby -g 1,2 -c 3 -o sum > CHM13.combined.v4.genes.SD-98.nr.protein_coding.body.accessible.tsv

bedtools intersect -wao -a CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.body.bed -b ../coordinates/combined_mask.bed | awk '{print $4"\t"$3-$2"\t"$8}' | bedtools groupby -g 1,2 -c 3 -o sum > CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.body.accessible.tsv
```

## 2. Depletion in T2T-CHM13

```bash
unique_auto_t2t="/share/dennislab/users/dcsoto/ms_hsd/wgs/depletion/regions/T2T/chm13.draft_v1.0.auto.no_SDs.no_cenAnnot.bed"
segdups_auto_t2t="/share/dennislab/users/dcsoto/ms_hsd/wgs/depletion/regions/T2T/chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed"
segdups98_auto_t2t="/share/dennislab/users/dcsoto/ms_hsd/wgs/depletion/regions/T2T/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed"
genome_auto_t2t="/share/dennislab/users/dcsoto/ms_hsd/wgs/depletion/regions/T2T/chm13.draft_v1.0.auto.bed"
censat_t2t="/share/dennislab/users/dcsoto/ms_hsd/wgs/depletion/regions/T2T/t2t_cenAnnotation.v2.021921.no_ct.bed"
```

### 2.1 Random regions

Generating random regions (within the same chromosomes) of the same size as SDs and SDs ≥98%:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/depletion

# SDs
mkdir -p random/random_segdup_t2t
for num in `seq 1 10000`; do bedtools shuffle -chrom -i ${segdups_auto_t2t} -g <(cut -f1,3 ${genome_auto_t2t}) -excl ${censat_t2t} -noOverlapping -maxTries 10000 -f 0.1 > random/random_segdup_t2t/random_segdup_${num}.bed; done

# SDs-98
mkdir -p random/random_segdup98_t2t
for num in `seq 1 10000`; do bedtools shuffle -chrom -i ${segdups98_auto_t2t} -g <(cut -f1,3 ${genome_auto_t2t}) -excl ${censat_t2t} -noOverlapping -maxTries 10000 -f 0.1 > random/random_segdup98_t2t/random_segdup98_${num}.bed; done
```

### 2.2 SNV depletion

First, we filtered only SNPs in 1KGP T2T-CHM13 VCF:
```bash
cd /share/dennislab/projects/hsd_genes/T2T/1kgp/

bcftools view --exclude-types indels 1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz > 1kgp.all.recalibrated.snp.pass.nogenos.labeled.sort.vcf
bgzip 1kgp.all.recalibrated.snp.pass.nogenos.labeled.sort.vcf
tabix 1kgp.all.recalibrated.snp.pass.nogenos.labeled.sort.vcf.gz

bcftools view --min-ac=1 --max-alleles 2 --exclude-types indels 1kgp.all.recalibrated.snp_indel.pass.nogenos.labeled.sort.vcf.gz > 1kgp.all.recalibrated.snp_biallelic.pass.nogenos.labeled.sort.vcf
bgzip 1kgp.all.recalibrated.snp_biallelic.pass.nogenos.labeled.sort.vcf
tabix 1kgp.all.recalibrated.snp_biallelic.pass.nogenos.labeled.sort.vcf.gz
```

Then we calculated the observed number of SNVs per region:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/depletion/1KGP_T2T

bedtools intersect -a /share/dennislab/projects/hsd_genes/T2T/1kgp/1kgp.all.recalibrated.snp.pass.nogenos.labeled.sort.vcf.gz -b ${segdups_auto_t2t} | cut -f1,2 | sort | uniq | wc -l # 4659742
bedtools intersect -a /share/dennislab/projects/hsd_genes/T2T/1kgp/1kgp.all.recalibrated.snp.pass.nogenos.labeled.sort.vcf.gz -b ${segdups98_auto_t2t} | cut -f1,2 | sort | uniq | wc -l # 1153189
```

Finally, we calculated the number of SNVs in 2000 random regions:
```bash
seq 1 2000 | xargs -n1 -P30 bash -c 'bedtools intersect -a /share/dennislab/projects/hsd_genes/T2T/1kgp/1kgp.all.recalibrated.snp.pass.nogenos.labeled.sort.vcf.gz -b ../random/random_segdup_t2t/random_segdup_$0.bed | cut -f1,2 | sort | uniq | wc -l' > random_segdup_t2t_1kgp_snps.txt
seq 1 2000 | xargs -n1 -P30 bash -c 'bedtools intersect -a /share/dennislab/projects/hsd_genes/T2T/1kgp/1kgp.all.recalibrated.snp.pass.nogenos.labeled.sort.vcf.gz -b ../random/random_segdup98_t2t/random_segdup98_$0.bed | cut -f1,2 | sort | uniq | wc -l' > random_segdup98_t2t_1kgp_snps.txt

awk 'BEGIN{OFMT="%f"}{sum+=$1}END{print sum/NR}' random_segdup98_t2t_1kgp_snps.txt # 3666857.3
```

* SD-98: 1153189/97797568*1000 = 11.79 variants/kbp
* SD: 4659742/179849378*1000 = 25.90 variants/kbp
* Random: 3666857/97797568*1000 = 37.49 variants/kbp

## 3. Benchmarking

### 3.1 Obtaining variant files

- PBCCS samples: HG002, HG003, HG004, HG005, HG006, HG007, HG01109, HG01243, HG01442, HG02055, HG02080, HG02109, HG02145, HG02723, HG03098, HG03492
- PBCCS samples in 1KGP: HG01109, HG01243, HG02055, HG02080, HG02145, HG02723, HG03098, HG03492

PBCCS calls were donwloaded from Globus to: `/share/dennislab/users/dcsoto/ms_hsd/wgs/variants/PBCCS`

We extracted biallelic SNPs from PacBio calls:
```bash
for smp in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do  
    bcftools view --max-alleles 2 --min-ac=1 --exclude-types indels ${smp}_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz | bgzip -c > ${smp}_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.biallelic_snps.vcf.gz
    tabix ${smp}_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.biallelic_snps.vcf.gz
done
```

Original Illumina joint calls from 1KGP are located in `/share/dennislab/projects/hsd_genes/T2T/1kgp`. 

We extracted common individuals from the main dataset:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/variants/Illumina
conda activate variants

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do
    (
    for smp in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do
        mkdir -p ${smp}   
        
        echo "Processing:" $smp $chr
        
        bcftools view --max-alleles 2 --min-ac=1 --exclude-types indels -s ${smp} /share/dennislab/projects/hsd_genes/T2T/1kgp/1kgp.${chr}.recalibrated.snp_indel.pass.vcf.gz > ${smp}/${smp}.${chr}.vcf     
    done
    ) &
done

for smp in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do
    vcf-concat $smp/${smp}.*.vcf > $smp.all.vcf
    bcftools sort $smp.all.vcf > $smp.all.srt.vcf
done
```

> Note: bcftools view updates INFO/AC and INFO/AN only.

### 3.2 PacBio and Illumina WGS concordance

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/concordance
```

First, we formatted T2T reference for the analysis:
```bash
conda activate variants # rtgtools
rtg format -o T2T_reference /share/dennislab/projects/hsd_genes/T2T/varcalling/reference/chm13.draft_v1.0.fasta
```
> Warning message: This looks like a genome reference, but is not recognized. You should enable sex-aware processing by manually installing an appropriate reference.txt containing chromosome metadata.  See the user manual for more information.

Using Illumina as "calls" and PBCCS as "benchmark", we generated an input file called `input.tsv` with 4 columns including: (1) filename, (2) variant calls, (3) bed benchmark, (4) vcf benchmark. (Note: VCF must be unzipped.)

For this analysis, we used as "all" all autosome coordinates (`/share/dennislab/users/dcsoto/ms_hsd/regions/T2T/chm13.draft_v1.1.auto.bed`)

Also, we used region files classified as non-SD, SD, and SD98. We generated symbolic links in regions folder:
```bash
cd regions
ln -s /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0.auto.no_SDs.no_cenAnnot.bed Non-SDs.bed
ln -s /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed SDs.bed
ln -s /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed SDs-98.bed

ln -s /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0.auto.no_SDs.no_cenAnnot.combined_mask.bed Non-SDs_mask.bed
ln -s /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs.merged.auto.combined_mask.bed SDs_mask.bed
ln -s /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.combined_mask.bed SDs-98_mask.bed
cd ..
```

Running evaluation:
```bash
conda activate vcfeval
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --snakefile snake-eval-squash.smk --config files=input.tsv reference=T2T_reference -p -j 10
```

> Note: This pipeline only compares SNPs; indels are removed from the VCF files prior to evaluation.

Generating evaluation summaries:
```bash
cd eval/

find * -name summary.txt | xargs -i bash -c "echo {}; cat {}" | grep -v '^-\|^Th' | grep '^eval\|None\|available' | awk '{c="\n"} NR%2 {c=" "} {printf("%s%s", $0, c)}' | tr -s ' ' > ../eval_summary.txt

cd ../eval_squash

find * -name summary.txt | xargs -i bash -c "echo {}; cat {}" | grep -v '^-\|^Th' | grep '^eval\|None\|available' | awk '{c="\n"} NR%2 {c=" "} {printf("%s%s", $0, c)}' | tr -s ' ' > ../eval_squash_summary.txt
```

### 3.3 PacBio and Illumina WGS SNP density

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/density/
```

We merged, LRS and SRS files for 8 common samples:
```bash
# LRS
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/variants/PBCCS/

for sample in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do
  bgzip ${sample}_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf
  tabix ${sample}_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz
dones

ls *vcf.gz | tr '\n' ' '
bcftools merge HG01109_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz HG01243_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz HG02055_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz HG02080_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz HG02145_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz HG02723_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz HG03098_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz HG03492_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf.gz > All_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf

awk '/^#/ || $7=="PASS" || $7=="."' All_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.vcf > All_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.pass.vcf

bcftools view --max-alleles 2 --min-ac=1 --exclude-types indels All_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.pass.vcf > All_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.pass.snps.vcf

# SRS
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/variants/Illumina

for sample in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do
  bgzip $sample.all.srt.vcf
  tabix $sample.all.srt.vcf.gz
done

ls *vcf.gz | tr '\n' ' '
bcftools merge HG01109.all.srt.vcf.gz HG01243.all.srt.vcf.gz HG02055.all.srt.vcf.gz HG02080.all.srt.vcf.gz HG02145.all.srt.vcf.gz HG02723.all.srt.vcf.gz HG03098.all.srt.vcf.gz HG03492.all.srt.vcf.gz > 1KGP.all.srt.vcf
```

Obtaining non-overlapping windows for each region of interest:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/density/windows

bedtools makewindows -w 1000 -s 1000 -b /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0.auto.no_SDs.no_cenAnnot.bed | awk '{if(($3-$2)==1000){print}}' > Non-SDs_windows_1000.bed
bedtools makewindows -w 1000 -s 1000 -b /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed | awk '{if(($3-$2)==1000){print}}' > SDs_windows_1000.bed
bedtools makewindows -w 1000 -s 1000 -b /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed | awk '{if(($3-$2)==1000){print}}' > SDs-98_windows_1000.bed

bedtools makewindows -w 1000 -s 1000 -b /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0.auto.no_SDs.no_cenAnnot.combined_mask.bed | awk '{if(($3-$2)==1000){print}}' > Non-SDs_windows_1000.combined_mask.bed
bedtools makewindows -w 1000 -s 1000 -b /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs.merged.auto.combined_mask.bed | awk '{if(($3-$2)==1000){print}}' > SDs_windows_1000.combined_mask.bed
bedtools makewindows -w 1000 -s 1000 -b /share/dennislab/users/dcsoto/ms_hsd/00_regions/T2T/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.combined_mask.bed | awk '{if(($3-$2)==1000){print}}' > SDs-98_windows_1000.combined_mask.bed
```

Obtaining SNP density per 1 kbp windows:
```bash
LRS_VCF="/share/dennislab/users/dcsoto/ms_hsd/wgs/variants/PBCCS/All_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.pass.snps.vcf"
SRS_VCF="/share/dennislab/users/dcsoto/ms_hsd/wgs/variants/Illumina/1KGP.all.srt.vcf"

bedtools coverage -a windows/Non-SDs_windows_1000.bed -b ${LRS_VCF} > Non-SDs_windows_1000.LRS_cov.bed &
bedtools coverage -a windows/Non-SDs_windows_1000.bed -b ${SRS_VCF} > Non-SDs_windows_1000.SRS_cov.bed &

bedtools coverage -a windows/SDs_windows_1000.bed -b ${LRS_VCF} > SDs_windows_1000.LRS_cov.bed &
bedtools coverage -a windows/SDs_windows_1000.bed -b ${SRS_VCF} > SDs_windows_1000.SRS_cov.bed &

bedtools coverage -a windows/SDs-98_windows_1000.bed -b ${LRS_VCF} > SDs-98_windows_1000.LRS_cov.bed &
bedtools coverage -a windows/SDs-98_windows_1000.bed -b ${SRS_VCF} > SDs-98_windows_1000.SRS_cov.bed &

bedtools coverage -a windows/Non-SDs_windows_1000.combined_mask.bed -b ${LRS_VCF} > Non-SDs_windows_1000.combined_mask.LRS_cov.bed &
bedtools coverage -a windows/Non-SDs_windows_1000.combined_mask.bed -b ${SRS_VCF} > Non-SDs_windows_1000.combined_mask.SRS_cov.bed &

bedtools coverage -a windows/SDs_windows_1000.combined_mask.bed -b ${LRS_VCF} > SDs_windows_1000.combined_mask.LRS_cov.bed &
bedtools coverage -a windows/SDs_windows_1000.combined_mask.bed -b ${SRS_VCF} > SDs_windows_1000.combined_mask.SRS_cov.bed &

bedtools coverage -a windows/SDs-98_windows_1000.combined_mask.bed -b ${LRS_VCF} > SDs-98_windows_1000.combined_mask.LRS_cov.bed &
bedtools coverage -a windows/SDs-98_windows_1000.combined_mask.bed -b ${SRS_VCF} > SDs-98_windows_1000.combined_mask.SRS_cov.bed &
```

## 4. Tajima's D

### 4.1 Tajima's D across SD98 genes

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/tajimasd/tajimasd_25kb
```

We first calculated outlier Tajima's D regions.

```bash
# Get Aarthi's D calculations
for pop in AFR EUR ASN SAS AMR; do
  for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22; do
    tail -n +2 /share/dennislab/projects/hsd_genes/coverage_stats/TajimaD_1KG/All_25kbTajimaD_Files/${pop}_1kg_${chr}_all_samples_25kb.Tajima.D
  done > ${pop}_1kg_all_samples_25kb.Tajima.D
done 

# Convert to bed preserving number of SNPs
awk '{print $1"\t"$2"\t"$2+25000"\t"$3"\t"$4}' AFR_1kg_all_samples_25kb.Tajima.D > AFR_1kg_all_samples_25kb.Tajima.D.bed
awk '{print $1"\t"$2"\t"$2+25000"\t"$3"\t"$4}' EUR_1kg_all_samples_25kb.Tajima.D > EUR_1kg_all_samples_25kb.Tajima.D.bed
awk '{print $1"\t"$2"\t"$2+25000"\t"$3"\t"$4}' ASN_1kg_all_samples_25kb.Tajima.D > ASN_1kg_all_samples_25kb.Tajima.D.bed
awk '{print $1"\t"$2"\t"$2+25000"\t"$3"\t"$4}' SAS_1kg_all_samples_25kb.Tajima.D > SAS_1kg_all_samples_25kb.Tajima.D.bed
awk '{print $1"\t"$2"\t"$2+25000"\t"$3"\t"$4}' AMR_1kg_all_samples_25kb.Tajima.D > AMR_1kg_all_samples_25kb.Tajima.D.bed

# Annotate SD98 regions and accessible regions
for pop in AFR EUR ASN SAS AMR; do
  bedtools intersect -wao -a ${pop}_1kg_all_samples_25kb.Tajima.D.bed -b /share/dennislab/users/dcsoto/ms_hsd/sd-98/coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed | cut -f1-5,9 | bedtools groupby -g 1,2,3,4,5 -c 6 -o sum | bedtools intersect -wao -a - -b /share/dennislab/users/dcsoto/ms_hsd/sd-98/coordinates/combined_mask.bed | cut -f1-6,10 | bedtools groupby -g 1,2,3,4,5,6 -c 7 -o sum > ${pop}_1kg_all_samples_25kb.Tajima.D.tsv
done

# Annotate SD98 genes
for pop in AFR EUR ASN SAS AMR; do
  bedtools intersect -wao -a ${pop}_1kg_all_samples_25kb.Tajima.D.tsv -b sd98_genes.txt | cut -f1-7,11 | bedtools groupby -g 1-7 -c 8 -o distinct > ${pop}_1kg_all_samples_25kb.Tajima.D.sd98_genes.tsv
done
```

Then, we selected only accessible regions (>=50%) with at least 5 SNPs to get the mean Tajima's D per SD-98 gene.

We uploaded these windows after pre-proccessing in R into `tajd_accessible_min5` folder.

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/tajimasd/tajimasd_sd98

cat <(cut -f1-4 /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed) <(cut -f1-4 /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes/CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed) > sd98_genes.bed

/share/dennislab/users/dcsoto/Miniconda3/envs/python3/bin/python3 genotype_cn_parallel.py \
--path tajd_accessible_min5 --genes sd98_genes.bed --output sd98_genes.tajimasd.tsv -t 16
```

### 4.2 Concordance in Tajima's D outliers

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/wgs/tajimasd
```

We then assessed the concordance between PacBio and Illumina calls on the 8 shared individuals for regions flagged as outliers in the Tajima's D analysis.

We generated bed files to calculate concordance per row:
```bash
split -l1 --numeric-suffixes=1 -a3 outliers.bed # stored in outliers folder
```

Running evaluation:
```bash
mkdir -p eval_default eval_squash eval_decompose

conda activate vcfeval

for bedfile in $(ls outliers); do
   for smp in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do
        rtg vcfeval \
        --baseline ../variants/PBCCS/"$smp"_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.biallelic_snps.vcf.gz \
        --bed-regions outliers/$bedfile \
        --calls ../variants/Illumina/$smp.all.srt.vcf.gz \
        -t ../concordance/T2T_reference \
        -o eval_default/${smp}_${bedfile}
   done
done
for bedfile in $(ls outliers); do
   for smp in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do
        rtg vcfeval \
        --baseline ../variants/PBCCS/"$smp"_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.biallelic_snps.vcf.gz \
        --bed-regions outliers/$bedfile \
        --calls ../variants/Illumina/$smp.all.srt.vcf.gz \
        -t ../concordance/T2T_reference \
        -o eval_squash/${smp}_${bedfile} \
        --squash-ploidy
   done
done
for bedfile in $(ls outliers); do
   for smp in HG01109 HG01243 HG02055 HG02080 HG02145 HG02723 HG03098 HG03492; do
        rtg vcfeval \
        --baseline ../variants/PBCCS/"$smp"_PacBio_HiFi_PEPPER_MARGIN_DEEPVARIANT_2_CHM13.biallelic_snps.vcf.gz \
        --bed-regions outliers/$bedfile \
        --calls ../variants/Illumina/$smp.all.srt.vcf.gz \
        -t ../concordance/T2T_reference \
        -o eval_decompose/${smp}_${bedfile} \
        --decompose
   done
done
```

> - By default, only considers records where FILTER is "." or "PASS"

Generating evaluation summaries:
```bash
find . -name summary.txt | xargs -i bash -c " cat {} | awk -v var={} '{{print var\" \"\$0}}' " | grep 'None\|available' | sed 's/0 total baseline variants, no summary statistics available/NA/g' | tr -s ' ' | tr ' ' '\t' | sed 's/\//\t/g;s/_/\t/g' | cut -f2-4,6- > eval_summary.txt
```

```r
library(tidyverse)
library(data.table)

df <- fread("eval_summary.txt", header = FALSE, fill = TRUE, 
              col.names = c("type", "sample", "region", "threshold",
                            "baseline_tp", "tp", "fp", "fn", "precision",
                            "sensitivity", "f_measure"))

df2 <- data %>% filter(type == "default") %>%
  select(region, sample, tp, fp, fn, precision, sensitivity) %>%
  group_by(region) %>% 
  summarize(sum(tp, na.rm=TRUE), sum(fp, na.rm=TRUE), 
            sum(fn, na.rm=TRUE), mean(precision, na.rm=TRUE), 
            mean(sensitivity, na.rm=TRUE))

fwrite(df2, "eval_summary_processed.txt")
```
