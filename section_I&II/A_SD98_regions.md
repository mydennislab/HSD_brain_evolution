# SD98 regions definition

## 1. Region coordinates

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/coordinates
```

We obtained three regions in T2T-CHM13v1.0: SDs, SDs with a sequence identity >=98%, and non SDs (excluding centromeres).

Obtaining T2T-CHM13v1.0 SegDups:
```bash
conda activate variants
wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/sedefSegDups/chm13.draft_v1.0_plus38Y.SDs.bed.bb
bigBedToBed chm13.draft_v1.0_plus38Y.SDs.bed.bb chm13.draft_v1.0_plus38Y.SDs.bed
cut -f1-3 chm13.draft_v1.0_plus38Y.SDs.bed | bedtools sort | bedtools merge > chm13.draft_v1.0_plus38Y.SDs.merged.bed
```

Selecting SDs with sequence identity over 98%:
```bash
awk '{if($24>=0.98){print}}' chm13.draft_v1.0_plus38Y.SDs.bed | bedtools sort | bedtools merge > chm13.draft_v1.0_plus38Y.SDs-98.merged.bed
```

Obtaining non-SD regions (without centromeres):
```bash
awk '{print $1"\t"0"\t"$2}' /share/dennislab/projects/t2t/assembly/v1.0/chm13.draft_v1.0.plusY.fasta.fai > chm13.draft_v1.0_plus38Y.genome.bed
bedtools subtract -a chm13.draft_v1.0_plus38Y.genome.bed -b chm13.draft_v1.0_plus38Y.SDs.merged.bed \
| bedtools sort | bedtools merge > chm13.draft_v1.0_plus38Y.non-SDs.merged.bed

# substracting CenSat file but keeping ct (pericentromeric)
awk '{if($4 !~ "^ct_"){print}}' t2t_cenAnnotation.v2.021921.bed > t2t_cenAnnotation.v2.021921.no_pericentromeric.bed
bedtools subtract -a chm13.draft_v1.0_plus38Y.non-SDs.merged.bed -b t2t_cenAnnotation.v2.021921.no_pericentromeric.bed \
| bedtools sort | bedtools merge > chm13.draft_v1.0_plus38Y.non-SDs_no-censat.merged.bed
```

Calculating sizes of each region:
```bash
grep -v "chrX\|chrY\|chrM" chm13.draft_v1.0_plus38Y.SDs.merged.bed > chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed
grep -v "chrX\|chrY\|chrM" chm13.draft_v1.0_plus38Y.SDs.merged.bed | awk '{sum+=$3-$2}END{print sum}' # 179849378

grep -v "chrX\|chrY\|chrM" chm13.draft_v1.0_plus38Y.SDs-98.merged.bed > chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed
grep -v "chrX\|chrY\|chrM" chm13.draft_v1.0_plus38Y.SDs-98.merged.bed | awk '{sum+=$3-$2}END{print sum}' # 97797568

grep -v "chrX\|chrY\|chrM" chm13.draft_v1.0_plus38Y.non-SDs_no-censat.merged.bed > chm13.draft_v1.0_plus38Y.non-SDs_no-censat.merged.auto.bed
grep -v "chrX\|chrY\|chrM" chm13.draft_v1.0_plus38Y.non-SDs_no-censat.merged.bed | awk '{sum+=$3-$2}END{print sum}' # 2550809796
```

SD-98 per chromosome:
```bash
awk '{print $1"\t"$3-$2}' chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed | bedtools groupby -g 1 -c 2 -o sum
```

## 2. Genes in SD-98

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/genes
```

We obtained genes with at least one exon fully encompassed in autosomal SD-98 region.

> Note: We noticed that the GFF3 file is missing "gene" annotations for somes transcripts. For example, SRGAP2A doesn't have a "gene" annotation associated. We realized this discrepancy was associated with the missanotation of some transcripts with a generic name "MSTRG". Therefore, we replaced MSTRG genes with their true gene names.

### 2.1 Removing redundant genes

We noticed that in some cases two different gene IDs pointed out to the same gene. We first look for these redundant genes by focusing on transcripts from protein coding genes and pseudogenes that were fully contained (90%) within another gene with a different gene ID.

Finding redundant gene names:
```bash 
conda activate variants

cat /share/dennislab/projects/t2t/annotation_v1.0/cat_v4/CHM13.combined.v4.gff3 \
| awk '{if($3=="transcript"){print}}' | ./parseGFF3tx.py | gff2bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10"\t0\t"$6}' \
| sed 's/gene_name=//g;s/source_gene=//g;s/biotype=//g;s/gene_id=//g' \
> CHM13.combined.v4.txs.bed

cat MSTRG.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0;/$1;/g" CHM13.combined.v4.txs.bed'

# obtaining tx coordinates of genes that have and exon fully covered in SD-98 regions
bedtools intersect -F 1 -s -wao -a CHM13.combined.v4.txs.bed -b CHM13.combined.v4.exons.SD-98.bed \
| awk '{if($4==$10){print}}' \
| cut -f1-6 | sort | uniq \
| grep "protein_coding\|unprocessed_pseudogene" \
> CHM13.combined.v4.txs.sd98.bed

# intersecting txs with itself
# selecting txs with different names but almost fully overlapping each other 
bedtools intersect -wao -s -f 0.9 -a CHM13.combined.v4.txs.sd98.bed -b CHM13.combined.v4.txs.sd98.bed \
| awk '{if($4!=$10){print}}' | cut -f4,10 \
| sort -k1,1 \
| bedtools groupby -g 1 -c 2 -o distinct \
| awk '{print $1","$2}' \
> overlapping_gene_names.bed

conda activate python3
./sortClusters.py overlapping_gene_names.bed | sort | uniq > overlapping_gene_names.2.bed
```

After manual curation, we landed on 71 genes that were redundant and removed them from downstream analyses.

> We kept only fusion genes that encoded an alternative protein.

### 2.2 Obtaining genes with exons overlapping SD-98 genes

```bash
# obtaining gene body annotations in CHM13v1.0
cat /share/dennislab/projects/t2t/annotation_v1.0/cat_v4/CHM13.combined.v4.gff3 \
| awk '{if($3=="gene"){print}}' | ./parseGFF3.py | gff2bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10"\t0\t"$6}' \
| sed 's/gene_name=//g;s/source_gene=//g;s/biotype=//g;s/gene_id=//g' \
> CHM13.combined.v4.genes.bed

cat MSTRG.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0;/$1;/g" CHM13.combined.v4.genes.bed'

# obtaining exon annotations in CHM13v1.1
cat /share/dennislab/projects/t2t/annotation_v1.0/cat_v4/CHM13.combined.v4.gff3 \
| awk '{if($3=="exon"){print}}' | ./parseGFF3_exons.py | gff2bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10"\t0\t"$6}' \
| sed 's/gene_name=//g;s/source_gene=//g;s/biotype=//g;s/gene_id=//g' \
| sort | uniq > CHM13.combined.v4.exons.bed

cat MSTRG.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0;/$1;/g" CHM13.combined.v4.exons.bed'

# obtaining exons fully contained in SD-98 regions (autosomes only)
bedtools intersect -wa -f 1 -a CHM13.combined.v4.exons.bed -b <(grep -v "chrY\|chrM" ../coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.bed) > CHM13.combined.v4.exons.SD-98.bed

# obtaining gene body coordinates of genes that have and exon fully covered in SD-98 regions
# adding gene body coordinate to name
# intersecting with SD-98 to get only duplicated part of the gene only
bedtools intersect -F 1 -s -wao -a CHM13.combined.v4.genes.bed -b CHM13.combined.v4.exons.SD-98.bed \
| awk '{if($4==$10){print}}' \
| cut -f1-6 | sort | uniq \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4";"$1":"$2"-"$3,$5,$6}' \
| bedtools intersect -a - -b <(grep -v "chrY\|chrM" ../coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.bed) \
> CHM13.combined.v4.genes.SD-98.bed

# removing redundant gene names
grep -Fvf genes_to_remove.txt CHM13.combined.v4.genes.SD-98.bed > CHM13.combined.v4.genes.SD-98.nr.bed

# selecting protein coding and unprocessed pseudogenes
grep "protein_coding" CHM13.combined.v4.genes.SD-98.nr.bed > CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed
grep "unprocessed_pseudogene" CHM13.combined.v4.genes.SD-98.nr.bed > CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed
```

We obtained gene-centric files, with SD-98 regions collapsed:
```bash
awk '{print $4"\t"$1":"$2"-"$3}' <(cat CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed) \
| bedtools groupby -g 1 -c 2 -o distinct,count_distinct > CHM13.combined.v4.genes.SD-98.nr.protein_pseudo.tsv
```

We also obtained a bed file with original gene body coordinates (not only SD-98 portion):
```bash
cut -f4 CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed \
| sort | uniq \
| sed 's/;chr/\tchr/g' | awk '{print $2"\t"$1}' \
| awk 'BEGIN{OFS="\t"}{gsub(":","\t",$1); gsub("-","\t",$1); print}' \
> CHM13.combined.v4.genes.SD-98.nr.protein_coding.body.bed

cut -f4 CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed \
| sort | uniq \
| sed 's/;chr/\tchr/g' | awk '{print $2"\t"$1}' \
| awk 'BEGIN{OFS="\t"}{gsub(":","\t",$1); gsub("-","\t",$1); print}' \
> CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.body.bed
```

```bash
cut -f4 CHM13.combined.v4.genes.SD-98.bed | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 5413
cut -f4 CHM13.combined.v4.genes.SD-98.nr.bed | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 5342
cut -f4 CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 825
cut -f4 CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 1120

cut -f4 CHM13.combined.v4.genes.SD-98.bed | grep -v "chrX" | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 5154
cut -f4 CHM13.combined.v4.genes.SD-98.nr.bed | grep -v "chrX" | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 5083
cut -f4 CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed | grep -v "chrX" | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 698
cut -f4 CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed | grep -v "chrX" | sed 's/;/\t/g' | cut -f3 | sort | uniq | wc -l # 1095

grep "^chrX" CHM13.combined.v4.genes.SD-98.bed | cut -f4 | sort | uniq | wc -l # 259
grep "^chrX" CHM13.combined.v4.genes.SD-98.bed | grep "protein_coding\|unprocessed_pseudogene" | cut -f4 | sort | uniq | wc -l # 152

# T2T-CHM13v2.0
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.0_to_v1.1_rdna_merged.chain

cat CHM13.combined.v4.genes.SD-98.nr.protein_coding.bed CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.bed \
| cut -f4 | sed 's/;/\t/g' | awk '{print $5"\t"$3}' | sed 's/:/\t/g' | sed 's/-/\t/g' > SD-98_gene_coords.t2t-chm13v1.0.bed

conda activate ucsc
liftOver SD-98_gene_coords.t2t-chm13v1.0.bed v1.0_to_v1.1_rdna_merged.chain SD-98_gene_coords.t2t-chm13v2.0.bed unMapped 
```

Full conversion table:
```bash
cat /share/dennislab/projects/t2t/annotation_v1.0/cat_v4/CHM13.combined.v4.gff3 \
| awk '{if($3=="transcript"){print}}' | ./parseGFF3tx_full.py | gff2bed \
| cut -f10 \
| sed 's/gene_name=//g;s/gene_id=//g;s/gene_ensembl=//g;s/transcript_id=//g;s/transcript_name=//g;s/transcript_ensembl=//g;s/biotype=//g' \
| sed 's/;/\t/g' \
>  geneid2txid.tsv
```

## 3. Gene families in SD-98

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/families
```

### 3.1 Clusters of gene with shared exons

We classified protein coding genes and pseudogenes based on shared exons (using recently duplicated exons only).

```bash
conda activate variants

# first, we obtained exon annotations in CHM13v1.0
# and added exon coordinates to names
cat /share/dennislab/projects/t2t/annotation_v1.0/cat_v4/CHM13.combined.v4.gff3 \
| awk '{if($3=="exon"){print}}' | ../genes/parseGFF3_exons.py | gff2bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10"\t0\t"$6}' \
| sed 's/gene_name=//g;s/source_gene=//g;s/biotype=//g;s/gene_id=//g' \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4";;"$1":"$2"-"$3,$5,$6}' \
| sort | uniq \
> CHM13.combined.v4.exons.bed

cat ../genes/MSTRG.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0;/$1;/g" CHM13.combined.v4.exons.bed'

# we selected exons from SD-98 protein coding and pseudogenes and removed known redundant gene annotations
bedtools intersect -wa -f 1 \
-a CHM13.combined.v4.exons.bed \
-b <(grep -v "chrX\|chrY\|chrM" ../coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.bed) \
| grep -Fvf ../genes/genes_to_remove.txt \
> CHM13.combined.v4.exons.SD98.bed
# | grep "protein_coding\|unprocessed_pseudogene" \

# we obtained nucleotide sequence of SD-98 protein coding and psuedogenes exons
bedtools getfasta -s -name -fi /share/dennislab/projects/t2t/assembly/v1.0/chm13.draft_v1.0.fasta -bed CHM13.combined.v4.exons.SD98.bed \
> CHM13.combined.v4.exons.SD98.fa

# ... and mapped exons back to reference genome using LiftOff minimap2 parameters
minimap2 -c --end-bonus 5 --eqx -N 50 -p 0.5 -t 64 /share/dennislab/projects/t2t/assembly/v1.0/chm13.draft_v1.0.plusY.fasta CHM13.combined.v4.exons.SD98.fa > CHM13.combined.v4.exons.SD98.paf

# then, we converted PAF file to bed format
cat CHM13.combined.v4.exons.SD98.paf \
| cut -f1-11 \
| sed 's/::chr[^\t]*//' \
| awk '{print $6"\t"$8"\t"$9"\t"$1"\t0\t"$5}' \
> CHM13.combined.v4.exons.SD98.paf.bed

# secondly, we intersected paf mappings with SD-98 exon annotations 
bedtools intersect -f 0.99 -s -wao \
-a CHM13.combined.v4.exons.SD98.bed \
-b CHM13.combined.v4.exons.SD98.paf.bed \
| awk '{if($7!="."){print}}' \
> CHM13.combined.v4.exons.SD98.bedpe

# removing exon coordinates, sorting each row alphabetically 
# and get unique pairwise comparisons
cut -f4,10 CHM13.combined.v4.exons.SD98.bedpe \
| sed 's/;;chr[^\t]*\t/\t/' \
| sed 's/;;chr.*//' \
| awk '{split($0, a, "\t"); asort(a); for(i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" );}' \
| awk '{if($1!=$2){print $1","$2}}' \
| sort | uniq \
| grep "protein_coding\|unprocessed_pseudogene" \
> CHM13.combined.v4.exons.SD98.clusters.tsv 
```

### 3.2 DupMasker SD units

We also obtained families using DupMasker annotations.

> From SegDup T2T paper: "We annotated all T2T-CHM13 SDs with DupMasker (26), which defines ancestral evolutionary units of duplication on the basis of mammalian outgroups and a repeat graph (27). Focusing on duplicons that carry genes or duplicated portions of genes, we identified 30 duplicons that show the greatest copy number change between T2T-CHM13 and GRCh38. These 30 genic SDs represent regions where gene annotation is most likely to change; all predicted differences favor an increase in copy number for the T2T-CHM13 assembly (Fig. 1D and table S2)."

DupMasker file was downloaded from: `https://github.com/mrvollger/SafFire/blob/main/datasets/T2T-CHM13v2_1/T2T-CHM13v2_1_dupmasker_colors.bed` and/or `https://github.com/mrvollger/SafFire/blob/main/datasets/CHM13_v1.1/CHM13_v1.1_dupmasker_colors.bed`

Since DupMasker files were in CHM13v2.0, we grouped genes using this coordinate system. LiftOver was not as efficient lifting all the SD ancestral units.

```bash
# obtaining gene body annotations in CHM13v2
gunzip -c /share/dennislab/projects/hsd_genes/T2T/annotation/vep.chm13.v2.0.gff3.gz \
| awk '{if($3=="gene"){print}}' | ./parseGFF3.py | gff2bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10"\t0\t"$6}' \
| sed 's/gene_name=//g;s/source_gene=//g;s/biotype=//g;s/gene_id=//g' \
> vep.chm13.v2.0.genes.bed

cat ../genes/MSTRG.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0;/$1;/g" vep.chm13.v2.0.genes.bed'

# intersecting with ancestral duplicons
bedtools intersect -wao -s -a vep.chm13.v2.0.genes.bed -b T2T-CHM13v2_1_dupmasker_colors.bed | cut -f4,10 | sort | uniq | bedtools groupby -g 1 -c 2 -o collapse > vep.chm13.v2.0.genes.duplicons.tsv

# selecting SD-genes only
grep -Ff SD98_genes_all.txt vep.chm13.v2.0.genes.duplicons.tsv > vep.chm13.v2.0.genes.duplicons.SD98.tsv
```

## 4. WSSD gene family CN

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/wssd
```

Preparing WSSD data:
```bash
mkdir -p data
cd data

# copying SGDP WSSD copy numbers into data folder
# important: we removed LP6005442-DNA_A08.CN.bed because it had outlier copy numbers
for file in /share/dennislab/databases/data/wssd-t2t/t2t_v1.0/sgdp_complete/*bed; do
  filename=$(basename -s _wssd.bed $file)
  tail -n +2 $file | cut -f1-3,10 > $filename.CN.bed
done

# copying NHP WSSD copy numbers into data folder
ls /share/dennislab/databases/data/wssd-t2t/t2t_v1.0/t2tdp/*CN.bed | grep -v "/HG\|/NA\|/CHM" | xargs -n1 -P1 bash -c 'cp $0 .'
```

Genotyping WSSD copy number of SD-98 protein coding genes and unprocessed pseudogenes:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/wssd

split -l 500 -d <(cat ../genes/CHM13.combined.v4.genes.SD-98.protein_coding.bed ../genes/CHM13.combined.v4.genes.SD-98.unprocessed_pseudogenes.bed | cut -f1-4) part_

conda activate python3
cp /share/dennislab/users/dcsoto/ms_san/03_SDs/quickmer2/genotype_cn_parallel.py . 

for num in 00 01 02 03; do
    python3 genotype_cn_parallel.py --path data --genes part_${num} --output part_${num}.wssd_cn.tsv -t 70
done
```

Merging loci:
```r
# module load R/4.0.1
files <- list.files(pattern = "_cn.tsv$")
list_df <- lapply(files, read.table, header=TRUE)
joined_df <- do.call("rbind", list_df)
write.table(joined_df, "wssd_sd98.t2t-chm13-v1.0.tsv", sep="\t", quote=FALSE, row.names = FALSE)
```

We genotyped GYPA unique space (t2t-chm13v1.0: chr4:147419057-147429057):
```bash
conda activate python3
python3 genotype_cn_parallel.py --path data --genes gypa_regions.bed --output gypa_regions.wssd_cn.tsv -t 120
```

## 5. Cancer somatic hotspots

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/somatic_hotspots/
```

We downloaded coordinates of pan-cancer somatic hotpots (somatic recurrent breakpoints or SRBs) in Hg19 coordinates, and lifted them over to T2T-CHM13v1.0.

```bash
# liftover chain from hg19 to hg38
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz

# liftover chain from hg38 to t2t-chm13v1.0
wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/hg38Lastz/hg38.t2t-chm13-v1.0.over.chain.gz
gunzip hg38.t2t-chm13-v1.0.over.chain.gz

# lifting over SRBs coordinates
conda activate variants
liftOver SRBs.hg19.bed hg19ToHg38.over.chain SRBs.hg38.bed SRBs.hg19.unmapped.bed
liftOver SRBs.hg38.bed hg38.t2t-chm13-v1.0.over.chain SRBs.t2t-chm13-v1.0.bed SRBs.hg38.unmapped.bed

# intersecting with SD-98 genes
bedtools intersect -wao -a <(cat ../genes/CHM13.combined.v4.genes.SD-98.protein_coding.gene_body.bed ../genes/CHM13.combined.v4.genes.SD-98.unprocessed_pseudogenes.gene_body.bed | sort | uniq) -b  SRBs.t2t-chm13-v1.0.bed | awk '{if($5!="."){print}}' \
> sd98_genes.somatic_hotspots.bed
```

Are pan-cancer somatic hotspots closer than expected by random chance to SDs and SD-98s?

Generating random regions of equal size as SD and SD98:
```bash
# SDs
mkdir -p random/random_SD
for num in `seq 1 10000`; do 
    bedtools shuffle -chrom -i ../coordinates/chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed -g <(cut -f1,2 /share/dennislab/projects/t2t/assembly/v1.0/chm13.draft_v1.0.plusY.fasta.fai | grep -v "^chrM\|^chrY\|^chrX") -noOverlapping -maxTries 10000 -f 0.1 | bedtools sort > random/random_SD/random_segdup_${num}.bed; 
done

# SDs-98
mkdir -p random/random_SD98
for num in `seq 1 10000`; do 
    bedtools shuffle -chrom -i ../coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed -g <(cut -f1,2 /share/dennislab/projects/t2t/assembly/v1.0/chm13.draft_v1.0.plusY.fasta.fai | grep -v "^chrM\|^chrY\|^chrX") -noOverlapping -maxTries 10000 -f 0.1 | bedtools sort > random/random_SD98/random_segdup98_${num}.bed; 
done
```

Obtaining distance of SRBs to random SD and SD98:
```bash
module load R/4.0.1

mkdir -p distance/distance_SD
for num in `seq 1 10000`; do 
    bedtools sort -i SRBs.t2t-chm13-v1.0.bed | bedtools closest -d -t first -a - -b random/random_SD/random_segdup_${num}.bed | cut -f8 > distance/distance_SD/random_segdup_${num}.tsv
    cat distance/distance_SD/random_segdup_${num}.tsv | Rscript median.R
done > SRBs.distance_SD.txt

mkdir -p distance/distance_SD98
for num in `seq 1 10000`; do 
    bedtools sort -i SRBs.t2t-chm13-v1.0.bed | bedtools closest -d -t first -a - -b random/random_SD98/random_segdup98_${num}.bed | cut -f8 > distance/distance_SD98/random_segdup98_${num}.tsv
    cat distance/distance_SD98/random_segdup98_${num}.tsv | Rscript median.R
done > SRBs.distance_SD98.txt
```

Calculating true observed distance:
```bash
module load R/4.0.1

bedtools sort -i SRBs.t2t-chm13-v1.0.bed \
| bedtools closest -d -t first -a - -b ../coordinates/chm13.draft_v1.0_plus38Y.SDs.merged.auto.bed \
| cut -f8 \
| Rscript median.R # 55806.5

bedtools sort -i SRBs.t2t-chm13-v1.0.bed \
| bedtools closest -d -t first -a - -b ../coordinates/chm13.draft_v1.0_plus38Y.SDs-98.merged.auto.bed \
| cut -f8 \
| Rscript median.R # 1880816
```

## 6. Genomic disorders

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/genomic_disorders
```

We obtained coordinates of autosomal genomic disorders in Hg19 from Satterstrom et al. 2020.

We lifted coordinates to T2T-CHM13v1.0:
```bash
conda activate variants
liftOver genomic_disorders.hg19.bed ../somatic/hg19ToHg38.over.chain genomic_disorders.hg38.bed genomic_disorders.hg19.bed_unmapped
liftOver genomic_disorders.hg38.bed ../somatic/hg38.t2t-chm13-v1.0.over.chain genomic_disorders.t2t-chm13-v1.0.bed genomic_disorders.hg38.bed_unmapped
```

From 51 regions, 48 were lifted over, and three were problematic.

We lifted over the borders of problematic regions:
```bash
conda activate variants

liftOver problematic_regions_1.hg19.bed ../somatic/hg19ToHg38.over.chain problematic_regions_1.hg38.bed problematic_regions_1.hg19.bed_unmapped
liftOver problematic_regions_1.hg38.bed ../somatic/hg38.t2t-chm13-v1.0.over.chain problematic_regions_1.t2t-chm13-v1.0.bed problematic_regions_1.hg38.bed_unmapped

liftOver problematic_regions_2.hg38.bed ../somatic/hg38.t2t-chm13-v1.0.over.chain problematic_regions_2.t2t-chm13-v1.0.bed problematic_regions_2.hg38.bed_unmapped
liftOver problematic_regions_3.hg38.bed ../somatic/hg38.t2t-chm13-v1.0.over.chain problematic_regions_3.t2t-chm13-v1.0.bed problematic_regions_3.hg38.bed_unmapped
```

We added problematic regions to coordinates file for genomic disorders and saved them as `genomic_disorders.t2t-chm13-v1.0.all.bed`. 

We also generated an expanded genomic disorders coordinate set, including additiona +- 500kbp saved as `sd98_genes.genomic_disorders_expanded.bed`.

Then, we intersected with SD-98 genes:
```bash
# all
bedtools intersect -wao -a <(cat ../genes/CHM13.combined.v4.genes.SD-98.protein_coding.gene_body.bed ../genes/CHM13.combined.v4.genes.SD-98.unprocessed_pseudogenes.gene_body.bed | sort | uniq) -b genomic_disorders.t2t-chm13-v1.0.all.bed | awk '{if($5!="."){print}}' \
> sd98_genes.genomic_disorders.bed

# all +- 500 kbp
bedtools intersect -wao -a <(cat ../genes/CHM13.combined.v4.genes.SD-98.nr.protein_coding.body.bed ../genes/CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.body.bed) -b genomic_disorders.t2t-chm13-v1.0.expanded.bed | awk '{if($5!="."){print}}' \
> sd98_genes.genomic_disorders_expanded.bed
```

And we intersect with all genes:
```bash
# all
bedtools intersect -wao -a genomic_disorders.t2t-chm13-v1.0.all.bed -b ../genes/CHM13.combined.v4.genes.bed | cut -f 1-4,8 | tr ';' '\t' | cut -f1-5 | bedtools groupby -g 1,2,3,4 -c 5 -o distinct

# all +- 500 kbp in bed
bedtools intersect -wao -a genomic_disorders.t2t-chm13-v1.0.expanded.bed -b ../genes/CHM13.combined.v4.genes.bed | cut -f 1-4,8 | tr ';' '\t' | cut -f1-5 | bedtools groupby -g 1,2,3,4 -c 5 -o distinct > genomic_disorders.t2t-chm13-v1.0.expanded.genes.bed

# all +- 500 kbp in tsv
bedtools intersect -wao -a genomic_disorders.t2t-chm13-v1.0.expanded.bed -b ../genes/CHM13.combined.v4.genes.bed | cut -f 1-4,8 | awk -F'\t' 'BEGIN{OFS="\t"}{gsub(";","\t",$5); print}' | awk -F'\t' '{OFS="\t"}{print $5"\t"$7"\t"$4}' > genomic_disorders_genes.tsv
```

## 7. Intersection w/Vollger table

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/vollger
```

```bash
conda activate variants

cut -f1-3,5 Genic_SD_expansions_in_T2T_relative_to_Clint.txt | tail -n +2 | bedtools sort | bedtools merge -c 4 -o collapse > Genic_SD_expansions_in_T2T_relative_to_Clint.bed

bedtools intersect -wao -a <(cat ../genes/CHM13.combined.v4.genes.SD-98.nr.protein_coding.body.bed ../genes/CHM13.combined.v4.genes.SD-98.nr.unprocessed_pseudogenes.body.bed | sed 's/;/\t/g' | cut -f1-4,6) -b Genic_SD_expansions_in_T2T_relative_to_Clint.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6":"$7"-"$8}' | sed 's/.:-1--1/./g' | awk '{if($6 != "."){print}}' | bedtools groupby -g 1,2,3,4,5 -c 6 -o collapse > SD98.Genic_SD_expansions_in_T2T_relative_to_Clint.tsv
```

## 8. Chromosome Y

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/sd-98/chrY

conda activate variants

wget https://hgdownload.soe.ucsc.edu/gbdb/hs1/sedefSegDups/sedefSegDups.bb
bigBedToBed sedefSegDups.bb sedefSegDups.bed
awk '{if($24>=0.98){print}}' sedefSegDups.bed | bedtools sort | bedtools merge > sedefSegDups.sd98.bed

# obtaining gene body annotations in CHM13v1.0
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3

cat chm13.draft_v2.0.gene_annotation.gff3 \
| awk '{if($3=="gene"){print}}' | ./parseGFF3.py | gff2bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10"\t0\t"$6}' \
| sed 's/gene_name=//g;s/source_gene=//g;s/biotype=//g;s/gene_id=//g' \
> chm13v2.0.genes.bed

cat MSTRG.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0;/$1;/g" chm13v2.0.genes.bed'

# obtaining exon annotations in CHM13v1.1
cat chm13.draft_v2.0.gene_annotation.gff3  \
| awk '{if($3=="exon"){print}}' | ./parseGFF3_exons.py | gff2bed \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$10"\t0\t"$6}' \
| sed 's/gene_name=//g;s/source_gene=//g;s/biotype=//g;s/gene_id=//g' \
| sort | uniq > chm13v2.0.exons.bed

cat MSTRG.tsv | xargs -n2 -P1 bash -c 'echo $0; sed -i "s/$0;/$1;/g" chm13v2.0.exons.bed'

# obtaining exons fully contained in SD-98 regions
bedtools intersect -wa -f 1 -a chm13v2.0.exons.bed -b <(grep -v "chrM" sedefSegDups.sd98.bed) > chm13v2.0.exons.sd98.bed

# obtaining gene body coordinates of genes that have an exon fully covered in SD-98 regions
# adding gene body coordinate to name
# intersecting with SD-98 to get only duplicated part of the gene
bedtools intersect -F 1 -s -wao -a chm13v2.0.genes.bed -b chm13v2.0.exons.sd98.bed \
| awk '{if($4==$10){print}}' \
| cut -f1-6 | sort | uniq \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4";"$1":"$2"-"$3,$5,$6}' \
| bedtools intersect -a - -b <(grep -v "chrM" sedefSegDups.sd98.bed) \
> chm13v2.0.genes.sd98.bed

# we obtained gene-centric files, with SD-98 regions collapsed:
awk '{print $4"\t"$1":"$2"-"$3}' chm13v2.0.genes.sd98.bed \
| bedtools groupby -g 1 -c 2 -o distinct,count_distinct \
> chm13v2.0.genes.sd98.tsv

# selecting genes only in chromosome Y
awk '$2 ~ /chrY:/' chm13v2.0.genes.sd98.tsv > chm13v2.0.genes.sd98.chrY.tsv

cut -f1 chm13v2.0.genes.sd98.chrY.tsv | sort | uniq | wc -l # 503
cut -f1 chm13v2.0.genes.sd98.chrY.tsv | sort | uniq | grep "protein_coding\|unprocessed_pseudogene\|unknown_likely_coding" | wc -l # 326
```

## 9. Novel/fixed SD98 genes in T2T-CHM13

```bash
bedtools intersect -wao -a sd98_genes.bed -b hg38_nonsyntenic.chm13_coords.bed \
| cut -f1-4,8 \
> sd98_genes.non_syntenic.bed

bedtools intersect -wao -a sd98_genes.bed -b hg38_collapsed_dups.chm13_coords.bed \
| cut -f1-4,8 \
| bedtools groupby -g 1,2,3,4 -c 5 -o sum \
> sd98_genes.collapsed_dups.bed

bedtools intersect -wao -a sd98_genes.bed -b hg38_false_dups.chm13_coords.bed \
| cut -f1-4,8 \
| bedtools groupby -g 1,2,3,4 -c 5 -o sum \
> sd98_genes.false_dups.bed
```
