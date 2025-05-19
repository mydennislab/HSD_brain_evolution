# Expression analysis of SD-98 genes

## 1. Download data

### 1.1. BrainSpan

Working directory:
```bash
cd /share/dennislab/databases/data/BrainSpan
```

The data was downloaded following the instructions in https://www.ncbi.nlm.nih.gov/sra/docs/sra-dbgap-download/#downloading-the-data.
Importantly, we needed to use SRA Toolkit > v2.10.2, since starting with SRA Toolkit version 2.10.2, there are several important changes:
- You no longer need to import the NGC file to the configuration
- The NGC file will need to be specified as part of the command line every time you run a tool
- For SRA Runs, you no longer have an option to create a cart, but will need to use a list of Run accessions

```bash
conda activate sratools

cat SraAccList.txt | xargs -n1 -P30 bash -c 'prefetch --ngc prj_10998_D37641.ngc $0'
cat SraAccList.txt | xargs -n1 -P30 bash -c 'fasterq-dump --ngc prj_10998_D37641.ngc $0'
ls *fastq | xargs -n1 -P30 bash -c 'gzip $0'
```

## 2. Remapping expression data

Working directory:
```bash
cd /share/dennislab/users/dcsoto/ms_hsd/rna-seq
```

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/rna-seq/annotations

gffread -w CHM13.v2.0.transcripts.fa -g /share/dennislab/databases/assemblies/T2T/CHM13_T2T_v2/CAT-Liftoff/chm13v2.0.fa /share/dennislab/databases/assemblies/T2T/CHM13_T2T_v2/CAT-Liftoff/CHM13.v2.0.sorted.gtf
```

Rotation student Cassie Olivas simplified the RNA-seq/salmon pipeline and remapped several previously published datasets to T2T-CHM13v2.0. Additionally, I modified the pipeline to provide raw counts instead of TPMs.

We re-mapped the CORTECON and BrainSpan datasets:
```bash
conda activate rna-seq

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_vandeLeemput.yaml -p -j 80
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_brainspan.yaml -p -j 120
```

Additionally, we remapped other datasets from human LCL and human developmental brain, including human apical and radial glia (Florio), SH-SY5Y neuroblastoma cells (Pezzini), germinal zones of human neocortex (Fietz).

LCL datasets where an individual was sequenced in different runs were concatenated as follows:
```bash
cat samples | xargs -n1 -P10 bash -c '
runs=$(ls $0* | tr "\n" " ")
cat $runs > concat/$0.fastq.gz'
```

```bash
conda activate rna-seq

/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_Pezzini.yaml -p -j 80
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_Florio.yaml -p -j 80
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_Fietz.yaml -p -j 80
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_Pickrell.yaml -p -j 80
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_GM12878.yaml -p -j 80 -n
/share/dennislab/users/dcsoto/Miniconda3/bin/snakemake --configfile config_blake.yaml -p -j 80 -n
```

### 3. Gene to transcript conversion

```bash
cd /share/dennislab/users/dcsoto/ms_hsd/rna-seq/annotations

/share/dennislab/users/dcsoto/ms_hsd/sd-98/genes

cat chm13.draft_v2.0.gene_annotation.gff3 \
| awk '{if($3=="transcript"){print}}' \
| ./parseGFF3tx_mod.py \
| gff2bed \
| cut -f 10 \
| sed 's/gene_name=//g; s/gene_id=//g; s/transcript_id=//g' \
| sed 's/;/\t/g' \
> geneid2txid.tsv
```
