# Running QuickMer2 (Dani's version)

QuickMer2 is a tool developed by the Kidd Lab that calculates paralog-specific copy-number based on unique k-mers.

QuickMer2 Snakemake pipeline takes as input a reference genome (in fasta format), a file containing genomic regions that do not exhibit copy-number polymorphism (or control file, in bed format), and a path to a folder containing sequencing data from samples (in cram, bam or fasta formats).

## Setting the working directory

The working directory needs to have the following structure:
```
├── quickmer2.smk
├── data
│   ├── sample1.cram
│   └── sample2.cram
│   └── ...
├── reference
│   ├── reference.fasta
│   └── control.bed
├── scripts
│   ├── fasta_len.awk
│   └── bedToBed9.py
└── output
```

Custom scripts can be found at the bottom of this Gist.

## Setting up software environment

The pipeline requires the following requirements: 
- Python3 installation with pandas and argparse
- Samtools
- UCSC's bedToBigBed
- Export path to [QuickMer2](https://github.com/KiddLab/QuicK-mer2) 

Using Conda:
```bash
conda create -n quickmer2 -c conda-forge -c bioconda python=3 pandas argparse samtools ucsc-bedtobigbed matplotlib
conda activate quickmer2
```

The latest version of QuickMer2 is already installed in the Dennis Lab cluster and can be exported as follows:
```bash
export PATH=/share/dennislab/databases/data/quickmer-t2t/QuicK-mer2/:$PATH
```

## Generating control file

The control file should contain regions that are not copy-number polymorphic. To obtain this files sex chromosomes, mithocondrial DNA, and SegDups should be excluded.

```bash
bedtools complement -i segdups.bed -g reference.fa.fai | grep -v "chrM\|chrX\|chrY" > control.bed
```

## Running QuickMer2 pipeline

```bash
snakemake \
--snakefile quickmer2.smk \
--config reference=reference/reference.fasta control=reference/control.bed input=data outdir=output distance=1 window=500 \
-p -j 20
```

## Notes

* If sample's raw data is in paired-end fastq files, concatenate into one fastq file prior to running the pipeline.
* QuicKmer2 does not process `fastq.gz` files correctly. Convert to `fasta` format before running the pipeline.
* 1KGP high-coverage genomes were mapped to the following reference: `/share/dennislab/programs/QuicK-mer2/references/GRCh38_full_analysis_set_plus_decoy_hla.fa`
