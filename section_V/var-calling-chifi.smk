# ------------------------------------------
# Snakefile to call variants in PacBio cHiFi
# ------------------------------------------

import pandas as pd

# CONFIGURATION

REFERENCE = config['reference']
REFBASENAME = os.path.splitext(REFERENCE)[0]
FILES = config['files'] # 2 cols: sample, path to BAM. BAM must be indexed.
ROI = config['roi']
MAPQ = [1, 2, 5, 10, 20]

# VARIABLE PARSING

df = pd.read_csv(FILES, sep="\t", names=['samples','paths'])
samples = df['samples'].to_list()
bams = dict(zip(df['samples'].to_list(), df['paths'].to_list())) 

# TARGETS

rule all:
  input:
    expand("results/snvs/{smp}.gatk_{mapq}.gvcf", smp=samples, mapq=MAPQ),
#    expand("results/svs/sniffles/{smp}.sniffles.vcf", smp=samples),
#    expand("results/svs/pbsv/{smp}.pbsv.vcf", smp=samples),
#    expand("results/svs/cutesv/{smp}.cutesv.vcf", smp=samples),
#    expand("results/svs/svim/{smp}.svim.vcf", smp=samples),
#    expand("results/svs/sniffles2/{smp}.sniffles2.snf", smp=samples)

# INDEXING

rule reference_faidx:
  input:
    REFERENCE
  output:
    REFERENCE+".fai"
  shell:
    """
    samtools faidx {input}
    """

rule reference_dict:
  input:
    REFERENCE
  output:
    REFBASENAME+".dict"
  shell:
    """
    picard CreateSequenceDictionary R={input} O={output}
    """

# MAPPING

rule mapping:
  input:
    lambda wildcards: bams.get(wildcards.smp)
  output:
    bam = temp("results/alignments/{smp}.dedup.pbmm2.bam"),
    bai = temp("results/alignments/{smp}.dedup.pbmm2.bam.bai")
  params:
    reference = REFERENCE
  threads: 5
  shell:
    """
    pbmm2 align -j {threads} --preset CCS --sort {params.reference} {input} {output.bam}
    """

#rule primary:
#  input:
#    "results/alignments/{smp}.dedup.pbmm2.bam"
#  output:
#    temp("results/alignments/{smp}.dedup.pbmm2.primary.bam")
#  threads: 2
#  shell:
#    """
#    samtools view -@ {threads} -b -F 260 {input} > {output}
#    """

rule read_groups:
  input:
    "results/alignments/{smp}.dedup.pbmm2.bam"
  output:
    temp("results/alignments/{smp}.dedup.pbmm2.rg.bam")
  shell:
    """
    picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I={input} O={output} RGID=1 RGLB=lib RGPL=PB RGPU=unit1 RGSM={wildcards.smp}
    """

rule md_tag:
  input:
    "results/alignments/{smp}.dedup.pbmm2.rg.bam"
  output:
    "results/alignments/{smp}.dedup.pbmm2.rg.md.bam"
  threads: 2
  shell:
    """
    samtools calmd -@ {threads} -b {input} {REFERENCE} > {output}
    """

rule bam_index:
  input:
    "results/alignments/{smp}.dedup.pbmm2.rg.md.bam"
  output:
    "results/alignments/{smp}.dedup.pbmm2.rg.md.bam.bai"
  threads: 2
  shell:
    """
    samtools index -@ {threads} {input}
    """

# SNV CALLING

rule gatk:
  input:
    bam = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam",
    index = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam.bai",
    dict = REFBASENAME+".dict",
    faidx = REFERENCE+".fai"
  output:
    vcf = "results/snvs/{smp}.gatk_{mapq}.gvcf"
  params:
    reference = REFERENCE,
    roi = ROI
  shell:
    """
    gatk HaplotypeCaller \
      --minimum-mapping-quality {wildcards.mapq} \
      --mapping-quality-threshold-for-genotyping {wildcards.mapq} \
      --ploidy 2 \
      -I {input.bam} -R {params.reference} -O {output.vcf} -L {params.roi} -ERC GVCF
    """

# SV CALLING

rule Sniffles:
  input:
    bam = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam",
    bai = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam.bai",
    ref = REFERENCE,
    fai = REFERENCE+".fai"
  output:
    vcf = "results/svs/sniffles/{smp}.sniffles.vcf",
    tmp = temp("results/svs/sniffles/{smp}.sniffles.tmp.vcf")
  threads: 5
  shell:
    """
    sniffles -t {threads} -m {input.bam} -v {output.tmp}
    awk 'BEGIN{{OFS="\\t"}}{{if($1~"^#CHROM"){{col=NF; $col="SAMPLE"; print}}else{{print}}}}' {output.tmp} > {output.vcf}
    """

rule cuteSV:
  input:
    bam = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam",
    bai = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam.bai",
    ref = REFERENCE,
    fai = REFERENCE+".fai"
  output:
    vcf = "results/svs/cutesv/{smp}.cutesv.vcf",
    tmp = temp("results/svs/cutesv/{smp}.cutesv.tmp.vcf"),
    dir = temp(directory("results/tmp_cuteSV/cutesv_{smp}"))
  shell:
    """
    mkdir -p {output.dir}
    cuteSV {input.bam} {input.ref} {output.tmp} {output.dir}
    awk 'BEGIN{{OFS="\\t"}}{{if($1~"^#CHROM"){{col=NF; $col="SAMPLE"; print}}else{{print}}}}' {output.tmp} > {output.vcf}
    """

rule SVIM:
  input:
    bam = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam",
    bai = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam.bai",
    ref = REFERENCE,
    fai = REFERENCE+".fai"
  output:
    vcf = "results/svs/svim/{smp}.svim.vcf",
    tmp = temp(directory("results/tmp_SVIM/svim_{smp}"))
  shell:
    """
    svim alignment {output.tmp} {input.bam} {input.ref}
    awk 'BEGIN{{OFS="\\t"}}{{if($1~"^#CHROM"){{col=NF; $col="SAMPLE"; print}}else{{print}}}}' {output.tmp}/variants.vcf > {output.vcf}
    """ 

rule pbsv_signatures:
  input:
    bam = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam",
    bai = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam.bai"
  output:
    sig = "results/svs/pbsv/{smp}.svsig.gz"
  shell:
    "pbsv discover {input.bam} {output.sig}"

rule pbsv_call:
  input:
    sig = "results/svs/pbsv/{smp}.svsig.gz",
    ref = REFERENCE,
    fai = REFERENCE+".fai"
  output:
    vcf = "results/svs/pbsv/{smp}.pbsv.vcf",
    tmp = temp("results/svs/pbsv/{smp}.pbsv.tmp.vcf")
  shell:
    """
    pbsv call --ccs {input.ref} {input.sig} {output.tmp}
    awk 'BEGIN{{OFS="\\t"}}{{if($1~"^#CHROM"){{col=NF; $col="SAMPLE"; print}}else{{print}}}}' {output.tmp} > {output.vcf}
    """

rule Sniffles2:
  input:
    bam = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam",
    bai = "results/alignments/{smp}.dedup.pbmm2.rg.md.bam.bai",
    ref = REFERENCE,
    fai = REFERENCE+".fai"
  output:
    snf = "results/svs/sniffles2/{smp}.sniffles2.snf"
  threads: 5
  shell:
    """
    sniffles --threads {threads} --input {input.bam} --snf {output.snf}
    """
 
