import pandas as pd


SAMPLES = config['samples']
REFERENCE = config['reference']
REGIONS = config['regions']


samples = pd.read_csv(SAMPLES, delimiter = '\t', index_col=0, names = ['sample', 'sample_name', 'hap1', 'hap2'])
samples_dict = samples.to_dict('index')
regions = pd.read_csv(REGIONS, delimiter = '\t', index_col=0, names = ['locus', 'coords'])
regions_dict = regions.to_dict('index') 


rule all:
  input:
    expand("02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt.bam", smp = samples_dict.keys(), hap = [1, 2], locus = regions_dict.keys()),
    expand("03_variants_pHSDs/{smp}.dip.{locus}.vcf.gz.tbi", smp = samples_dict.keys(), locus = regions_dict.keys()), 
    #expand("03_callable_pHSDs/{smp}.1.{locus}.bed", smp = samples_dict.keys(), locus = regions_dict.keys())    


def get_hap1(wildcards):
  return samples_dict[wildcards.smp]['hap1']
def get_hap2(wildcards):
  return samples_dict[wildcards.smp]['hap2']


rule mapping:
  input:
    h1 = get_hap1,
    h2 = get_hap2
  output:
    h1 = "01_alignments/{smp}.1.bam",
    h2 = "01_alignments/{smp}.2.bam"
  threads: 8
  shell:
    '''
    minimap2 -t {threads} -a --eqx --cs -x asm5 --secondary=no -s 25000 -K 8G {REFERENCE} {input.h1} \
    | samtools view -F 4 -b - > {output.h1}

    minimap2 -t {threads} -a --eqx --cs -x asm5 --secondary=no -s 25000 -K 8G {REFERENCE} {input.h2} \
    | samtools view -F 4 -b - > {output.h2}
    '''


rule bam_sorting:
  input:
    "01_alignments/{smp}.{hap}.bam"
  output:
    "01_alignments/{smp}.{hap}.sort.bam"
  threads: 4
  shell:
    '''
    samtools sort -@ {threads} {input} > {output} 
    '''


rule bam_index:
  input:
    "01_alignments/{smp}.{hap}.sort.bam"
  output:
    "01_alignments/{smp}.{hap}.sort.bam.bai"
  threads: 1
  shell:
    '''
    samtools index {input}
    '''


def get_coords(wildcards):
  return regions_dict[wildcards.locus]['coords']


rule parse_alignments:
  input:
    bam = "01_alignments/{smp}.{hap}.sort.bam",
    index = "01_alignments/{smp}.{hap}.sort.bam.bai"
  output:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.bam"
  params:
    coords = get_coords
  threads: 1
  shell:
    '''
    samtools view -Sb -F 256 {input.bam} {params.coords} > {output}
    '''


rule convert_to_paf:
  input:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.bam"
  output:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf"
  threads: 1
  shell:
    '''
    samtools view -h {input} | paftools.js sam2paf - > {output}
    '''


rule select_longest:
  input:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf"
  output:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt"
  threads: 1
  shell:
    """
    max=$(awk 'BEGIN{{a=0}}{{if ($11>0+a) a=$11}} END{{print a}}' {input})
    awk -v var=$max '{{if ($11==var) {{print $0}} }}' {input} > {output}
    """


def get_fasta(wildcards):
  if wildcards.hap == "1":
     return samples_dict[wildcards.smp]['hap1']
  elif wildcards.hap == "2":
     return samples_dict[wildcards.smp]['hap2']  


rule convert_to_sam:
  input:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt"
  output:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt.sam"
  params:
    fasta = get_fasta
  shell:
    '''
    rb paf-to-sam -f {params.fasta} {input} > {output}
    '''

 
rule convert_to_bam:
  input:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt.sam"
  output:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt.bam"
  shell:
    '''
    samtools view -Sb {input} > {output}
    '''


rule get_callable_haps:
  input:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt"
  output:
    "02_alignments_pHSDs/{smp}.{hap}.{locus}.paf.filt.var"
  threads: 1
  shell:
    '''
    cat {input} | sort -k6,6 -k8,8n | paftools.js call -q 0 -L 0 -l 0 - > {output}
    '''


rule get_callable_diploid:
  input:
    h1 =  "02_alignments_pHSDs/{smp}.1.{locus}.paf.filt.var",
    h2 =  "02_alignments_pHSDs/{smp}.2.{locus}.paf.filt.var"
  output:
    h1 = "03_callable_pHSDs/{smp}.1.{locus}.bed",
    h2 = "03_callable_pHSDs/{smp}.2.{locus}.bed"
  threads: 1
  shell:
    '''
    cat {input.h1} | grep ^R | cut -f2- > {output.h1}
    cat {input.h2} | grep ^R | cut -f2- > {output.h2}
    '''


rule call_variants:
  input:
    h1 = "02_alignments_pHSDs/{smp}.1.{locus}.paf.filt.bam",
    h2 = "02_alignments_pHSDs/{smp}.2.{locus}.paf.filt.bam"
  output:
    "03_variants_pHSDs/{smp}.pair.{locus}.vcf"
  shell:
    '''
    htsbox pileup -q 0 -evcf {REFERENCE} {input.h1} {input.h2} > {output}
    '''


rule get_diploid_calls:
  input:
    "03_variants_pHSDs/{smp}.pair.{locus}.vcf"
  output:
    "03_variants_pHSDs/{smp}.dip.{locus}.vcf"
  shell:
    '''
    dipcall-aux.js vcfpair -s {wildcards.smp} -a {input} \
    | bcftools norm -Ov -m-any \
    | bcftools norm -Ov -d exact \
    | bcftools norm -Ov -m-any --fasta-ref {REFERENCE} --check-ref w \
    | bcftools sort -m 4G > {output} 
    '''


rule compress_calls:
  input:
    "03_variants_pHSDs/{smp}.dip.{locus}.vcf"
  output:
    "03_variants_pHSDs/{smp}.dip.{locus}.vcf.gz"
  shell:
    '''
    bgzip -c {input} > {output}
    '''


rule index_calls:
  input:
    "03_variants_pHSDs/{smp}.dip.{locus}.vcf.gz"
  output:
    "03_variants_pHSDs/{smp}.dip.{locus}.vcf.gz.tbi" 
  shell:
    '''
    tabix {input}
    '''


