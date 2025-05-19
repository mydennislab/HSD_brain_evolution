# --------------------
# QuickMer-2 Snakefile
# Version: dcsoto
# --------------------

import re

# CONFIGURATION

REFERENCE = config['reference']
CONTROL = config['control']
INPUT = config["input"]
OUTDIR = config["outdir"]
DISTANCE = config["distance"]
WINDOW = config["window"]

# PARSE FILENAMES

FILES, = glob_wildcards(INPUT+"/{files}")

def parse_filename(file):
  filename = re.split("\.", file)
  sample = filename[0]
  extension = filename[-2] if filename[-1]=="gz" else filename[-1]
  return sample, (extension, file)

sample_info = dict([parse_filename(x) for x in FILES])
samples = list(sample_info.keys())

# TARGETS

rule all:
  input:
     expand(OUTDIR+"/CN/{smp}.CN.bed", smp = samples),
     expand(OUTDIR+"/bigBed/{smp}.CN.bb", smp = samples)

# QUICKMER2

rule search:
  input: 
    ref = REFERENCE, 
    control = CONTROL
  output: 
    OUTDIR+"/index/ref.fasta",
    OUTDIR+"/index/ref.fasta.qgc",
    OUTDIR+"/index/ref.fasta.qm",
    OUTDIR+"/index/ref.fasta.bed"
  params:
    ref = OUTDIR+"/index/ref.fasta"
  threads: 20
  shell:
    """
    cp {REFERENCE} {params.ref}
    quicKmer2 search -t {threads} -c {input.control} {params.ref} -e {DISTANCE} -w {WINDOW}
    """

rule count:
  input:
    reads = lambda wildcards: INPUT+"/"+sample_info.get(wildcards.smp)[1],
    qgc = OUTDIR+"/index/ref.fasta.qgc",
    qm = OUTDIR+"/index/ref.fasta.qm",
    bed = OUTDIR+"/index/ref.fasta.bed"
  output:
    OUTDIR+"/bins/{smp}.count.bin"
  params:
    ref = OUTDIR+"/index/ref.fasta",
    prefix = OUTDIR+"/bins/{smp}.count"
  threads: 6
  run:
    if sample_info.get(wildcards.smp)[0] == "cram":
      shell("""samtools view -F 3840 {input.reads} -T /share/dennislab/programs/QuicK-mer2/references/GRCh38_full_analysis_set_plus_decoy_hla.fa --input-fmt-option required_fields=0x202 | awk '{{print ">\\n"$10}}' | quicKmer2 count -t {threads} {params.ref} /dev/fd/0 {params.prefix}""")
    elif sample_info.get(wildcards.smp)[0] == "bam":
      shell("""samtools view -F 3840 {input.reads} | awk '{{print ">\\n"$10}}' | quicKmer2 count -t {threads} {params.ref} /dev/fd/0 {params.prefix}""") 
    else:
      shell("quicKmer2 count -t {threads} {params.ref} {input.reads} {params.prefix}")

rule estimate:
  input:
    bin = OUTDIR+"/bins/{smp}.count.bin",
    qgc = OUTDIR+"/index/ref.fasta.qgc",
    qm = OUTDIR+"/index/ref.fasta.qm",
    bed = OUTDIR+"/index/ref.fasta.bed"
  output:
    OUTDIR+"/CN/{smp}.CN.bed"
  params:
    ref = OUTDIR+"/index/ref.fasta",
    prefix = OUTDIR+"/bins/{smp}.count"
  shell:
    """
    quicKmer2 est {params.ref} {params.prefix} {output}
    """

# VISUALIZATION

rule chromsizes:
  input:
    OUTDIR+"/index/ref.fasta"
  output:
    OUTDIR+"/index/ref.fasta.chromsizes"
  shell:
   """
   scripts/fasta_len.awk {input} > {output}
   """

rule bed2bigBed:
  input:
    bedGraph = OUTDIR+"/CN/{smp}.CN.bed",
    chromsizes = OUTDIR+"/index/ref.fasta.chromsizes"
  output:
    bed9 = temp(OUTDIR+"/bigBed/{smp}.CN.bed9"),
    sorted = temp(OUTDIR+"/bigBed/{smp}.CN.srt.bed9"),
    bigbed = OUTDIR+"/bigBed/{smp}.CN.bb"
  shell:
    '''
    python3 scripts/bedToBed9.py {input.bedGraph} {output.bed9}
    LC_COLLATE=C sort -k1,1 -k2,2n {output.bed9} > {output.sorted}
    bedToBigBed -type=bed9 {output.sorted} {input.chromsizes} {output.bigbed}
    '''

