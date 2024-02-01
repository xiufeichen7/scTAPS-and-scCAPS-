"""
Workflow for lowinput taps/caps (tagmentation)

samples include:
    100_cells_taps
    10_cells_taps
    10k_cells_taps
    1k_cells_taps

    100_cells_caps
    10_cells_caps
    10k_cells_caps
    1k_cells_caps

taps: read length 110 /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS
caps: read length 151 /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_CAPS

run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake --use-envmodules --max-status-checks-per-second 0.01  --snakefile code/lowinput_analysis.smk --cluster "sbatch -p long --cpus-per-task 3 " -j 8 -np

# split bam -> bam convert -> convert merge -> merge read -> convert call
for i in `ls 1*filter.meth.sta.txt.gz`;do zcat $i |grep 144hmC|
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["TAG_SMPS"]
REF = config["TAG_SMPS_REF"]
print(SAMPLES)

rule all:
    input:
        expand("stats/{sample}.align.txt", sample = SAMPLES),
        # expand("align/{sample}.bwa.bam", sample = SAMPLES),
        # expand("align/{sample}.md.bam", sample = SAMPLES),
        # expand("stats/{sample}.depth.txt", sample=SAMPLES),
        # expand("stats/{sample}.mapping.txt", sample = SAMPLES),
        # expand("{sample}.call_temp/{sample}.1.bam", sample=SAMPLES)
        # expand("align/{sample}.spikeins.bam", sample = SAMPLES),
        # expand("align/{sample}.spikeins.read.meth.txt.gz", sample = SAMPLES)
        # expand("align/{sample}.filtered.read.meth.txt.gz", sample=SAMPLES),



################################ PREPROCESS #################################

rule fastp: 
    input:
        expand("fastq/{{sample}}_{readDirection}.fastq.gz",readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    params:
        prefix="fastq/{sample}.fastp",
        fastp=config["fastp"]
    shell:
        """
        {params.fastp} \
            -i {input[0]} \
            -I {input[1]} \
            -o {output[0]} \
            -O {output[1]} \
            -j {params.prefix}.json \
            -h {params.prefix}.html 
        """

rule fastqc: 
    input:
        raw=expand("fastq/{{sample}}_{readDirection}.fastq.gz",readDirection=['1','2']),
        clp=expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}_{readDirection}_fastqc.html",readDirection=['1','2'])
    params:
        fastqc=config["fastqc"]
    shell:
        """
        module purge
        module load {params.fastqc}
        fastqc {input.raw[0]}
        fastqc {input.raw[1]}
        fastqc {input.clp[0]}
        fastqc {input.clp[1]}
        """

rule align_fastp: 
    input:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    output:
        "align/{sample}.bwa.bam"
    log:
        "logs/{sample}.bwa.log"
    params:
        ref=REF,
        tmp="{sample}",
        bwa=config["bwa"],
        samtools=config["samtools"]
    threads: 2
    shell:
        """
        module purge
        module load {params.bwa}
        module load {params.samtools}
        (bwa mem -t {threads} {params.ref} {input} |\
        samtools sort -@ 8 -O BAM -T {params.tmp} >{output}) 1>{log} 2>&1
        """
################################ MM9 #################################
rule markdup: 
    input:
        "align/{sample}.bwa.bam"
    output:
        mdbam="align/{sample}.md.bam",
        matrix="align/{sample}.md_report.txt"
    log:
        "logs/{sample}.picard.log"
    params:
        mdtmp=temp("markdup"),
        picard=config["picard"]
    shell:
        """
        module purge
        module load {params.picard}
        (
            java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar MarkDuplicates \
            I={input} \
            O={output.mdbam} \
            M={output.matrix} \
            TMP_DIR={params.mdtmp}) 1>{log} 2>&1
        """

rule genome_depth:
    input:
        bam="align/{sample}.md.bam"
    output:
        "stats/{sample}.depth.txt"
    params:
        sample="{sample}",
        samtools=config["samtools"]
    shell:
        """
        module purge
        module load {params.samtools}
        depth=`samtools depth -a {input}|awk '{{sum+=$3}}END{{print sum/NR}}'`
        echo -e "{params.sample}\\t$depth" >{output}
        """

rule mapping: 
    input:
        bam="align/{sample}.md.bam",
        dup="align/{sample}.md_report.txt"
    output:
        "stats/{sample}.mapping.txt"
    params:
        sample="{sample}",
        samtools=config["samtools"]
    shell:
        """
        module purge
        module load {params.samtools}
        mkdir {params.sample}_temp
        q1_nmap=`samtools view -q 1 {input.bam}|cut -f1|sort -u -T {params.sample}_temp |wc -l `
        q10_nmap=`samtools view -q 10 {input.bam}|cut -f1|sort -u -T {params.sample}_temp|wc -l `
        dup_rate=`grep ^LI {input.dup} -A1|tail -n +2|cut -f9`
        proper_nmap=`samtools view -q 10 {input.bam} |awk '$2==99||$2==163'|awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{sum+=$9}}END{{print NR,sum/NR}}'`
        echo -e "sample\\tq1_nmap\\tq10_nmap\\tproper_nmap\\tmean_isize\\tdup_rate" >{output}
        echo -e "{params.sample}\\t$q1_nmap\\t$q10_nmap\\t$proper_nmap\\t$dup_rate" >>{output}
        rm {params.sample}_temp -rf
        """

rule align_sta: 
    input:
        mapping="stats/{sample}.mapping.txt",
        rawfqc="fastq/{sample}_1_fastqc.html",
        cleanfqc="fastq/{sample}_fastp_1_fastqc.html"
    output:
        "stats/{sample}.align.txt"
    params:
        sample="{sample}",
    shell:
        """
        nraw_reads=`sed 's/<[^>]*>/\\n/g' {input.rawfqc} |grep -i total -A3|sed -n 3p`
        nclean_reads=`sed 's/<[^>]*>/\\n/g' {input.cleanfqc} |grep -i total -A3|sed -n 3p`
        echo -e "sample\\traw_reads\\tclean_reads\\tq1_nmap\\tq10_nmap\\tproper_nmap\\tmean_isize\\tdup_rate" >{output}
        echo -e "{params.sample}\\t$nraw_reads\\t$nclean_reads\\t`tail -n +2 {input.mapping}|cut -f2-`" >>{output}
        """

######################## with split bam into smaller pieces when bam are too big ########################
# rule split_bam: 
#     input:
#         bam="align/{sample}.md.bam"
#     output:
#         "{sample}.call_temp/{sample}.1.bam"
#     log: "logs/{sample}.split_bam.log"
#     shell:
#         """
#         (
#             module purge
#             module load Python/3.7.4-GCCcore-8.3.0
#             python3 code/nd_taps_split_bam.py -b {input.bam} -n 100000
#         ) 1>{log} 2>&1
#         """
# use lowinput_analysis_parallel.smk to call meth for each splitted bam files
#rule merge_read: 
#    input:
#        bam="align/{sample}.md.bam"
#    output:
#        "align/{sample}.read.meth.txt.gz"
#    log: "logs/{sample}.merge_call.log"
#    shell:
#        """
#        (
#            module purge
#            module load Python/3.7.4-GCCcore-8.3.0
#            python3 code/nd_taps_merge_read.py -b {input.bam}
#        ) 1>{log} 2>&1
#        """
#
######################## spike-ins ########################
rule extract_spikein: 
    input:
        "align/{sample}.bwa.bam"
    output:
        "align/{sample}.spikeins.bam",
    params:
        ref=REF,
        samtools=config["samtools"]
    shell:
        """
        module purge
        module load {params.samtools}
        samtools view -bS {input} -L <(grep -v "chr" {params.ref}.fai|awk 'BEGIN{{OFS="\\t"}}{{print $1,"0",$2}}') >{output}
        """

rule spikein_calls: 
    input:
        "align/{sample}.spikeins.bam"
    output:
        "align/{sample}.spikeins.read.meth.txt.gz"
    params:
        ref=REF,
        python=config["python"]
    shell:
        """
        module purge
        module load {params.python}
        python3 code/nd_taps_extract_parallel.py -b {input} -t 5 -n 100000
        """
# for i in `ls align/1*caps*md.read.meth.txt.gz`;do python3 code/nd_taps_convert_call.py -r $i -s 10 -e 140 -c True;done
# for i in `ls align/1*caps*spikeins.read.meth.txt.gz`;do python3 code/nd_taps_convert_call.py -r $i -s 10 -e 140 -c True;done
# for i in `ls align/1*taps*md.read.meth.txt.gz`;do python3 code/nd_taps_convert_call.py -r $i -s 10 -e 100 -c True;done
# for i in `ls align/1*taps*spikeins.read.meth.txt.gz`;do python3 code/nd_taps_convert_call.py -r $i -s 10 -e 100 -c True;done
