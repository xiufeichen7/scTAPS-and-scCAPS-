"""
Workflow for single-cell caps

samples include:
    mESC

read length: 128 bp

rawdata: /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/P230072

run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake --use-envmodules --max-status-checks-per-second 0.01 --snakefile code/mESC_sccaps.smk --cluster "sbatch -p short " -j 96 -np

dir=/gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/P230072/WTCHG_981766
paste <(ls ${dir}/*gz) <(ls ${dir}/*gz|sed 's/WTCHG_981766/WTCHG_981767/g') <(ls ${dir}/*gz|awk -F'/' '{print $NF}'|sed 's,^,>fastq/,g;s/[ATCG+]*_R//g')|sed 's/^/cat /g;' >P230072_rename_fastq.sh
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["SCCAPS_SMPS"] #"C203_N702_i7_10_rev_N502_i5_1_rev" # 
REF = config["SCTAPS_SMPS_REF"]
print(SAMPLES)

rule all:
    input:
        expand("fastq/{sample}_{readDirection}_fastqc.html",sample = SAMPLES, readDirection=['1','2']),
        expand("align/{sample}.bwa.bam", sample = SAMPLES),
        expand("align/{sample}.md.bam", sample = SAMPLES),
        expand("stats/{sample}.mapping.txt", sample = SAMPLES),
        expand("align/{sample}.md.read.meth.txt.gz", sample=SAMPLES),
        expand("align/{sample}.spikeins.bam", sample = SAMPLES),
        expand("align/{sample}.spikeins.read.meth.txt.gz", sample = SAMPLES),
        expand("align/{sample}.md.filter.meth.sta.txt.gz", sample = SAMPLES),
        expand("align/{sample}.spikeins.filter.meth.sta.txt.gz", sample = SAMPLES)



################################ PREPROCESS #################################
rule mergefastq: 
    input:
        "code/P230072_rename_fastq.sh"
    output:
        expand("fastq/{{sample}}_{readDirection}.fastq.gz", readDirection=['1','2'])
    params:
        "{sample}"
    shell:
        """
        cmdline=$(grep {params} {input}|tr '\\n' ';')
        echo $cmdline|bash -
        """

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

rule extract_meth: 
    input:
        bam="align/{sample}.md.bam"
    output:
        "align/{sample}.md.read.meth.txt.gz"
    params:
        python=config["python"]
    shell:
        """
        module purge
        module load {params.python}
        python3 code/nd_taps_extract_parallel.py -b {input.bam} -t 20 -n 100000
        """

rule meth_trim: 
    input:
        "align/{sample}.md.read.meth.txt.gz"
    output:
        "align/{sample}.md.filter.meth.sta.txt.gz"
    params:
        python=config["python"]
    shell:
        """
        module purge
        module load {params.python}
        python3 code/nd_taps_convert_call.py -r {input} -s 10 -e 118 -c True
        """

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

rule spikein_trim: 
    input:
        "align/{sample}.spikeins.read.meth.txt.gz"
    output:
        "align/{sample}.spikeins.filter.meth.sta.txt.gz"
    params:
        python=config["python"]
    shell:
        """
        module purge
        module load {params.python}
        python3 code/nd_taps_convert_call.py -r {input} -s 10 -e 118 -c True
        """

