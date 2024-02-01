# ln -s /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS/C4_S2_L001_R1_001.fastq.gz cd4_tcells_R1.fastq.gz
# ln -s /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS/C4_S2_L001_R2_001.fastq.gz cd4_tcells_R2.fastq.gz
# ln -s /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS/C5_S3_L001_R1_001.fastq.gz cd8_tcells_R1.fastq.gz
# ln -s /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS/C5_S3_L001_R2_001.fastq.gz cd8_tcells_R2.fastq.gz
# ln -s /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS/C8_S4_L001_R1_001.fastq.gz bcells_R1.fastq.gz
# ln -s /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS/C8_S4_L001_R2_001.fastq.gz bcells_R2.fastq.gz
"""
Workflow for standard taps
samples include:
bcells
cd8_tcells
cd4_tcells


taps: read length 110 /gpfs3/well/ludwig/users/cfo155/scTAPS_CAPS/mESC/rawdata/mESC_TAPS
samtools faidx ../../mESC/resource/caps_mm9_lambda.fa -r <(grep -v chr ../../mESC/resource/caps_mm9_lambda.fa.fai|cut -f1) |\
cat - /gpfs3/well/ludwig/users/cfo155/cfDNA/012020_cfDNA/resource/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  >GRCh38_spike_ins.fasta

run:
module purge
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake --use-envmodules --max-status-checks-per-second 0.01 --snakefile code/standard_taps.smk --cluster "sbatch -p long " -j 3 -np

"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["STD_SMPS"]
REF = config["STD_SMPS_REF"]
print(SAMPLES)

rule all:
    input:
        expand("fastq/{sample}_{readDirection}_fastqc.html",sample = SAMPLES, readDirection=['1','2']),
        expand("align/{sample}.md.bam", sample = SAMPLES),
        expand("stats/{sample}.depth.txt", sample=SAMPLES),
        expand("stats/{sample}.mapping.txt", sample = SAMPLES),
        expand("meth/{sample}.mbias.{seq}.txt", sample = SAMPLES, seq=['chr1','J02459.1']),
        expand("meth/{sample}.mbias.pdf", sample = SAMPLES),
        expand("meth/{sample}_merge_CpG.bedGraph.gz", sample=SAMPLES)



################################ PREPROCESS #################################

rule fastp: 
    input:
        expand("fastq/{{sample}}_{readDirection}.fastq.gz",readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    params:
        prefix="fastq/{sample}.fastp",
        fastp=config["fastp"]
    threads: 1
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
################################ SPIKE-INS #################################
rule extract_spikein: 
    input:
        "align/{sample}.bwa.bam"
    output:
        unmodified_2kb="align/{sample}.unmodified_2kb.bam",
        ncnn="align/{sample}.ncnn.bam"
    params:
        ref=REF,
        samtools=config["samtools"]
    shell:
        """
        module purge
        module load {params.samtools}
        samtools view -bS {input} -L <(grep unmodified_2kb {params.ref}.fai|awk 'BEGIN{{OFS="\\t"}}{{print $1,"0",$2}}') >{output.unmodified_2kb}
        samtools view -bS {input} -L <(grep 237mer {params.ref}.fai|awk 'BEGIN{{OFS="\\t"}}{{print $1,"0",$2}}') >{output.ncnn}
        """

################################ GENOME #################################
rule markdup: 
    input:
        "align/{sample}.bwa.bam"
    output:
        mdbam="align/{sample}.md.bam",
        matrix="align/{sample}.md.matrix.txt"
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
        dup="align/{sample}.md.matrix.txt"
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

rule meth_mbias:
    input:
        "align/{sample}.md.bam"
    output: 
        mbias=expand("meth/{{sample}}.mbias.{seq}.txt", seq=['chr1','J02459.1'])
    params:
        pars="-q 10 -p 13 -t 2 ",
        prefix="meth/{sample}",
        ref=REF,
        methydackel=config["methydackel"],
        matplotlib=config["matplotlib"]
    log: "meth/{sample}.suggested_mbias.txt" 
    shell: 
        """
        module purge
        module load {params.matplotlib}
        (
            {params.methydackel} mbias -r chr1 \
            {params.ref} {input} {params.prefix} --txt >{output.mbias[0]}
            {params.methydackel} mbias -r J02459.1 \
            {params.ref} {input} {params.prefix} --txt  >{output.mbias[1]}
        ) >{log} 2>&1
        """

rule mbias_plot:
    input:
        t1="meth/{sample}.mbias.chr1.txt",
        t2="meth/{sample}.mbias.J02459.1.txt"
    output: 
        pdf="meth/{sample}.mbias.pdf"
    params:
        r=config["r"],
        mbiasplot=config["mbiasplot"]
    shell: 
        """   
        module purge
        module load {params.r}
        Rscript {params.mbiasplot} -t1 {input.t1} -t2 {input.t2} -o {output.pdf}
        """

rule meth_call:
    input:
        "align/{sample}.md.bam"
    output: 
        "meth/{sample}_merge_CpG.bedGraph.gz"
    params:
        pars="-q 10 -p 13 -t 2 --mergeContext",
        trim="$(cat meth/{sample}.suggested_mbias.txt|head -1 |cut -d ':' -f2)",
        prefix="meth/{sample}_call",
        ref=REF,
        methydackel=config["methydackel"],
        matplotlib=config["matplotlib"]
    shell: 
        """
        module purge
        module load {params.matplotlib}
            {params.methydackel} extract -o {params.prefix} \
            {params.ref} {input} \
            {params.pars} \
            {params.trim}
        cat {params.prefix}_CpG.bedGraph | awk 'BEGIN{{OFS="\\t"}}{{if(NR>1)print $1,$2,$3,100-$4,$6,$5+$6}}'|\
            cat <(echo -e "#chr\\tstart\\tend\\tratio\\tmC\ttotal") -|gzip - >{output}
        #rm {params.prefix}_CpG.bedGraph -rf
        """

