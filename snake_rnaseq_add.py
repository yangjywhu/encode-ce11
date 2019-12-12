import os
import pandas as pd

configfile: "../config_rnaseq.yaml"

base_dir = config["BASE_DIR"]
anno_dir = config["ANNO_DIR"]
info1 = pd.read_csv("ce_data_info_add.csv")
samples = list(info1['id'])
read1s = list(info1['read1'])
read2s = list(info1['read2'])
matchs = list(info1['match'])
targets = list(info1['target'])
times = list(info1['time'])

# samples = ["ce1"]
# read1s = ["ENCFF759EVC"]
# read2s = ["ENCFF334LYW"]
# matchs = ["ENCFF759EVC_ENCFF334LYW"]

os.chdir(base_dir)

rule all:
    input:
        expand("result_treat/fastqc/{read}", read=read1s+read2s),
        expand("result_treat/mapping_rRNA/{match}/Unmapped.out.mate1", match=matchs),
        expand("result_treat/mapping_geno/{match}/Uniq_Sort.bam", match=matchs),
        expand("result_treat/featurecounts/get/{match}.txt", match=matchs),

rule fastqc:
    input:
        fq = "data/adddata/{read}.fastq.gz",
    output:
        qc = directory("result_treat/fastqc/{read}"),
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        mkdir -p {output.qc}
        fastqc -o {output.qc} {input.fq}
        set +u; conda deactivate; set -u
        """

rule mapping_rRNA:
    input:
        index_dir = anno_dir + "/ce11_star_rRNA",
        fq1 = "data/adddata/{read1}.fastq.gz",
        fq2 = "data/adddata/{read2}.fastq.gz",
    output:
        un1 = "result_treat/mapping_rRNA/{read1}_{read2}/Unmapped.out.mate1",
        un2 = "result_treat/mapping_rRNA/{read1}_{read2}/Unmapped.out.mate2",
    params:
        bam_dir = "result_treat/mapping_rRNA/{read1}_{read2}/",
    log: "result_treat/mapping_rRNA/{read1}_{read2}/mapping_rRNA.log",
    threads: 10
    shell:"""
        set +u; source ~/miniconda3/bin/activate star253; set -u
        STAR --genomeDir {input.index_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand gunzip -c \
            --runThreadN {threads} \
            --outFileNamePrefix {params.bam_dir} \
            --outReadsUnmapped Fastx \
            --outSAMtype None \
            --outFilterMatchNmin 40 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            > {log} 2>&1
        set +u; conda deactivate; set -u
        """

rule mapping_geno:
    input:
        index_dir = anno_dir + "/ce11_star_geno",
        fq1 = rules.mapping_rRNA.output.un1,
        fq2 = rules.mapping_rRNA.output.un2,
    output:
        bam = temp("result_treat/mapping_geno/{read1}_{read2}/Aligned.out.bam"),
        uniq = temp("result_treat/mapping_geno/{read1}_{read2}/Uniq.bam"),
        sort = "result_treat/mapping_geno/{read1}_{read2}/Uniq_Sort.bam",
        bai = "result_treat/mapping_geno/{read1}_{read2}/Uniq_Sort.bam.bai",
    params:
        bam_dir = "result_treat/mapping_geno/{read1}_{read2}/",
    log: "result_treat/mapping_geno/{read1}_{read2}/mapping_geno.log",
    threads: 5
    shell:"""
        set +u; source ~/miniconda3/bin/activate star253; set -u
        STAR --genomeDir {input.index_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.bam_dir} \
            --outSAMtype BAM Unsorted \
            --outFilterMatchNmin 40 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --limitBAMsortRAM 5000000000 \
            > {log} 2>&1
        samtools view -@ {threads} -h -bq 255 \
            {output.bam} > {output.uniq}
        samtools sort -@ {threads} {output.uniq} -o {output.sort}
        samtools index -@ {threads} {output.sort}
        set +u; conda deactivate; set -u
        """

rule featurecounts:
    input:
        bam = "result_treat/mapping_geno/{match}/Uniq_Sort.bam",
        gtf = anno_dir + "/ce11_anno/c_elegans.PRJNA13758.WS258.canonical_geneset.gtf",
    output:
        fc_res = "result_treat/featurecounts/get/{match}.txt",
    params:
        ss = 0,
        id = lambda wildcards: samples[matchs.index(wildcards.match)],
    threads: 1
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        mkdir -p result_treat/featurecounts/{params.id}
        featureCounts -a {input.gtf} -s {params.ss} -T {threads} -p \
            -o result_treat/featurecounts/{params.id}/{params.id}_featureCounts.txt \
            {input.bam} \
            > result_treat/featurecounts/{params.id}/{params.id}_featureCounts.log 2>&1
        echo "{params.id}:{wildcards.match} finished" > {output.fc_res}
        set +u; conda deactivate; set -u
        """
