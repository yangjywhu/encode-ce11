import os
import pandas as pd

configfile: "config_rnaseq.yaml"

base_dir = config["BASE_DIR"]
anno_dir = config["ANNO_DIR"]
info = pd.read_csv("deseq_raw_vs_ctrl_noours.csv")
compares = list(info['compare'])

os.chdir(base_dir)

rule all:
    input:
        # expand("result/degenes_0.05_0.585/{compare}", compare=compares),
        "/data4/zhoulab/yangjiayi/project/yxh_ce/result/degene_0.05_0.585_table",

rule deseq2_raw_vs_ctrl:
    input:
        c1 = "result_ctrl/featurecounts/{c1}/{c1}_featureCounts.txt",
        c2 = "result_ctrl/featurecounts/{c2}/{c2}_featureCounts.txt",
        t1 = "result_treat/featurecounts/{t1}/{t1}_featureCounts.txt",
        t2 = "result_treat/featurecounts/{t2}/{t2}_featureCounts.txt",
    output:
        degene = directory("result/degenes_0.05_0.585/{c1}-{c2}_vs_{t1}_{t2}"),
    params:
        fc = 0.585,
        padj = 0.05,
    shell:"""
        set +u; source ~/miniconda3/bin/activate deseq2; set -u
        Rscript script/deseq2_data.R \
            {input.c1} {input.c2} {input.t1} {input.t2} \
            {output.degene} {params.fc} {params.padj}
        set +u; conda deactivate; set -u
        """

rule get_overlap:
    input:
        base_dir = "/data4/zhoulab/yangjiayi/project/yxh_ce/script",
    output:
        out_dir = directory("/data4/zhoulab/yangjiayi/project/yxh_ce/result/degene_0.05_0.585_table"),
    shell:"""
        mkdir -p {output.out_dir}
        python script/get_deseq2_overlap_upanddown.py \
            {input.base_dir} {output.out_dir}
        """