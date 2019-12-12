import os
import pandas as pd

configfile: "config_rnaseq.yaml"

base_dir = config["BASE_DIR"]
anno_dir = config["ANNO_DIR"]
info = pd.read_csv("deseq_raw_vs_ctrl_noours.csv")
compares = list(info['compare'])

target_info = pd.read_csv("deseq_raw_vs_ctrl_sort.csv")
ntof_dict = target_info.set_index("id").to_dict()["compare"]

overlap_info = pd.read_csv("deseq2_overlap_data.csv")
overlap_dict = overlap_info.set_index("target").to_dict()["list"]

gene_info = pd.read_csv("ce11_genename.csv")
gene_names = list(target_info["overlap"])

os.chdir(base_dir)

rule all:
    input:
        expand("result/degene_overlap/{target}", target=overlap_dict.keys())

rule diff_res_table:
    input:
    output:
        table = directory("result/degene_overlap/{target}"),
    shell: 
        target = overlap_dict[wildcards.target]
        target_times = target.split(' ')
        up_table = output.table + "/%s_upgene_stat.tsv"%(target_times[0])
        down_table = output.table + "/%s_downgene_stat.tsv"%(target_times[0])
        with open(up_table, 'a') as fo:
            fo.write('\t'.join(["gene_name"] + target_times + ["count"]) + '\n')
            for gene_name in gene_names[:100]:
                count = 0
                #upgene
                fo.write(gene_name + '\t')
                for target_time in target_times:
                    file_name = ntof_dict[target_time]
                    uplist = "result/degenes_noours/%s/upgenelist.txt"%(file_name)
                    downlist = "result/degenes_noours/%s/downgenelist.txt"%(file_name)
                    with open(uplist, 'r') as fi:
                        upgenes = fi.read().splitlines()
                        if gene_name in upgenes:
                            fo.write('1' + '\t')
                            count += 1
                        else:
                            fo.write('0' + '\t')
                fo.write(str(count) + '\n')