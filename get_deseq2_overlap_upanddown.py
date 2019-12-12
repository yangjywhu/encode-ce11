import os
import sys
import pandas as pd

# base_dir = "D:/MyDoc/OneDrive/Doc/Lab/analysis/yxh_ce/metadata"
base_dir = sys.argv[1]
out_dir = sys.argv[2]
os.chdir(base_dir)

os.chdir
with open("ce11_genename.csv", 'r') as f:
    gene_names = f.read().splitlines()
with open("deseq2_overlap_data2.csv", 'r') as f:
    targets = f.read().splitlines()

# create dict, from id to file name
name_to_file = pd.read_csv("deseq_raw_vs_ctrl_sort.csv")
ntof_dict = name_to_file.set_index("id").to_dict()["compare"]

def list_to_res(out_table, diff):
    with open(out_table, 'a') as fo:
        fo.write('\t'.join(["gene_name"] + target_times + ["count"]) + '\n')
        for gene_name in gene_names:
            count = 0
            fo.write(gene_name + '\t')
            for target_time in target_times:
                file_name = ntof_dict[target_time]
                if diff == "up":
                    degene = "../result/degenes_0.05_0.585/%s/upgenelist.txt"%(file_name)
                else:
                    degene = "../result/degenes_0.05_0.585/%s/downgenelist.txt"%(file_name)
                with open(degene, 'r') as f:
                    upgenes = f.read().splitlines()
                    if gene_name in upgenes:
                        fo.write('1' + '\t')
                        count += 1
                    else:
                        fo.write('0' + '\t')
            fo.write(str(count) + '\n')

# targets1 = [targets[1]]
for target in targets:
    target_times = target.split(',')
    up_table = out_dir + "/%s_upgene_stat.tsv"%(target_times[0])
    down_table = out_dir + "/%s_downgene_stat.tsv"%(target_times[0])
    list_to_res(up_table, "up")
    list_to_res(down_table, "down")
