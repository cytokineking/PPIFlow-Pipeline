import os, shutil, subprocess
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import PDB
import time
import glob
import re
import argparse
import ast


def get_score_map_auto(rese_file, inter_res_a, inter_res_b, chain_id_1, chain_id_2):
    # 创建映射字典，根据链ID选择对应的数字
    chain_map = {"L": "76", "R": "82"}
    # 读取 Rosetta 输出文件
    temp = pd.read_csv(rese_file, delim_whitespace=True)

    # 对 Res1 进行分割并提取链和残基编号
    split1 = temp['Res1'].apply(lambda x: pd.Series(x.split('_')))
    # print(split1)
    split1['part_a_chain'], split1['part_a_idx'] = split1[1].str[:2], split1[1].str[2:]
    # split1['part_a_chain'].replace(["66", "67"], ["B", "C"], inplace=True)
    split1['part_a_chain'].replace(chain_map, inplace=True)

    # 对 Res2 进行分割并提取链和残基编号
    split2 = temp['Res2'].apply(lambda x: pd.Series(x.split('_')))
    split2['part_b_chain'], split2['part_b_idx'] = split2[1].str[:2], split2[1].str[2:]
    # split2['part_b_chain'].replace(["66", "67"], ["B", "C"], inplace=True)
    split2['part_b_chain'].replace(chain_map, inplace=True)

    # 将结果合并到 temp1
    temp1 = pd.concat([temp, split1['part_a_chain'], split1['part_a_idx'], split2['part_b_chain'], split2['part_b_idx']], axis=1)
    # print(temp1)

    # 创建一个初始值为零的矩阵 score_map
    score_map = np.zeros((len(inter_res_a), len(inter_res_b)))

    # 遍历所有行并填充 score_map
    for idx, row in temp1.iterrows():
        try:
            res1 = int(row['part_a_idx'])
            res2 = int(row['part_b_idx'])
        except:
            print('weights row are NaN')
            continue  # 跳过异常行

        # 检查残基编号是否在界面残基列表中
        # if (row['part_a_chain'] == 'B') and (row['part_b_chain'] == 'C'):
        if (row['part_a_chain'] == chain_map.get(chain_id_1)) and (row['part_b_chain'] == chain_map.get(chain_id_2)):
            idx_a = np.where(np.array(inter_res_a) == res1)[0]
            idx_b = np.where(np.array(inter_res_b) == res2)[0]

            # 检查是否找到匹配的索引并填入 score_map
            if len(idx_a) > 0 and len(idx_b) > 0:
                try:
                    score_map[idx_a[0], idx_b[0]] = float(row['total'])  # 确保总能量为浮点数
                except ValueError:
                    print(f"Invalid total value at row {idx}")
            else:
                print(f"Residue {res1} or {res2} not found in inter_res_a or inter_res_b")

    return score_map, temp1


def plot_score_map(plot_path,
                    score_map, 
                   inter_res_a, 
                   inter_res_b, 
                   chain_a, 
                   chain_b,
                   title,
                   x_fontsize=8,
                   y_fontsize=8,
                   title_fontsize=10,
                   vmin=-5,
                   vmax=1):
    # 去除得分全为零的行和列
    non_zero_rows = np.any(score_map != 0, axis=1)
    non_zero_cols = np.any(score_map != 0, axis=0)
    filtered_score_map = score_map[non_zero_rows][:, non_zero_cols]
    filtered_inter_res_a = np.array(inter_res_a)[non_zero_rows]
    filtered_inter_res_b = np.array(inter_res_b)[non_zero_cols]
    # 不去除得分全为零的行和列
    # filtered_score_map = score_map
    # filtered_inter_res_a = np.array(inter_res_a)
    # filtered_inter_res_b = np.array(inter_res_b)

    row_labels = [str(x) for x in filtered_inter_res_a]
    column_labels = [str(x) for x in filtered_inter_res_b]
    # 只显示大于0的数值，否则显示为空白
    annot = np.where(filtered_score_map < 0, filtered_score_map, 0).T
    
    # fig, ax = plt.subplots(figsize=(len(filtered_inter_res_a)//4, len(filtered_inter_res_b)//4))
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(filtered_score_map.T, annot=False, linewidths=.5, xticklabels=row_labels, yticklabels=column_labels,
    vmin=vmin, vmax=vmax,  cbar=True)
    # sns.heatmap(filtered_score_map.T, annot=False, linewidths=.5, xticklabels=row_labels, yticklabels=column_labels,
    #  cbar=True)

    # 调整标题和轴标签的字体大小
    ax.set_xlabel(f'chain {chain_a}', fontsize=x_fontsize)
    ax.set_ylabel(f'chain {chain_b}', fontsize=y_fontsize)
    ax.set_title(f'{title}', fontsize=title_fontsize)

    ax.invert_yaxis()
    fig.savefig(f'{plot_path}/{title}-chain_{chain_a}_chain_{chain_b}.png', dpi=200)
    # plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot residue-residue interface energy heatmap.")
    parser.add_argument("--log_path", required=True, help="Path to the Rosetta log file")
    parser.add_argument("--pdb_name", required=True, help="PDB name (used for title)")
    parser.add_argument("--plot_path", required=True, help="Output path to save the plot")
    parser.add_argument("--chain1", required=True, help="Chain ID 1")
    parser.add_argument("--chain2", required=True, help="Chain ID 2")
    parser.add_argument("--residues1", required=True, help="Comma-separated list of residue IDs for chain 1")
    parser.add_argument("--residues2", required=True, help="Comma-separated list of residue IDs for chain 2")
    args = parser.parse_args()

    log_path = args.log_path
    pdb_name = args.pdb_name
    plot_path = args.plot_path
    chain_id_1 = args.chain1
    chain_id_2 = args.chain2
    inter_res_a = ast.literal_eval(args.residues1)
    inter_res_b = ast.literal_eval(args.residues2)
    # inter_res_a = [24,29,30,31,32,33,34,35,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,93]
    # inter_res_b = [5,6,7,8,9,10,11,12,13,14,15,17,18,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,74]


    score_map, temp = get_score_map_auto(log_path, inter_res_a, inter_res_b, chain_id_1, chain_id_2)

    # plot
    plot_score_map(plot_path, score_map, inter_res_a, inter_res_b, chain_id_1, chain_id_2, title=pdb_name, 
                    x_fontsize=20, y_fontsize=20, title_fontsize=24)