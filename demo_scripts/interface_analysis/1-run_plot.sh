#!/bin/bash

# 设置项目和路径
project='20250324_20_test'
input_dir='pdbs1'
output_dir='jobs'

# 切换到输出目录
cd $project


# 获取上层目录中 .pdb 文件列表
pdb_files=(../$input_dir/*.pdb)

# 遍历所有 .pdb 文件
for file in "${pdb_files[@]}"; do
    folder_name=$(basename "$file" .pdb)

    # 进入子文件夹
    cd "$output_dir/$folder_name"

    # 打印当前所在目录
    echo "📌 当前目录: $(pwd)"

    out_dir="out"
    out_path="$out_dir/${folder_name}.out"
    log_path="$out_dir/${folder_name}.log"
    plot_dir="plots"
    

    echo "📄 PDB Path: $file"
    echo "📁 Output Path: $out_path"
    echo "🔤 PDB Name: $folder_name"
    echo "📝 Log Path: $log_path"
    echo "plot path: $plot_dir"

    # extract res energy from out file
    grep "ResResE" "$out_path" > "$log_path"
    echo "✅ 已将匹配行写入: $log_path"

    # # *****1.计算interface residue******
    OUTPUT=$(python /lustre/grp/cmclab/wuxy/process/interface_analysis/codes/get_interface_residues.py -i "${folder_name}.pdb" -c1 L -c2 R -d 15.0)

    # 解析输出
    IFS=' ' read -r -a CHAIN_A_RES <<< "$(echo "$OUTPUT" | grep "Interface residues in chain L" | awk -F: '{print $2}' )"
    IFS=' ' read -r -a CHAIN_B_RES <<< "$(echo "$OUTPUT" | grep "Interface residues in chain R" | awk -F: '{print $2}' )"


    # 格式化输出
    CHAIN_A_RES_LIST="$(echo "${CHAIN_A_RES[@]}")"
    CHAIN_B_RES_LIST="$(echo "${CHAIN_B_RES[@]}")"

    echo $CHAIN_A_RES_LIST
    echo $CHAIN_B_RES_LIST
    # *****计算interface residue******


    # 绘制能量的heatmap图

    # sbatch /lustre/grp/cmclab/wuxy/process/interface_analysis/codes/run_plot.sh \
    # "$log_path" "$folder_name" "$plot_dir" "L" "R" "$CHAIN_A_RES_LIST" "$CHAIN_B_RES_LIST"
    # cmd="sbatch /lustre/grp/cmclab/wuxy/process/interface_analysis/codes/run_plot.sh \
    # \"$log_path\" \"$folder_name\" \"$plot_dir\" \"L\" \"R\" \"$CHAIN_A_RES_LIST\" \"$CHAIN_B_RES_LIST\""
    python /lustre/grp/cmclab/share/wuxy/process/interface_analysis/codes/plot_heatmap.py --log_path "$log_path" --pdb_name "$folder_name" --plot_path "$plot_dir" --chain1 'L' --chain2 'R' --residues1 "$CHAIN_A_RES_LIST" --residues2 "$CHAIN_B_RES_LIST"

    # echo "$cmd"
    # eval $cmd



    
    # 回到项目根目录
    cd ../../
done
