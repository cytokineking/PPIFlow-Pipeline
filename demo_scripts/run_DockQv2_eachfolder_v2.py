import os
import argparse

def generate_shell_script(subfolder, reference_dir, output_dir, shell_dir):
    subfolder_name = os.path.basename(subfolder)

    # 遍历参考目录
    reference_files = os.listdir(reference_dir)
    reference_pdb_path = None

    for f in reference_files:
        if f.lower() == f"{subfolder_name.lower()}.pdb":
            # 直接匹配到 pdb 文件（可以是 link）
            reference_pdb_path = os.path.join(reference_dir, f)
            break
        elif f.lower() == subfolder_name.lower():
            # 匹配到文件夹或符号链接目录
            subdir = os.path.join(reference_dir, f)
            pdb_candidates = [
                os.path.join(subdir, x)
                for x in os.listdir(subdir)
                if x.endswith(".pdb")
            ]
            if pdb_candidates:
                reference_pdb_path = pdb_candidates[0]  # 取第一个 pdb
                break

    if reference_pdb_path is None:
        print(f"❌ Reference PDB not found for {subfolder_name}")
        return None
    else:
        print(f"✅ Found reference: {reference_pdb_path}")

    # 输出目录
    sub_output_dir = os.path.join(output_dir, subfolder_name)
    os.makedirs(sub_output_dir, exist_ok=True)

    shell_script_path = os.path.join(shell_dir, f"run_dockq_{subfolder_name}.sh")

    # 写入 shell 脚本
    script_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=dockq_{subfolder_name}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=1",
        "#SBATCH --time=12:00:00",
        "#SBATCH -p cpu_short,cpu1,cpu2,fat,hygon,dcu",
        f"#SBATCH --output={sub_output_dir}/dockq_{subfolder_name}_%j.log",
        "",
        "set -e",  # 出现错误立即停止
        "source /lustre/grp/cmclab/qinxy/new_miniconda3/etc/profile.d/conda.sh",
        "conda activate py39",
    ]

    # 为每个 PDB 文件生成 DockQ 命令
    for pdb_file in os.listdir(subfolder):
        if not pdb_file.endswith(".pdb"):
            continue
        model_pdb_path = os.path.join(subfolder, pdb_file)
        output_file = os.path.join(
            sub_output_dir, f"{os.path.splitext(pdb_file)[0]}_dockq_score"
        )
        # 如果输出已存在，直接跳过
        if os.path.exists(output_file):
            continue
        cmd = f"DockQ --allowed_mismatches 10 {model_pdb_path} {reference_pdb_path} --short > {output_file}"
        script_lines.append(cmd)

    # 写入 shell 文件
    with open(shell_script_path, "w") as f:
        f.write("\n".join(script_lines) + "\n")

    os.chmod(shell_script_path, 0o755)
    return shell_script_path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_dir',
        default='/lustre/grp/cmclab/qinxy/antibody_ppiflow/iptm_dockq_cor/nanobody/af3/pdb',
        help='Directory containing subfolders with PDB files'
    )
    parser.add_argument(
        '--reference_dir',
        default='/lustre/grp/cmclab/qinxy/antibody_ppiflow/iptm_dockq_cor/nanobody/pdb_renamed_chain',
        help='Directory containing reference PDB files'
    )
    parser.add_argument(
        '--output_dir',
        default='/lustre/grp/cmclab/qinxy/antibody_ppiflow/iptm_dockq_cor/nanobody/dockq2',
        help='Directory to save DockQ results'
    )
    parser.add_argument(
        '--shell_dir',
        default='/lustre/grp/cmclab/qinxy/antibody_ppiflow/iptm_dockq_cor/nanobody/shell2',
        help='Directory to save generated shell scripts'
    )

    args = parser.parse_args()

    os.makedirs(args.shell_dir, exist_ok=True)

    generated_scripts = []
    for subfolder in os.listdir(args.input_dir):
        subfolder_path = os.path.join(args.input_dir, subfolder)
        if not os.path.isdir(subfolder_path):
            continue
        script_path = generate_shell_script(subfolder_path, args.reference_dir, args.output_dir, args.shell_dir)
        if script_path:
            generated_scripts.append(script_path)

    # 生成统一提交脚本
    submit_all_path = os.path.join(args.shell_dir, "submit_all_jobs.sh")
    with open(submit_all_path, "w") as f:
        f.write("#!/bin/bash\n\n")
        for script in generated_scripts:
            f.write(f"sbatch {script}\n")
    os.chmod(submit_all_path, 0o755)
    print(f"Generated {len(generated_scripts)} shell scripts in {args.shell_dir}")

if __name__ == "__main__":
    main()