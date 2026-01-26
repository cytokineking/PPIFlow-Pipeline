import os,sys
import math
import subprocess

def batch_submit_jobs(batch_size, input_dir, output_dir, shell_script, cache_dir):
    """
    批量提交作业，将 input_dir 中的 PDB 文件分批处理，并提交到 SLURM 作业管理器。
    如果 output_dir 中已存在对应 relax 开头的文件，则跳过该 PDB 文件。

    参数：
    batch_size (int): 每个批次处理的文件数
    input_dir (str): 包含目标 PDB 文件的目录
    output_dir (str): 输出目录
    shell_script (str): 提交到 SLURM 的 shell 脚本路径
    cache_dir (str): 存放 batch 文件夹的目录
    """
    # 获取所有PDB文件，并过滤掉已处理的
    pdb_files = sorted([
        f for f in os.listdir(input_dir) if f.endswith('.pdb')
        and not os.path.exists(os.path.join(output_dir, f"relax_{os.path.splitext(f)[0]}.pdb"))
    ])
    total_files = len(pdb_files)
    num_batches = math.ceil(total_files / batch_size)
    print("需要处理的文件个数：",total_files)
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)

    for batch_idx in range(num_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, total_files)
        batch_files = pdb_files[batch_start:batch_end]

        if not batch_files:
            continue  # 这个批次没有文件要处理

        # 创建当前批次的目录
        batch_dir = os.path.join(cache_dir, f"batch_{batch_idx}")
        os.makedirs(batch_dir, exist_ok=True)

        for pdb_file in batch_files:
            src = os.path.realpath(os.path.join(input_dir, pdb_file))  # 获取真实路径
            dest = os.path.join(batch_dir, pdb_file)
            if not os.path.exists(dest):
                os.symlink(src, dest)


        # 构建 sbatch 命令并提交
        cmd = f"sbatch {shell_script} {batch_dir} {output_dir} {batch_idx}"
        status, jobnum = subprocess.getstatusoutput(cmd)

        # 打印提交信息
        print(cmd)
        print(status, jobnum)


batch_size = 1

input_dir = sys.argv[1]
cache_dir  = f"{input_dir}/cache/rlx-comp"
output_dir = f"{input_dir}/rst_itf_nofix"

shell_script = "demo_scripts/submit_relax_complex.sh"
batch_submit_jobs(batch_size, input_dir, output_dir, shell_script, cache_dir)