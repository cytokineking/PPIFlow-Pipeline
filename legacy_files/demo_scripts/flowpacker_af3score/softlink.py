import os
import glob
import pandas as pd
from multiprocessing import Pool, cpu_count

softlink_dir = "/lustre/grp/cmclab/share/wuxy/process/flowpacker_af3score/data/qilin/softlink/"
pdb_dir = "/lustre/grp/cmclab/share/yuql/design/IFNR2/partial_flow_motif/result"
pdb_paths = glob.glob(rf"{pdb_dir}/**/*.pdb", recursive=True)

os.makedirs(softlink_dir, exist_ok=True)

# ========== 创建任务函数 ==========
def create_symlink(pdb_file):
    try:
        pdb_name = "_".join(pdb_file.split("/")[-2].split("_")[1:-1]) + "_" + os.path.basename(pdb_file)[6] + ".pdb"
        link_name = os.path.join(softlink_dir, pdb_name)

        # 若已有旧链接，先删除
        if os.path.exists(link_name):
            os.remove(link_name)

        os.symlink(pdb_file, link_name)

        return [pdb_file, link_name]

    except Exception as e:
        print(f"❌ 错误: {pdb_file} → {e}")
        return None

# ========== 并行执行 ==========
if __name__ == "__main__":
    num_workers = max(1, cpu_count() - 2)  # 留出 CPU 给系统
    print(f"🚀 正在使用 {num_workers} 个进程并行创建 symlinks...")

    with Pool(processes=num_workers) as pool:
        results = pool.map(create_symlink, pdb_paths)

    # 过滤掉 None 的结果
    link_correspondence = [r for r in results if r is not None]

    # 写入 CSV
    df = pd.DataFrame(link_correspondence, columns=["pdb_file", "link_name"])
    df.to_csv(os.path.join(softlink_dir, "link_correspondence.csv"), index=False)

    print(f"✅ 创建完毕，共 {len(link_correspondence)} 个软链接")
