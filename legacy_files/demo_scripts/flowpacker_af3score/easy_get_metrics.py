#!/usr/bin/env python3
import os
import orjson
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed


def process_single_summary(summary_path, description):
    if not summary_path.exists():
        return {"description": description, "iptm": None, "ptm_A": None, "status": "missing"}

    try:
        with open(summary_path, "rb") as f:
            summary = orjson.loads(f.read())

        iptm = summary.get("iptm", None)
        chain_ptm = summary.get("chain_ptm", [])
        if not isinstance(chain_ptm, list):
            chain_ptm = []

        return {
            "description": description,
            "iptm": iptm,
            "ptm_A": chain_ptm[0] if chain_ptm else None,
            "status": "ok"
        }
    except Exception as e:
        return {"description": description, "iptm": None, "ptm_A": None, "status": f"error: {e}"}


def extract_summary_metrics(base_dir, sample_dir="seed-10_sample-0",
                            num_workers=64, batch_size=10000,
                            save_csv="results.csv", resume=False):
    descriptions = [d for d in os.listdir(base_dir) if (Path(base_dir) / d).is_dir()]
    paths = [(Path(base_dir) / d / sample_dir / "summary_confidences.json", d) for d in descriptions]

    processed = set()
    if resume and Path(save_csv).exists():
        try:
            df_existing = pd.read_csv(save_csv, usecols=["description"])
            processed = set(df_existing["description"].astype(str))
            print(f"🔄 检测到已有 {len(processed)} 条记录，将跳过这些任务")
        except Exception:
            print("⚠️ 已存在的 CSV 读取失败，无法断点续跑，将重新开始")
            os.remove(save_csv)

    elif not resume and Path(save_csv).exists():
        print("⚠️ 未启用 --resume，旧结果将被覆盖")
        os.remove(save_csv)

    total = len(paths)
    for i in range(0, total, batch_size):
        batch = paths[i:i+batch_size]
        batch = [(p, d) for p, d in batch if d not in processed]
        if not batch:
            print(f"✅ 跳过批次 {i//batch_size+1}（已处理过）")
            continue

        print(f"\n🚀 开始处理批次 {i//batch_size+1}，包含 {len(batch)} 个文件...")

        results = []
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_single_summary, p, d): d for p, d in batch}
            for j, fut in enumerate(tqdm(as_completed(futures), total=len(batch), desc=f"Batch {i//batch_size+1}")):
                res = fut.result()
                results.append(res)
                if (j + 1) % 1000 == 0:
                    print(f"   🔹 已完成 {j+1}/{len(batch)} 个文件")

        df_batch = pd.DataFrame(results)
        df_batch.to_csv(save_csv, mode="a", header=(not Path(save_csv).exists()), index=False)

        print(f"✅ 批次 {i//batch_size+1} 处理完成，本批保存 {len(df_batch)} 条记录")

    print(f"\n🎉 全部完成，共 {total} 项任务，结果保存到 {save_csv}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="批量提取 AF3 summary_confidences.json 的 iptm / ptm 信息")
    parser.add_argument("af3score_output_dir", help="AF3 结果的主目录")
    parser.add_argument("save_csv", help="输出 CSV 路径")
    parser.add_argument("--sample_dir", type=str, default="seed-10_sample-0", help="子目录名，例如 seed-10_sample-0")
    parser.add_argument("--num_workers", type=int, default=64, help="并行线程数")
    parser.add_argument("--batch_size", type=int, default=5000, help="每批次处理的文件数量")
    parser.add_argument("--resume", action="store_true", help="是否启用断点续跑")
    args = parser.parse_args()

    extract_summary_metrics(
        args.af3score_output_dir,
        sample_dir=args.sample_dir,
        num_workers=args.num_workers,
        batch_size=args.batch_size,
        save_csv=args.save_csv,
        resume=args.resume
    )
