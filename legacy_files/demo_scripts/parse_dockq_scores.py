import os,sys
import csv
from pathlib import Path
from multiprocessing import Pool, cpu_count
import time

# ----------------------------------------------------------------------
# Worker Function: 提取单个文件夹中的 DockQ 分数
# ----------------------------------------------------------------------
def extract_dockq_scores(sub_folder_path):
    """
    遍历指定子文件夹下的所有 *\_dockq\_score 文件，提取所需的 DockQ 值。
    
    Args:
        sub_folder_path (str): 要扫描的子文件夹路径。

    Returns:
        list: 包含该子文件夹中所有提取结果的列表 (不包含表头)。
    """
    # 使用 Pathlib 处理路径，更健壮
    base_path = Path(sub_folder_path)
    # 仅存储数据行，不包含表头
    data_rows = []

    # 检查文件夹是否存在
    if not base_path.is_dir():
        # 如果路径无效，返回空列表，避免进程崩溃
        print(f"  [ERROR] 路径不是一个有效的文件夹，跳过: {sub_folder_path}")
        return []
        
    print(f"  [Worker {os.getpid()}] 开始处理文件夹: {base_path.name}")

    # 请确保文件名匹配您的实际文件后缀。这里保留了 *_dockq_score 以与提供的代码一致。
    for score_file in base_path.glob("*_dockq_score"):
        overall_avg_dockq = None
        a_dockq_values = []
        
        try:
            # 使用 'r' 模式，默认编码通常是 UTF-8
            with open(score_file, 'r') as f:
                lines = f.readlines()
                
                if not lines:
                    print(f"    [WARN] 文件为空，跳过: {score_file.name}")
                    continue

                # 1. 改进的总体平均 DockQ 读取 (第一行)
                first_line = lines[0].strip()
                try:
                    if ':' in first_line:
                        after_colon = first_line.split(':', 1)[1].strip()
                        parts = after_colon.split()
                        
                        for part in parts:
                            try:
                                val = float(part)
                                overall_avg_dockq = val
                                break 
                            except ValueError:
                                continue
                    
                    if overall_avg_dockq is None:
                         print(f"    [WARN] 无法从第一行中提取总体平均 DockQ: {first_line[:40]}...")

                except Exception as e:
                    print(f"    [ERROR] 处理第一行时发生错误: {e}")
                    
                # 2. 读取后续行，提取 A 参与的 DockQ
                for line in lines[1:]:
                    line = line.strip()
                    if not line.startswith("DockQ"):
                        continue
                        
                    # 检查 ' A ' 是否存在于 mapping 字段 (注意两侧空格确保是链名 A)
                    if ' mapping ' in line and ' A ' in line:
                        try:
                            dockq_str = line.split()[1] # DockQ 值是每行的第二个元素
                            dockq_val = float(dockq_str)
                            a_dockq_values.append(dockq_val)
                        except (ValueError, IndexError) as e:
                            print(f"    [WARN] 无法解析 DockQ 行: {line[:50]}... 错误: {e}")

            # 3. 计算 A 的平均 DockQ
            a_avg_dockq = sum(a_dockq_values) / len(a_dockq_values) if a_dockq_values else 0.0
            
            # 4. 汇总结果 (仅数据)
            data_rows.append([
                base_path.name,  # 当前文件夹的名字
                score_file.name, # 文件名
                f"{overall_avg_dockq:.3f}" if overall_avg_dockq is not None else "N/A",
                f"{a_avg_dockq:.3f}" if a_dockq_values else "N/A"
            ])
            
        except FileNotFoundError:
            print(f"    [ERROR] 文件未找到: {score_file}")
        except Exception as e:
            print(f"    [ERROR] 处理文件 {score_file.name} 时发生意外错误: {e}")

    return data_rows

# ----------------------------------------------------------------------
# Main Execution: 主程序 (并行调用)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    start_time = time.time()
    
    # --- 配置 ---
    target_folder = sys.argv[1] 
    output_csv_name = "summary_dockq_scores.csv"
    output_csv_path = Path(target_folder) / output_csv_name

    # 2. 预设 CSV 文件的头部
    header = ["FolderName", "FileName", "Overall_Avg_DockQ", "A_Participating_Avg_DockQ"]
    final_summary_results = []
    final_summary_results.append(header)

    print(f"开始遍历父目录: {target_folder}")

    # 3. 收集所有子文件夹路径
    sub_folder_paths = []
    try:
        for sub_folder_name in os.listdir(target_folder):
            sub_folder_path = os.path.join(target_folder, sub_folder_name)
            if os.path.isdir(sub_folder_path):
                sub_folder_paths.append(sub_folder_path)
    except FileNotFoundError:
        print(f"[FATAL] 目标父目录未找到: {target_folder}")
        exit(1)

    if not sub_folder_paths:
        print(f"[WARN] 在 {target_folder} 中未找到任何子文件夹进行处理。")
        exit(0)

    # 4. 使用多进程池进行并行处理
    # 默认使用所有可用的 CPU 核心
    num_processes = cpu_count()
    print(f"检测到 {num_processes} 个 CPU 核心，将使用 {num_processes} 个进程并行处理 {len(sub_folder_paths)} 个文件夹。")

    # 使用 with Pool 确保进程池被正确关闭
    with Pool(processes=num_processes) as pool:
        # map() 函数将 extract_dockq_scores 应用于 sub_folder_paths 列表中的每个元素
        # results 是一个包含所有子文件夹返回的数据行列表的列表
        results_from_pool = pool.map(extract_dockq_scores, sub_folder_paths)

    # 5. 累积所有结果
    # results_from_pool 是 [[data_rows_folder1], [data_rows_folder2], ...]
    for result_list in results_from_pool:
        final_summary_results.extend(result_list)
        
    # 6. 写入总的 CSV 文件
    if len(final_summary_results) > 1:
        try:
            with open(output_csv_path, "w", newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerows(final_summary_results)
            
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"\n--- 任务完成 ---")
            print(f"总计处理了 {len(sub_folder_paths)} 个文件夹，提取了 {len(final_summary_results) - 1} 条记录。")
            print(f"所有结果已汇总写入: {output_csv_path.resolve()}")
            print(f"总耗时: {elapsed_time:.2f} 秒")

        except Exception as e:
            print(f"\n--- 写入总 CSV 文件时发生错误: {e} ---")
    else:
        print("\n--- 未找到任何有效的 DockQ 结果来汇总。---")
