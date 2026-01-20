import os, glob, time
import sys
import pandas as pd
from multiprocessing import Pool
import pyrosetta
from pyrosetta import *
from pyrosetta.teaching import *

#Protocol Includes
#from rosetta.protocols import relax as rel
#from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
#from rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

# PyRosetta 初始化
init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res \
     -ignore_zero_occupancy false -load_PDB_components false  -relax:default_repeats 2 -no_fconfig')

def bool_type(bool_str: str):
    bool_str_lower = bool_str.lower()
    if bool_str_lower in ('false', 'f', 'no', 'n', '0'):
        return False
    elif bool_str_lower in ('true', 't', 'yes', 'y', '1'):
        return True
    else:
        raise ValueError(f'Cannot interpret {bool_str} as bool')

# file_name 函数在 main 中未使用，保持不变
def file_name(file_dir):
    L = []
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if os.path.splitext(file)[1] == '.pdb':
                L.append(os.path.join(root , file))
        return L

def main(args):
    
    output = []
    
    # --- 修复 Bug 4: 将 ScoreFunction 定义移到循环外部以提高效率 ---
    # 定义用于 Relax 和总分的 ScoreFunction (通常是默认的 ref2015 或其他)
    scorefxn = get_score_function() 
    
    # 定义用于单项评分的 ScoreFunction
    scorefxn_fa_atr = ScoreFunction(); scorefxn_fa_atr.set_weight(fa_atr, 1.0)
    scorefxn_fa_rep = ScoreFunction(); scorefxn_fa_rep.set_weight(fa_rep, 1.0)
    scorefxn_fa_intra_rep = ScoreFunction(); scorefxn_fa_intra_rep.set_weight(fa_intra_rep, 1.0)
    scorefxn_fa_sol = ScoreFunction(); scorefxn_fa_sol.set_weight(fa_sol, 1.0)
    scorefxn_lk_ball_wtd = ScoreFunction(); scorefxn_lk_ball_wtd.set_weight(lk_ball_wtd, 1.0)
    scorefxn_fa_intra_sol = ScoreFunction(); scorefxn_fa_intra_sol.set_weight(fa_intra_sol, 1.0)
    scorefxn_fa_elec = ScoreFunction(); scorefxn_fa_elec.set_weight(fa_elec, 1.0)
    scorefxn_hbond_lr_bb = ScoreFunction(); scorefxn_hbond_lr_bb.set_weight(hbond_lr_bb, 1.0)
    scorefxn_hbond_sr_bb = ScoreFunction(); scorefxn_hbond_sr_bb.set_weight(hbond_sr_bb, 1.0)
    scorefxn_hbond_bb_sc = ScoreFunction(); scorefxn_hbond_bb_sc.set_weight(hbond_bb_sc, 1.0)
    scorefxn_hbond_sc = ScoreFunction(); scorefxn_hbond_sc.set_weight(hbond_sc, 1.0)
    scorefxn_dslf_fa13 = ScoreFunction(); scorefxn_dslf_fa13.set_weight(dslf_fa13, 1.0)
    scorefxn_rama_prepro = ScoreFunction(); scorefxn_rama_prepro.set_weight(rama_prepro, 1.0)
    scorefxn_p_aa_pp = ScoreFunction(); scorefxn_p_aa_pp.set_weight(p_aa_pp, 1.0)
    scorefxn_fa_dun = ScoreFunction(); scorefxn_fa_dun.set_weight(fa_dun, 1.0)
    scorefxn_omega = ScoreFunction(); scorefxn_omega.set_weight(omega, 1.0)
    scorefxn_pro_close = ScoreFunction(); scorefxn_pro_close.set_weight(pro_close, 1.0)
    scorefxn_yhh_planarity = ScoreFunction(); scorefxn_yhh_planarity.set_weight(yhh_planarity, 1.0)
    scorefxn_ref = ScoreFunction(); scorefxn_ref.set_weight(ref, 1.0)
    scorefxn_rg = ScoreFunction(); scorefxn_rg.set_weight(rg , 1 )
    
    # 预处理固定链信息
    fixed_chain = args.fixed_chain.split("_")
    
    # --- 修复 Bug 1 & 2: 分离两种模式，解决 NameError 和迭代源错误 ---

    if args.csv_path != "":
        # --- 模式 1: 使用 CSV 文件 (包含额外信息) ---
        df = pd.read_csv(args.csv_path)
        
        if len(df) == 0:
            print(f"CSV file at {args.csv_path} is empty.")
            return

        print(f"Processing {len(df)} pdb from CSV")
        # 循环遍历 DataFrame 的每一行，获取 PDB 路径和界面信息
        for idx, row in df.iterrows():
            pdb = row["pdb"]
            pdb_path = pdb # 保持变量名与原始代码一致

            # --------------------------- 核心计算代码 (保持不变) ---------------------------
            start = time.time()
            pdb_name = os.path.basename(pdb_path).replace(".pdb","")
            pose = pose_from_pdb(pdb_path)
            original_pose = pose.clone()
            
            # FastRelax 设置和应用
            if args.relax:
                fr = FastRelax()
                fr.set_scorefxn(scorefxn)
                fr.max_iter(int(args.max_iter))
                movemap = MoveMap()
                movemap.set_bb(True)
                movemap.set_chi(True)
                
                if args.fixbb:
                    for i in range(1, pose.total_residue() + 1):
                        chain = pose.pdb_info().chain(i)
                        if chain in fixed_chain:
                            movemap.set_bb(i, False)
                        else:
                            movemap.set_bb(i, True)
                fr.set_movemap(movemap)
            
            # 应用 Relax
            if (not os.getenv("DEBUG")) and (args.relax):
                # 确保在 fr 存在时才调用 apply
                if 'fr' in locals():
                     fr.apply(pose)
            
            # interface analysis and energy calculation
            ia = InterfaceAnalyzerMover()
            # 在 CSV 模式下，可以安全地从 row 中获取 interface 信息
            if ("ligand" in df.columns) and ("receptor" in df.columns):
                interface = row["ligand"].replace(",","")+"_"+row["receptor"].replace(",","")
                ia.set_interface(interface)
            else:
                 # ⚠️ 警告：如果 CSV 缺少这两列，InterfaceAnalyzerMover 将使用默认界面（所有链）
                 print(f"Warning: 'ligand' or 'receptor' column missing for {pdb_name}. Using default interface.")

            ia.set_skip_reporting(True)
            ia.apply(pose)

            time_consumed = time.time()-start

            # 结果收集 (保持不变)
            temp_dict = {'relaxed':scorefxn(pose),'interface_score':ia.get_interface_dG(), 
            'original':scorefxn(original_pose),'delta':scorefxn(pose) - scorefxn(original_pose),
            'fa_atr':scorefxn_fa_atr(pose),'fa_rep':scorefxn_fa_rep(pose),'fa_intra_rep':scorefxn_fa_intra_rep(pose),
            'fa_sol':scorefxn_fa_sol(pose),'lk_ball_wtd':scorefxn_lk_ball_wtd(pose),'fa_intra_sol':scorefxn_fa_intra_sol(pose),
            'fa_elec':scorefxn_fa_elec(pose),'hbond_lr_bb':scorefxn_hbond_lr_bb(pose),'hbond_sr_bb(pose)':scorefxn_hbond_sr_bb(pose),
            'hbond_bb_sc':scorefxn_hbond_bb_sc(pose),'hbond_sc':scorefxn_hbond_sc(pose),'dslf_fa13':scorefxn_dslf_fa13(pose),
            'rama_prepro':scorefxn_rama_prepro(pose),'p_aa_pp':scorefxn_p_aa_pp(pose),'fa_dun':scorefxn_fa_dun(pose),
            'omega':scorefxn_omega(pose),'pro_close':scorefxn_pro_close(pose),'yhh_planarity':scorefxn_yhh_planarity(pose),
            'ref':scorefxn_ref(pose), 'get_complexed_sasa':ia.get_complexed_sasa(),
            'get_interface_delta_sasa':ia.get_interface_delta_sasa(), 'time_consumed': time_consumed}

            if args.dump_pdb: pose.dump_pdb(os.path.join(args.output_dir, 'relax_' + os.path.basename(pdb_path)))

            output.append([pdb_name]+[v for k, v in temp_dict.items()])
            # --------------------------- 核心计算代码结束 ---------------------------

    else:
        # --- 模式 2: 使用 PDB 文件夹 (只获取文件路径) ---
        all_pdb_name = sorted(glob.glob(os.path.join(args.pdb_dir, '*.pdb')))
        
        if len(all_pdb_name) > 0:
            print(f"Processing {len(all_pdb_name)} pdb in {args.pdb_dir}")
        else:
            print(f"No pdb available in {args.pdb_dir}")
            return
            
        # 循环遍历 PDB 文件的路径列表
        for pdb_path in all_pdb_name: 

            # --------------------------- 核心计算代码 (保持不变) ---------------------------
            start = time.time()
            pdb_name = os.path.basename(pdb_path).replace(".pdb","")
            pose = pose_from_pdb(pdb_path)
            original_pose = pose.clone()
            
            # FastRelax 设置和应用
            if args.relax:
                fr = FastRelax()
                fr.set_scorefxn(scorefxn)
                fr.max_iter(int(args.max_iter))
                movemap = MoveMap()
                movemap.set_bb(True)
                movemap.set_chi(True)
                
                if args.fixbb:
                    for i in range(1, pose.total_residue() + 1):
                        chain = pose.pdb_info().chain(i)
                        if chain in fixed_chain:
                            movemap.set_bb(i, False)
                        else:
                            movemap.set_bb(i, True)
                fr.set_movemap(movemap)
            
            # 应用 Relax
            if (not os.getenv("DEBUG")) and (args.relax):
                 if 'fr' in locals():
                     fr.apply(pose)

            # interface analysis and energy calculation
            ia = InterfaceAnalyzerMover()
            # ⚠️ 注意：在此模式下，没有 interface 信息，InterfaceAnalyzerMover 将使用默认界面
            print(f"Warning: No interface information provided for {pdb_name}. Using default interface.")
            ia.set_skip_reporting(True)
            ia.apply(pose)

            time_consumed = time.time()-start

            # 结果收集 (保持不变)
            temp_dict = {'relaxed':scorefxn(pose),'interface_score':ia.get_interface_dG(), 
            'original':scorefxn(original_pose),'delta':scorefxn(pose) - scorefxn(original_pose),
            'fa_atr':scorefxn_fa_atr(pose),'fa_rep':scorefxn_fa_rep(pose),'fa_intra_rep':scorefxn_fa_intra_rep(pose),
            'fa_sol':scorefxn_fa_sol(pose),'lk_ball_wtd':scorefxn_lk_ball_wtd(pose),'fa_intra_sol':scorefxn_fa_intra_sol(pose),
            'fa_elec':scorefxn_fa_elec(pose),'hbond_lr_bb':scorefxn_hbond_lr_bb(pose),'hbond_sr_bb(pose)':scorefxn_hbond_sr_bb(pose),
            'hbond_bb_sc':scorefxn_hbond_bb_sc(pose),'hbond_sc':scorefxn_hbond_sc(pose),'dslf_fa13':scorefxn_dslf_fa13(pose),
            'rama_prepro':scorefxn_rama_prepro(pose),'p_aa_pp':scorefxn_p_aa_pp(pose),'fa_dun':scorefxn_fa_dun(pose),
            'omega':scorefxn_omega(pose),'pro_close':scorefxn_pro_close(pose),'yhh_planarity':scorefxn_yhh_planarity(pose),
            'ref':scorefxn_ref(pose), 'get_complexed_sasa':ia.get_complexed_sasa(),
            'get_interface_delta_sasa':ia.get_interface_delta_sasa(), 'time_consumed': time_consumed}

            if args.dump_pdb: pose.dump_pdb(os.path.join(args.output_dir, 'relax_' + os.path.basename(pdb_path)))

            output.append([pdb_name]+[v for k, v in temp_dict.items()])
            # --------------------------- 核心计算代码结束 ---------------------------

    # 循环结束后，统一将 output 写入 CSV
    if len(output) > 0:
        # 使用其中一个 temp_dict 结构来确定列名
        score_df = pd.DataFrame(output, columns=['pdb_name']+list(temp_dict.keys()))
        score_df.to_csv(os.path.join(args.output_dir, f'rosetta_complex_{args.batch_idx}.csv'),index=False)
    else:
        print("No results to write. Output list is empty.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_dir',type=str, default="", help="input pdb path.")
    parser.add_argument('--csv_path',type=str, default="", help="save absolute path of pdb file in the 'pdb' column")
    parser.add_argument('--output_dir',type=str, help="output relaxed pdb path", required=True)
    parser.add_argument('--dump_pdb', default= False, type=bool_type, help="dump pdb or not")
    parser.add_argument('--batch_idx', type=int, default= False, help="batch index for output filename") # 修复 default 类型
    parser.add_argument('--relax', default= False, type=bool_type, help="run relax or not") # 修复 help 文本
    parser.add_argument('--fixbb', default= False, type=bool_type, help="fix backbone or not") # 修复 help 文本
    parser.add_argument('--fixed_chain', default="", type=str, help="chains whose backbone should be fixed (e.g., 'A_B')") # 修复 help 文本
    parser.add_argument('--max_iter', default=1, type=int, help="max iteration for FastRelax") # 修复 help 文本
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    main(args)