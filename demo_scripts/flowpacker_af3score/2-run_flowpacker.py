import os
import yaml
import glob
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--base_yaml_path', type=str, default='/lustre/grp/cmclab/share/wuxy/flowpacker-main/config/inference/base.yaml')
parser.add_argument('--base_pdb_path', type=str, default='/lustre/grp/cmclab/share/wuxy/process/flowpacker_af3score/data/input_pdb_batch')
parser.add_argument('--output_yaml_dir', type=str, default='/lustre/grp/cmclab/share/wuxy/process/flowpacker_af3score/data/yaml_batch')
args = parser.parse_args()

base_yaml_path = args.base_yaml_path
base_pdb_path = args.base_pdb_path
output_yaml_dir = args.output_yaml_dir

with open(base_yaml_path, 'r') as f:
    base_config = yaml.safe_load(f)

batch_dirs = glob.glob(os.path.join(base_pdb_path, 'batch_*'))

for sub in batch_dirs:
    batch_name = os.path.basename(sub) 
    
    batch_yaml_path = f"{output_yaml_dir}/{batch_name}.yml"

    if os.path.exists(batch_yaml_path):
        print(f"🗑️ 删除已存在的 YML 文件: {batch_yaml_path}")
        os.remove(batch_yaml_path)


    config = base_config.copy()
    config["data"]["test_path"] = sub


    with open(batch_yaml_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    print(f"✔ Yaml Generated: {batch_yaml_path}")
