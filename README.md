# PPIFlow

![](./assets/model.png)

PPIFlow is a flow-matching framework and unified design pipeline for de novo
protein binders. This repo supports binder, antibody, and VHH protocols with a
single resumable CLI pipeline. Legacy scripts and notebooks remain available for
reference.

## Quickstart (Unified Pipeline)

### 1) Install

Recommended: use the install script to set up the environment and external tool
paths in one place.

```bash
./install_ppiflow.sh \
  --af3-weights-path /path/to/af3.bin.zst \
  --ppiflow-checkpoints-path /path/to/ppiflow_checkpoints \
  --install-flowpacker \
  --install-proteinmpnn \
  --install-af3score \
  --write-env
```

Optional: pass `--rosetta-db-path /path/to/rosetta/database` to override the
default Rosetta database bundled with the installer.

Manual option (if you want to manage tools yourself):

```bash
conda env create -f environment.yml
conda activate ppiflow
```

You will still need paths to external tools (FlowPacker, ProteinMPNN or
AbMPNN, AF3Score/AF3) and PPIFlow checkpoints. Rosetta CLI (`rosetta_scripts`)
is installed into the main `ppiflow` conda environment by the installer.

### 2) Run with a YAML input

Create `design.yaml` (see schema below), then run:

```bash
python ppiflow.py pipeline \
  --protocol binder \
  --preset fast \
  --input /path/to/design.yaml \
  --output /path/to/out_dir
```

### 3) Run with CLI-only inputs (no YAML)

```bash
python ppiflow.py pipeline \
  --protocol binder \
  --name il7ra_001 \
  --target_pdb /path/to/target.pdb \
  --target_chains B \
  --hotspots B119,B141,B200 \
  --binder_length 75-90 \
  --samples_per_target 100 \
  --ppiflow_ckpt /path/to/binder.ckpt \
  --output /path/to/out_dir
```

The CLI will write `config/pipeline_input.yaml` and `pipeline_input.json` in the
output directory and then run the pipeline.

## Pipeline CLI

The unified CLI lives in `ppiflow.py` and exposes four commands:

- `pipeline`: configure + execute
- `configure`: write per-step configs + `steps.yaml`
- `execute`: run from `steps.yaml`
- `rank`: run the rank step only

Common flags:

- `--protocol`: `binder|antibody|vhh` (required unless input YAML already has it)
- `--preset`: `fast|full|custom` (default: `full`)
- `--input`: path to `design.yaml` (optional)
- `--output`: output directory (required)
- `--steps`: `all|gen,seq,score,rosetta,partial,rank`
- `--reuse`: skip outputs that already exist
- `--partial_start_t`: partial flow start_t (default 0.6)
- `--af3score1_iptm_min`: round-1 ipTM threshold (default 0.2)
- `--af3score2_iptm_min`: round-2 ipTM threshold (default 0.5)
- `--af3score2_ptm_min`: round-2 pTM threshold (default 0.8)
- `--interface_energy_min`: interface energy cutoff (default -5.0)
- `--interface_distance`: interface distance cutoff (default 12.0)
- `--relax_max_iter`: Rosetta relax max iterations (default 170)
- `--relax_fixbb`: fix backbone for chains listed in `--fixed_chains`
- `--fixed_chains`: chains to fix (underscore or comma-separated, e.g. `A_B` or `A,B`)
- `--dockq_min`: DockQ threshold (default 0.49)
- `--rank_top_k`: number of finalists (default 30)

CLI-only input flags (when `--input` is omitted):

- `--name`
- `--target_pdb`, `--target_chains`, `--hotspots`
- Binder: `--binder_length`
- Antibody/VHH: `--framework_pdb`, `--heavy_chain`, `--light_chain` (omit for VHH), `--cdr_length`
- Tools: `--ppiflow_ckpt`, `--mpnn_ckpt`, `--abmpnn_ckpt`, `--af3score_repo`,
  `--flowpacker_repo`, `--af3_weights`, `--mpnn_repo`, `--rosetta_bin`,
  `--rosetta_db`, `--abmpnn_repo`, `--mpnn_run`, `--abmpnn_run`, `--dockq_bin`
- Sequence design knobs: `--seq1_num_per_backbone`, `--seq1_temp`,
  `--seq1_bias_large_residues`, `--seq1_bias_num`, `--seq2_num_per_backbone`,
  `--seq2_temp`, `--seq2_use_soluble_ckpt`

## Input Schema (`design.yaml`)

Minimal binder example:

```yaml
protocol: binder
name: il7ra_001

target:
  pdb: /path/to/target.pdb
  chains: ["B"]
  hotspots: ["B119", "B141", "B200"]

binder:
  length: "75-90"

sampling:
  samples_per_target: 100

tools:
  ppiflow_ckpt: /path/to/binder.ckpt
  mpnn_ckpt: /path/to/proteinmpnn
  mpnn_run: /path/to/protein_mpnn_run.py
  af3score_repo: /path/to/af3score
  flowpacker_repo: /path/to/flowpacker
  # rosetta_bin optional if installer wires assets/tools/rosetta_scripts
  rosetta_bin: /path/to/rosetta_scripts
  # rosetta_db is optional; installer wires it if available.
  rosetta_db: /path/to/rosetta/database
  af3_weights: /path/to/af3.bin.zst
sequence_design:
  round1:
    num_seq_per_backbone: 16
    sampling_temp: 0.2
    bias_large_residues: true
    bias_num: 8
  round2:
    num_seq_per_backbone: 8
    sampling_temp: 0.1
partial:
  start_t: 0.6
  samples_per_target: 8
filters:
  af3score:
    round1:
      iptm_min: 0.2
    round2:
      iptm_min: 0.5
      ptm_min: 0.8
  rosetta:
    interface_energy_min: -5.0
    interface_distance: 12.0
    relax_max_iter: 170
    relax_fixbb: false
    fixed_chains: ""   # e.g. "A_B" to fix A and B backbones
    # Defaults above mirror the legacy VHH notebook relax settings.
  dockq:
    min: 0.49
ranking:
  top_k: 30
  composite_score: "iptm*100 - interface_score"

# Set dockq.min to null to disable the DockQ step.
```

Antibody or VHH example:

```yaml
protocol: antibody  # or vhh
name: il13_001

target:
  pdb: /path/to/antigen.pdb
  chains: ["C"]
  hotspots: ["C11", "C14", "C15"]

framework:
  pdb: /path/to/framework.pdb
  heavy_chain: A
  light_chain: B        # omit for vhh
  cdr_length: "CDRH1,8-8,CDRH2,8-8,CDRH3,10-20,CDRL1,6-9,CDRL2,3-3,CDRL3,9-11"

sampling:
  samples_per_target: 100

tools:
  ppiflow_ckpt: /path/to/antibody.ckpt
  abmpnn_ckpt: /path/to/abmpnn
  abmpnn_run: /path/to/protein_mpnn_run.py
  af3score_repo: /path/to/af3score
  flowpacker_repo: /path/to/flowpacker
  rosetta_db: /path/to/rosetta/database
  af3_weights: /path/to/af3.bin.zst
```

Validation highlights:

- `protocol`, `name`, `target.pdb`, `target.chains` are required.
- Binder: `binder.length` required; no `framework` block.
- Antibody: `framework.pdb`, `framework.heavy_chain`, `framework.light_chain`, `framework.cdr_length` required.
- VHH: `framework.pdb`, `framework.heavy_chain`, `framework.cdr_length` required; omit `framework.light_chain`.

## Multi-Chain Targets, Concatenation, and Chain Conventions

Supported hotspot specs (ColabDesign/BindCraft-style):

- `C` → whole chain C is hotspot
- `C62-73` → range on chain C
- `C62` → single residue
- Mixed list: `A,C1-10,C45`

Behavior (pipeline runs):

- **All protocols** accept multi-chain targets via `target.chains`. During
  `configure`, all target chains are concatenated into a single **chain B** PDB
  written under `output/inputs/`.
- **Chain conventions**:
  - Binder/VHH: binder is **A**, target is **B**
  - Antibody: heavy **A**, light **C**, target **B**
  Framework PDBs are rewritten under `output/inputs/` to enforce this, and the
  binder chain is forced to **A** in pipeline inputs.
- **Hotspots** are specified on the original target PDB and then mapped onto
  concatenated chain **B**. The concatenated target uses fresh integer residue
  numbering; mapping is recorded in `output/inputs/target_chain_map.json`.

See `docs/target_chain_concatenation_spec.md` for details.

## External Tool Wiring

The unified pipeline orchestrates multiple external tools. The pipeline engine
runs PPIFlow generation, Rosetta interface/relax, and partial-flow steps
directly. Sequence design can run automatically if `mpnn_run` or `abmpnn_run`
is provided; FlowPacker, AF3Score, and DockQ are executed through external
commands.

After `configure`, edit these per-step configs to provide commands:

- `config/step_seq1.yaml`
- `config/step_flowpacker1.yaml`
- `config/step_af3score1.yaml`
- `config/step_rosetta_interface.yaml` (Rosetta CLI runs automatically)
- `config/step_interface_enrich.yaml` (uses Rosetta interface energies)
- `config/step_seq2.yaml`
- `config/step_flowpacker2.yaml`
- `config/step_af3score2.yaml`
- `config/step_relax.yaml` (Rosetta CLI runs automatically)
- `config/step_dockq.yaml` (optional)
- `config/step_af3_refold.yaml` (optional)

Each step config supports a `command` field that can be a string or list, and
will be executed as a subprocess during `execute`.

## Output Layout and Resume

The pipeline writes a deterministic layout under your output directory:

- `config/`: per-step YAMLs generated by `configure`
- `steps.yaml`: manifest of pipeline steps
- `pipeline_input.json`: normalized input (resume identity)
- `pipeline_state.json`: run metadata + tool stamps
- `runs/run_XXX/`: per-run artifacts
- `runs/run_XXX/interface_enrich/`: fixed positions from Rosetta interface energies
- `runs/run_XXX/flowpacker_round2/`, `runs/run_XXX/af3score_round2/`: round-2 maturation
- `runs/run_XXX/dockq/`: optional DockQ results
- `results/` and `results_vN/`: ranked outputs (rank step snapshots)

Resume behavior:

- `--reuse` skips outputs that already exist
- missing-index fill supports changing `WORLD_SIZE` across resumes
- mismatched input hash or tool stamps requires a new output directory

Heartbeat/status files (optional):

- `status.json` and `status_rank*.json` written to the output root
- set `PPIFLOW_HEARTBEAT=0` to disable

## Troubleshooting

- `steps.yaml not found`: run `configure` (or `pipeline`) first.
- `pipeline_state.json ... mismatch`: change output directory or ensure tool paths
  and inputs match the original run.
- External tools failing: verify paths in `tools` or CLI flags, and confirm
  per-step `command` entries are correct.

## Appendix: Legacy Scripts (Short)

Legacy demos and shell pipelines from the original PPIFlow repo live under `legacy/`, and the original entrypoint scripts remain in `src/entrypoints/`. These are not the recommended path for end-to-end usage but are kept for reference.

## License

The original PPIFlow code is licensed under an Attribution-NonCommercial-ShareAlike 4.0 International license. This pipeline is licensed under the same license to respect the upstream repository.
