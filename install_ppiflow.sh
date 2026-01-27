#!/usr/bin/env bash
set -euo pipefail

#############################################
# PPIFlow unified pipeline install script
# - Resumable / idempotent
# - Flag-driven
#############################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_PREFIX="$SCRIPT_DIR"

pkg_manager="conda"
prefix="$DEFAULT_PREFIX"
install_flowpacker=false
install_af3score=false
install_dockq=true
skip_gpu_check=false
install_os_deps=false
write_env=false
no_af3_data_pipeline=true
install_conda=false
conda_prefix=""
af3score_env_name="ppiflow-af3score"
rosetta_env_name="ppiflow-rosetta"

af3_weights_path=""
af3_db_path=""
ppiflow_ckpt_path=""
rosetta_db_path=""
rosetta_bin_path=""
mpnn_weights_dir=""
mpnn_weights_soluble_dir=""
abmpnn_weights_path=""
abmpnn_weights_dir=""
abmpnn_weights_file=""

usage() {
  cat <<EOF
Usage: $0 [options]

Required:
  --af3-weights-path <path>        Path to AF3 weights (af3.bin.zst)
  --ppiflow-checkpoints-path <dir> Dir containing binder.ckpt/antibody.ckpt/nanobody.ckpt/monomer.ckpt

Options:
  --prefix <dir>                   Install root (default: repo root)
  --pkg-manager <conda|mamba>      Package manager (default: conda)
  --install-flowpacker             Clone FlowPacker repo + LFS weights
  --no-install-dockq               Skip DockQ install (default: install)
  --abmpnn-weights-path <dir>      AbMPNN weights dir (override; default: assets/weights/abmpnn)
  --install-af3score               Clone + build AF3Score in a separate env
  --af3score-env-name <name>       AF3Score conda env name (default: ppiflow-af3score)
  --rosetta-env-name <name>        Rosetta conda env name (default: ppiflow-rosetta)
  --no-af3-data-pipeline           Skip AF3 database download/build (default: true)
  --af3-db-path <dir>              AlphaFold3 database directory (required for AF3Score inference)
  --rosetta-db-path <path>         Optional Rosetta database dir (defaults to ppiflow env database)
  --install-os-deps                Install OS deps (apt-get)
  --install-conda                  Install Miniforge if conda is missing
  --conda-prefix <dir>             Miniforge install prefix (default: <prefix>/.miniforge3)
  --skip-gpu-check                 Skip nvidia-smi check
  --write-env                      Write env.sh with paths/exports
  -h, --help                       Show this help
EOF
}

log() { echo -e "[ppiflow-install] $*"; }
die() { echo -e "[ppiflow-install] ERROR: $*" >&2; exit 1; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
}

env_exists() {
  conda env list | awk '{print $1}' | grep -qx "$1"
}

maybe_git_lfs_pull() {
  local repo="$1"
  if [[ ! -d "$repo/.git" ]]; then
    return
  fi
  if command -v git-lfs >/dev/null 2>&1; then
    (cd "$repo" && git lfs install && git lfs pull) || true
  else
    log "Warning: git-lfs not found; LFS files may be missing in $repo"
  fi
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --prefix) prefix="$2"; shift 2 ;;
      --pkg-manager) pkg_manager="$2"; shift 2 ;;
      --install-flowpacker) install_flowpacker=true; shift ;;
      --no-install-dockq) install_dockq=false; shift ;;
      --install-af3score) install_af3score=true; shift ;;
      --no-af3-data-pipeline) no_af3_data_pipeline=true; shift ;;
      --af3-db-path) af3_db_path="$2"; shift 2 ;;
      --af3score-env-name) af3score_env_name="$2"; shift 2 ;;
      --rosetta-env-name) rosetta_env_name="$2"; shift 2 ;;
      --abmpnn-weights-path) abmpnn_weights_path="$2"; shift 2 ;;
      --install-os-deps) install_os_deps=true; shift ;;
      --install-conda) install_conda=true; shift ;;
      --conda-prefix) conda_prefix="$2"; shift 2 ;;
      --skip-gpu-check) skip_gpu_check=true; shift ;;
      --write-env) write_env=true; shift ;;
      --af3-weights-path) af3_weights_path="$2"; shift 2 ;;
      --ppiflow-checkpoints-path) ppiflow_ckpt_path="$2"; shift 2 ;;
      --rosetta-db-path) rosetta_db_path="$2"; shift 2 ;;
      -h|--help) usage; exit 0 ;;
      *) die "Unknown option: $1" ;;
    esac
  done
}

verify_required_flags() {
  [[ -n "$af3_weights_path" ]] || die "--af3-weights-path is required"
  [[ -n "$ppiflow_ckpt_path" ]] || die "--ppiflow-checkpoints-path is required"
  [[ -f "$af3_weights_path" ]] || die "AF3 weights not found: $af3_weights_path"
  [[ -d "$ppiflow_ckpt_path" ]] || die "PPIFlow checkpoints dir not found: $ppiflow_ckpt_path"
  if [[ -n "$af3_db_path" ]]; then
    [[ -d "$af3_db_path" ]] || die "AF3 DB not found: $af3_db_path"
  fi
  if [[ -n "$rosetta_db_path" ]]; then
    [[ -d "$rosetta_db_path" ]] || die "Rosetta DB not found: $rosetta_db_path"
  fi
}

maybe_install_os_deps() {
  if ! $install_os_deps; then
    log "Skipping OS deps install"
    return
  fi
  need_cmd apt-get
  log "Installing OS deps (apt-get)..."
  sudo apt-get update -y
  sudo apt-get install -y \
    git git-lfs wget bzip2 ca-certificates build-essential \
    zlib1g-dev
}

check_gpu() {
  if $skip_gpu_check; then
    log "Skipping GPU check"
    return
  fi
  if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi || log "Warning: nvidia-smi failed"
  else
    log "Warning: nvidia-smi not found"
  fi
}

ensure_conda() {
  if command -v conda >/dev/null 2>&1; then
    if [[ "$pkg_manager" == "mamba" ]]; then
      need_cmd mamba
    fi
    return
  fi
  if ! $install_conda; then
    die "Missing required command: conda (rerun with --install-conda to bootstrap Miniforge)"
  fi
  install_miniforge
  if [[ "$pkg_manager" == "mamba" ]]; then
    need_cmd mamba
  fi
}

install_miniforge() {
  local target="${conda_prefix:-$prefix/.miniforge3}"
  if [[ -d "$target" ]]; then
    export PATH="$target/bin:$PATH"
    return
  fi
  log "Installing Miniforge to $target"
  local installer="/tmp/Miniforge3.sh"
  curl -L -o "$installer" https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash "$installer" -b -p "$target"
  rm -f "$installer"
  export PATH="$target/bin:$PATH"
}

ensure_dirs() {
  mkdir -p "$prefix"/assets/{external,weights,checkpoints,tools}
}


install_ppiflow_env() {
  log "Setting up ppiflow env from environment.yml"
  if env_exists "ppiflow"; then
    log "Conda env 'ppiflow' already exists; skipping"
    return
  fi
  "$pkg_manager" env create -f "$prefix/environment.yml"
}

install_dockq_env() {
  if ! $install_dockq; then
    log "Skipping DockQ install"
    return
  fi
  if [[ ! -d "$prefix/assets/external/DockQ" ]]; then
    log "DockQ repo missing; skipping DockQ install"
    return
  fi
  if ! env_exists "ppiflow"; then
    log "PPIFlow env missing; skipping DockQ install"
    return
  fi
  log "Installing DockQ into 'ppiflow' env..."
  "$pkg_manager" run -n ppiflow pip install -e "$prefix/assets/external/DockQ"
}

install_af3score_env() {
  if ! $install_af3score; then
    log "Skipping AF3Score install"
    return
  fi
  if ! env_exists "$af3score_env_name"; then
    log "Creating AF3Score env '$af3score_env_name' (python 3.11)"
    "$pkg_manager" create -n "$af3score_env_name" -y python=3.11 pip
  fi
  log "Installing AF3Score deps into '$af3score_env_name' env..."
  "$pkg_manager" run -n "$af3score_env_name" pip install -r "$prefix/assets/external/AF3Score/dev-requirements.txt"
  "$pkg_manager" run -n "$af3score_env_name" "$pkg_manager" install -y cmake ninja pybind11 scikit-build-core zlib
  "$pkg_manager" run -n "$af3score_env_name" "$pkg_manager" install -y -c conda-forge biopython h5py pandas
  "$pkg_manager" run -n "$af3score_env_name" pip install --no-build-isolation --no-deps -e "$prefix/assets/external/AF3Score"
  # Ensure CUDA toolkit in env provides a recent ptxas for JAX compilation.
  "$pkg_manager" run -n "$af3score_env_name" "$pkg_manager" install -y -c nvidia cuda-nvcc=12.4 || true
  log "Running AF3Score build_data (ccd/chemical_component_sets pickles)"
  "$pkg_manager" run -n "$af3score_env_name" build_data
}

clone_repo() {
  local url="$1"
  local dest="$2"
  local branch="${3:-}"
  if [[ -d "$dest/.git" ]]; then
    log "Repo already exists: $dest"
    return
  fi
  log "Cloning $url -> $dest"
  if [[ -n "$branch" ]]; then
    git clone --branch "$branch" --depth 1 "$url" "$dest"
  else
    git clone "$url" "$dest"
  fi
}

install_repos() {
  clone_repo "https://github.com/dauparas/ProteinMPNN.git" "$prefix/assets/external/ProteinMPNN"
  maybe_git_lfs_pull "$prefix/assets/external/ProteinMPNN"
  if $install_flowpacker; then
    clone_repo "https://gitlab.com/mjslee0921/flowpacker.git" "$prefix/assets/external/flowpacker"
    if command -v git-lfs >/dev/null 2>&1; then
      (cd "$prefix/assets/external/flowpacker" && git lfs install && git lfs pull) || true
    else
      log "Warning: git-lfs not found; FlowPacker weights may be missing"
    fi
  fi
  if $install_af3score; then
    clone_repo "https://github.com/Mingchenchen/AF3Score.git" "$prefix/assets/external/AF3Score" "v2.0.0"
  fi
  if $install_dockq; then
    clone_repo "https://github.com/bjornwallner/DockQ.git" "$prefix/assets/external/DockQ"
  fi
}

install_rosetta_cli() {
  log "Installing Rosetta CLI (rosetta_scripts) into '$rosetta_env_name' env"
  "$pkg_manager" config --add channels https://conda.rosettacommons.org
  "$pkg_manager" config --add channels defaults
  if ! env_exists "$rosetta_env_name"; then
    log "Creating Rosetta env '$rosetta_env_name' (python 3.12)"
    "$pkg_manager" create -n "$rosetta_env_name" -y python=3.12
  fi
  if ! "$pkg_manager" install -n "$rosetta_env_name" -y rosetta; then
    log "Warning: Rosetta install failed for env '$rosetta_env_name'."
    log "You can retry manually or set ROSETTA_BIN/ROSETTA_DB to an existing install."
    return 0
  fi
  local conda_base
  conda_base="$("$pkg_manager" info --base)"
  rosetta_bin_path="$conda_base/envs/$rosetta_env_name/bin/rosetta_scripts"
  local rosetta_db_default="$conda_base/envs/$rosetta_env_name/database"
  [[ -x "$rosetta_bin_path" ]] || die "rosetta_scripts not found at $rosetta_bin_path"
  if [[ -z "$rosetta_db_path" ]]; then
    rosetta_db_path="$rosetta_db_default"
  fi
  [[ -d "$rosetta_db_path" ]] || die "Rosetta DB not found: $rosetta_db_path"
  ln -s "$rosetta_bin_path" "$prefix/assets/tools/rosetta_scripts" || true
  ln -s "$rosetta_db_path" "$prefix/assets/tools/rosetta_db" || true
}

place_weights_and_ckpts() {
  log "Placing PPIFlow checkpoints"
  local ckpt_dir="$prefix/assets/checkpoints"
  mkdir -p "$ckpt_dir"
  for f in binder.ckpt antibody.ckpt nanobody.ckpt monomer.ckpt; do
    if [[ -f "$ckpt_dir/$f" ]]; then
      log "Checkpoint already present: $ckpt_dir/$f"
    else
      [[ -f "$ppiflow_ckpt_path/$f" ]] || die "Missing checkpoint: $ppiflow_ckpt_path/$f"
      ln -s "$ppiflow_ckpt_path/$f" "$ckpt_dir/$f"
    fi
  done

  log "Placing AF3 weights"
  local af3_dir="$prefix/assets/weights/af3"
  mkdir -p "$af3_dir"
  local af3_target="$af3_dir/af3.bin.zst"
  local af3_bin="$af3_dir/af3.bin"
  if [[ -f "$af3_target" ]]; then
    log "AF3 weights already present: $af3_target"
  else
    ln -s "$af3_weights_path" "$af3_target"
  fi
  if [[ -f "$af3_bin" ]]; then
    log "AF3 weights already decompressed: $af3_bin"
  else
    if command -v zstd >/dev/null 2>&1; then
      log "Decompressing AF3 weights -> $af3_bin"
      zstd -d -f "$af3_target" -o "$af3_bin"
    else
      log "Warning: zstd not found; AF3 weights remain compressed at $af3_target"
    fi
  fi
  if [[ -n "$af3_db_path" ]]; then
    mkdir -p "$prefix/assets/tools"
    ln -s "$af3_db_path" "$prefix/assets/tools/af3_db" || true
  elif $no_af3_data_pipeline; then
    # Create a stub AF3 database directory to satisfy AF3Score file checks.
    local stub_db="$prefix/assets/tools/af3_db_stub"
    mkdir -p "$stub_db/mmcif_files"
    touch "$stub_db/bfd-first_non_consensus_sequences.fasta"
    touch "$stub_db/mgy_clusters_2022_05.fa"
    touch "$stub_db/uniprot_all_2021_04.fa"
    touch "$stub_db/uniref90_2022_05.fa"
    touch "$stub_db/nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta"
    touch "$stub_db/rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta"
    touch "$stub_db/rnacentral_active_seq_id_90_cov_80_linclust.fasta"
    touch "$stub_db/pdb_seqres_2022_09_28.fasta"
    af3_db_path="$stub_db"
    log "Created stub AF3 DB dir at $stub_db (set --af3-db-path to override)."
  fi
}

place_mpnn_weights() {
  # ProteinMPNN ships its own weights via git-lfs.
  local mpnn_repo="$prefix/assets/external/ProteinMPNN"
  local vanilla_dir="$mpnn_repo/vanilla_model_weights"
  local soluble_dir="$mpnn_repo/soluble_model_weights"
  local ca_dir="$mpnn_repo/ca_model_weights"
  local legacy_dir="$mpnn_repo/model_weights"

  if [[ ! -d "$vanilla_dir" && ! -d "$soluble_dir" && -d "$mpnn_repo" ]]; then
    maybe_git_lfs_pull "$mpnn_repo"
  fi

  mkdir -p "$prefix/assets/weights/mpnn"
  if [[ -d "$vanilla_dir" ]]; then
    mpnn_weights_dir="$vanilla_dir"
    ln -s "$vanilla_dir" "$prefix/assets/weights/mpnn/vanilla" || true
    ln -s "$vanilla_dir" "$prefix/assets/weights/mpnn/weights" || true
    log "Using ProteinMPNN vanilla weights -> $vanilla_dir"
  elif [[ -d "$legacy_dir" ]]; then
    mpnn_weights_dir="$legacy_dir"
    ln -s "$legacy_dir" "$prefix/assets/weights/mpnn/weights" || true
    log "Using ProteinMPNN legacy weights -> $legacy_dir"
  else
    log "Warning: ProteinMPNN vanilla_model_weights not found at $vanilla_dir (binder sequence design may fail)"
  fi

  if [[ -d "$soluble_dir" ]]; then
    mpnn_weights_soluble_dir="$soluble_dir"
    ln -s "$soluble_dir" "$prefix/assets/weights/mpnn/soluble" || true
    log "Using ProteinMPNN soluble weights -> $soluble_dir"
  else
    log "Warning: ProteinMPNN soluble_model_weights not found at $soluble_dir"
  fi

  if [[ -d "$ca_dir" ]]; then
    ln -s "$ca_dir" "$prefix/assets/weights/mpnn/ca" || true
  fi

  mkdir -p "$prefix/assets/weights/abmpnn"
  if [[ -n "$abmpnn_weights_path" ]]; then
    [[ -d "$abmpnn_weights_path" ]] || die "AbMPNN weights dir not found: $abmpnn_weights_path"
    ln -s "$abmpnn_weights_path/abmpnn.pt" "$prefix/assets/weights/abmpnn/abmpnn.pt" || true
    abmpnn_weights_dir="$prefix/assets/weights/abmpnn"
    abmpnn_weights_file="$prefix/assets/weights/abmpnn/abmpnn.pt"
    log "Linked AbMPNN torch weights -> $abmpnn_weights_file"
  else
    # default bundled location
    if [[ -f "$prefix/assets/weights/abmpnn/abmpnn.pt" ]]; then
      abmpnn_weights_dir="$prefix/assets/weights/abmpnn"
      abmpnn_weights_file="$prefix/assets/weights/abmpnn/abmpnn.pt"
      log "Using bundled AbMPNN weights -> $abmpnn_weights_file"
    else
      maybe_git_lfs_pull "$prefix"
      log "Warning: AbMPNN weights not found at $prefix/assets/weights/abmpnn/abmpnn.pt"
    fi
  fi

  # If ProteinMPNN repo exists, place a symlink in its model_weights/ for convenience
  if [[ -d "$prefix/assets/external/ProteinMPNN" && -f "$prefix/assets/weights/abmpnn/abmpnn.pt" ]]; then
    mkdir -p "$prefix/assets/external/ProteinMPNN/model_weights"
    ln -s "$prefix/assets/weights/abmpnn/abmpnn.pt" "$prefix/assets/external/ProteinMPNN/model_weights/abmpnn.pt" || true
    log "Symlinked AbMPNN weights into ProteinMPNN model_weights/"
  fi
}

write_env_file() {
  if ! $write_env; then
    return
  fi
  local env_path="$prefix/env.sh"
  local conda_base
  conda_base="$("$pkg_manager" info --base)"
  local dockq_bin_path=""
  if [[ -x "$conda_base/envs/ppiflow/bin/DockQ" ]]; then
    dockq_bin_path="$conda_base/envs/ppiflow/bin/DockQ"
  elif [[ -f "$prefix/assets/external/DockQ/DockQ.py" ]]; then
    dockq_bin_path="$prefix/assets/external/DockQ/DockQ.py"
  fi
  cat > "$env_path" <<EOF
export PPIFLOW_ROOT="$prefix"
export ROSETTA_BIN="${rosetta_bin_path}"
export ROSETTA_DB="${rosetta_db_path}"
export ROSETTA_ENV="${rosetta_env_name}"
export AF3_WEIGHTS="$prefix/assets/weights/af3"
export AF3_DB_DIR="${af3_db_path}"
export AF3SCORE_REPO="$prefix/assets/external/AF3Score"
export AF3SCORE_ENV="${af3score_env_name}"
export AF3SCORE_PYTHON="$conda_base/envs/${af3score_env_name}/bin/python"
export AF3SCORE_CUDA_HOME="$conda_base/envs/${af3score_env_name}"
export FLOWPACKER_REPO="$prefix/assets/external/flowpacker"
export PROTEINMPNN_REPO="$prefix/assets/external/ProteinMPNN"
export MPNN_WEIGHTS="${mpnn_weights_dir}"
export MPNN_SOLUBLE_WEIGHTS="${mpnn_weights_soluble_dir}"
export ABMPNN_WEIGHTS_DIR="${abmpnn_weights_dir}"
export ABMPNN_WEIGHTS_FILE="${abmpnn_weights_file}"
export ABMPNN_WEIGHTS="${abmpnn_weights_file}"
export DOCKQ_BIN="${dockq_bin_path}"
export PATH="$prefix/assets/tools:\$PATH"
EOF
  log "Wrote env file: $env_path"
}

main() {
  parse_args "$@"
  verify_required_flags
  ensure_conda
  maybe_install_os_deps
  check_gpu
  ensure_dirs
  install_repos
  install_ppiflow_env
  install_dockq_env
  install_af3score_env
  install_rosetta_cli
  place_weights_and_ckpts
  place_mpnn_weights
  write_env_file
  log "Install complete."
}

main "$@"
