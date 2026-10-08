#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# IsoForge — install from source
# ─────────────────────────────────────────────────────────────
#
#   bash install.sh                    ← conda environment 'isoforge'
#   bash install.sh myenv              ← another environment name
#   bash install.sh /path/to/env       ← an environment prefix (shared clusters)
#   bash install.sh --use-current-env  ← build into the active environment
#
#   CONDA_BASE=/path/to/conda bash install.sh   ← conda outside PATH
#
# With conda, the script installs Python, Rust and the C build tools itself.
# ─────────────────────────────────────────────────────────────
set -eo pipefail

PYTHON_VER="3.12"
RUST_MIN="1.85"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Arguments ────────────────────────────────────────────────
USE_CURRENT_ENV=""
_args=()
for _a in "$@"; do
    case "${_a}" in
        --use-current-env) USE_CURRENT_ENV=1 ;;
        -h|--help)
            sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) _args+=("${_a}") ;;
    esac
done
# A name is a named environment; a path is an environment prefix.
ENV_NAME="${_args[0]:-isoforge}"

# ── Messages ─────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; NC='\033[0m'
info()  { echo -e "${BLUE}[INFO]${NC}  $*"; }
ok()    { echo -e "${GREEN}[ OK ]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
die()   { echo -e "${RED}[FAIL]${NC}  $*" >&2; exit 1; }

# ── 1. Conda ─────────────────────────────────────────────────
# Prints the conda base directory: CONDA_BASE, the conda on PATH, the active
# environment's base, or a usual install location.
find_conda_base() {
    if [[ -n "${CONDA_BASE:-}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
        echo "${CONDA_BASE}"; return 0
    fi
    if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
        echo "$(dirname "$(dirname "${CONDA_EXE}")")"; return 0
    fi
    local exe; exe="$(command -v conda 2>/dev/null || true)"
    if [[ -n "${exe}" ]]; then
        echo "$(dirname "$(dirname "${exe}")")"; return 0
    fi
    if [[ -n "${CONDA_PREFIX:-}" ]]; then
        [[ -f "${CONDA_PREFIX}/etc/profile.d/conda.sh" ]] && { echo "${CONDA_PREFIX}"; return 0; }
        local parent="${CONDA_PREFIX%/envs/*}"
        [[ -f "${parent}/etc/profile.d/conda.sh" ]] && { echo "${parent}"; return 0; }
    fi
    local p
    for p in "${HOME}/miniforge3" "${HOME}/mambaforge" "${HOME}/miniconda3" \
             "${HOME}/anaconda3" /opt/conda /usr/local/conda; do
        [[ -f "${p}/etc/profile.d/conda.sh" ]] && { echo "${p}"; return 0; }
    done
    return 1
}

# Runs a conda command; if it fails, runs it once more with the classic solver.
conda_do() {
    conda "$@" && return 0
    warn "conda $1 failed; trying again with the classic solver"
    conda "$@" --solver=classic
}

# Conda's activation scripts read unset variables.
set +u
CONDA_FOUND=""
if _base="$(find_conda_base)"; then
    # shellcheck source=/dev/null
    source "${_base}/etc/profile.d/conda.sh"
    CONDA_FOUND=1
    ok "conda $(conda --version 2>&1 | awk '{print $2}')  (base: ${_base})"
elif [[ -n "${USE_CURRENT_ENV}" || -n "${CONDA_PREFIX:-}${VIRTUAL_ENV:-}" ]]; then
    USE_CURRENT_ENV=1
    warn "conda not found — building into the active interpreter instead"
else
    die "conda not found.
       If it lives outside PATH:  CONDA_BASE=/path/to/conda bash install.sh ${ENV_NAME}
       Or activate an environment yourself and re-run:
                                  bash install.sh --use-current-env"
fi

# ── 2. Environment ───────────────────────────────────────────
if [[ -n "${USE_CURRENT_ENV}" ]]; then
    [[ -n "${CONDA_PREFIX:-}${VIRTUAL_ENV:-}" ]] \
        || die "--use-current-env given but no environment is active. Activate one first."
    info "Using the active environment: ${CONDA_PREFIX:-${VIRTUAL_ENV}}"
elif [[ "${ENV_NAME}" == */* ]]; then
    if [[ -d "${ENV_NAME}" ]]; then
        info "Environment prefix '${ENV_NAME}' exists — reusing it."
    else
        info "Creating conda environment at '${ENV_NAME}' (Python ${PYTHON_VER})..."
        conda_do create -y -p "${ENV_NAME}" -c conda-forge python="${PYTHON_VER}" -q
    fi
    conda activate "${ENV_NAME}"
else
    if conda env list 2>/dev/null | awk '!/^#/{print $1}' | grep -qx "${ENV_NAME}"; then
        info "Environment '${ENV_NAME}' already exists — reusing it."
    else
        info "Creating conda environment '${ENV_NAME}' (Python ${PYTHON_VER})..."
        conda_do create -y -n "${ENV_NAME}" -c conda-forge python="${PYTHON_VER}" -q
    fi
    conda activate "${ENV_NAME}"
fi
ok "Python $(python --version 2>&1 | awk '{print $2}')  (env: ${CONDA_PREFIX:-${VIRTUAL_ENV}})"

# ── 3. Build tools ───────────────────────────────────────────
if [[ -n "${CONDA_FOUND}" ]]; then
    info "Installing build tools from conda-forge..."
    conda_do install -y -c conda-forge -q \
        compilers make cmake pkg-config \
        zlib bzip2 xz
    ok "Build tools ready"
else
    warn "Without conda, provide a C compiler, make, cmake and pkg-config"
    warn "yourself (e.g. via module load)."
fi

# ── 4. Rust ──────────────────────────────────────────────────
rust_version() { rustc --version 2>/dev/null | awk '{print $2}'; }
rust_ok() {
    command -v rustc &>/dev/null && command -v cargo &>/dev/null || return 1
    [[ "$(printf '%s\n%s\n' "${RUST_MIN}" "$(rust_version)" | sort -V | head -1)" == "${RUST_MIN}" ]]
}
if rust_ok; then
    ok "Rust $(rust_version)"
else
    command -v rustc &>/dev/null && warn "Rust $(rust_version) is older than ${RUST_MIN}; installing a newer toolchain"
    if [[ -n "${CONDA_FOUND}" ]]; then
        info "Installing Rust into the environment from conda-forge..."
        conda_do install -y -c conda-forge -q "rust>=${RUST_MIN}"
    else
        info "Installing Rust with rustup (no root needed)..."
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
            | sh -s -- -y --default-toolchain stable 2>&1
        # shellcheck source=/dev/null
        source "${CARGO_HOME:-$HOME/.cargo}/env"
    fi
    hash -r
    rust_ok || die "Rust >= ${RUST_MIN} is required. Install it with: conda install -c conda-forge rust (or https://rustup.rs)"
    ok "Rust $(rust_version)"
fi
set -u

# ── 5. Build and install ─────────────────────────────────────
pip install -q "maturin>=1.7,<2.0"
ok "maturin $(python -c 'import importlib.metadata; print(importlib.metadata.version("maturin"))' 2>/dev/null)"

cd "$SCRIPT_DIR"

# One Rust job per 4 GB of memory, at most one per core, at least two.
if [[ -f /proc/meminfo ]]; then
    _total_gb=$(awk '/MemTotal/{printf "%d", $2/1024/1024}' /proc/meminfo)
else
    _total_gb=$(sysctl -n hw.memsize 2>/dev/null | awk '{printf "%d", $1/1024/1024/1024}' || echo 16)
fi
_max_by_ram=$(( _total_gb / 4 ))
_nproc=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
_jobs=$(( _max_by_ram < _nproc ? _max_by_ram : _nproc ))
_jobs=$(( _jobs > 2 ? _jobs : 2 ))
export CARGO_BUILD_JOBS="$_jobs"
info "Rust parallel jobs: ${_jobs}  (${_total_gb} GB RAM, ${_nproc} cores)"

# htslib headers of the environment are moved aside during the build and put
# back when the script ends.
set +u
_htslib_inc="${CONDA_PREFIX:-}/include/htslib"
_hid=""
restore_htslib() {
    if [[ -n "${_hid}" && -d "${_htslib_inc}._hidden" ]]; then
        mv "${_htslib_inc}._hidden" "${_htslib_inc}"
        _hid=""
    fi
}
trap restore_htslib EXIT
if [[ -n "${CONDA_PREFIX:-}" && -d "${_htslib_inc}" ]]; then
    warn "Moving the environment's htslib headers aside during the build"
    mv "${_htslib_inc}" "${_htslib_inc}._hidden"
    _hid=1
fi
set -u

_build_log="$(mktemp -t isoforge_build.XXXXXX.log)"
_build_rc=0
info "Compiling isoforge_core (Rust); this takes a few minutes"
pip install . -v > "${_build_log}" 2>&1 &
_build_pid=$!
_t0=$SECONDS
# Progress: seconds elapsed and crates compiled so far.
while kill -0 "${_build_pid}" 2>/dev/null; do
    sleep 5
    if [[ -t 1 ]]; then
        printf '\r    %3ds  %s crates compiled' \
            "$(( SECONDS - _t0 ))" "$(grep -c '^ *Compiling ' "${_build_log}" 2>/dev/null || echo 0)"
    fi
done
wait "${_build_pid}" || _build_rc=$?
if [[ -t 1 ]]; then printf '\r%*s\r' 40 ''; fi
info "Build finished in $(( SECONDS - _t0 ))s"

restore_htslib

if [[ ${_build_rc} -ne 0 ]]; then
    tail -40 "${_build_log}" >&2
    die "Build failed; full log: ${_build_log}"
fi
grep -E '^Successfully (built|installed)' "${_build_log}" | tail -2
rm -f "${_build_log}"

ok "IsoForge installed"

# ── 6. Check ─────────────────────────────────────────────────
echo ""
info "Checking the installation..."
_fail=0
if command -v isoforge &>/dev/null; then
    ok "isoforge"
    for sub in annotate index post export; do
        if isoforge "$sub" --help &>/dev/null; then
            ok "isoforge $sub"
        else
            echo -e "${RED}[MISS]${NC}  isoforge $sub" >&2
            _fail=1
        fi
    done
else
    echo -e "${RED}[MISS]${NC}  isoforge" >&2
    _fail=1
fi

python -c "import isoforge_core; print(f'  isoforge_core {isoforge_core.version()}')" 2>/dev/null \
    && ok "isoforge_core importable" \
    || { echo -e "${RED}[MISS]${NC}  isoforge_core import failed" >&2; _fail=1; }

echo ""
if [[ $_fail -eq 0 ]]; then
    echo -e "${GREEN}══════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}  IsoForge v$(python -c 'import isoforge_core; print(isoforge_core.version())' 2>/dev/null) installed successfully!${NC}"
    echo -e "${GREEN}══════════════════════════════════════════════════════${NC}"
    echo ""
    if [[ -z "${USE_CURRENT_ENV}" ]]; then
        echo -e "  Activate:  ${YELLOW}conda activate ${ENV_NAME}${NC}"
        echo ""
    fi
    echo "  Quick start:"
    echo "    isoforge annotate --reads sample1.fastq.gz sample2.fastq.gz \\"
    echo "        --genome hg38.fa --annotation annotation.gtf --o results/project"
    echo ""
else
    die "Some components failed to install — see errors above."
fi
