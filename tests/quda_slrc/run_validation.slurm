#!/bin/bash
#SBATCH --job-name=quda_slrc_validate
#SBATCH --output=quda_slrc_validate_%j.out
#SBATCH --error=quda_slrc_validate_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --mem=32G

# ================================================================
# QUDA SLRC Solver Validation Job
# ================================================================
#
# This script runs the QUDA SLRC validation test which:
#   1. Computes a reference propagator using BiCGStab (double precision)
#   2. Computes the same propagator using QUDA SLRC direct solver (double precision)
#   3. Computes hadron correlators from both propagators
#   4. Compares the correlator outputs
#
# BEFORE RUNNING:
#   1. Edit the variables below to match your environment
#   2. Ensure the gauge configuration file exists
#   3. Update validate_quda_slrc.xml with your lattice parameters
#
# ================================================================

# ---- User Configuration ----

# Path to chroma executable
CHROMA_BIN="/path/to/chroma"

# Path to gauge configuration
CONFIG_FILE="/path/to/your/config.lime"

# Lattice geometry (must match nrow in XML)
LATT_SIZE="32 32 32 64"

# MPI geometry (must multiply to SLURM ntasks)
GEOM="1 1 1 1"

# QUDA tuning cache directory (avoids re-autotuning)
QUDA_TUNECACHE="${HOME}/quda_tunecache"

# ---- End User Configuration ----

# ================================================================
# Environment Setup
# ================================================================

# Load required modules (adapt to your cluster)
# module load cuda/12.0
# module load openmpi/4.1
# module load chroma

# QUDA environment
export QUDA_RESOURCE_PATH="${QUDA_TUNECACHE}"
mkdir -p "${QUDA_RESOURCE_PATH}"
export QUDA_ENABLE_TUNING=1
export QUDA_TUNE_VERSION_CHECK=0

# QDP-JIT GPU memory pool (leave room for QUDA)
export QDP_JIT_GPU_POOL_SIZE=4000

# Work directory
WORKDIR="${SLURM_SUBMIT_DIR}"
cd "${WORKDIR}"

# ================================================================
# Patch XML with actual config file path
# ================================================================

INPUT_XML="validate_quda_slrc.xml"
RUNTIME_XML="validate_quda_slrc_runtime.xml"

sed "s|/path/to/your/config.lime|${CONFIG_FILE}|g" "${INPUT_XML}" > "${RUNTIME_XML}"

echo "================================================================"
echo "QUDA SLRC Validation Test"
echo "================================================================"
echo "Date:        $(date)"
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Chroma:      ${CHROMA_BIN}"
echo "Config:      ${CONFIG_FILE}"
echo "Lattice:     ${LATT_SIZE}"
echo "MPI Geom:    ${GEOM}"
echo "GPU:         ${CUDA_VISIBLE_DEVICES:-auto}"
echo "QUDA cache:  ${QUDA_RESOURCE_PATH}"
echo "================================================================"

# ================================================================
# Run Chroma
# ================================================================

echo ""
echo "Starting Chroma at $(date)..."
echo ""

srun ${CHROMA_BIN} -i "${RUNTIME_XML}" -geom ${GEOM}

CHROMA_EXIT=$?

echo ""
echo "Chroma finished at $(date) with exit code ${CHROMA_EXIT}"
echo ""

if [ ${CHROMA_EXIT} -ne 0 ]; then
    echo "ERROR: Chroma exited with non-zero status. Check error log."
    exit ${CHROMA_EXIT}
fi

# ================================================================
# Extract and Check Chroma-Side Residual (Level 1)
# ================================================================

echo "================================================================"
echo "Level 1: Chroma-Side Residual Check"
echo "================================================================"

# Extract the Chroma-side residual from the log
grep "QUDA_SLRC_SOLVER.*Chroma.*resid" quda_slrc_validate_${SLURM_JOB_ID}.out || \
    echo "WARNING: Could not find Chroma-side residual in log output"

echo ""

# ================================================================
# Compare Correlators (Level 3)
# ================================================================

echo "================================================================"
echo "Level 3: Correlator Comparison"
echo "================================================================"

if [ -f "ref_hadspec.xml" ] && [ -f "quda_hadspec.xml" ]; then
    echo "Running correlator comparison..."
    echo ""

    # Use the companion Python script
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    python3 "${SCRIPT_DIR}/compare_correlators.py" ref_hadspec.xml quda_hadspec.xml

    COMPARE_EXIT=$?

    if [ ${COMPARE_EXIT} -eq 0 ]; then
        echo ""
        echo "Correlator comparison PASSED."
    else
        echo ""
        echo "Correlator comparison FAILED or had warnings (exit code ${COMPARE_EXIT})."
    fi
else
    echo "WARNING: Correlator output files not found."
    echo "  Expected: ref_hadspec.xml and quda_hadspec.xml"
    echo "  Check Chroma log for errors in HADRON_SPECTRUM measurements."
fi

echo ""
echo "================================================================"
echo "Validation Complete at $(date)"
echo "================================================================"

# Cleanup runtime XML
rm -f "${RUNTIME_XML}"
