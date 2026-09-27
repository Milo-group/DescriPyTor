#!/bin/bash
#$ -N esp_grid
#$ -S /bin/bash
#$ -cwd
#$ -q fairshare.q,k32.q
#$ -l h_vmem=16G
#$ -t 1-1
#$ -o logs/$JOB_NAME.$TASK_ID.out
#$ -e logs/$JOB_NAME.$TASK_ID.err

# Cylindrical ESP sampling on a list of structures, one array task per structure.
#
# The whole reason this runs in a private directory is esp.f:
#
#     grid_file = get_namespace('esp_coord')
#     inquire(file=grid_file,exist=ex)
#     if(.not.ex) call surfac(grid_file,n,xyz,at)
#
# surfac WRITES esp_coord when it is absent and esp.f READS it when present --
# same filename, both directions. A leftover file from a previous structure is
# therefore silently reused as the grid for the next one, with no error and no
# warning. That is almost certainly what grimme-lab/xtb#501 reported as
# order-dependent failures on ~10% of QM9, and on a shared scratch directory it
# does not crash: it returns the potential of this molecule sampled on some
# other molecule's surface.
#
# Two independent guards below: a per-task directory, and an explicit delete.
# The point-count check in the Python driver is the third.
#
# Usage:
#   mkdir -p logs
#   qsub -t 1-4477 M1_pre_calculations/esp_grid_job.sh structures.list

set -euo pipefail

STRUCTURE_LIST="${1:?usage: qsub esp_grid_job.sh <structure list>}"
PYTHON=/gpfs0/gaus/projects/miniconda3/envs/molfeatures/bin/python
XTB=/gpfs0/gaus/projects/xtb-6.4.1/bin/xtb
RESULTS_DIR="${PWD}/esp_results"

XYZ=$(sed -n "${SGE_TASK_ID}p" "${STRUCTURE_LIST}")
if [[ -z "${XYZ}" ]]; then
    echo "no entry at line ${SGE_TASK_ID} of ${STRUCTURE_LIST}" >&2
    exit 1
fi
XYZ=$(readlink -f "${XYZ}")
NAME=$(basename "${XYZ}" .xyz)

# Scratch if the node has it, task-private directory either way. The task id is
# in the path so two array tasks can never share an esp_coord.
SCRATCH_BASE="${TMPDIR:-/scratch/${USER}}"
[[ -d "${SCRATCH_BASE}" ]] || SCRATCH_BASE="${PWD}/scratch"
WORKDIR="${SCRATCH_BASE}/esp_${JOB_ID}_${SGE_TASK_ID}"
mkdir -p "${WORKDIR}" "${RESULTS_DIR}"
trap 'rm -rf "${WORKDIR}"' EXIT

cd "${WORKDIR}"
cp "${XYZ}" ./input.xyz

# Guard two: even in a fresh directory, refuse to inherit a grid we did not write.
rm -f esp_coord xtb_esp.dat xtb_esp_profile.dat xtb_esp.cosmo

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export XTBHOME=$(dirname $(dirname "${XTB}"))

"${PYTHON}" -m M1_pre_calculations.esp_grid_driver write \
    --xyz input.xyz --grid grid.pkl --coord esp_coord --framed-xyz framed.xyz

# xtb is given the SAME frame the grid was built in. The potential is a scalar
# field so the rotation changes nothing physical, and it removes the only way
# the grid and the structure could end up in different axes.
#
# esp_coord now exists, so esp.f reads it instead of calling surfac, and the
# potential is evaluated at exactly these points.
"${XTB}" framed.xyz --esp --namespace "${NAME}" > "${NAME}.xtb.out" 2>&1

"${PYTHON}" -m M1_pre_calculations.esp_grid_driver read \
    --grid grid.pkl --esp xtb_esp.dat --out "${NAME}.esp.csv" \
    --moments "${NAME}.moments.csv"

cp "${NAME}.esp.csv" "${NAME}.moments.csv" "${NAME}.xtb.out" "${RESULTS_DIR}/"
echo "done ${NAME}"
