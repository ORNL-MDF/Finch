#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd "${script_dir}/../.." && pwd)

case_name=${1:-full_physics}
mpi_ranks=${2:-${FINCH_MPI_RANKS:-1}}
application=${3:-${FINCH_EXECUTABLE:-${repository_root}/build/install/bin/finch}}
input_file="${script_dir}/inputs_${case_name}.json"
scan_path_file="${script_dir}/scan_path.txt"
run_root=${FINCH_EXAMPLE_OUTPUT:-${script_dir}/output}
run_dir="${run_root}/${case_name}"

mkdir -p "${run_dir}"
cp "${input_file}" "${run_dir}/inputs.json"
cp "${scan_path_file}" "${run_dir}/$(basename "${scan_path_file}")"
cp "${script_dir}/tabulated_profile.txt" "${run_dir}/tabulated_profile.txt"

read -r -a mpi_extra_args <<< "${FINCH_MPI_ARGS:-}"

echo "Running ${case_name} with ${mpi_ranks} MPI rank(s) in ${run_dir}"
cd "${run_dir}"
"${MPIEXEC:-mpirun}" "${mpi_extra_args[@]}" -np "${mpi_ranks}" \
    "${application}" -i inputs.json

if [[ "${FINCH_COMBINE_SOLIDIFICATION:-0}" == "1" ]]; then
    "${repository_root}/utilities/combine_solidification_data.sh" \
        -i solidification -o solidification.csv
fi
