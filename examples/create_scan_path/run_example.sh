#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd "${script_dir}/../.." && pwd)
application=${1:-${FINCH_SCAN_PATH_EXECUTABLE:-${repository_root}/build/install/bin/create_scan_paths}}
run_dir=${FINCH_EXAMPLE_OUTPUT:-${script_dir}/output}

mkdir -p "${run_dir}"
cp "${script_dir}/inputs.json" "${run_dir}/inputs.json"

echo "Writing scan paths under ${run_dir}"
cd "${run_dir}"
"${application}" -i inputs.json
