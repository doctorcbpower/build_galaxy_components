#!/usr/bin/env bash
# Build ICs: a GalIC snapshot (halo + bulge + stellar disc) plus an analytic gas disc.
# The command-line version of build_gas_disc_ics.ipynb; both use galic_m1.toml.
#
#   ./build.sh GALIC_SNAPSHOT.hdf5 [OUT.hdf5]
#
# Environment: PARAMS (default galic_m1.toml), NDISC (override the particle number),
#              STRICT=1 (fail if the snapshot disagrees with the analytic model),
#              BACKGROUND_GRID=1 (add an AREPO background grid), GC (the command; default
#              "galaxy-components", or e.g. "python -m galaxy_components" without installing).
set -euo pipefail
cd "$(dirname "$0")"

SNAPSHOT=${1:?usage: ./build.sh GALIC_SNAPSHOT.hdf5 [OUT.hdf5]}
OUT=${2:-ic_with_gas.hdf5}
PARAMS=${PARAMS:-galic_m1.toml}
GC=${GC:-galaxy-components}
EXTRA=()
[[ -n "${NDISC:-}" ]] && EXTRA+=(--ndisc "$NDISC")

echo "== 1. model, rotation curve and vertical resolution"
$GC check "$PARAMS" ${EXTRA[@]+"${EXTRA[@]}"}

echo; echo "== 2. snapshot check, gas disc, merge"
FLAGS=()
[[ -n "${STRICT:-}" ]] && FLAGS+=(--strict)
[[ -n "${BACKGROUND_GRID:-}" ]] && FLAGS+=(--background-grid)
$GC add-gas "$PARAMS" --snapshot "$SNAPSHOT" --out "$OUT" ${EXTRA[@]+"${EXTRA[@]}"} ${FLAGS[@]+"${FLAGS[@]}"}
