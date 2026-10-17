#!/bin/bash
# Demo runs, then the figures. Results go to demo/runs/, figures to demo/figures/.
# About 2 minutes with 3 parallel jobs (JOBS=n to change); the slowest run is f = 0.15.
#  - one disorder realization, L = 32, disorder periodic in u with period M = 64, at several
#    drives f = F/F_c, until the line has advanced by 3M: the dynamics becomes a periodic
#    cycle of metastable states, from which the dominant configuration is taken.
#  - statistics of events in steady state: L = 64, f = 0.6, M large enough that the disorder
#    never repeats (M = 16384 > advance = 6000), three disorder realizations.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
make -s -C "$ROOT"
BIN="$ROOT"/build/krgk_creep
R="$ROOT"/demo/runs
mkdir -p "$R"

{
	for f in 0.15 0.25 0.4 0.6 0.9; do
		echo "--rng mt --L 32 --M 64 --seed 7 --f $f --advance 192 --ground --out $R/single_f$f"
	done
	for s in 1 2 3; do
		echo "--rng mt --L 64 --M 16384 --seed $s --f 0.6 --advance 6000 --out $R/stats_s$s"
	done
} | xargs -P "${JOBS:-3}" -I{} bash -c 'set -- {}; mkdir -p "${!#}"; "$0" --events 1000000 --save-events "$@"' "$BIN"

python3 "$ROOT"/demo/plot_demo.py "$R" "$ROOT"/demo/figures
