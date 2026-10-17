#!/bin/bash
# Regression test: the clean code (src/krgk_creep.cpp) must give the same critical force and the
# same creep events as the original code (original/vmc, creep driver -o, with the missing-return
# fix of original/patches/vmc-return.patch) for the same disorder seed.
# Needs Blitz++ and FFTW3 for the original. Usage: tests/regression.sh [seconds for the original, default 60]
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SECS=${1:-60}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT

# original: copy, patch, build without graphics
cp "$ROOT"/original/vmc/*.cpp "$ROOT"/original/vmc/*.hpp "$W"/
(cd "$W" && patch -s -p1 --directory=. < <(sed 's#[ab]/vmc/#./#' "$ROOT"/original/patches/vmc-return.patch))
g++ -O2 -DOPTIMIZE2 -DOPTIMIZE3 -DEXACT -DOPTIMIZE5 -DINSTANTONIC_WALK \
	"$W"/actions.cpp "$W"/mainalbe.cpp "$W"/minimize.cpp "$W"/mt.cpp "$W"/system.cpp "$W"/vmc.cpp \
	-lblitz -lfftw3 -lm -o "$W"/VMC 2>/dev/null
make -s -C "$ROOT" >/dev/null

# the original reads "M Mcut L disorder F/Fc Ndis c beta shift movie" on stdin and the seed from idum.dat
mkdir "$W"/orig "$W"/new
echo 12345 > "$W"/orig/idum.dat
(cd "$W"/orig && timeout "$SECS" ../VMC -o <<< "64 1 16 1.0 0.08 1 2.0 -1.0 0.0 0" > /dev/null || true)
n=$(($(wc -l < "$W"/orig/events.dat) - 1))
[ "$n" -ge 5 ] || { echo "FAIL: the original produced only $n events in $SECS s"; exit 1; }
"$ROOT"/build/krgk_creep --rng mt --L 16 --M 64 --mcut 1 --disorder 1.0 --c 2.0 --f 0.08 --seed 12345 --events "$n" --out "$W"/new

# original events.dat: L_rel area zlow zhigh barrier L_opt (6 significant digits)
# new events.dat:      step barrier L_opt L_rel area mean_u E expanded
python3 - "$W" "$n" <<'EOF'
import sys
W, n = sys.argv[1], int(sys.argv[2])
fc_o = float(open(f"{W}/orig/fc.dat").read().split()[0])
fc_n = float([l.split()[1] for l in open(f"{W}/new/params.txt") if l.startswith("Fc ")][0])
o = [l.split() for l in open(f"{W}/orig/events.dat").readlines()[1:] if l.strip()]
m = [l.split() for l in open(f"{W}/new/events.dat") if not l.startswith("#")]
bad = 0
for k, (a, b) in enumerate(zip(o, m), 1):
    same = abs(float(a[4]) - float(b[1])) <= 1e-5 * max(1.0, abs(float(a[4]))) and (a[5], a[0], a[1]) == (b[2], b[3], b[4])
    if not same:
        bad += 1
        print(f"event {k}: original barrier {a[4]} L_opt {a[5]} L_rel {a[0]} area {a[1]} | new {b[1]} {b[2]} {b[3]} {b[4]}")
ok = bad == 0 and abs(fc_o - fc_n) <= 1e-5 * fc_o
print(f"F_c original {fc_o} new {fc_n}; {len(o)} events compared, {bad} different: {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
EOF
