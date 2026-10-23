#!/bin/bash
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
TMPD=$(mktemp -d); trap 'rm -rf "$TMPD"' EXIT

# The binary reads canon maps from ./canon_maps, like the shell script does.
mkdir -p "$TMPD/canon_maps"
cp "$ROOT/test/fixture/canon_list5.txt" "$ROOT/test/fixture/orbit_map5.txt" "$TMPD/canon_maps/"
cd "$TMPD"

FIX="$ROOT/test/fixture/fixture.el"
run() { "$ROOT/blant-clusters" -D "$ROOT/test/fixture" "$@" ./blant 5 0.5 "$FIX"; }

run -m3 -r42 > out1.txt
"$ROOT/test/check" "$FIX" out1.txt 0.5 0.5 0.5

# determinism: same seed, same bytes
run -m3 -r42 > out2.txt
diff out1.txt out2.txt

# golden file
diff out1.txt "$ROOT/test/fixture/expected.txt"

# -nm drops the nodeSets
run -m3 -r42 -nm > out3.txt
if grep -q nodeSet out3.txt; then
    echo "FAIL: -nm still printed nodeSet" >&2; exit 1
fi

# a missing -D file is a named error and a non-zero exit
if "$ROOT/blant-clusters" -D /nonexistent -m3 ./blant 5 0.5 "$FIX" \
        > /dev/null 2> err.txt; then
    echo "FAIL: missing -D file exited 0" >&2; exit 1
fi
grep -q "blant5.out" err.txt

echo "e2e ok"
