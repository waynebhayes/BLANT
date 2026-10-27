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

# without -D the tool really runs BLANT: the command line, quoting, stdout
# redirect and exit status. The space in the path must survive the shell.
cp "$FIX" "my net.el"
"$ROOT/blant-clusters" -m3 -r42 "$ROOT/test/fixture/fake-blant" 5 0.5 "my net.el" > out4.txt
diff out4.txt "$ROOT/test/fixture/expected.txt"

# a failing BLANT is reported and the tool exits non-zero
if "$ROOT/blant-clusters" -m3 "$ROOT/test/fixture/fake-blant-fail" 5 0.5 "$FIX" \
        > /dev/null 2> err2.txt; then
    echo "FAIL: BLANT exit 3 gave exit 0" >&2; exit 1
fi
grep -q "BLANT exited 3" err2.txt

# a weighted run can be checked too
mkdir wdir && cp "$ROOT/test/fixture/blant5.weighted.out" wdir/blant5.out
"$ROOT/blant-clusters" -D wdir -w -m3 -r42 ./blant 5 0.5 "$ROOT/test/fixture/fixture.elw" > out5.txt
"$ROOT/test/check" "$ROOT/test/fixture/fixture.elw" out5.txt 0.5 0.5 0.5

# the checker reads a cluster line longer than 64 KB: a 300-clique with
# 250-character names
awk 'BEGIN{n=300; pad=sprintf("%240s",""); gsub(/ /,"x",pad)
    for(i=1;i<=n;i++) name[i]=pad i
    for(i=1;i<=n;i++) for(j=i+1;j<=n;j++) print name[i], name[j] > "big.el"
    printf "%d nodes, %d/%d edges, 100%% density from k5, nodeSet {", n, n*(n-1)/2, n*(n-1)/2
    for(i=1;i<=n;i++) printf " %s", name[i]; print " }"}' > big.txt
[ "$(wc -c < big.txt)" -gt 65536 ]
"$ROOT/test/check" big.el big.txt 1 0.5 0.5

echo "e2e ok"
