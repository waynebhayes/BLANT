#!/bin/bash
for i in canon_maps/alpha_list_MCMC8.txt canon_maps/alpha_list_NBE8.txt canon_maps/alpha_list_EBE8.txt canon_maps/canon_list8.txt canon_maps/canon_map8.bin canon_maps/canon-ordinal-to-signature8.txt canon_maps/orbit_map8.txt canon_maps/perm_map8.bin; do
    [ ! -f "$i" ] && echo "exiting since '$i' is missing" >&2 && exit 0
done

[ "$CI" = true ] && echo "assuming continuous integration based on CI = $CI" >&2 && exit 0
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
./blant -k8 -lDEG2 -sINDEX -mr -a1 networks/syeast0/syeast0.el >$TMPDIR/syeast0.index 2>/dev/null
./dedup.sh $TMPDIR/syeast0.index

echo Checking for index correctness
if cmp <(sort $TMPDIR/syeast0.index) <(sort seed_mining/examples/syeast0.index); then
    echo "Matches Patrick's original"
elif cmp <(sort $TMPDIR/syeast0.index) <(sort seed_mining/examples/syeast0.index2); then
    echo "Matches the HACK 2nd one Wayne added after Patrick's started failing..."
else
    echo "ERROR: generated index doesn't match either Patrick's or Wayne's" >&2
    trap "" 0 1 2 3 15
    exit 1
fi
