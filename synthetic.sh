#!/bin/bash -x
# Near the start of a run, we need smaller samples (say 1,000) and shorter stagnation (eg, 100).
# As things improve, we can up the samples to 100,000 and stagnation to eg 10,000.
# Basic idea seems to be this: we initially want VERY few samples, so we can see which graphlets
# are most crazy frequent, so that we initially will be optimizing adding some of those crazy
# frequent graphlets (with so few samples, all the ones with fewer frequencies will just be
# zero and thus ignored by the objective function).  Then we slowly increase the number of
# samples in order to start seeing the less frequent graphlets in the target network, and
# slowly start optimizing for those as well.  As we SLOWLY increase the sample size, it's
# probably good enough to get to 1 million samples.  And you KEEP the same sample size
# until you stop seeing improvement in the objective at that sample size.
die() { echo "$@" >&2; exit 1
}
TMPDIR=/tmp/synthetic$$
#trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15
mkdir $TMPDIR

k=$1
INPUT="$2"
baseI=`basename $INPUT`
OUTPUT="$3"
baseO=`basename $OUTPUT`
STAG=1024
sampleSize=102400
[ -f "$OUTPUT" ] || die "output graph $OUTPUT must already exist, and will be modified"
BLANTFILES_IN=''
BLANTFILES_OUT=''
for j in `seq 3 $k`; do
    BLANTFILES_IN="$BLANTFILES_IN $TMPDIR/blant.$baseI.k$j"
    BLANTFILES_OUT="$BLANTFILES_OUT $TMPDIR/blant.$baseO.k$j"
    ./blant -s NBE -k $j -mi -n $sampleSize $INPUT > $TMPDIR/blant.$baseI.k$j
    ./blant -s NBE -k $j -mi -n $sampleSize $OUTPUT > $TMPDIR/blant.$baseO.k$j
done
score=`for j in $(seq 3 $k); do paste <(awk '{print $1}' $TMPDIR/blant.$baseI.k$j | sort -n | uniq -c) <(awk '{print $1}' $TMPDIR/blant.$baseO.k$j | sort -n | uniq -c); done | awk '{sum2+=($3-$1)^2}END{print int(sqrt(sum2))}'`
while [ $STAG -le 1048576 ]; do # 1 Mebisample
    make synthetic
    md5sum $OUTPUT
    ./synthetic -s $STAG -k $j $INPUT $OUTPUT $BLANTFILES_IN $BLANTFILES_OUT > $TMPDIR/x || exit
    mv $TMPDIR/x $OUTPUT
    BLANTFILES_IN=''
    BLANTFILES_OUT=''
    for j in `seq 3 $k`; do
    BLANTFILES_IN="$BLANTFILES_IN $TMPDIR/blant.$baseI.k$j"
    BLANTFILES_OUT="$BLANTFILES_OUT $TMPDIR/blant.$baseO.k$j"
    ./blant -s NBE -k $j -mi -n $sampleSize $OUTPUT > $TMPDIR/blant.$baseO.k$j
    ./blant -s NBE -k $j -mi -n $sampleSize $INPUT > $TMPDIR/blant.$baseI.k$j
    done
    newScore=`paste <(awk '{print $1}' $TMPDIR/blant.$baseI.k? | sort -n | uniq -c) <(awk '{print $1}' $TMPDIR/blant.$baseO.k? | sort -n | uniq -c) | awk '{sum2+=($3-$1)^2}END{print int(sqrt(sum2))}'`
    if [ $newScore -gt $score ]; then # improvement has stopped at this sample size.
    STAG=`expr $STAG \* 4`
    echo increasing STAG to $STAG
    score=$(($newScore*4))
    for j in `seq 3 $k`; do
        ./blant -s NBE -k $j -mi -n $sampleSize $OUTPUT > $TMPDIR/blant.$baseO.k$j
        ./blant -s NBE -k $j -mi -n $sampleSize $INPUT > $TMPDIR/blant.$baseI.k$j
    done
    else
    score="$newScore"
    fi
done
