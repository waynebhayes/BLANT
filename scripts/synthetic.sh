#!/bin/bash
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
USAGE="USAGE $0 k targetNetwork.el currentSynthetic.el
PURPOSE: given a target network and a starting synthetic network, incrementally move the synthetic network towards
the target network, attempting to duplicate its distribution of graphlets from size 3 up to k.
    k - the maximum size of graphlets to use in the effort to duplicate the graphlet distribution
    target - the network who's graphlet distribution we're trying to replicate
    currentSynthetic - the network we're doing to modify to move it towards the target.
    
EXAMPLE USAGE: some example starting networks for syeast are in the directory SyntheticInit. DO NOT OVERWRITE THESE.
    Instead, make a copy:
	cp SyntheticInit/syeast.random.good syeast.synth # syeast.synth is now random, and will be modified below
    and then run the script on the copy:
	./synthetic.sh 5 networks/syeast.el syeast.synth # iteratively improve syeast.synth towards real syeast."

die() { (echo "$USAGE"; echo "FATAL ERROR: $@" ) >&2; exit 1
}
TMPDIR=/tmp/synthetic$$
#trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15
mkdir $TMPDIR
ulimit -s unlimited

[ $# -eq 3 ] || die "expecting exactly 3 argumets"
k=$1; [ $k -ge 3 -a $k -le 8 ] || die "first argument '$k' must be between 3 and 8"
TARGET="$2";  [ -f "$TARGET" ] || die "second argument '$TARGET' is the target network; it will not be touched"
baseT=`basename $TARGET`
OUTPUT="$2"; [ -f "$OUTPUT" ] || die "output graph $OUTPUT must already exist, and will be modified"
baseO=`basename $OUTPUT`

# The stagnation threshold, STAG, tells us how stringent to be when deciding when to move to the next meta-iteration.
# We start with a very lax threshold (1024), since the initial guess network is assumed to be random and thus far from
# the target, so there's no point in trying too hard. As we approach the target network, we make the stagnation threshold
# higher and higher, to push harder and harder towards the "correct" distribution.
STAG=1024

sampleSize=102400 # how many samples BLANT will take.


BLANTFILES_IN=''
BLANTFILES_OUT=''
echo -n "Creating distribution of target '$TARGET' for graphlets of size"
for j in `seq 3 $k`; do
    echo -n " $j"
    BLANTFILES_IN="$BLANTFILES_IN $TMPDIR/blant.$baseT.k$j"
    BLANTFILES_OUT="$BLANTFILES_OUT $TMPDIR/blant.$baseO.k$j"
    ./blant -s NBE -k $j -mi -n $sampleSize $TARGET > $TMPDIR/blant.$baseT.k$j
    ./blant -s NBE -k $j -mi -n $sampleSize $OUTPUT > $TMPDIR/blant.$baseO.k$j
done
echo " ...done"

score=`for j in $(seq 3 $k); do paste <(awk '{print $1}' $TMPDIR/blant.$baseT.k$j | sort -n | uniq -c) <(awk '{print $1}' $TMPDIR/blant.$baseO.k$j | sort -n | uniq -c); done | awk '{sum2+=($3-$1)^2}END{print int(sqrt(sum2))}'`
while [ $STAG -le 1048576 ]; do # 2^20 (~1 million) samples
    #make synthetic
    echo -n "Currently: STAG $STAG score $score md5sum "
    md5sum $OUTPUT
    ./synthetic -s $STAG -k $j $TARGET $OUTPUT $BLANTFILES_IN $BLANTFILES_OUT > $TMPDIR/x || exit
    mv $TMPDIR/x $OUTPUT
    #ount.el networks/syeast.el $OUTPUT; for j in 3 4 5 6 7 8; do echo $j: $(paste <(./blant -mf -s MCMC -n 800000 -k $j networks/syeast.el) <(./blant -mf -s MCMC -n 800000 -k $j $OUTPUT)|hawk '$2==$4{print $2,$1,$3}' | grep -v ' 0 0$' | awk '$3!=0&&$2!=0{print $3/$2}' | stats -g); done; pushkar-graph-compare.sh networks/syeast.el $OUTPUT
    BLANTFILES_IN=''
    BLANTFILES_OUT=''
    echo -n "Checking distribution of synthetic '$OUTPUT' for graphlets of size"
    for j in `seq 3 $k`; do
	echo -n " $j"
	BLANTFILES_IN="$BLANTFILES_IN $TMPDIR/blant.$baseT.k$j"
	BLANTFILES_OUT="$BLANTFILES_OUT $TMPDIR/blant.$baseO.k$j"
	./blant -s NBE -k $j -mi -n $sampleSize $OUTPUT > $TMPDIR/blant.$baseO.k$j
	./blant -s NBE -k $j -mi -n $sampleSize $TARGET > $TMPDIR/blant.$baseT.k$j
    done
    echo " ...done"
    newScore=`paste <(awk '{print $1}' $TMPDIR/blant.$baseT.k? | sort -n | uniq -c) <(awk '{print $1}' $TMPDIR/blant.$baseO.k? | sort -n | uniq -c) | awk '{sum2+=($3-$1)^2}END{print int(sqrt(sum2))}'`
    echo -n "New: STAG $STAG newScore $score"
    if [ $newScore -gt $score ]; then
	echo -n .... improvement has stopped at this sample size.
	STAG=`expr $STAG \* 2`
	echo increasing STAG to $STAG
	score=$(($newScore*2))
	for j in `seq 3 $k`; do
	    ./blant -s NBE -k $j -mi -n $sampleSize $OUTPUT > $TMPDIR/blant.$baseO.k$j
	    ./blant -s NBE -k $j -mi -n $sampleSize $TARGET > $TMPDIR/blant.$baseT.k$j
	done
    else
	echo .... improvement has occurred, keep STAG at $STAG
	score="$newScore"
    fi
done
