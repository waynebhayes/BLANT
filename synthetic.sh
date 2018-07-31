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
trap "/bin/rm -f x$$; exit" 0 1 2 3 15
k=5
STAG=1024
sampleSize=1024
./blant -k $k -mi -s $sampleSize syeast.random > blant.ryeast.k$k
./blant -k $k -mi -s $sampleSize networks/syeast.el > blant.syeast.k$k
score=`paste <(awk '{print $1}' blant.syeast.k$k | sort -n | uniq -c) <(awk '{print $1}' blant.ryeast.k$k | sort -n | uniq -c) | hawk '{sum2+=($3-$1)^2}END{print int(sqrt(sum2))}'`
while [ $sampleSize -le 1048576 ]; do # 1 Mebisample
    make synthetic
    md5sum syeast.random
    ./synthetic -s $STAG -k $k networks/syeast.el syeast.random blant.syeast.k$k blant.ryeast.k$k > x$$
    mv x$$ syeast.random
    ./blant -k $k -mi -s $sampleSize syeast.random > blant.ryeast.k$k
    ./blant -k $k -mi -s $sampleSize networks/syeast.el > blant.syeast.k$k
    newScore=`paste <(awk '{print $1}' blant.syeast.k$k | sort -n | uniq -c) <(awk '{print $1}' blant.ryeast.k$k | sort -n | uniq -c) | hawk '{sum2+=($3-$1)^2}END{print int(sqrt(sum2))}'`
    if [ $newScore -gt $score ]; then # improvement has stopped at this sample size.
	sampleSize=`expr $sampleSize \* 4`
	echo increasing sampleSize to $sampleSize
	score=`parse "$newScore*4"`
	./blant -k $k -mi -s $sampleSize syeast.random > blant.ryeast.k$k
	./blant -k $k -mi -s $sampleSize networks/syeast.el > blant.syeast.k$k
    else
	score="$newScore"
    fi
done
