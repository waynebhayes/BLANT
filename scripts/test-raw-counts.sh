#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NUL=" "; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME BLANT-exec SampleMeth N k correct.txt [correct files for 3e9 samples in regression-tests/sanity]
PURPOSE: test the BLANT-exec for reasonable counts when run with SampleMeth with N samples for given k.
    Here 'reasonable' means testing the raw frequency of graphlets is within 6 sigma of the 'theoretical' Poisson frequency
    based on 3e9 samples, as stored in the directory regression-tests/sanity. We assume that, when NOT corrected for
    the sampling method's bias (ie., using the '-R' option which gives 'raw' sample counts), the sampled frequency of
    a graphlet is a Poisson process with PDF given by the samples in the stored 3e9 sample file.

    Output columns are:
	standard-deviations-off-expected
	expected-count-from-N-samples
	BLANT-exec-raw-count
	BLANT-exec-graphletID
	3e9-correct-count
	3e9-graphletID"

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" == skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

[ $# -ge 5 ] || die "expecting 5 required args and 1 optional arg (NUM_THREADS)"
BLANT=$1
S=$2
N=$3
k=$4
correct=$5
NUM_THREADS=${6:-1}

[ -x "$BLANT" ] || die "'$BLANT' is not executable or doesn't exist"
[ -f "$correct" ] || die "can't open correct file '$correct'"

"$BLANT" -t $NUM_THREADS -R -s $S -mf -Fc -n $N -k $k networks/syeast.el |
    paste - "$correct" |
    awk "BEGIN{N=$N;}"'
	function ABS(x){return x<0?-x:x}
	$2!=$4{print "ERROR: graphlet IDs must agree in columns 2 and 4:\n", $0 >"/dev/stderr"; exit 1}
	$1>30 && $3>30 { # rule of thumb: minumum number of samples for graphlet count to be reliable
	    p=$3/3e9; # estimated PDF of this graphlet according to 3e9 samples
	    l=N*p; # expected count from BLANT for N samples
	    B=$1; # actual count from BLANT for N samples
	    d=ABS(l-B);
	    if(l==0) {
	      if(d) printf "%g %g\t%s\n", d, l, $0
	    } else {
		sigmas=d/sqrt(l);
		# MCMC has high variance of on high ordinals
		if("'"$S"'"=="MCMC") sigmas /= 10; # MCMC has higher variance than Poisson... fudge factor
		if(sigmas>0) printf "%g %g\t%s\n", sigmas, l, $0
	    }
	}'
