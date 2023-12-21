#!/bin/bash
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ $# -eq 5 ] || die "need 5 args, BLANT-exec SampleMeth N k correct.txt [correct files for 3e9 samples in regression-tests/sanity]"
BLANT=$1
S=$2
N=$3
k=$4
correct=$5
"$BLANT" -R -s $S -mf -n $N -k $k networks/syeast.el |
    paste "$correct" - | 
    awk "BEGIN{N=$N}"'
	function ABS(x){return x<0?-x:x}
	$1>30 && $3>30 { # rule of thum: minumum number of samples for graphlet count to be reliable
	    p=$1/3e9; # estimated PDF of this graphlet according to 3e9 samples
	    l=N*p; # expected count from BLANT for N samples
	    c=$3; # actual count from BLANT for N samples
	    d=ABS(l-c);
	    if(l==0) {
	      if(d) printf "%g\t%g %g %s\n", d, l, c, $0
	    } else if(d/sqrt(l)>0) printf "%g\t%g %g %s\n", d/sqrt(l), l, c, $0
	}'
