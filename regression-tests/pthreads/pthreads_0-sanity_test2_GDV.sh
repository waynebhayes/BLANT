#!/bin/bash
echo 'testing Graphlet (not orbit) Degree Vectors [MULTITHREADED]'

# "correct" values came from 3e9 samples... but we no longer need to divide by 1000 since it's already absolute estimates

N=9000000

export k=1
echo "Skipping MCMC and SEC sampling since it cannot run in threads"
for S in NBE EBE; do
    case $S in
    MCMC) TOL=0.006; exp=2;;
    SEC)  TOL=1.1e-4; exp=3;;
    NBE)  TOL=1.1e-4; exp=3;;
    EBE)  TOL=1.5e-4; exp=3;;
    esac

    for k in 3 4 5 6 7 8
    do
	CORRECT=regression-tests/0-sanity/syeast.$S.gdv.abs.3e9.k$k.txt.xz
    for t in $CORES; do
        if [ -f canon_maps/canon_map$k.bin -a -f $CORRECT ]; then
            /bin/echo -n "$S:$k:$t "
            ./blant -q -R -s $S -mg -n $N -k $k -t $t networks/syeast.el |
            sort -n | cut -d' ' -f2- |
            paste - <(unxz < $CORRECT) |
            awk '{  cols=NF/2;
                for(c1=1;c1<=cols;c1++){
                    c2=cols+c1; if($c1&&$c2) print $c1/$c2
                }
                }' |
            $LIBWAYNE_HOME/bin/stats -g | # the -g option means "geometric mean"
            sed -e 's/#/num/' -e 's/[	 ][	 ]*/ /g' |
            $LIBWAYNE_HOME/bin/named-next-col '
                BEGIN{k='$k'; numBad=0}
                {
		    diff=ABS(1-mean)/(k*stdDev)^'$exp';
		    if(diff > '"$TOL"') {
			printf "BEYOND TOLERANCE: %g\n%s\n", diff, $0;
			numBad++;
		    } else
			printf "diff %.4e\t%s\n", diff, $0;
                }
		END{exit(numBad)}' || exit 1
        fi
    done || exit 1
    done || exit 1
done || exit 1
