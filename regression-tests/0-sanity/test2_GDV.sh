#!/bin/bash
echo 'testing Graphlet (not orbit) Degree Vectors'

# "correct" values came from 3e9 samples, so we divide by 1000 to correspond to testing with 3e6 samples

N=3000000
S=SEC
export k=1
for k in 3 4 5 6 7 8; do
    echo "trying $k..."
    if [ -f canon_maps/canon_map$k.bin ]; then
	/bin/echo -n "$k: "
	./blant -R -s $S -mg -n $N -k $k networks/syeast.el |
	    sort -n | cut -d' ' -f2- |
	    paste - <(unxz < regression-tests/0-sanity/syeast.$S.gdv.3e9.k$k.txt.xz) |
	    awk 'function MIN(a,b){return (a<b)?a:b}
		 function MAX(a,b){return (a>b)?a:b}
		 {
		    cols=NF/2;
		    for(i=1;i<=cols;i++) if($i>1000 && $(cols+i)/1000>1000)
			printf "%.9f\n", 1-MIN($i,$(cols+i)/1000)/MAX($i,$(cols+i)/1000)
		}' |
	    $LIBWAYNE_HOME/bin/stats | sed -e 's/#/num/' -e 's/var.*//' |
	    $LIBWAYNE_HOME/bin/named-next-col '
		{
		    if(num<900 || mean>.005*'$k' || max>0.2 || stdDev>0.005*'$k'){
			printf "BEYOND TOLERANCE:\n%s\n", $0;
			exit 1
		    } else
			print $0
		}' || break
    fi
done
[ $k -eq 8 ] || exit 1
