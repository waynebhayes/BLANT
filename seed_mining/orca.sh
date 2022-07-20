#!/bin/sh
USAGE="$0 k net.el"
die() { echo "$@" >&2; exit 1
}
echo "$1" | egrep '^[45]$' || die "k must be 4 or 5; USAGE: $USAGE"
[ $# -eq 2 ] || die "USAGE: $USAGE"
TMPDIR=/tmp/orca.$$
trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15
mkdir $TMPDIR
k=$1; shift
for NET
do
    # FUCKING Orca needs node IDs to be fucking integers 0 through n-1. Morons!
    # But internally to awk we need to use IDs 1 through n so that id[0] is not zero, which is interpreted as false by awk
    awk '!id[$1]{name[++n]=$1;id[$1]=n}!id[$2]{name[++n]=$2;id[$2]=n}{m++;edge[$1][$2]=1}
	END{for(i=1;i<=n;i++)print name[i] >"'$TMPDIR'/names.txt";
	print n,m;
	for(u in edge)for(v in edge[u])print id[u]-1,id[v]-1}' "$NET" > $TMPDIR/orca.in

    cat -v $TMPDIR/names.txt

    orca.g $k $TMPDIR/orca.in | awk 'ARGIND==1{name[FNR]=$1}ARGIND==2{print name[FNR], $0}' $TMPDIR/names.txt -
done
