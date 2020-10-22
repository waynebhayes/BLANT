#!/bin/bash
die() { echo "FATAL ERROR: $@" >&2; exit 1
}
warn() { echo "WARNING: $@" >&2;
}
PATH=`pwd`:`pwd`/scripts:$PATH
export PATH

EXE=./blant
MAKE=false
while [ $# -gt -0 ]; do
    case "$1" in
    -make) MAKE=true; shift;;
    *) [ -x "$1" -o "$MAKE" = true ] || die "unknown argument '$1'; valid arguments are '-make', and an optional executable filename"
	EXE="$1"; shift;;
    esac
done

export EXE
export CORES=${CORES:=`cpus 2>/dev/null | awk '{c2=int($1/2); print c2?c2:1}'`}
[ "$CORES" -gt 0 ] || die "can't figure out how many cores this machine has"
echo "Using $CORES as number of cores"
if $MAKE ; then
	make realclean
	make all || die "failed to make"
fi

[ -x "$EXE" ] || die "no executable '$EXE' exists to test!"

NUM_FAILS=0
for REG_DIR in regression-tests/*; do
    if [ -d "$REG_DIR" ]; then
	export REG_DIR
	echo --- in directory $REG_DIR ---
	for r in $REG_DIR/*.sh; do
	    echo --- running test $r ---
	    if "$r"; then
		:
	    else
		(( NUM_FAILS+=$? ))
	    fi
	done
    fi
done
echo Number of failures: $NUM_FAILS
exit $NUM_FAILS
