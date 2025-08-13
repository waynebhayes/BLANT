#!/bin/bash 

export LANG=C # to ensure sort is not fucked up by stupid POSIX standards

if [ "$RUNNING_UNDER_TIME" != true ]; then
    export RUNNING_UNDER_TIME=true
    #echo "restarting with timing" >&2
    time $0 "$@"
fi

case "$1" in
-use-git-at)
    if [ -f git-at ] && [ `wc -l < git-at` -eq 2 -a `git log -1 --format=%at` -eq `tail -1 git-at` ]; then
	echo -n "Repo unchanged; returning same status code as "
	tail -1 git-at | xargs -I{} date -d @{} +%Y-%m-%d-%H:%M:%S
	exit `head -1 git-at`
    fi
    shift
    ;;
esac

if [ "X$GCC_VER" = "X" ] && hostname | egrep Waynes-Air; then
    echo "Assuming this is Wayne's MacBook Air, needing gcc-14"
    export GCC_VER=-14
fi

USAGE="USAGE: $0 [ -make ] [ -x BLANT_EXE ][ list of tests to run, defaults to regression-tests/*/*.sh ]"


# Bash Functions
die(){ (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $@}" </dev/null; }
cpus() {
    TMP=/tmp/cpus.$$
    trap "/bin/rm -f $TMP; exit" 0 1 2 3 15

    # Most Linux machines:
    lscpu >$TMP 2>/dev/null && awk '/^CPU[(s)]*:/{cpus=$NF}END{if(cpus)print cpus; else exit 1}' $TMP && return

    # MacOS:
    (uname -a | egrep 'Darwin|arm64' >/dev/null) && sysctl -n hw.ncpu && return

    # Cygwin:
    case "`uname -a`" in
    *CYGWIN*) grep -c '^processor[ 	]*:' /proc/cpuinfo; return ;;
    *) if [ -d /dev/cpu -a ! -f /dev/cpu/microcode ]; then
	ls -F /dev/cpu | fgrep -c
	return
       fi
	;;
    esac

    # Oops
    echo "couldn't figure out number of CPUs" >&2; exit 1
}

BLANT_HOME=`/bin/pwd`
LIBWAYNE_HOME="$BLANT_HOME/libwayne"
PATH="$BLANT_HOME:$BLANT_HOME/scripts:$BLANT_HOME/libwayne/bin:$PATH"
export PATH BLANT_HOME LIBWAYNE_HOME

if [ ! -f libwayne/Makefile ]; then
    echo "you need the submodule libwayne; trying to get it now" >&2
    (git submodule init libwayne && git submodule update libwayne && cd libwayne && git checkout master && git pull) || die "failed to get libwayne"
    [ -f libwayne/Makefile ] || die "Still can't find libwayne"
fi

EXE=$BLANT_HOME/blant
MAKE=false
while [ $# -gt -0 ]; do
    case "$1" in
    -make) MAKE=true; shift;;
    -x) EXE="$2"; shift 2;;
    -*) die "unknown option '$1";;
    *) break;;
    esac
done
[ -x "$EXE" -o "$MAKE" = true ] || die "Executable '$EXE' must exist or you must specify -make"

NUM_THREADS=$(getconf _NPROCESSORS_ONLN)
MAX_THREADS=4  # we don't want to hog all the threads when unit testing
CORES=$(( NUM_THREADS > MAX_THREADS ? MAX_THREADS : NUM_THREADS ))
#CORES=${CORES:=`cpus 2>/dev/null | awk '{c2=int($1/2); if(c2>0)print c2; else print 1}'`}
[ "$CORES" -gt 0 ] || die "can't figure out how many cores this machine has"
MAKE_CORES=1 # for BLANT, we don't want or need paralellism during make

WHAT=most
if $MAKE ; then
    make pristine
    WHAT=all
    #WHAT='"DEBUG=1" all'
fi

export EIGHT=8
if [ "$CI" = true ]; then # continuous integration needs to run faster
    echo '$CI'" variable is '$CI', so assuming we are doing continuous integration" >&2
    unset EIGHT
    export NO8=1
    export PAUSE=0
else
    : #echo '$CI'" = '$CI'; NOT doing continuous integration" >&2
fi

echo "Using $MAKE_CORES cores to make and $CORES cores for regression tests"
export EXE CORES MAKE_CORES

if [ "$NO8" != "" ]; then unset EIGHT; fi
make -j$MAKE_CORES $WHAT || die "failed to make"
# The gzip below is now done in the Makefile
#F=canon_maps/canon_map8.txt; [ -f $F ] && nice -19 gzip -9 $F &

[ -x "$EXE" ] || die "no executable '$EXE' exists to test!"

NUM_FAILS=0
STDBUF=''
if which stdbuf >/dev/null; then
    STDBUF='stdbuf -oL -eL'
fi
if [ $# -eq 0 ]; then
    set regression-tests/*/*.sh
fi
for r
do
    REG_DIR=`dirname "$r"`
    NEW_FAILS=0
    export REG_DIR
    echo --- running test $r ---
    if eval time $STDBUF "$r"; then # force output and error to be line buffered
	:
    else
	NEW_FAILS=$?
	(( NUM_FAILS+=$NEW_FAILS ))
	if echo "$r" | grep sanity; then
	    echo ""
	    echo "***********************************************************************************************"
	    echo "****** This test really needs to work... stopping early... please contact whayes@uci.edu"
	    echo "****** Please provide a screenshot of the above (or cut-and-paste into the message)"
	    echo "***********************************************************************************************"
	    echo ""
	    break
	fi
    fi
    echo --- test $r incurred $NEW_FAILS failures, cumulative failures is $NUM_FAILS ---
done
echo Total number of failures: $NUM_FAILS
(echo $NUM_FAILS; git log -1 --format=%at) > git-at
exit $NUM_FAILS
