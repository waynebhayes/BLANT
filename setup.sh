#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
USAGE="USAGE: source $0
PURPOSE: set up PATH and other things necessary for BLANT and its scripts to work seamlessly."

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }
which(){ echo "$PATH" | tr : "$NL" | awk '!seen[$0]{print}{++seen[$0]}' | while read d; do eval /bin/ls $d/$N; done 2>/dev/null | newlines; }

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
BIGTMP=`for i in /scratch/preserve/RaidZ3/tmp /var/tmp /scratch/preserve /var/tmp /tmp; do mkdir -p "$i/wayne" && (df $i | awk 'NR==1{for(av=1;av<=NF;av++)if(match($av,"[Aa]vail"))break;}NR>1{print $av,"'"$i"'"}'); done 2>/dev/null | sort -nr | awk 'NR==1{print $2}'`
[ "$MYTMP" ] || MYTMP="$BIGTMP/wayne"
TMPDIR=`mktemp -d $MYTMP/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

# if SETUP_USAGE have set, exit
if [ "`echo $0 | sed 's,.*/,,'`" = setup.sh ]; then
    die "You've run this as a script; source it instead by typing:
    source $0"
fi

MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
case $# in
1) export BLANT_HOME="$1";;
*) export BLANT_HOME="$MYDIR" ;;
esac
if [ "$BLANT_HOME" = . ]; then BLANT_HOME=`/bin/pwd`; fi

[ -d libwayne ] || die "need libwayne directory, please clone libwayne while in $BLANT_HOME"
[ -d libwayne/bin ] || die "need libwayne/bin directory, please clone libwayne while in $BLANT_HOME"
[ -x libwayne/bin/hawk ] || die "no hawk in libwayne/bin???"

export PATH="$BLANT_HOME/libwayne/bin:$BLANT_HOME/scripts:$PATH"
