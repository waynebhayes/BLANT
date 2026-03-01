#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
USAGE="USAGE: source $0
PURPOSE: set up PATH and other things necessary for BLANT and its scripts to work seamlessly."

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }
export TAB='	'; export NL='
'
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

[ -d $BLANT_HOME/libwayne ] || die "need libwayne directory, please clone libwayne while in $BLANT_HOME"
[ -d $BLANT_HOME/libwayne/bin ] || die "need libwayne/bin directory, please clone libwayne while in $BLANT_HOME"
[ -x $BLANT_HOME/libwayne/bin/hawk ] || die "no hawk in libwayne/bin???"

#tty -s && echo Setting PATH appropriately...
# the back-quote command removes duplicate directories in the PATH
export PATH=`echo "$BLANT_HOME/libwayne/bin:$BLANT_HOME/scripts:$PATH" | tr : "$NL" | awk '!seen[$0]{print}{++seen[$0]}' | tr "$NL" :`
