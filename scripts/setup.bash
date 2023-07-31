#!/bin/bash
SETUP_USAGE="USAGE:

    source setup.bash [BLANT Repo Directory]

This script does not change any files, but sets the environment variables BLANT_HOME and
LIBWAYNE_HOME appropriately--the former which is usually just the same directory this
file resides. If you give it an argument, it'll use that directory instead. Otherwise it'll
assume that it's own directory is the SpArcFiRe directory, and all it'll do is echo that
directory.
Note: the 'source' command assumes you're using bash as your shell; if not, you'll need
to figure out how to set the environment variable BLANT_HOME yourself.
"
fail() { echo "$@" >&2; return 1
}

if echo $0 | fgrep setup.; then
    fail "$SETUP_USAGE
    SETUP ERROR: You've run this script; source it instead by typing:
    source $0"
    exit 1
fi

[ "$BLANT_HOME" != "" ] && echo "already set" && return

NL='
' # newline as a variable, for ease later

# Source of following line: https://stackoverflow.com/questions/59895/how-to-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && cd .. && pwd )"
case $# in
1) export BLANT_HOME="$1";;
*) export BLANT_HOME="$MYDIR" ;;
esac
if [ "$BLANT_HOME" = . ]; then BLANT_HOME=`/bin/pwd`; fi

LIBWAYNE_HOME=$BLANT_HOME/libwayne
PATH="$LIBWAYNE_HOME/bin:$BLANT_HOME/scripts:$PATH"
export BLANT_HOME LIBWAYNE_HOME PATH

echo "finished setting BLANT_HOME LIBWAYNE_HOME PATH, they are"
echo "$BLANT_HOME $LIBWAYNE_HOME $PATH"
