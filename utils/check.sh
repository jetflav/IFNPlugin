#!/bin/bash
#
# usage:
#   ./check.sh <program> <datafile> [args]
#
# this script checks if "<program>" exists, runs it using <datafile>
# as an input, and compares the output with "<program>.ref" (it fails if
# "<program>.ref" does not exist),
#
# Note that <program> is actually run using "./<program> [args] < <datafile>"

# check if tput is present before setting colors, just to be on the safe side
if [ x`which tput` != "x" ]; then
   GREEN=$(tput setaf 2)
   RED=$(tput setaf 1)
   NORMAL=$(tput sgr0)
fi

function print_status_and_exit {
  if [ x"$2" == x"Success" ]; then
     col=$GREEN
     exit_code=0
  else
     col=$RED
     exit_code=1
  fi      
  OUTPUT="${col}$2${NORMAL}"

  if [[ x"$FJCONTRIB_SUBMAKE" == x ]]; then
    # if this is a standalone make call then just give the output and 
    # exit with the relevant exit code;
    echo "$OUTPUT"  
    exit $exit_code
  else
    # if this is is a subdirectory make, then discard the exit code,
    # and record success/failure in a file that will be output from 
    # the parent makefile
    printf "  %-20s %-25s %s\n" `pwd | sed 's/.*\///g'` "$1" "$OUTPUT" >> ../test_summary.tmp
    exit
  fi
}

echo
echo -n "In directory" `pwd | sed 's/.*\///'`
echo "/ checking the output of '$1' using the input file"
echo "$2 and reference file $1.ref"

# check that all the necessary files are in place
#
# Note that, for the executable, the Makefile should in principle
# rebuild the program if it is absent or not executable
test -e ./$1 || { echo "ERROR: the executable $1 cannot be found."; print_status_and_exit "$1" "Failed (executable not found)";}
test -x ./$1 || { echo "ERROR: $1 is not executable."; print_status_and_exit "$1" "Failed (program not executable)";}
test -e ./$1.ref || { echo "ERROR: the expected output $1.ref cannot be found."; print_status_and_exit "$1" "Failed (reference output not found)";}
test -e ./$2 || { echo "ERROR: the datafile $2 cannot be found."; print_status_and_exit "$1" "Failed (datafile not found)";}

# make sure that the reference output file is not empty (after removal of
# comments and empty lines)
[ -z "$(cat ./$1.ref | grep -v '^#' | grep -v '^$')" ] && {
    echo "ERROR: the reference output, $1.ref"
    echo "should contain more than comments and empty lines"
    echo
    print_status_and_exit "$1" "Failed (no valid reference)"
}

# run the example
#./$1 < $2 2>/dev/null | grep -v "^#" > $1.tmp_ref
command=$1
infile=$2
shift;shift
./$command $* < $infile 2>$command.tmp_err | grep -v "^#" > $command.tmp_ref

DIFF=`cat $command.ref | grep -v "^#" | diff $command.tmp_ref -`
if [[ -n $DIFF ]]; then 
    cat $command.ref | grep -v "^#" | diff $command.tmp_ref - > $command.diff
    echo "ERROR: Outputs differ (output, stderr and diff available in "
    echo "       $command.tmp_ref , $command.tmp_err and $command.diff)"
    echo
    #rm $command.tmp_ref
    print_status_and_exit "$command" "Failed (see difference file)"
fi

rm $command.tmp_err
rm $command.tmp_ref

print_status_and_exit "$command" "Success"
