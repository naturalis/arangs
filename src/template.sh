#!/bin/bash

# all shell scripts (and perl, ruby, python, etc) get $0, which
# stores the quoted path to the executable being run.
# use $0 with a list of arguments required.  You should change this to
# fit your scripts argument needs.  Change arg1 to be something more
# explanatory, e.g. path_to_reference_file
usage=$0" arg1 arg2"

# shell scripts give you commandline arguments in numeric order from 1 onward
arg1=$1

# -z tests for zero bytes in a string
if [ -z $arg1 ]
then
  # print usage to STDERR
  echo $usage 1>&2
  exit 1
fi

# you may only need one argument, or you may need more than 2.
# remove the next 7 lines, or copy them and increment the number
# to fit your scripts needs
arg2=$2
if [ -z $arg2 ]
then
  # print usage to STDERR
  echo $usage 1>&2
  exit 1
fi

# prints feedback to the user on STDERR
date "+running on %Y-%m-%d at %H:%M:%S" 1>&2

# this might come in handy later on. It creates a variable
# that holds a timestamp that could be used for, say, a file or
# directory name.  Change as needed %Y is YYYY, %m is MM, %d is DD, %H is 24hr, %M is minutes, %S is seconds
running_date=`date "+%Y-%m-%d_%H_%M_%S"`
echo $running_date 1>&2

# do things with arg1 (and arg2, etc)
# for each command that you run, check its exit status to make sure it ran ok
if [ 0 != $? ]
then
  echo "Problem encountered" 1>&2
  exit 1
fi

# exit with a zero exit status, to signal shell, or wrapper scripts that it ran to completion
exit
