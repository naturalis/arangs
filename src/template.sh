#!/bin/bash

# all shell scripts (and perl, ruby, python, etc) get $0, which
# stores the quoted path to the executable being run
usage=$0" arg1 arg2"

# shell scripts give you commandline arguments in numeric order from 1 onward
arg1=$1

# -z tests for zero bytes in a string
if [ -z $arg1 ]
then
  echo $usage
  exit 1
fi

arg2=$2
if [ -z $arg2 ]
then
  echo $usage
  exit 1
fi

date "+running on %Y-%m-%d at %H:%M:%S"

# this might come in handy later on
running_date=`date "+%Y-%m-%d_%H_%M_%S"`
echo $running_date

# do things with arg1 (and arg2, etc)
# for each command that you run, check its exit status to make sure it ran ok
if [ 0 != $? ]
then
  echo "WARNING ${?}"
  exit 1
fi
