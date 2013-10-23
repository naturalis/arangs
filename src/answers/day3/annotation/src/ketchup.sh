#!/bin/bash

usage=$0" <query> <database>"
QUERY=$1
DB=$2

if [ -z $QUERY ]
then
  echo "Usage: bash "$usage
  exit 1
fi

if [ -z $DB ]
then
  echo "Usage: bash "$usage
  exit 1
fi

RESULT="${QUERY}_out"

if [ ! -e $DB ]
	then
	echo "No database found, please add the reference .fasta file to the current directory."
	exit 1
fi

if [ ! -e $QUERY ]
	then
	echo "The query file does not exist. Mamma mia!"
	exit 1
fi

if [ ! -e "${DB}.nsd" ]
	then 
	echo "be patient: smashing tomatoes..."
	formatdb -i $DB -p F -o
fi

blastn -db $DB -query $QUERY -evalue 0.01 -out $RESULT

python src/ketchup.py $RESULT


