#!/bin/bash
###submit with './roary.sh prokkaDir outputDir'

set -e

date
echo "Running roary..."
echo "roary -p 20 -e -n *.gff"

roary -p 20 -g 150000 -n -f $2 $1/*.gff

