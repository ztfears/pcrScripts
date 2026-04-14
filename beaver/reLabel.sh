#!/bin/bash


# To relabel a portion of your genomes to include more species under your "spec" group
# such as when dealing with complexes
# Usage: ./reLabel.sh <dir> <species>

dir=$1
species=$2


for file in $dir/gen*; do
        [ -f "$file" ] || continue;
        if head -n 1 "$file" | grep -q "$species"; then
                mv "$file" "${file/gen/spec}";
        fi;
done
