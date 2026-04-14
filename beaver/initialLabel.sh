#!/bin/bash


# To add either a gen or spec prefix to the beginning of your FNA files
# Usage: ./initialLabel.sh <dir> <species>

dir=$1
species=$2


for f in $dir/*.fna; do
	base=$(basename $f)
        if [[ "$(head -n 1 $f)" == *$species* ]]; then
                mv "$f" "$dir/spec$base";
        else
                mv "$f" "$dir/gen$base";
        fi;
done
