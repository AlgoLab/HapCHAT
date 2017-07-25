#!/bin/bash
#
# extract a sample (set of line numbers) from a file.  Sample was
# shuffled according to seed and downsampled according to maxcov
file=$1
sample=$2
seed=$3
maxcov=$4

# extract the lines
awk '{printf "%.20d %s\n", NR, $0}' "$file" | join - \
    <(awk '{printf "%.20d\n", $1}' "$sample" | sort) | \
    sed 's/^[0-9]* //' \
    > "$file.downs.s$seed.m$maxcov.wif"
