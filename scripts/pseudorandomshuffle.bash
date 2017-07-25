#!/bin/bash
#
# do a seeded psuedorandom shuffle of the lines (numbers) of a file
file="$1"
seed="$2"

# obtain the line (numbers)
awk '{print NR}' "$file" > "$file.lines"

# shuffle the lines according to seed (cf. gnu.org)
shuf "$file.lines" --random-source=<(openssl enc -aes-256-ctr \
    -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null) \
    > "$file.lines.shuf$seed"
