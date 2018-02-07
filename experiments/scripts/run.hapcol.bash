#!/bin/bash
#
usage="usage: bash run.hapcol.bash version input output"
#
# run version of hapcol in all-heterozygous mode with decreasing
# alphas, starting with 0.01 (the default) and using decreasing orders
# of magnitude until a solution exists
#
default=2 # 1e-2 = 0.01 (the default)
#
# note: more precisely, we decrease alpha up to a limit (an extreme
# limit that no run should ever need/want), because having infinite
# loops is never a good idea in principle
#
limit=100 # 1e-100 (the limit)
#----------------------------------------------------------------------
#
version=$1
input=$2
output=$3

# run hapcol with decreasing alphas
echo "running hapcol version $version on $input"
for a in $(seq $default $limit) ; do

    echo " --- with alpha = 1e-$a ..."
    $version \
	-i $input \
	-o $output \
	-A -a 1e-$a

    # break should the above exit with 0, i.e., a solution exists
    if [[ $? == 0 ]] ; then
	break
    else
	rm -f $output
	echo " --- no solution at alpha = 1e-$a"
    fi

done
