#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage: "$0" input_dir output_dir nsite"
    exit
fi

input_dir=$1
output_dir=$2
nsite=$3

echo "input_dir= "$input_dir
echo "output_dir= "$output_dir
echo "nsite= "$nsite

files=$(ls $input_dir/*.chk)

for input in $files
do
    output=$(echo $input | sed -e "s:$input_dir:$output_dir:g")
    ./change_nsite.out $nsite < $input > $output
done
