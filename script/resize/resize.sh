#!/bin/sh

nmol=$1
dim=$2
input_dir=$3
output_dir=$4

for input in `ls $input_dir/*chk`
do
    output=`echo $input | sed -e "s:$input_dir:$output_dir:g"`
    ./resize.out $nmol $dim $input $output
done
