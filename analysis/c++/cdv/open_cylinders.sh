#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage: $0 prefix input_dir output_dir"
    exit 0
fi

prefix=$1
input_dir=$2
output_dir=$3

files=$(ls ${input_dir}/${prefix}*.cdv)

for file in $files;
do
    output=$(echo $file | sed -e "s:$input_dir:$output_dir:g")
    echo $output
    ./define_inner.out  < $file > ${output}
    output=$(echo $output | sed -e "s:$prefix:open:g")
    ./open_cylinder.out < $file > ${output}
done
