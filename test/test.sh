#!/bin/sh

echo "This test was confirmed using GeForce GTX 980 Ti."

prefix=`pwd`
input_dir=$prefix/input
output_dir=$prefix/output

#echo $input_dir
#echo $output_dir

mkdir -p $output_dir/cdv
sed -e s:INPUT_DIR:$input_dir/:g $input_dir/gpu.inp | sed -e s:OUTPUT_DIR:$output_dir/:g > tmp.inp

#sed -e s:CURRENT_DIR:$prefix:g correct > tmp

#../bin/md.out < tmp.inp > output/log
../src/md.out < tmp.inp > output/log

#diff output/log tmp

grep rep output/log > tmp
diff correct tmp

if [ $? -eq 0 ];
then
    echo "Test successfully finished!"
else
    echo "Test failed..."
fi

rm tmp tmp.inp
