#!/bin/sh

if [ $# != 6 ]
then
    echo "usage: $0 root_dir src_dir D_min dD  D_max node_name"
    exit 1
fi

ROOTDIR=$1
SRCDIR=$2
DIAMETERS=`seq $3 $4 $5`
NODE=$6

prev_dir=$SRCDIR

for d in $DIAMETERS
do
    radious=`echo $d | awk '{printf "%lf",0.5*$1}'`
    curr_dir=$ROOTDIR/$d/runeq
    echo " --- making $curr_dir ---"
    copy_working_dir.sh $prev_dir $curr_dir
    input=$curr_dir/input/gpu.inp
    replace_setting.sh wall_length $radious $input
    replace_setting.sh input_prefix $curr_dir/input/ $input
    replace_setting.sh output_prefix $curr_dir/output/ $input
    ssh $NODE md.out < $curr_dir/input/gpu.inp > $curr_dir/output/log
    prev_dir=$curr_dir
done
