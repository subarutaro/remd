#!/bin/sh

function replace_setting(){
    local ITEM=$1
    local VALUE=$2
    local FILE=$3

    local tmpfile=tmp`date +"%s"`
    sed -e "s:$ITEM[ \t]\+[a-zA-Z0-9\/_.]\+:$ITEM $VALUE:g" $FILE > $tmpfile
    mv $tmpfile $FILE
}

function copy_dir(){
    local SRCDIR=$1
    local DSTDIR=$2

    mkdir -p $DSTDIR/output/cdv $DSTDIR/input
    cp $SRCDIR/input/gpu.inp $SRCDIR/input/moltype.txt $SRCDIR/output/*chk $SRCDIR/output/phys_value.dat $DSTDIR/input/.
}

if [ $# -ne 1 ];then
    echo "usage: $0 [prefix] [server]"
    exit 1
fi

prefix=$1
current_run=`ls | grep ${prefix} | tail -n1 | sed -e "s:$prefix::"`
current_dir=`pwd`/$prefix$current_run

#echo $current_run $current_dir

nstep=`grep ninterval ${current_dir}/input/gpu.inp | sed 's/[\t ]\+/\t/g' | cut -f2`

model=`grep -v "#" $current_dir/input/moltype.txt | grep TIP | sed -e 's/[\t ]/\t/g' | cut -f2`
model=`grep -v \# $current_dir/input/moltype.txt | grep TIP | sed -e 's/[\t ]/\t/g' | cut -f2`
nmol=`grep nmol ${current_dir}/input/gpu.inp | sed 's/[\t ]\+/\t/g' | cut -f2`
diameter=`grep wall_length ${current_dir}/input/gpu.inp | sed 's/[\t ]\+/\t/g' | cut -f2`
diameter=`echo $diameter | awk '{printf("%3.2f",$1*2.0)}'`
press=`grep press_min ${current_dir}/input/gpu.inp | sed 's/[\t ]\+/\t/g' | cut -f2`
press=`echo $press | awk '{printf("%.1f",$1/1000000)}'`
gpu=`grep gpu_offset ${current_dir}/input/gpu.inp | sed 's/[\t ]\+/\t/g' | cut -f2`

echo "model="$model "nmol="$nmol "wall_length="$diameter "pressure="$press

nstep=`expr $nstep \* 2`
next_run=`expr $current_run + 1`
next_dir=`pwd`/${prefix}`printf %02d ${next_run}`

echo "Next directory is $next_dir, ninterval= $nstep"
copy_dir $current_dir $next_dir

input=$next_dir/input/gpu.inp
replace_setting input_prefix  $next_dir/input/  $input
replace_setting output_prefix $next_dir/output/ $input
replace_setting ninterval     $nstep            $input

cp `which md.out` $next_dir/${model}_N${nmol}_D${diameter}_P${press}_GPU$gpu.out
#$next_dir/.out < $next_dir/input/gpu.inp > $next_dir/output/log

