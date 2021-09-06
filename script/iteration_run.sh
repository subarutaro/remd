#!/bin/sh

if [ $# != 6 ]
then
    echo "usage: $0 src_dir D_min dD D_max ninterval_min hosts_file"
    exit 1
fi

DIAMETERS=`seq $2 $3 $4`
nstep=$5
NODES=(`cat $6`)
echo "assigned nodes are ${NODES[@]}"

TOTAL=0
for node in ${NODES[@]}
do
    ngpus=`ssh $node nvidia-smi -L | wc -l`
    TOTAL=`expr $TOTAL + $ngpus`
done

ROOTDIR=`pwd`

inode=0
node=${NODES[$inode]}
ngpus=`ssh $node nvidia-smi -L | wc -l`
igpu=0
ipid=0
npid=0

equilibration_run.sh $ROOTDIR $1 $2 $3 $4 ${NODES[0]}

curr=0
prev=eq

while true
do
    for d in $DIAMETERS
    do
	echo "==== making $ROOTDIR/$d/run${curr} for nstep=$nstep ===="
	echo " --- gpu$igpu of $node ---"
	if [ $npid -ge $TOTAL ]
	then
	    ps ${pid[$ipid]} > /dev/null
	    if [ $? = 0 ]
	    then
		echo "waiting ${pid[$ipid]}"
		wait ${pid[$ipid]}
	    fi
	    ipid=`expr $ipid + 1`
	fi
	curr_dir=$ROOTDIR/$d/run$curr
	prev_dir=$ROOTDIR/$d/run$prev
	copy_working_dir.sh $prev_dir $curr_dir

	input=$curr_dir/input/gpu.inp
	replace_setting.sh input_prefix $curr_dir/input/ $input
	replace_setting.sh output_prefix $curr_dir/output/ $input
	replace_setting.sh ninterval $nstep $input
	replace_setting.sh gpu_offset $igpu $input

	ssh $node md.out < $curr_dir/input/gpu.inp > $curr_dir/output/log &

	pid+=($!)
	npid=`expr $npid + 1`
	echo "pid is $!"
	igpu=`expr $igpu + 1`
	if [ $igpu = $ngpus ]
	then
	    igpu=0
	    inode=`expr $inode + 1`
	    if [ $inode = ${#NODES[@]} ]
	    then
		inode=0
	    fi
	    node=${NODES[inode]}
	    ngpus=`ssh $node nvidia-smi -L | wc -l`
	fi
    done
    prev=$curr
    curr=`expr $curr + 1`
    nstep=`expr $nstep \* 2`
done
