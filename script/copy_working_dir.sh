#!/bin/sh

if [ $# != 2 ]
then
    echo "usage: $0 [src_dir] [dist_dir]"
    exit 1
fi

SRCDIR=$1
DSTDIR=$2

mkdir -p $DSTDIR/output/cdv $DSTDIR/input
cp $SRCDIR/input/moltype.txt $SRCDIR/output/*chk $SRCDIR/output/phys_value.dat $DSTDIR/input/.
sed -e "s:$SRCDIR:$DSTDIR:g" $SRCDIR/input/gpu.inp > $DSTDIR/input/gpu.inp
