#!/bin/sh

if [ $# != 3 ]
then
    echo "usage: $0 [item] [value] [src_file]"
    exit 1
fi

ITEM=$1
VALUE=$2
FILE=$3

tmpfile=tmp`date +"%s"`
sed -e "s:$ITEM[ \t]\+[a-zA-Z0-9\/_.]\+:$ITEM $VALUE:g" $FILE > $tmpfile
mv $tmpfile $FILE
 
