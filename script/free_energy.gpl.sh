#!/bin/sh

P=$1
T=$2
N=$3
input=$4
output=$5

#cat << EOF > tmp.gpl
gnuplot << EOF
set term post enh color eps 20

n0 = 1e38
nmin(n) = (n < n0) ? n0 = n : 0

fe(n,e,v) = (n > n0) ? (-n + (e+$P*v)*$N/$T) : 1e38
g0 = 1e38
gmin(g) = (g < g0) ? g0 = g : 0

set table "/dev/null"
sp "$input" u 1:2:(nmin(\$3))
n0 = n0 + 1
sp "$input" u 1:2:(gmin(fe(\$3,\$2,\$1)))
unset table

gmax = g0 + 12
gmin = g0

set autoscale x
set autoscale y
set zrange [gmin:gmax]

set xtics 10
set ytics 10
set grid

set xlabel "V / N"
set ylabel "E / N"


set contour
set cntrparam levels incremental gmin,4,gmax

set output "$output"
sp "$input" u 1:2:(fe(\$3,\$2,\$1)) w l not

EOF

