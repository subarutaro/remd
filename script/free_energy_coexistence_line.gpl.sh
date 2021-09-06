#!/bin/sh

N=$1
input_dir=$2
output_dir=$3

PD=$input_dir/phys_value.dat
DOS=$input_dir/dens_state.dat

echo "physvalue file is "$PD
echo "dos file is "$DOS

for i in `seq 0 9`
do
    gnuplot << EOF
set term post enh color eps 20

set table "/dev/null"

c0=-1e38
T=0
P=0
tmax(t,c)= (c > c0) ? T=t : 0
pmax(p,c)= (c > c0) ? P=p : 0
cmax(c)  = (c > c0) ? c0=c : 0
sp "$PD" i $i u (pmax(\$1,\$4)):(tmax(\$2,\$4)):(cmax(\$4))

n0 = 1e38
nmin(n) = (n < n0) ? n0 = n : 0

fe(n,e,v) = (n > n0) ? (-n + (e+P*v)*$N/T) : 1e38
g0 = 1e38
gmin(g) = (g < g0) ? g0 = g : 0

sp "$DOS" u 1:2:(nmin(\$3))
n0 = n0 + 1
sp "$DOS" u 1:2:(gmin(fe(\$3,\$2,\$1)))

unset table

height=10

gmax = g0 + height
gmin = g0

set autoscale x
set autoscale y
set zrange [gmin:gmax]

#set xtics 10
#set ytics 10
set grid

#set xrange [15:30]
#set yrange [-50:-30]

set xlabel "V / N"
set ylabel "E / N"

set contour
set cntrparam levels incremental gmin,2,gmax

pcoef=0.1/0.000060
tcoef=200/1.6621

set output sprintf("${output_dir}/P%.2fT%.2f.eps",P*pcoef,T*tcoef)
sp "$DOS" u 1:2:(fe(\$3,\$2,\$1)) ti sprintf("%.2f MPa, %.2f K",P*pcoef,T*tcoef) w l

EOF
done
