#!/bin/sh

src=$1
output_dir=$2
ymax=$3
ymin=$4
end=`expr $5 - 1`
nmol=$6

for i in `seq 0 $end`
do
    output=${output_dir}/`printf "enthalpy%03g.eps" $i`
    tmp=$output_dir/tmp
    grep "rep $i " $src | sed -e "s:rep\ $i::g" > $tmp

    gnuplot << EOF
set term post enh color eps

set yrange [${ymin}:${ymax}]

set format x "%4.0f"

set ylabel "Enthalpy [kJ/mol]"
set xlabel "Time [ns]"

set output "$output"
p "$tmp" every 100 u (\$1/1000):((\$2+\$3*\$6+\$4+\$7+\$8)/${nmol}) w l not
EOF

    rm $tmp
done


