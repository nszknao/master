#!/bin/bash
##### Handling if "plot.plt" exists. #############
plot="plot.plt"

if [ -e ${plot} ]
then
rm ${plot}
touch ${plot}
else
touch ${plot}
fi

# 変位応答分布（対数）
cat <<EOF >${plot}
set terminal postscript eps enhanced color
set style line 1 lt 1 lw 2
set size 0.45,0.55
unset key
EOF

for (( i = 0 ; i < 120 ; i++ ))
do
cat <<EOF >>${plot}
set output "first_cross_prob_${i}.eps"
set xlabel "threshold {/Symbol x}"
set ylabel "up-crossing rate"
set xtics 2
set ytics 0.1
p [0:8][0:0.6] "results/l=0.4,b2=14.2857142/dat/dat_a=0.1/firstcross_${i}.dat" with lines linetype 1 lw 4 lc rgb "red"
EOF
done
for (( i = 0 ; i < 120 ; i++ ))
do
cat <<EOF >>${plot}
set output "log_pdf_of_y1_${i}.eps"
set xlabel "displacement {/Italic-Times y_1}"
set ylabel "probability density {/Italic-Times p_{Y_1}}"
set xtics 2
set ytics 0.1
set logscale y
set format y "10^{%L}"
p [-6:6][0.00001:1.0] "results/l=0.07,b2=14.2857142/dat/dat_a=0.1/gsay1pdf_${i}.dat" with lines linetype 1 lw 4 lc rgb "red"
EOF
done
gnuplot ${plot}
