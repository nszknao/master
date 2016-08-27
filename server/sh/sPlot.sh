#!/bin/bash
RESULT_PATH="/usr/local/src/master/sResults"

##### Handling whether argument exists. #############
if [ $# -ne 2 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには2個の引数が必要です。" 1>&2
  exit 1
fi

lambda=$1
alpha=$2

##### Declare variables ################
params="l=${lambda}"
dat="dat_a=${alpha}"
eps="eps_a=${alpha}"

##### Handling whether directory exists. #############
if [ ! -e ${RESULT_PATH}/${params}/${dat} ]; then
	echo "指定されたアルファの結果が存在しません。" 1>&2
	exit 1
fi

##### Handlings if a directry "eps" exists. ######
if [ -e ${RESULT_PATH}/${params}/${eps} ]; then
	rm -rf ${RESULT_PATH}/${params}/${eps}/*
else
	mkdir ${RESULT_PATH}/${params}/${eps}
fi

##### Handling if "plot.plt" exists. #############
plot="plot.plt"

if [ -e ${plot} ]; then
	cp /dev/null ${plot}
else
	touch ${plot}
fi

# 初期設定
cat <<EOF >${plot}
set terminal postscript eps enhanced color
set xlabel font "Times New Roman"
set ylabel font "Times New Roman"
unset key
EOF

# plot exitation
if [ -e ${RESULT_PATH}/${params}/${dat}/"t_force.dat" ]; then
cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/excitation.eps"
set size 0.7,0.3
set xlabel "Time {/Italic-Times t}"
set xtics 100
set ytics 150
set ylabel "input"
p [0:500][-300:300] "${RESULT_PATH}/${params}/${dat}/t_force.dat" with lines lt 1 notitle
EOF
fi

# plot response distribution
if [ -e ${RESULT_PATH}/${params}/${dat}/"t_x1.dat" ]; then
cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/response.eps"
set size 0.7,0.3
set xlabel "Time {/Italic-Times t}"
set ylabel "displacement {/Italic-Times y_1}"
set xtics 100
set ytics 3
p [0:500][-6:6] "${RESULT_PATH}/${params}/${dat}/t_x1.dat" with lines lt 1 lw 3 lc rgb "red"
EOF
fi

# plot distribution pdf
if [ -e ${RESULT_PATH}/${params}/${dat}/"y1_pdf.dat" ]; then
cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/pdf_of_y1.eps"
set size 0.45,0.55
set xlabel "displacement {/Italic-Times y_1}"
set ylabel "probability density {/Italic-Times p_{Y_1}}"
set xtics 2
set ytics 0.2
p [-6:6][0.:1.0] "${RESULT_PATH}/${params}/${dat}/y1_pdf.dat" with p pt 6 ps 1
EOF

cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/log_pdf_of_y1.eps"
set size 0.45,0.55
set xlabel "displacement {/Italic-Times y_1}"
set ylabel "probability density {/Italic-Times p_{Y_1}}"
set xtics 2
set ytics 0.1
set logscale y
set format y "10^{%L}"
p [-6:6][0.00001:1.0] "${RESULT_PATH}/${params}/${dat}/y1_pdf.dat" with p pt 6 ps 1
EOF
fi

gnuplot ${plot}

rm ${plot}