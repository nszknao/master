#!/bin/bash
######## plot.sh ########
# datファイルからグラフをプロットする
# （シミュレーション）入力，変位応答，変位・速度のPDF，
# （解析）変位・速度のPDF，閾値通過率
#
# 実行例：/usr/local/src/master/sh/plot.sh 0.4 0.3 -1
# 【引数】1番目：ラムダ，2番目：アルファ，3番目：プロットする分布番号（-1のときは最小二乗法）
#########################

SRC_PATH="/usr/local/src/master"
RESULT_PATH="/usr/local/src/master/results"

##### Handling whether argument exists. #############
if [ $# -ne 3 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには3個の引数が必要です。" 1>&2
  exit 1
fi

lambda=$1
alpha=$2
if [ $3 -eq -1 ]; then
	num=""
else
	num=_$3
fi

##### Declare variables ################
params="l=${lambda}"
dat="dat_a=${alpha}"
eps="eps_a=${alpha}"

##### Handling whether directory exists. #############
if [ ! -e ${RESULT_PATH}/${params}/${dat} ]; then
	echo "指定されたラムダorアルファの結果が存在しません。" 1>&2
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
if [ -e "${RESULT_PATH}/${params}/${dat}/t_force.dat" ]; then
cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/excitation.eps"
set size 0.6,0.4
set xlabel "time {/Symbol-Oblique t}"
set ylabel "input"
p [0:50][-50:50] "${RESULT_PATH}/${params}/${dat}/t_force.dat" with lines lt 1 notitle
EOF
fi

# plot response distribution
if [ -e "${RESULT_PATH}/${params}/${dat}/t_x1.dat" ]; then
cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/response.eps"
set size 0.6,0.4
set xlabel "time {/Italic-Times t}"
set ylabel "displacement {/Italic-Times y_1}"
set ytics 1
p [0:1300][-5:5] "${RESULT_PATH}/${params}/${dat}/t_x1.dat" with lines lt 1 lw 3 lc rgb "red"
EOF
fi

# plot first-cross pdf
if [ -e "${RESULT_PATH}/${params}/${dat}/firstcross${num}.dat" ]; then
cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/firstcross${num}.eps"
set size 0.45,0.55
set xlabel "threshold {/Symbol x}"
set ylabel "up-clossing rate"
set xtics 2
set ytics 0.2
p [0:8][0:0.6] "${RESULT_PATH}/${params}/${dat}/firstcross${num}.dat" with l lt 1 lw 4 lc rgb "red"
EOF
fi

# plot distribution pdf
if [ -e "${RESULT_PATH}/${params}/${dat}/gsay1pdf${num}.dat" ]; then
cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/gsay1pdf${num}.eps"
set size 0.45,0.55
set xlabel "displacement {/Italic-Times y_1}"
set ylabel "probability density {/Italic-Times p_{Y_1}}"
set xtics 2
set ytics 0.2
p [-6:6][0.:1.0] "${RESULT_PATH}/${params}/${dat}/gsay1pdf${num}.dat" with l lt 1 lw 4 lc rgb "blue",\
					"${RESULT_PATH}/${params}/${dat}/y1_pdf.dat" with p pt 6 ps 1 lc rgb "red",\
					"${SRC_PATH}/dat/y1_Gpdf.dat" with lines lt 2 lw 4
EOF

cat <<EOF >>${plot}
set output "${RESULT_PATH}/${params}/${eps}/log_gsay1pdf${num}.eps"
set size 0.45,0.55
set xlabel "displacement {/Italic-Times y_1}"
set ylabel "probability density {/Italic-Times p_{Y_1}}"
set xtics 2
set ytics 0.1
set logscale y
set format y "10^{%L}"
p [-6:6][0.00001:1.0] "${RESULT_PATH}/${params}/${dat}/gsay1pdf${num}.dat" with l lt 1 lw 4 lc rgb "blue",\
						"${RESULT_PATH}/${params}/${dat}/y1_pdf.dat" with p pt 6 ps 1 lc rgb "red",\
						"${SRC_PATH}/dat/y1_Gpdf.dat" with lines lt 2 lw 4
EOF
fi

gnuplot ${plot}

rm ${plot}