#!/bin/bash
######## run.sh ########
# シミュレーションと解析のプログラムを実行する
#
# 実行例：/usr/local/src/master/sh/run.sh 1
# 【引数】実行コード（0:シミュレーション＋解析，1:シミュレーション，2:解析）
#########################

SRC_PATH="/usr/local/src/master"
RESULT_PATH="${SRC_PATH}/results"
DAT_PATH="${SRC_PATH}/dat"
SH_PATH="${SRC_PATH}/sh"

##### Handling whether argument exists. #############
if [ $# -ne 1 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには1個の引数が必要です。" 1>&2
  exit 1
fi
code=$1

############ Compile-force #############################################
if [ ${code} -eq 0 ]; then
	g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/common.cpp ${SRC_PATH}/src/simulation.cpp -o ${SRC_PATH}/simulation.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread
	g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/analysis.cpp ${SRC_PATH}/src/expfit.cpp ${SRC_PATH}/src/nsga2.cpp ${SRC_PATH}/src/common.cpp -o ${SRC_PATH}/analysis.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread
elif [ ${code} -eq 1 ]; then
	g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/common.cpp ${SRC_PATH}/src/simulation.cpp -o ${SRC_PATH}/simulation.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread
elif [ ${code} -eq 2 ]; then
	g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/analysis.cpp ${SRC_PATH}/src/expfit.cpp ${SRC_PATH}/src/nsga2.cpp ${SRC_PATH}/src/common.cpp -o ${SRC_PATH}/analysis.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread
fi

############ Reading init_value.txt ####################################
initfile=init_value_disp.txt
cat "${SRC_PATH}/sh/${initfile}" | grep -v ^# | while read line; do
	lambda=`echo ${line} | cut -d ',' -f1`
	alpha=`echo ${line} | cut -d ',' -f2`
	lininit=`echo ${line} | cut -d ',' -f3`
	loginit=`echo ${line} | cut -d ',' -f4`

	beta2=`echo "scale=7; 1 / ${lambda}" | bc`
	echo "lambda = ${lambda}, a=${alpha} is processing."

##### Declare variables ################
	params="l=${lambda}"
	dat="dat_a=${alpha}"
	a="a=${alpha}"

#### Run #####################################
	if [ ${code} -eq 0 ]; then
		${SH_PATH}/run_sim.sh ${lambda} ${alpha}
		${SH_PATH}/run_ana.sh ${lambda} ${alpha} ${loginit}
	elif [ ${code} -eq 1 ]; then
		${SH_PATH}/run_sim.sh ${lambda} ${alpha}
	elif [ ${code} -eq 2 ]; then
		${SH_PATH}/run_ana.sh ${lambda} ${alpha} ${loginit}
	fi

	echo "lambda = ${lambda}, a=${alpha} was processed."
done