#!/bin/bash
SRC_PATH="/usr/local/src/master"
RESULT_PATH="/usr/local/src/master/aResults"

############ Compile-force #############################################
g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/aMain.cpp ${SRC_PATH}/src/analysis.cpp ${SRC_PATH}/src/expfit.cpp ${SRC_PATH}/src/nsga2.cpp ${SRC_PATH}/src/common.cpp -o analysis.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread

############ Reading init_value.txt ####################################
initfile=init_value_disp.txt
for line in `cat "${SRC_PATH}/sh/${initfile}" | grep -v ^#`
do
	lambda=`echo ${line} | cut -d ',' -f1`
	alpha=`echo ${line} | cut -d ',' -f2`
	lininit=`echo ${line} | cut -d ',' -f3`
	loginit=`echo ${line} | cut -d ',' -f4`

	beta2=`echo "scale=7; 1 / ${lambda}" | bc`
	echo "lambda = ${lambda}, beta2 = ${beta2}"

##### Declare variables ################
	params="l=${lambda}"
	dat="dat_a=${alpha}"

#### Run #####################################
	./analysis.exe ${lambda} ${beta2} ${alpha} ${loginit}

##### Make directory. #############
	if [ ! -e ${RESULT_PATH} ]; then
		mkdir ${RESULT_PATH}
		mkdir ${RESULT_PATH}/${params}
		mkdir ${RESULT_PATH}/${params}/${dat}
	else
		rm -rf ${RESULT_PATH}/*
		mkdir ${RESULT_PATH}/${params}/
		mkdir ${RESULT_PATH}/${params}/${dat}
	fi

############ Move directory ##############################
	mv ./*.dat ${RESULT_PATH}/${params}/${dat}/

	echo -e "a=${alpha}, mu1_width=${loginit} was processed.\n\n"
done
