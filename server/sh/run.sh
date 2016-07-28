#!/bin/bash
SRC_PATH="/usr/local/src/master"
RESULT_PATH="/usr/local/src/master/results"

############ Compile-force #############################################
g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/sMain.cpp ${SRC_PATH}/src/common.cpp ${SRC_PATH}/src/research.cpp -o ${SRC_PATH}/simulation.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread
g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/aMain.cpp ${SRC_PATH}/src/analysis.cpp ${SRC_PATH}/src/expfit.cpp ${SRC_PATH}/src/nsga2.cpp ${SRC_PATH}/src/common.cpp -o ${SRC_PATH}/analysis.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread
# g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ${SRC_PATH}/src/main.cpp ${SRC_PATH}/src/analysis.cpp ${SRC_PATH}/src/expfit.cpp ${SRC_PATH}/src/nsga2.cpp ${SRC_PATH}/src/common.cpp ${SRC_PATH}/src/research.cpp -o ${SRC_PATH}/both -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread

##### Make directory. #############
if [ ! -e ${RESULT_PATH} ]; then
	mkdir ${RESULT_PATH}
else
	rm -rf ${RESULT_PATH}/*
fi

############ Reading init_value.txt ####################################
initfile=init_value_disp.txt
count=0
cat "${SRC_PATH}/sh/${initfile}" | grep -v ^# | while read line; do
	lambda=`echo ${line} | cut -d ',' -f1`
	alpha=`echo ${line} | cut -d ',' -f2`
	lininit=`echo ${line} | cut -d ',' -f3`
	loginit=`echo ${line} | cut -d ',' -f4`

	beta2=`echo "scale=7; 1 / ${lambda}" | bc`
	echo "lambda = ${lambda}, a=${alpha}"

##### Declare variables ################
	params="l=${lambda}"
	dat="dat_a=${alpha}"

#### Run #####################################
	${SRC_PATH}/simulation.exe ${lambda} ${beta2} ${alpha}
	${SRC_PATH}/analysis.exe ${lambda} ${beta2} ${alpha} ${loginit}
	# ${SRC_PATH}/both ${lambda} ${beta2} ${alpha} ${loginit}

##### Make directory. #############
	if [ ! -e ${RESULT_PATH}/${params} ]; then
		mkdir ${RESULT_PATH}/${params}
		mkdir ${RESULT_PATH}/${params}/${dat}
	else
		if [ ${count} -eq 0 ]; then
			rm -rf ${RESULT_PATH}/${params}/*
		fi
		mkdir ${RESULT_PATH}/${params}/${dat}
	fi

############ Move directory ##############################
	mv ./*.dat ${RESULT_PATH}/${params}/${dat}/

	count=$(( count + 1 ))
done
