#!/bin/bash

############ Compile-force #############################################
g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ./sMain.cpp ../src/common.cpp ../src/research.cpp -o simulation.exe -std=c++11 -lm -lgsl -lgslcblas -lshark -lpthread

############ Reading init_value.txt ####################################
initfile=init_value_disp.txt
for line in `cat ${initfile} | grep -v ^#`
do
 lambda=`echo ${line} | cut -d ',' -f1`
 alpha=`echo ${line} | cut -d ',' -f2`

beta2=`echo "scale=7; 1 / ${lambda}" | bc`

echo "lambda = " $lambda "beta2 = " $beta2

##### Handling if "results" exists. #############
results="sResults"

if [ ! -e ./${results}/ ]
then
mkdir ./${results}/
fi

##### Handling if "results" exists. #############
params="l=${lambda},b2=${beta2}"

X=`echo "${alpha} * 10" | bc`
XX=${X%.*}
if [[ $XX -eq 1 ]]
then
if [ -e ./${results}/${params}/ ]
then
rm -Rf ./${results}/${params}/*
else
mkdir ./${results}/${params}/
fi
fi

##### Handlings if a directry "dat" exists. #####
dat="dat_a=${alpha}"

if [ ! -e ${dat} ]
then
mkdir ${dat}
fi

##### Handlings if a directry "eps" exists. ######
eps="eps_a=${alpha}"

if [ ! -e ${eps} ] 
then
mkdir ${eps}
mkdir ${eps}/"displacement(log)"
mkdir ${eps}/"displacement(linear)"
mkdir ${eps}/"velocity(log)"
mkdir ${eps}/"velocity(linear)"
fi

#### Run #####################################
./simulation.exe ${lambda} ${beta2} ${alpha}

##### Handling if "eps" exists in "results". #############
if [ ! -e ./${results}/${params}/"eps" ]
then
mkdir ./${results}/${params}/"eps"
fi

##### Handling if "dat" exists in "results". #############
if [ ! -e ./${results}/${params}/"dat" ]
then
mkdir ./${results}/${params}/"dat"
fi

############ Move directory ##############################
mv ./${eps} ./${results}/${params}/"eps"
mv ./*.dat ${dat}/
mv ./${dat} ./${results}/${params}/"dat"


echo -e "a=${alpha}, mu1_width=${loginit} was processed.\n\n"
done
