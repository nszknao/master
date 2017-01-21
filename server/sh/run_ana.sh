#!/bin/bash
######## run_ana.sh ########
# 解析のプログラムを実行する
# バイナリはすでに生成済み
#
# 実行例：/usr/local/src/master/sh/run_ana.sh 0.05 0.5 0.43124
# 【引数】到着率，入力強度，初期値
#########################

SRC_PATH="/usr/local/src/master"
RESULT_PATH="${SRC_PATH}/results"

##### Handling whether argument exists. #############
if [ $# -ne 2 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには2個の引数が必要です。" 1>&2
  exit 1
fi
lambda=$1
alpha=$2
beta2=`echo "scale=7; 1 / ${lambda}" | bc`
initValue=$3

##### Declare variables ################
params="l=${lambda}"
dat="dat_a=${alpha}"
fig="fig_a=${alpha}"

#### Run #####################################
${SRC_PATH}/exe/analysis.exe ${lambda} ${beta2} ${alpha} ${initValue}

##### Make directory. #############
if [ ! -e ${RESULT_PATH}/${params} ]; then
    mkdir ${RESULT_PATH}/${params}
    mkdir ${RESULT_PATH}/${params}/${dat}
    mkdir ${RESULT_PATH}/${params}/${fig}
elif [ ! -e ${RESULT_PATH}/${params}/${dat} ]; then
    mkdir ${RESULT_PATH}/${params}/${dat}
    mkdir ${RESULT_PATH}/${params}/${fig}
else
    rm -rf ${RESULT_PATH}/${params}/${dat}/*
    rm -rf ${RESULT_PATH}/${params}/${fig}/*
fi

############ Move directory ##############################
mv ./ana_* ${RESULT_PATH}/${params}/${dat}/
