#!/bin/bash
######## run_sim.sh ########
# シミュレーションのプログラムを実行する
# バイナリはすでに生成済み
#
# 実行例：/usr/local/src/master/sh/run_sim.sh 0.05 0.5
# 【引数】到着率，入力強度
#########################

SRC_PATH="/usr/local/src/master"
DAT_PATH="/usr/local/src/master/dat"

##### Handling whether argument exists. #############
if [ $# -ne 2 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには2個の引数が必要です。" 1>&2
  exit 1
fi
lambda=$1
alpha=$2
beta2=`echo "scale=7; 1 / ${lambda}" | bc`

##### Declare variables ################
params="l=${lambda}"
a="a=${alpha}"

#### Run #####################################
${SRC_PATH}/exe/simulation.exe ${lambda} ${beta2} ${alpha}

##### Make directory. #############
if [ ! -e ${DAT_PATH}/${params} ]; then
    mkdir ${DAT_PATH}/${params}
    mkdir ${DAT_PATH}/${params}/${a}
elif [ ! -e ${DAT_PATH}/${params}/${a} ]; then
    mkdir ${DAT_PATH}/${params}/${a}
else
    rm -rf ${DAT_PATH}/${params}/${a}/*
fi

############ Move directory ##############################
mv ./sim_* ${DAT_PATH}/${params}/${a}/
