#!/bin/bash
SRC_PATH="/usr/local/src/master"

##### Handling whether argument exists. #############
if [ $# -ne 2 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには2個の引数が必要です。" 1>&2
  exit 1
fi

lambda=$1
alpha=$2
num=$3

##### Plot all. #############
for (( i = 0; i < ${num}; i++ )); do
	sh ${SRC_PATH}/sh/plot.sh ${lambda} ${alpha} ${i}
	echo "sh ${SRC_PATH}/sh/plot.sh 0.2 0.1 ${i}"
	if [ "$?" -eq 1 ]; then
		echo "エラーが発生しました。"
		exit 1
	fi
done
