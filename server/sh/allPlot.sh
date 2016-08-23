#!/bin/bash
######## allPlot.sh ########
# plot.shを呼んで複数の分布を同時にプロットする
# （シミュレーション）入力，変位応答，変位・速度のPDF，
# （解析）変位・速度のPDF，閾値通過率
#
# 実行例：/usr/local/src/master/sh/allPlot.sh 0.05 0.7 120
# 【引数】1番目：ラムダ，2番目：アルファ，3番目：プロットする分布の数
#########################
SRC_PATH="/usr/local/src/master"

##### Handling whether argument exists. #############
if [ $# -ne 3 ]; then
  echo "指定された引数は$#個です。" 1>&2
  echo "実行するには3個の引数が必要です。" 1>&2
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
