#!/bin/bash
SRC_PATH="/usr/local/src/master"

for (( i = 0; i < 120; i++ )); do
	sh ${SRC_PATH}/sh/plot.sh 0.2 0.1 ${i}
	echo "sh ${SRC_PATH}/sh/plot.sh 0.2 0.1 ${i}"
	if [ "$?" -eq 1 ]; then
		echo "エラーが発生しました。"
		exit 1
	fi
done
