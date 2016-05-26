#!/bin/bash
# lab_autoDeploy.sh
# ソースファイルが置いてあるディレクトリで実行する
export TARGET_HOST=192.168.33.11
export SRCPATH=/usr/local/src/master
export PEMPATH=/c/Users/N.Nishizaka/.ssh/id_rsa
export DESTPATH=/c/Users/N.Nishizaka/git/Master

# サーバ上のデータを全消し
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo rm -rf ${SRCPATH}/*
# 必要なものだけをアップロード
scp -i ${PEMPATH} -r  ${DESTPATH}/ga vagrant@${TARGET_HOST}:${SRCPATH}/
scp -i ${PEMPATH} -r  ${DESTPATH}/ggd vagrant@${TARGET_HOST}:${SRCPATH}/
scp -i ${PEMPATH} -r  ${DESTPATH}/include vagrant@${TARGET_HOST}:${SRCPATH}/
scp -i ${PEMPATH} -r  ${DESTPATH}/siman vagrant@${TARGET_HOST}:${SRCPATH}/
scp -i ${PEMPATH} -r  ${DESTPATH}/src vagrant@${TARGET_HOST}:${SRCPATH}/
scp -i ${PEMPATH} -r  ${DESTPATH}/test vagrant@${TARGET_HOST}:${SRCPATH}/
#rsync -rlptvu -e "ssh -i ${PEMPATH}" --delete ${DESTPATH}/* vagrant@${TARGET_HOST}:${SRCPATH}/
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo chmod 777 -R ${SRCPATH}/