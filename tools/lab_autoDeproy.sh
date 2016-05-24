#!/bin/bash
# lab_autoDeploy.sh
# ソースファイルが置いてあるディレクトリで実行する
export TARGET_HOST=192.168.33.11
export SRCPATH=/usr/local/src/master
export PEMPATH=/c/Users/N.Nishizaka/.ssh/id_rsa
export DESTPATH=/c/Users/N.Nishizaka/git/Master

ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo rm -rf ${SRCPATH}/*
scp -i ${PEMPATH} -r  ${DESTPATH}/* vagrant@${TARGET_HOST}:${SRCPATH}/
#rsync -rlptvu -e "ssh -i ${PEMPATH}" --delete ${DESTPATH}/* vagrant@${TARGET_HOST}:${SRCPATH}/
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo chmod 777 -R ${SRCPATH}/