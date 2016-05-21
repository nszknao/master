#!/bin/bash
# mac_autoDeploy.sh
# ソースファイルが置いてあるディレクトリで実行する
export TARGET_HOST=192.168.33.10
export SRCPATH=/usr/local/src/master
export PEMPATH=/c/Users/N.Nishizaka/.ssh/vagrant
export DESTPATH=/c/Users/N.Nishizaka/git/Master

rsync -rlptvu -e "ssh -i ${PEMPATH}" --delete ${DESTPATH}/* vagrant@${TARGET_HOST}:${SRCPATH}/
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo chmod 777 -R ${SRCPATH}