#!/bin/bash
# mac_autoDeploy.sh
# ソースファイルが置いてあるディレクトリで実行する
export TARGET_HOST=192.168.33.10
export SRCPATH=/usr/local/src/master
export DESTPATH=/Users/Naoto/git/master
export PEMPATH=~/.ssh/chrp_l

# 同期する
rsync -rlptvu -e "ssh -i ${PEMPATH}" --delete ${DESTPATH}/* vagrant@${TARGET_HOST}:${SRCPATH}/
# 権限を与える
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo chmod 777 -R ${SRCPATH}/
# シンボリックリンクを貼る
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo ln -s ${SRCPATH}/web/common ${SRCPATH}/web/api/common