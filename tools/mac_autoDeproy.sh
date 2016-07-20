#!/bin/bash
# mac_autoDeploy.sh
# ソースファイルが置いてあるディレクトリで実行する
TARGET_HOST=192.168.33.10
SRCPATH=/usr/local/src/master
DESTPATH=/Users/Naoto/git/master
PEMPATH=~/.ssh/chrp_l

# 同期する
rsync -rlptvu -e "ssh -i ${PEMPATH}" --delete ${DESTPATH}/server/ vagrant@${TARGET_HOST}:${SRCPATH}/
# 権限を与える
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo chmod 777 -R ${SRCPATH}/
# シンボリックリンクを貼る
ssh -i ${PEMPATH} vagrant@${TARGET_HOST} sudo ln -s ${SRCPATH}/web/common ${SRCPATH}/web/api/common