#!/bin/bash
# mac_autoDeploy.sh
# ソースファイルが置いてあるディレクトリで実行する
# 本番環境(scdp_c1とscdp_c2)にbatch関係をアップロード
export TARGET_HOST=192.168.33.10
export SRCPATH=/usr/local/src/master
export DESTPATH=/Users/Naoto/git/Master

rsync -rlptvu --delete ${DESTPATH}/ga/* vagrant@${TARGET_HOST}:${SRCPATH}/ga