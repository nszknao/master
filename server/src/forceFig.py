# -*- coding: utf-8 -*-
u"""入力の図を作成します．"""
import matplotlib.pyplot as plt
import numpy as np
import sys

DAT_PATH	= '/usr/local/src/master/dat'
OUTPUT_PATH	= '/usr/local/src/master/fig'

def drange(begin, end, step):
	u"""小数刻みでrangeできるように拡張
	【引数】begin: 開始の値，end: 終了の値，step: ステップ幅"""
	n = begin
	while n+step < end:
		yield n
		n += step

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print('ERROR: Wrong number of arguments.')
		sys.exit()
	arg_lambda	= sys.argv[1]
	arg_alpha	= sys.argv[2]

	# ファイル読み込み
	forceData	= [[],[]]
	for line_a in open(DAT_PATH + '/l=' + str(arg_lambda) + '/a=' + str(arg_alpha) + '/sim_force.dat'):
		col	= line_a.strip().split(' ')
		forceData[0].append(float(col[0]))
		forceData[1].append(float(col[1]))
	# プロット
	plt.figure(figsize=(16, 5))
	plt.xlabel("Time", fontsize=18)
	plt.ylabel("Input", fontsize=18)
	plt.plot(forceData[0], forceData[1], lw=1, color='r')
	plt.xlim([0, 100])
	plt.ylim([-500, 500])
	plt.grid(which='major',color='black',linestyle='--')
	filename	= "force.png"
	plt.savefig(OUTPUT_PATH+"/"+filename)
	plt.clf()

	# ファイル読み込み
	forceGaussData	= [[],[]]
	for line_a in open(DAT_PATH + '/l=' + str(arg_lambda) + '/a=' + str(arg_alpha) + '/sim_force_gaussian.dat'):
		col	= line_a.strip().split(' ')
		forceGaussData[0].append(float(col[0]))
		forceGaussData[1].append(float(col[1]))
	# プロット
	plt.figure(figsize=(16, 5))
	plt.xlabel("Time", fontsize=18)
	plt.ylabel("Input_Gauusian", fontsize=18)
	plt.plot(forceGaussData[0], forceGaussData[1], lw=1, color='r')
	plt.xlim([0, 50])
	plt.ylim([-100, 100])
	plt.grid(which='major',color='black',linestyle='--')
	filename	= "force_gaussian.png"
	plt.savefig(OUTPUT_PATH+"/"+filename)
	plt.clf()

	# ファイル読み込み
	forcePulseData	= [[],[]]
	for line_a in open(DAT_PATH + '/l=' + str(arg_lambda) + '/a=' + str(arg_alpha) + '/sim_force_pulse.dat'):
		col	= line_a.strip().split(' ')
		forcePulseData[0].append(float(col[0]))
		forcePulseData[1].append(float(col[1]))
	# プロット
	plt.figure(figsize=(16, 5))
	plt.xlabel("Time", fontsize=18)
	plt.ylabel("Input_Pulse", fontsize=18)
	plt.plot(forcePulseData[0], forcePulseData[1], lw=1, color='r')
	plt.xlim([0, 100])
	plt.ylim([-300, 300])
	plt.grid(which='major',color='black',linestyle='--')
	filename	= "force_pulse.png"
	plt.savefig(OUTPUT_PATH+"/"+filename)
	plt.clf()
