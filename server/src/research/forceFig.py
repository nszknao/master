# -*- coding: utf-8 -*-
u"""入力の図を作成します．"""
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import sys

DAT_PATH    = '/usr/local/src/master/dat'
OUTPUT_PATH = '/usr/local/src/master/fig'

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
    arg_lambda  = sys.argv[1]
    arg_alpha   = sys.argv[2]

    # ファイル読み込み
    forceData   = [[],[]]
    for line_a in open(DAT_PATH + '/l=' + str(arg_lambda) + '/a=' + str(arg_alpha) + '/sim_force.dat'):
        col = line_a.strip().split(' ')
        forceData[0].append(float(col[0]))
        forceData[1].append(float(col[1]))
    # プロット
    plt.figure(figsize=(16, 5))
    plt.xlabel("Time", fontsize=30)
    plt.ylabel("Input", fontsize=30)
    plt.plot(forceData[0], forceData[1], lw=1, color='r')
    plt.xlim([0, 100])
    plt.ylim([-500, 500])
    plt.subplots_adjust(left=0.1, bottom=0.2)
    plt.tick_params(labelsize=24)
    plt.gca().yaxis.set_major_locator(tick.MultipleLocator(250))
    plt.grid(which='major',color='black',linestyle='--')
    filename    = "force.eps"
    plt.savefig(OUTPUT_PATH+"/"+filename, format="eps", dpi=300)
    plt.clf()

    # ファイル読み込み
    forceGaussData  = [[],[]]
    for line_a in open(DAT_PATH + '/l=' + str(arg_lambda) + '/a=' + str(arg_alpha) + '/sim_force_gaussian.dat'):
        col = line_a.strip().split(' ')
        forceGaussData[0].append(float(col[0]))
        forceGaussData[1].append(float(col[1]))
    # プロット
    plt.figure(figsize=(16, 5))
    plt.xlabel("Time", fontsize=30)
    plt.ylabel("Input_Gauusian", fontsize=30)
    plt.plot(forceGaussData[0], forceGaussData[1], lw=1, color='r')
    plt.xlim([0, 50])
    plt.ylim([-100, 100])
    plt.subplots_adjust(left=0.1, bottom=0.2)
    plt.tick_params(labelsize=24)
    plt.grid(which='major',color='black',linestyle='--')
    filename    = "force_gaussian.eps"
    plt.savefig(OUTPUT_PATH+"/"+filename, format="eps", dpi=300)
    plt.clf()

    # ファイル読み込み
    forcePulseData  = [[],[]]
    for line_a in open(DAT_PATH + '/l=' + str(arg_lambda) + '/a=' + str(arg_alpha) + '/sim_force_pulse.dat'):
        col = line_a.strip().split(' ')
        forcePulseData[0].append(float(col[0]))
        forcePulseData[1].append(float(col[1]))
    # プロット
    plt.figure(figsize=(16, 5))
    plt.xlabel("Time", fontsize=30)
    plt.ylabel("Input_Pulse", fontsize=30)
    plt.plot(forcePulseData[0], forcePulseData[1], lw=1, color='r')
    plt.xlim([0, 100])
    plt.ylim([-300, 300])
    plt.subplots_adjust(left=0.1, bottom=0.2)
    plt.grid(which='major',color='black',linestyle='--')
    filename    = "force_pulse.eps"
    plt.savefig(OUTPUT_PATH+"/"+filename, format="eps", dpi=300)
    plt.clf()
    plt.tick_params(labelsize=24)
