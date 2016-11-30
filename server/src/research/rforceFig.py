# -*- coding: utf-8 -*-
u"""復元力特性の図を作成します．"""
import matplotlib.pyplot as plt

HIGH_EPSILON    = 0.3
LOW_EPSILON = 0.1

def drange(begin, end, step):
    u"""小数刻みでrangeできるように拡張
    【引数】begin: 開始の値，end: 終了の値，step: ステップ幅"""
    n = begin
    while n+step < end:
        yield n
        n += step

def duffingRestoringForce(x, epsilon):
    u"""Duffing系の復元力の値を返す
    【引数】x: xの値
    【戻り値】y: 復元力の値，epsilon: 非線形強度"""
    return x + epsilon*x**3


if __name__ == "__main__":
    # xは0.01刻みで-8~+8
    x   = []
    for i in drange(-10,10,0.01):
        x.append(i)
    # Duffing系の復元力
    high_duffing    = []
    for i in range(0, len(x)):
        rf  = duffingRestoringForce(x[i], HIGH_EPSILON)
        high_duffing.append(rf)
    low_duffing = []
    for i in range(0, len(x)):
        rf  = duffingRestoringForce(x[i], LOW_EPSILON)
        low_duffing.append(rf)
    # プロット
    outputPath  = "/usr/local/src/master/fig"
    plt.xlabel("displacement", fontsize=18)
    plt.ylabel("Restoring Force", fontsize=18)
    plt.plot(x, x, lw=3)
    plt.plot(x, high_duffing,  lw=3)
    plt.plot(x, low_duffing,  lw=3)
#   plt.legend(['Linear', 'Duffing($\\varepsilon$='+str(HIGH_EPSILON)+')', 'Duffing($\\varepsilon$='+str(LOW_EPSILON)+')'], loc='upper left')
    plt.xlim([-10, 10])
    plt.ylim([-50, 50])
    plt.grid(which='major',color='black',linestyle='--')
    plt.tick_params(labelsize=18)
    plt.tight_layout()
    filename    = "rf.eps"
    plt.savefig(outputPath+"/"+filename, format='eps', dpi=300) # 解像度高めに
    plt.clf()
