# -*- coding: utf-8 -*-
u"""個体群に関するいろいろな図を作成します．"""
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import common as cmn
import math
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fp = FontProperties(fname='/usr/share/fonts/TakaoGothic.ttf')

MOMENT_ARRAY = ["$E[y_1^2]$", "$E[y_1y_2]$", "$E[y_2^2]$", "$E[y_1^4]$", "$E[y_1^3y_2$", "$E[y_1^2y_2^2]$", "$E[y_1y_2^3]$", "$E[y_2^4]$", "$E[y_1^6]$", "$E[y_1^5y_2]$", "$E[y_1^4y_2^2]$", "$E[y_1^3y_2^3]$", "$E[y_1^2y_2^4]$", "$E[y_1y_2^5]$", "$E[y_2^6]$", "$E[y_1^8]$", "$E[y_1^7y_2]$", "$E[y_1^6y_2^2]$", "$E[y_1^5y_2^3]$", "$E[y_1^4y_2^4]$", "$E[y_1^3y_2^5]$"]

def plotDispPdf_Ana_Sim_Gauss(ind, sim, gwn, savePath):
    u"""個体の変位応答分布をプロット（シミュレーション解，ホワイトノイズ解と合わせて）
    【引数】ind: 個体，sim: シミュレーション解，gwn: ホワイトノイズ解，savePath: 保存先のパス"""
    int_dPdf    = cmn.culcIntegralPdf(ind.dPdf)
    if (0.95 > int_dPdf) or (int_dPdf > 1.05):
        print("PDF integral violation: " + str(int_dPdf))
        return
    plt.xlim([-6, 6])
    plt.ylim([0.00001, 1])
    plt.xlabel("displacement", fontsize=18)
    plt.ylabel("Probability density", fontsize=18)
    plt.grid(which='major',color='black',linestyle='--')
    plt.plot(ind.disp, ind.dPdf)
    plt.plot(sim.dispx, sim.dispy, 'o')
    plt.plot(gwn.dispx, gwn.dispy)
    plt.legend(['Analytical solution', 'Simulation', 'Gauss-White'], loc='lower center')
    plt.yscale('log')
    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()

def plotVelPdf_Ana_Sim(ind, sim, savePath):
    u"""個体の速度応答分布をプロット（シミュレーション解，ホワイトノイズ解と合わせて）
    【引数】ind: 個体，sim: シミュレーション解，gwn: ホワイトノイズ解，savePath: 保存先のパス"""
    int_vPdf    = cmn.culcIntegralPdf(ind.vPdf)
    if (0.95 > int_vPdf) or (int_vPdf > 1.05):
        print("PDF integral violation: " + str(int_vPdf))
        return
    plt.xlim([-12, 12])
    plt.ylim([0.00001, 1])
    plt.xlabel("velocity", fontsize=18)
    plt.ylabel("Probability density", fontsize=18)
    plt.grid(which='major',color='black',linestyle='--')
    plt.plot(ind.vel, ind.vPdf)
    plt.plot(sim.velx, sim.vely, 'o')
#    plt.plot(gwn.velx, gwn.vely)
#    plt.legend(['Analytical solution', 'Simulation', 'Gauss-White'], loc='lower center')
    plt.legend(['Analytical solution', 'Simulation'], loc='lower center')
    plt.yscale('log')
    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()

def plot3DDispPdf_Ana(ind, sim, savePath):
    fig = plt.figure()
    ax = Axes3D(fig)
    x = sim.dispx[::1]
    y = sim.dispy[::1]
    X, Y = np.meshgrid(x, y)
    Z = [cmn.create3DPdf(x[i], y[i], ind.detailPrm) for i in range(len(x))]
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)

    ax.set_xlabel(r'$y_1$', fontsize=20)
    ax.set_ylabel(r'$y_2$', fontsize=20)
    ax.set_zlabel(r'$p_y (y_1, y_2)$', fontsize=20)
#    plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))  # 指数表記
    plt.savefig(savePath, format="eps", dpi=300)

def plot3DDispPdf_Sim(sim, savePath):
    fig = plt.figure()
    ax = Axes3D(fig)
    x = sim.disp3d[::1]
    y = sim.vel3d[::1]
    Z = sim.pdf3d[0::1, 0::1]
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    ax.set_xlabel(r'$y_1$', fontsize=20)
    ax.set_ylabel(r'$y_2$', fontsize=20)
    ax.set_zlabel(r'$p_y (y_1, y_2)$', fontsize=20)
    plt.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))  # 指数表記
    plt.savefig(savePath, format="eps", dpi=300)

def plotRelation_KLDivergence_Objective(pop, sim, savePath):
    u"""KLダイバージェンスと目的関数値の二乗誤差を散布図にプロット
    【引数】pop: 個体群，sim: シミュレーション解，savePath: 保存先のパス"""
    meanSdList = cmn.getStandardDeviationList(pop)
    squareObjList = []
    dKLList = []
    for i in range(len(pop)):
        squareValue = 0.
        for ii in range(cmn.getConstValue("NUM_OF_MOMENTEQ")):
#            squareValue += (pop[i].o[ii]/meanSdList[1][ii])**2
            squareValue += abs(pop[i].o[ii]/meanSdList[1][ii])
        squareObjList.append(squareValue)
        tmp_dPdf = [cmn.createDispPdf(sim.dispx[ii], pop[i].detailPrm) for ii in range(len(sim.dispx))]
        dKLList.append(cmn.klDivergence(tmp_dPdf, sim.dispy))
    
    plt.scatter(squareObjList, dKLList)
#    plt.xlim([0, 5])
#    plt.ylim([0.0001, 50])
    plt.yscale('log')
    plt.gca().xaxis.set_major_locator(tick.MultipleLocator(1))
    plt.title(u'KLダイバージェンスと目的関数値の二乗誤差の関係', fontproperties=fp)
    plt.xlabel('Square Objective Value')
    plt.ylabel('KL-Divergence')
    plt.grid(which='major',color='black',linestyle='--')

    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()

def plotSpecificObjectValue(pop, key, savePath):
    u"""個体群の指定した番号の目的関数値をプロット
    【引数】pop: 個体群，key: 目的関数の番号，savePath: 保存先のパス"""
    plt.scatter(range(len(pop)), [pop[i].o[key] for i in range(len(pop))])
#    plt.xlim([0, 10])
#    plt.ylim([0.0001, 50])
    plt.title(u'特定の目的関数値', fontproperties=fp)
    plt.ylabel('Specific Objective Value')
    plt.xlabel('Number')
    plt.grid(which='major',color='black',linestyle='--')

    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()

def plotIndObjectiveValue(ind, meanSdList, savePath):
    u"""個体の目的関数値をプロット
    【引数】ind: 個体，meanSdList: 平均と不偏分散のリスト，savePath: 保存先のパス"""
    plt.title(u'個体の目的関数値', fontproperties=fp)
    plt.plot(range(len(ind.o)), [ind.o[i]/meanSdList[1][i] for i in range(cmn.getConstValue("NUM_OF_MOMENTEQ"))])
    plt.xlabel('Objective Number')
    plt.ylabel('Objective Value')
#    plt.ylim([-0.4, 0.4])
    plt.gca().yaxis.set_major_locator(tick.MultipleLocator(0.2))
    plt.grid(which='major',color='black',linestyle='--')

    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()

def plotPopObjectiveValue(pop, meanSdList, num, savePath):
    u"""上位個体群の目的関数値をプロット
    【引数】pop: 個体群，meanSdList: 平均と不偏分散のリスト，num: 個体数，savePath: 保存先のパス"""
    plt.title(u'個体の目的関数値', fontproperties=fp)
    for i in range(num):
        plt.plot(range(len(pop[i].o)), [pop[i].o[ii]/meanSdList[1][ii] for ii in range(cmn.getConstValue("NUM_OF_MOMENTEQ"))], label="No."+str(i))
    plt.xlabel('Objective Number')
    plt.ylabel('Objective Value')
    plt.ylim([-0.4, 0.4])
    plt.gca().yaxis.set_major_locator(tick.MultipleLocator(0.2))
    plt.grid(which='major',color='black',linestyle='--')
    plt.legend()

    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()

def plotMomentValue(ind, savePath):
    u"""個体のモーメント値をプロット
    【引数】ind: 個体，savePath: 保存先のパス"""
    plt.title(u'個体のモーメント値', fontproperties=fp)
    plt.plot(range(len(ind.m)), ind.m)
    plt.xlabel('Moment')
    plt.ylabel('Moment Value')
    plt.xticks(range(len(ind.m)), MOMENT_ARRAY, rotation=45)
#    plt.ylim([-1500, 1500])
#    plt.gca().yaxis.set_major_locator(tick.MultipleLocator(300))
    plt.grid(which='major',color='black',linestyle='--')

    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()
