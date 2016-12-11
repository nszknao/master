# -*- coding: utf-8 -*-
u"""個体群に関するいろいろな図を作成します．"""
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import common as cmn
import math
from matplotlib.font_manager import FontProperties

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
    filename    = "KL=" + str(cmn.klDivergence(cmn.createDispPdf(sim.dispx, ind.detailPrm), sim.dispy)) + ".eps"
    plt.savefig(savePath + "/" + filename, format="eps", dpi=300)
    plt.clf()

def plotRelation_KLDivergence_Objective(pop, sim, savePath):
    u"""KLダイバージェンスと目的関数値の二乗誤差を散布図にプロット
    【引数】pop: 個体群，sim: シミュレーション解，savePath: 保存先のパス"""
    meanSdList = cmn.getStandardDeviationList(pop)
    squareObjList = []
    dKLList = []
    for i in range(len(pop)):
        squareValue = 0.
        for ii in range(cmn.getConstValue("NUM_OF_MOMENTEQ")):
            squareValue += ((pop[i].o[ii]-meanSdList[0][ii])/meanSdList[1][ii])**2
        squareObjList.append(squareValue)
        dKLList.append(cmn.klDivergence(cmn.createDispPdf(sim.dispx, pop[i].detailPrm), sim.dispy))
    
    plt.scatter(squareObjList, dKLList)
    plt.xlim([0, 20])
    plt.ylim([0.0001, 100])
    plt.yscale('log')
    plt.gca().xaxis.set_major_locator(tick.MultipleLocator(1))
    plt.title(u'KLダイバージェンスと目的関数値の二乗誤差の関係', fontproperties=fp)
    plt.xlabel('Square Objective Value')
    plt.ylabel('KL-Divergence')
    plt.grid(which='major',color='black',linestyle='--')

    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()

def plotObjectiveValue(pop, savePath):
    u"""個体の目的関数値をプロット
    【引数】ind: 個体，savePath: 保存先のパス"""
    meanSdList = cmn.getStandardDeviationList(pop)
    plt.title(u'個体の目的関数値', fontproperties=fp)
    #plt.plot(range(len(pop[0].o)), [(pop[100].o[i]-meanSdList[0][i])/meanSdList[1][i] for i in range(cmn.getConstValue("NUM_OF_MOMENTEQ"))])
    for ii in range(1):
        plt.plot(range(len(pop[0].o)), [(pop[ii].o[i]-meanSdList[0][i])/meanSdList[1][i] for i in range(cmn.getConstValue("NUM_OF_MOMENTEQ"))])
    plt.xlabel('Objective Number')
    plt.ylabel('Objective Value')
    plt.ylim([-0.2, 1.6])
    plt.gca().yaxis.set_major_locator(tick.MultipleLocator(0.2))
    plt.grid(which='major',color='black',linestyle='--')

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
    plt.ylim([-1500, 1500])
    plt.gca().yaxis.set_major_locator(tick.MultipleLocator(300))
    plt.grid(which='major',color='black',linestyle='--')

    plt.savefig(savePath, format="eps", dpi=300)
    plt.clf()
