# -*- coding: utf-8 -*-
u"""個体群に関するいろいろな図を作成します．"""
import matplotlib.pyplot as plt
import common as cmn
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
    filename    = "KL=" + str(cmn.klDivergence(cmn.createDispPdf(sim.dispx, ind.detailPrm), sim.dispy)) + ", E[y1^2]=" + str(ind.m[0]) + ", E[y1^4]=" + str(ind.m[3]) + ", E[y1^6]=" + str(ind.m[8]) + ".png"
    plt.savefig(savePath + "/" + filename)
    plt.clf()

def plotRelation_KLDivergence_Objective(pop, sim, savePath):
    u"""KLダイバージェンスと目的関数値の二乗誤差を散布図にプロット
    【引数】pop: 個体群，sim: シミュレーション解，savePath: 保存先のパス"""
    squareObjList = []
    dKLList = []
    for i in range(len(pop)):
        sumSquareObj = 0.
        for ii in range(cmn.getConstValue("NUM_OF_MOMENTEQ")):
            sumSquareObj += pop[i].o[ii]**2
        squareObjList.append(sumSquareObj)
        dKLList.append(cmn.klDivergence(cmn.createDispPdf(sim.dispx, pop[i].detailPrm), sim.dispy))
    
    plt.scatter(squareObjList, dKLList)
    plt.title(u'KLダイバージェンスと目的関数値の二乗誤差の関係', fontproperties=fp)
    plt.xlabel('KL-Divergence')
    plt.ylabel('Square Objective Value')

    figname = "Relation_of_KLDivergence-SquareObjectiveValue"
    plt.savefig(savePath + "/" + figname)
    plt.clf()


