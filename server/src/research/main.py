# -*- coding: utf-8 -*-
u"""メインの処理を実行します．"""
import common as cmn
import analysis as ana
import plotFig as fig
import sys

if __name__ == "__main__":
    # コマンドライン引数読み込み
    arg_l   = sys.argv[1]
    arg_a   = sys.argv[2]
    
    sim = cmn.getSimulationFromFile(arg_l, arg_a)
    gwn = cmn.getGWhiteNoiseFromFile()
    
    pop = cmn.getPopFromFile(arg_l, arg_a)
    pop = sorted(pop, key=lambda x: cmn.klDivergence(cmn.createDispPdf(sim.dispx, x.detailPrm), sim.dispy))
    i = 0
    while True:
        kl = cmn.klDivergence(cmn.createDispPdf(sim.dispx, pop[i].detailPrm), sim.dispy)
        if kl < 0:
            del pop[i]
        else:
            i += 1
        if i == len(pop):
            break

    minSquareInd = ana.getMinSquareObjectInd(pop)
    fig.plotDispPdf_Ana_Sim_Gauss(minSquareInd, sim, gwn, "/usr/local/src/master/fig")
    
#    selectedPop = ana.getSpecifiedKLPop(pop, sim, 1.)
#    print("Selected pop number(using for PCA) is " + str(len(selectedPop)))
#    ana.getCorMatrix(selectedPop)

#    fig.plotDispPdf_Ana_Sim_Gauss(pop[1], sim, gwn, "/usr/local/src/master/fig")
    fig.plotRelation_KLDivergence_Objective(pop, sim, "/usr/local/src/master/fig")
    fig.plotObjectiveValue(pop, "/usr/local/src/master/fig")
