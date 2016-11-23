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
    
    selectedPop = ana.getSpecifiedKLPop(pop, sim, 5)
    print("Selected pop number(using for PCA) is " + str(len(selectedPop)))
    ana.getCorMatrix(selectedPop)

#    fig.plotDispPdf_Ana_Sim_Gauss(pop[0], sim, gwn, "/usr/local/src/master/fig")
    fig.plotRelation_KLDivergence_Objective(pop, sim, "/usr/local/src/master/fig")
