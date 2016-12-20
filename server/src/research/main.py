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
#    _pop = ana.getSpecifiedKLPop(pop, sim, 4)
    _pop = ana.getAroundSpecifiedSquareObjectValue(pop, 0, 4.)
#    for i in range(cmn.getConstValue('NUM_OF_MOMENTEQ')):
#        fig.plotSpecificObjectValue(_pop, i, '/usr/local/src/master/fig/a' + str(i) + '.eps')
#    for i in range(len(_pop)):
#        filename    = "_Obj=" + str(round(cmn.getSquareObjectiveValue(_pop[i], cmn.getStandardDeviationList(pop)), 5)) + ".eps"
#        fig.plotDispPdf_Ana_Sim_Gauss(_pop[i], sim, gwn, '/usr/local/src/master/fig/' + filename)
    meanSdList = cmn.getStandardDeviationList(pop)
    minSquareInd = ana.getMinSquareObjectInd(pop)
    fig.plotDispPdf_Ana_Sim_Gauss(minSquareInd, sim, gwn, "/usr/local/src/master/fig/aaa.eps")
    
#    selectedPop = ana.getSpecifiedKLPop(pop, sim, 1.)
#    print("Selected pop number(using for PCA) is " + str(len(selectedPop)))
#    ana.getCorMatrix(selectedPop)

    print(pop[0].detailPrm)
    fig.plotDispPdf_Ana_Sim_Gauss(_pop[0], sim, gwn, "/usr/local/src/master/fig/bbb.eps")
    fig.plotRelation_KLDivergence_Objective(pop, sim, "/usr/local/src/master/fig/_Relation_of_KLDivergence-SquareObjectiveValue_" + str(arg_l) + "_" + str(arg_a) + ".eps")

    fig.plotObjectiveValue(pop, "/usr/local/src/master/fig/ObjectiveValue_" + str(arg_l) + "_" + str(arg_a) + ".eps")
#    fig.plotMomentValue(pop[0], "/usr/local/src/master/fig/_MomentValue_" + str(arg_l) + "_" + str(arg_a) + ".eps")
