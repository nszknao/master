# -*- coding: utf-8 -*-
u"""メインの処理を実行します．"""
import common as cmn
import numpy as np
import analysis as ana
import plotFig as fig
import sys

if __name__ == "__main__":
    # コマンドライン引数読み込み
    arg_l   = float(sys.argv[1])
    arg_a   = float(sys.argv[2])
    
    sim = cmn.getSimulationFromFile(arg_l, arg_a)
    gwn = cmn.getGWhiteNoiseFromFile()
    pop = cmn.getPopFromFile(arg_l, arg_a)

#    pop = sorted(pop, key=lambda x: cmn.klDivergence([cmn.createDispPdf(sim.dispx[i], x.detailPrm) for i in range(len(sim.dispx))], sim.dispy))
    pop = sorted(pop, key=lambda x: cmn.klDivergence(np.reshape(cmn.create3DPdf(np.meshgrid(sim.disp3d, sim.vel3d)[0], np.meshgrid(sim.disp3d, sim.vel3d)[1], x.detailPrm), len(sim.disp3d)*len(sim.vel3d)), np.reshape(sim.pdf3d, len(sim.disp3d)*len(sim.vel3d))))
    
#    meanSdList = cmn.getStandardDeviationList(pop)
#    minSquareInd = ana.getMinSquareObjectInd(pop)
#    fig.plotDispPdf_Ana_Sim_Gauss(minSquareInd, sim, gwn, "/usr/local/src/master/fig/aaa.eps")
    
    selectedInd = ana.culcEquivalentLinearizationMethod(pop, arg_l, arg_a)
    print(selectedInd.detailPrm)    
    # 3D系
    fig.plot3DDispPdf_Ana(selectedInd, sim, "/usr/local/src/master/fig/JointPdf_Ana_" + str(arg_l) + "_" + str(arg_a) + ".eps")
    fig.plot3DDispPdf_Sim(sim, "/usr/local/src/master/fig/JointPdf_Sim_" + str(arg_l) + "_" + str(arg_a) + ".eps")
    fig.plotRelation_KLDivergence_Objective(pop, sim, "/usr/local/src/master/fig/Relation_of_KLDivergence-SquareObjectiveValue_" + str(arg_l) + "_" + str(arg_a) + ".eps")
    # 2D系
    fig.plotDispPdf_Ana_Sim_Gauss(selectedInd, sim, gwn, "/usr/local/src/master/fig/bbb.eps")
    fig.plotVelPdf_Ana_Sim(selectedInd, sim, "/usr/local/src/master/fig/ddd.eps")
#    fig.plotIndObjectiveValue(pop[0], meanSdList, "/usr/local/src/master/fig/ObjectiveValue_" + str(arg_l) + "_" + str(arg_a) + ".eps")
#    fig.plotMomentValue(pop[0], "/usr/local/src/master/fig/MomentValue_" + str(arg_l) + "_" + str(arg_a) + ".eps")
#    fig.plotPopObjectiveValue(pop, meanSdList, 10, "/usr/local/src/master/fig/10_ObjectiveValue_" + str(arg_l) + "_" + str(arg_a) + ".eps")
