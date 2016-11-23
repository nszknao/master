# -*- coding: utf-8 -*-
u"""個体群の変位応答分布を作成します．"""
import matplotlib.pyplot as plt
import common as cmn
import numpy as np
import math
from sklearn.decomposition import PCA

# ---- public ----
def getCorMatrix(pop):
    u"""相関係数行列を求める．
    【引数】pop: 個体群
    【戻り値】CorMatrix: 相関係数行列"""
    # 平均を引いた行列を作成
    popObjList = np.zeros((cmn.getConstValue("NUM_OF_MOMENTEQ"), len(pop)))
    for i in range(popObjList.shape[0]):
        sumObj = 0.
        for ii in range(popObjList.shape[1]):
            sumObj += pop[ii].o[i]
        muObj = sumObj/popObjList.shape[1]
        for ii in range(popObjList.shape[1]):
            popObjList[i, ii] = pop[ii].o[i] - muObj

    # 分散共分散行列の作成
    VaCovMatrix = np.zeros((cmn.getConstValue("NUM_OF_MOMENTEQ"), cmn.getConstValue("NUM_OF_MOMENTEQ")))
    for i in range(popObjList.shape[0]):
        for ii in range(popObjList.shape[0]):
            covValue = 0.
            for iii in range(popObjList.shape[1]):
                covValue += (popObjList[i, iii]*popObjList[ii, iii])
            VaCovMatrix[i, ii] = covValue/(float)(popObjList.shape[1]-1.)

    # 相関係数行列を作成
    CorMatrix = np.zeros((cmn.getConstValue("NUM_OF_MOMENTEQ"), cmn.getConstValue("NUM_OF_MOMENTEQ")))
    for i in range(cmn.getConstValue("NUM_OF_MOMENTEQ")):
        for ii in range(cmn.getConstValue("NUM_OF_MOMENTEQ")):
            CorMatrix[i, ii] = VaCovMatrix[i, ii] / math.sqrt(VaCovMatrix[i, i]*VaCovMatrix[ii, ii])

    # 固有値解析
    pca = PCA(n_components=cmn.getConstValue("NUM_OF_MOMENTEQ"))
    pca.fit(CorMatrix)
#    print('---- variance_ratio ----')
#    print(pca.explained_variance_ratio_)
    eigenValueRatio = pca.explained_variance_ratio_
#    print('---- components ----')
#    print(pca.components_)
    eigenVector = pca.components_

    # 次元削減
    selectedObj = []
    sumContRatio = 0.
    num = 0
    while sumContRatio < 0.997:
        allNegativeFlg = True
        allPositiveFlg = True
        # 固有ベクトルの各要素の符号チェック
        for i in range(cmn.getConstValue("NUM_OF_MOMENTEQ")):
            if eigenVector[num, i] > 0:
                allNegativeFlg = False
            if eigenVector[num, i] < 0:
                allPositiveFlg = False

        tmpEigenVector = eigenVector[num]
        if allNegativeFlg == True:
            tmpEigenVector = sorted(tmpEigenVector)
            selectedObj.append(_getKeyInList(eigenVector[num], tmpEigenVector[0]))
            selectedObj.append(_getKeyInList(eigenVector[num], tmpEigenVector[1]))
        elif allPositiveFlg == True:
            tmpEigenVector = sorted(tmpEigenVector, reverse=True)
            selectedObj.append(_getKeyInList(eigenVector[num], tmpEigenVector[0]))
            selectedObj.append(_getKeyInList(eigenVector[num], tmpEigenVector[1]))
        else:
            tmpEigenVector = sorted(tmpEigenVector, key=abs, reverse=True)
            selectedObj.append(_getKeyInList(eigenVector[num], tmpEigenVector[0]))
#            if tmpEigenVector[0] > 0:
#                for i in range(len(eigenVector[num])):
#                    if eigenVector[num, i] < 0:
#                        selectedObj.append(i)
#            if tmpEigenVector[0] < 0:
#                for i in range(len(eigenVector[num])):
#                    if eigenVector[num, i] > 0:
#                        selectedObj.append(i)


        sumContRatio += eigenValueRatio[num]
        num += 1

    selectedObj = sorted(set(selectedObj), key=selectedObj.index)
    print(selectedObj)

def getSpecifiedKLPop(pop, sim, kl):
    selectedPop = []
    for i in range(len(pop)):
        dKL = cmn.klDivergence(cmn.createDispPdf(sim.dispx, pop[i].detailPrm), sim.dispy)
        if dKL < 0:
            continue
        if dKL < kl:
            selectedPop.append(pop[i])
    return selectedPop

# ---- private ----
def _getKeyInList(targetList, value):
    for i in range(len(targetList)):
        if targetList[i] == value:
            return i
    print('value: ' + str(value) + 'is not found in the list')
    sys.exit()
