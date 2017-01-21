# -*- coding: utf-8 -*- 
import numpy as np
import sys
import os
import csv
from scipy import integrate

NUM_OF_GAUSS    = 3     # 足し合わせるガウス分布の数
NUM_OF_MOMENT   = 21
NUM_OF_MOMENTEQ = 12
NUM_OF_PARAM    = 7

dx  = 0.01  # pdfpdfのX軸刻み幅

class Individual:
    u"""個体情報を格納
    メンバ: m[], o[], disp[], vel[], dPdf[], vPdf[], detailPrm[]
    Pythonでは，代入すれはメンバが生成される．元からメンバを明記しておくと識別子が同一として扱われてしまう．"""
    pass

class Simulation:
    u"""シミュレーション解の情報を格納
    メンバ: dispx[], dispy[], velx[], vely[], disp3d[], vel3d[], pdf3d[]"""

class GWhiteNoise:
    u"""ガウス性ホワイトノイズ解の情報を格納
    メンバ: dispx[], dispy[]"""

#---- public ----
def culcIntegralJointPdf(prm):
    u"""結合確率密度関数の積分値を計算します．
    【引数】ind: 個体
    【戻り値】integration: 積分値"""
    integration = 0.
    for i in range(NUM_OF_GAUSS):
        integration += prm['a'][i]*integrate.dblquad(lambda x, y: (1./(2.*np.pi*prm['sigma1'][i]*prm['sigma2'][i])) * np.exp(-1.*((x-prm['mu1'][i])**2/prm['sigma1'][i]**2 + (y-prm['mu2'][i])**2/prm['sigma2'][i]**2)/2.), -np.inf, np.inf, lambda y: -np.inf, lambda y: np.inf)[0]
    return integration

def createLevelCrossingRate(xi, prm):
    u"""閾値通過率を計算します．
    【引数】xi: 閾値，prm: 詳細のパラメータ
    【戻り値】prob: 閾値通過率"""
    prob = 0.
    for i in range(NUM_OF_GAUSS):
        pp_c = prm['rho'][i]*prm['sigma1'][i]*prm['sigma2'][i]/prm['sigma1'][i]/prm['sigma2'][i]
        pp_g = prm['mu2'][i] + pp_c*prm['sigma2'][i]*(xi - prm['mu1'][i])/prm['sigma1'][i]
        pp_sigma    = prm['sigma2'][i]*np.sqrt(1. - pp_c**2)
        # 閾値通過率
        prob += prm['a'][i]*np.exp(-1.*(pp_xi - prm['mu1'][i])**2/2./prm['sigma1'][i]**2)/2./np.pi/prm['sigma1'][i]/prm['sigma2'][i]/np.sqrt(1. - pp_c**2)* (pp_sigma**2*np.exp(-pp_g**2/2./pp_sigma**2)
                        + np.sqrt(np.pi/2.)*pp_g*pp_sigma*(1. + np.erf(pp_g/np.sqrt(2.)/pp_sigma)));
    return prob;

def getConstValue(name):
    if name == "NUM_OF_MOMENT":
        return NUM_OF_MOMENT
    elif name == "NUM_OF_MOMENTEQ":
        return NUM_OF_MOMENTEQ
    elif name == "NUM_OF_PARAM":
        return NUM_OF_PARAM
    elif name == "NUM_OF_GAUSS":
        return NUM_OF_GAUSS
    else:
        print("ERROR: 指定した定数は存在しません．")
        sys.exit()

def createDispPdf(x, prm):
    u"""パラメータから変位応答を計算します．
    【引数】x: はシミュレーションから取得，prm: 詳細のパラメータ
    【戻り値】pdf: 変位応答分布の値リスト"""
    pdf = 0.
    for i in range(0, NUM_OF_GAUSS):
        pdf += prm['a'][i]*(1./np.sqrt(2.*np.pi)/prm['sigma1'][i])*np.exp(-1.*(x-prm['mu1'][i])**2/2./prm['sigma1'][i]**2)
    return pdf

def createVelPdf(x, prm):
    u"""パラメータから速度応答を計算します．
    【引数】x: はシミュレーションから取得，prm: 詳細のパラメータ
    【戻り値】pdf: 速度応答分布の値リスト"""
    pdf = 0.
    for i in range(0, NUM_OF_GAUSS):
        pdf += prm['a'][i]*(1./np.sqrt(2.*np.pi)/prm['sigma2'][i])*np.exp(-1.*(x-prm['mu2'][i])**2/2./prm['sigma2'][i]**2)
    return pdf

def create3DPdf(x, y, prm):
    u"""パラメータから3次元応答分布を計算します．
    【引数】x: はシミュレーションから取得，prm: 詳細のパラメータ
    【戻り値】pdf: 応答分布の値リスト"""
    pdf = 0
    for i in range(0, NUM_OF_GAUSS):
#        pdf += prm['a'][i]*(1./(2.*np.pi*prm['sigma1'][i]*prm['sigma2'][i]*np.sqrt(1.-prm['rho'][i]**2))) * np.exp(-1.*((x-prm['mu1'][i])**2/prm['sigma1'][i]**2 - 2.*prm['rho'][i]*(x-prm['mu1'][i])*(y-prm['mu2'][i])/prm['sigma1'][i]/prm['sigma2'][i] + (y-prm['mu2'][i])**2/prm['sigma2'][i]**2)/(2.*(1.-prm['rho'][i]**2)))
        pdf += prm['a'][i]*(1./(2.*np.pi*prm['sigma1'][i]*prm['sigma2'][i])) * np.exp(-1.*((x-prm['mu1'][i])**2/prm['sigma1'][i]**2 + (y-prm['mu2'][i])**2/prm['sigma2'][i]**2)/2.)
    return pdf

def klDivergence(comp, true):
    u"""KLダイバージェンスを計算する．
    【引数】comp: 比較対象の確率分布の値リスト，true: 真値の確率分布の値リスト
    【戻り値】ret: KLダイバージェンスの値"""
    if len(comp) != len(true):
        print("comp size: " + str(len(comp)) + "true size: " + str(len(true)))
        print("KLダイバージェンス計算時のlistのサイズが異なります．")
        sys.exit()
    ret = 0.
    for i in range(len(true)):
        if (comp[i] == 0) or (true[i] == 0):
            continue
        else:
            ret += comp[i]*np.log(comp[i]/true[i])
    return ret

def getSimulationFromFile(arg_l, arg_a):
    u"""シミュレーションファイル読み込み
    【引数】arg_l: ラムダ，arg_a: アルファ
    【戻り値】sim: プロット用のシミュレーション解"""
    sim = Simulation()
    # 変位
    sim.dispx = []
    sim.dispy = []
    for line_s in open('/usr/local/src/master/dat/l='+str(arg_l)+'/a='+str(arg_a)+'/sim_y1pdf.dat'):
        col_s   = line_s.strip().split(' ')
        sim.dispx.append(float(col_s[0]))
        sim.dispy.append(float(col_s[1]))

    # 速度
    sim.velx = []
    sim.vely = []
    for line_s in open('/usr/local/src/master/dat/l='+str(arg_l)+'/a='+str(arg_a)+'/sim_y2pdf.dat'):
        col_s   = line_s.strip().split(' ')
        sim.velx.append(float(col_s[0]))
        sim.vely.append(float(col_s[1]))

    # 結合
    sim.disp3d = []
    sim.vel3d = []
    sim.pdf3d = []
    f = open('/usr/local/src/master/dat/l='+str(arg_l)+'/a='+str(arg_a)+'/sim_jointpdf.dat', 'r')
    dataReader = csv.reader(f)
    lineCount = 0    
    for row in dataReader:
        if lineCount == 0:
            del row[0]
            sim.disp3d = [float(row[i]) for i in range(len(row))]            
        else :
            row = [float(row[i]) for i in range(len(row))]
            sim.vel3d.append(row.pop(0))
            if len(sim.pdf3d) == 0:
                sim.pdf3d = np.array([row])
            else :
                sim.pdf3d = np.append(sim.pdf3d, np.array([row]), axis=0)
        lineCount += 1
    
    return sim

def getGWhiteNoiseFromFile():
    gwn = GWhiteNoise()
    # ホワイトノイズのみの厳密解ファイルの読み込み
    gwn.dispx = []
    gwn.dispy = []
    for line_s in open('/usr/local/src/master/dat/sim_y1Gpdf.dat'):
        col_s = line_s.strip().split(' ')
        gwn.dispx.append(float(col_s[0]))
        gwn.dispy.append(float(col_s[1]))
    return gwn

def getPopFromFile(arg_l, arg_a):
    # プロット用解析解X軸
    # for文の中で１回づつ読み込むのは無駄だから上で宣言しておく
    anaDispX = np.arange(-6., 6., dx)
    anaVelX   = np.arange(-12., 12., dx)

    # 解析ファイル読み込み
    fileNum = _getFileNum('/usr/local/src/master/results/l='+str(arg_l)+'/dat_a='+str(arg_a)+'/')
    pop = []
    for i in range(0, fileNum):
        print('ana_gsay1pdf_' + str(i) + '.pdf is reading')
        for line_a in open('/usr/local/src/master/results/l='+str(arg_l)+'/dat_a='+str(arg_a)+'/ana_gsay1pdf_'+str(i)+'.dat'):
            col_a = line_a.strip().split(' ')
            ind = Individual()
            #--- モーメント値 ---
            ind.m = [float(col_a[i]) for i in range(NUM_OF_MOMENT)]
            #--- 目的関数値 ---
            ind.o = [float(col_a[i]) for i in range(NUM_OF_MOMENT, NUM_OF_MOMENT + NUM_OF_MOMENTEQ)]
            #--- パラメータ ---
            simplePrm = [float(col_a[i]) for i in range(NUM_OF_MOMENT + NUM_OF_MOMENTEQ, NUM_OF_MOMENT + NUM_OF_MOMENTEQ + NUM_OF_PARAM)]
            ind.detailPrm = _getDetailParameterFromSimpleNotation(simplePrm)
            #--- 変位 ---
            ind.disp = anaDispX
            #--- 速度 ---
            ind.vel = anaVelX
            #--- 変位応答分布 ---
            ind.dPdf = createDispPdf(anaDispX, ind.detailPrm)
            #--- 速度応答分布 ---
            ind.vPdf = createVelPdf(anaVelX, ind.detailPrm)

            # PDFの積分がおかしいものを除去（引っかかる個体はなさそうなのでカット）
#            sumJointPdf = culcIntegralJointPdf(ind.detailPrm)
#            print(sumJointPdf)
#            if sumJointPdf > 1.005 or sumJointPdf < 0.995:
#                continue

            pop.append(ind)
            del ind
    
    return pop

def getStandardDeviationList(pop):
    u"""個体群の不偏標準偏差を求める
    【引数】pop: 個体群
    【戻り値】meanSdList: 2次元配列で，0が平均・1が不偏標準偏差のリスト"""
    meanSdList = [[0 for i in range(NUM_OF_MOMENTEQ)] for ii in range(2)]
    # 平均
    for i in range(NUM_OF_MOMENTEQ):
        sumObj = 0.
        for ii in range(len(pop)):
            sumObj += pop[ii].o[i]
        meanSdList[0][i] = sumObj / len(pop)
    # 分散
    for i in range(NUM_OF_MOMENTEQ):
        sumSquareObj = 0.
        for ii in range(len(pop)):
            sumSquareObj += (pop[ii].o[i] - meanSdList[0][i])**2
        if sumSquareObj == 0:
            meanSdList[1][i] = 1.
        else:
            meanSdList[1][i] = np.sqrt(sumSquareObj / (len(pop)))

    return meanSdList

def getSquareObjectiveValue(ind, meanSdList):
    u"""個体の目的関数の二乗誤差を取得
    【引数】ind: 個体，meanSdList: 平均と不偏標準偏差のリスト
    【戻り値】squareValue: 目的関数の二乗誤差"""
    squareValue = 0.
    for i in range(NUM_OF_MOMENTEQ):
        if meanSdList[1][i] == 0:
            squareValue += ind.o[i]**2
        else:
            squareValue += (ind.o[i]/meanSdList[1][i])**2
    return squareValue

def getKurtosisValue(x, y, dt):
    if len(x) != len(y):
        print("ERROR: Wrong size of x and y")
        sys.exit()
    mean = np.dot(x, y)*dt
    
    center = x - np.array([mean]*len(x))
    th2Moment = np.dot(np.power(center, 2), y)*dt
    th4Moment = np.dot(np.power(center, 4), y)*dt

    kurtosis = th4Moment / th2Moment**2 - 3.
    return

#---- private ----
def _getFileNum(path):
    u"""指定したディレクトリ配下のファイル数をカウント
    【引数】path: ディレクトリへのパス
    【戻り値】counter: ファイル数"""
    ch = os.listdir(path)
    counter = 0
    for i in ch:
        if (os.path.isdir(path + i)):
            checkFileNum(path + i + "/")
        else:
            counter += 1
    return counter

def _getDetailParameterFromSimpleNotation(simple):
    u"""パラメータの形式を変換します。
    【引数】simple: 簡易のパラメータ
    【戻り値】detail: 詳細のパラメータ"""
    detail  = {}
    # 重み
    detail['a'] = [0]*NUM_OF_GAUSS
    detail['a'][0]  = (1. - simple[0]) / 2.
    detail['a'][1]  = simple[0]
    detail['a'][2]  = (1. - simple[0]) / 2.
    # 変位
    detail['mu1']   = [0]*NUM_OF_GAUSS
    detail['mu1'][0]    = simple[1]
    detail['mu1'][1]    = 0.
    detail['mu1'][2]    = -1.*simple[1]
    detail['sigma1']    = [0]*NUM_OF_GAUSS
    detail['sigma1'][0] = simple[3]
    detail['sigma1'][1] = simple[5]
    detail['sigma1'][2] = simple[3]
    # 速度
    detail['mu2']   = [0]*NUM_OF_GAUSS
    detail['mu2'][0]    = simple[2]
    detail['mu2'][1]    = 0.
    detail['mu2'][2]    = -1.*simple[2]
    detail['sigma2']    = [0]*NUM_OF_GAUSS
    detail['sigma2'][0] = simple[4]
    detail['sigma2'][1] = simple[6]
    detail['sigma2'][2] = simple[4]
    # 共分散係数
#    detail['rho'] = [0]*NUM_OF_GAUSS
#    detail['rho'][0]  = simple[7]
#    detail['rho'][1]  = simple[8]
#    detail['rho'][2]  = simple[9]
    return detail

