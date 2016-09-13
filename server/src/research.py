# -*- coding: utf-8 -*- 
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import math
import sys

NUM_OF_GAUSS	= 3		# 足し合わせるガウス分布の数
NUM_OF_MOMENT	= 21
NUM_OF_MOMENTEQ	= 15
NUM_OF_PARAM	= 10

class Individual:
	u"""個体情報を格納
	メンバ: m[], o[], disp[], vel[], dPdf[], vPdf[], dKL, vKL
	Pythonでは，代入すれはメンバが生成される．元からメンバを明記しておくと識別子が同一として扱われてしまう．"""
	pass

def createDispPdf(x, prm):
	u"""パラメータから変位応答を計算します．
	【引数】x: はシミュレーションから取得，prm: 詳細のパラメータ
	【戻り値】pdf: 変位応答分布の値リスト"""
	pdf	= []
	integrantion	= 0.
	for i in range(0, len(x)):
		tmp_pdf	= 0.
		for ii in range(0, NUM_OF_GAUSS):
			tmp_pdf += prm['a'][ii]*(1./math.sqrt(2.*math.pi)/prm['sigma1'][ii])*math.exp(-1.*(x[i]-prm['mu1'][ii])**2/2./prm['sigma1'][ii]**2)
		pdf.append(tmp_pdf)
		integrantion	+= tmp_pdf*0.01
	print("disp pdf integration: " + str(integrantion))
	return pdf

def createVelPdf(x, prm):
	u"""パラメータから変位応答を計算します．
	【引数】x: はシミュレーションから取得，prm: 詳細のパラメータ
	【戻り値】pdf: 速度応答分布の値リスト"""
	pdf	= []
	integrantion	= 0.
	for i in range(0, len(x)):
		tmp_pdf	= 0.
		for ii in range(0, NUM_OF_GAUSS):
			tmp_pdf += prm['a'][ii]*(1./math.sqrt(2.*math.pi)/prm['sigma2'][ii])*math.exp(-1.*(x[i]-prm['mu2'][ii])**2/2./prm['sigma2'][ii]**2)
		pdf.append(tmp_pdf)
		integrantion	+= tmp_pdf*0.01
	print("vel pdf integration: " + str(integrantion))
	return pdf

def createLevelCrossingRate(xi, prm):
	u"""閾値通過率を計算します．
	【引数】xi: 閾値，prm: 詳細のパラメータ
	【戻り値】prob: 閾値通過率"""
	prob	= 0.
	for i in range(0, NUM_OF_GAUSS):
		pp_c	= prm['kappa'][i]/prm['sigma1'][i]/prm['sigma2'][i]
		pp_g	= prm['mu2'][i] + pp_c*prm['sigma2'][i]*(xi - prm['mu1'][i])/prm['sigma1'][i]
		pp_sigma	= prm['sigma2'][i]*math.sqrt(1. - pp_c**2)
		# 閾値通過率
		prob	+= prm['a'][i]*math.exp(-1.*(pp_xi - prm['mu1'][i])**2/2./prm['sigma1'][i]**2)/2./math.pi/prm['sigma1'][i]/prm['sigma2'][i]/math.sqrt(1. - pp_c**2)* (pp_sigma**2*math.exp(-pp_g**2/2./pp_sigma**2)
						+ math.sqrt(math.pi/2.)*pp_g*pp_sigma*(1. + math.erf(pp_g/math.sqrt(2.)/pp_sigma)));
	return prob;


def getDetailParameterFromSimpleNotation(simple):
	u"""パラメータの形式を変換します。
	【引数】simple: 簡易のパラメータ
	【戻り値】detail: 詳細のパラメータ"""
	detail	= {}
	# 重み
	detail['a']	= [0]*NUM_OF_GAUSS
	detail['a'][0]	= (1. - simple[0]) / 2.
	detail['a'][1]	= simple[0]
	detail['a'][2]	= (1. - simple[0]) / 2.
	# 変位
	detail['mu1']	= [0]*NUM_OF_GAUSS
	detail['mu1'][0]	= simple[1]
	detail['mu1'][1]	= 0.
	detail['mu1'][2]	= -1.*simple[1]
	detail['sigma1']	= [0]*NUM_OF_GAUSS
	detail['sigma1'][0]	= simple[3]
	detail['sigma1'][1]	= simple[5]
	detail['sigma1'][2]	= simple[3]
	# 速度
	detail['mu2']	= [0]*NUM_OF_GAUSS
	detail['mu2'][0]	= simple[2]
	detail['mu2'][1]	= 0.
	detail['mu2'][2]	= -1.*simple[2]
	detail['sigma2']	= [0]*NUM_OF_GAUSS
	detail['sigma2'][0]	= simple[4]
	detail['sigma2'][1]	= simple[6]
	detail['sigma2'][2]	= simple[4]
	# 共分散
	detail['kappa']	= [0]*NUM_OF_GAUSS
	detail['kappa'][0]	= simple[7]
	detail['kappa'][1]	= simple[8]
	detail['kappa'][2]	= simple[9]
	return detail

def klDivergence(comp, true):
	u"""KLダイバージェンスを計算する．
	【引数】comp: 比較対象の確率分布の値リスト，true: 真値の確率分布の値リスト
	【戻り値】ret: KLダイバージェンスの値"""
	if len(comp) != len(true):
		print("KLダイバージェンス計算時のlistのサイズが異なります．")
		sys.exit()
	ret	= 0.
	for i in range(0, len(true)):
		if (comp[i] == 0) or (true[i] == 0):
			continue
		else:
			ret	+= comp[i]*math.log(comp[i]/true[i])
	return ret

if __name__ == "__main__":
	# コマンドライン引数読み込み
	arg_l	= sys.argv[1]
	arg_a	= sys.argv[2]

	# シミュレーションファイル読み込み
	simX = []
	simY = []
	for line_s in open('/usr/local/src/master/dat/l='+str(arg_l)+'/a='+str(arg_a)+'/sim_y1pdf.dat'):
		col_s	= line_s.strip().split(' ')
		simX.append(float(col_s[0]))
		simY.append(float(col_s[1]))

	# 解析ファイル読み込み
	pop	= []
	for line_a in open('/usr/local/src/master/results/l='+str(arg_l)+'/dat_a='+str(arg_a)+'/ana_gsay1pdf.dat'):
		col_a	= line_a.strip().split(' ')
		ind	= Individual()
		elm	= 0
		# モーメント値
		tmp_m	= []
		for i in range(0, NUM_OF_MOMENT):
			tmp_m.append(float(col_a[i]))
			elm	+= 1
		ind.m	= tmp_m
		tmp_elm	= elm
		# 目的関数値
		tmp_o	= []
		for i in range(elm, tmp_elm+NUM_OF_MOMENTEQ):
			tmp_o.append(float(col_a[i]))
			elm	+= 1
		ind.o	= tmp_o
		tmp_elm	= elm
		# パラメータ
		simplePrm	= []
		for i in range(elm, tmp_elm+NUM_OF_PARAM):
			simplePrm.append(float(col_a[i]))
			elm	+= 1
		detailPrm	= getDetailParameterFromSimpleNotation(simplePrm)
		# 変位
		ind.disp	= simX
		# 変位応答分布
		ind.dPdf	= createDispPdf(simX, detailPrm)
		# 変位応答分布のKLダイバージェンス
		ind.dKL	= klDivergence(ind.dPdf, simY)

		pop.append(ind)
		del ind
		del simplePrm

	# キーを指定してソート
	pop	= sorted(pop, key=lambda x: x.m[0])

	# プロットする用のシミュレーションデータを作成
	tmp_simX	= []
	tmp_simY	= []
	for k, v in enumerate(simX):
		if k % 10 == 0:
			tmp_simX.append(simX[k])
			tmp_simY.append(simY[k])

	# プロット
	# KLダイバージェンスは8辺りを越えるとやばい
	outputPath	= "/usr/local/src/master/results"
	for i in range(0, len(pop)):
		if pop[i].dKL > 8:
			continue
		plt.plot(pop[i].disp, pop[i].dPdf)
		plt.plot(tmp_simX, tmp_simY, 'o')
		plt.yscale('log')
		filename	= "E[y1^2]=" + str(pop[i].m[0]) + ".png"
		plt.savefig(outputPath+"/l="+str(arg_l)+"/fig_a="+str(arg_a)+"/"+filename)
		plt.clf()