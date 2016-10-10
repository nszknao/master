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
MOMENT_ARRAY	= ["$E[y_1^2]$", "$E[y_1y_2]$", "$E[y_2^2]$", "$E[y_1^4]$", "$E[y_1^3y_2$", "$E[y_1^2y_2^2]$", "$E[y_1y_2^3]$", "$E[y_2^4]$", "$E[y_1^6]$", "$E[y_1^5y_2]$", "$E[y_1^4y_2^2]$", "$E[y_1^3y_2^3]$", "$E[y_1^2y_2^4]$", "$E[y_1y_2^5]$", "$E[y_2^6]$", "$E[y_1^8]$", "$E[y_1^7y_2]$", "$E[y_1^6y_2^2]$", "$E[y_1^5y_2^3]$", "$E[y_1^4y_2^4]$", "$E[y_1^3y_2^5]$"]
DISP_OBJ_NUMBER	= [1, 4 ,9]
# DISP_OBJ_NUMBER	= [15, 16, 17, 18, 19, 20]

dx	= 0.01	# pdfpdfのX軸刻み幅

class Individual:
	u"""個体情報を格納
	メンバ: m[], o[], disp[], vel[], dPdf[], vPdf[], dKL, vKL
	Pythonでは，代入すれはメンバが生成される．元からメンバを明記しておくと識別子が同一として扱われてしまう．"""
	pass

def drange(begin, end, step):
	u"""小数刻みでrangeできるように拡張
	【引数】begin: 開始の値，end: 終了の値，step: ステップ幅"""
	n = begin
	while n+step < end:
		yield n
		n += step

def culcIntegralPdf(pdf):
	u"""PDFの積分値を計算します．
	【引数】pdf: 確率分布の値リスト
	【戻り値】integration: 積分値"""
	integrantion	= 0.
	for i in range(0, len(pdf)):
		integrantion	+= pdf[i]*dx
	return integrantion

def createDispPdf(x, prm):
	u"""パラメータから変位応答を計算します．
	【引数】x: はシミュレーションから取得，prm: 詳細のパラメータ
	【戻り値】pdf: 変位応答分布の値リスト"""
	pdf	= []
	for i in range(0, len(x)):
		tmp_pdf	= 0.
		for ii in range(0, NUM_OF_GAUSS):
			tmp_pdf += prm['a'][ii]*(1./math.sqrt(2.*math.pi)/prm['sigma1'][ii])*math.exp(-1.*(x[i]-prm['mu1'][ii])**2/2./prm['sigma1'][ii]**2)
		pdf.append(tmp_pdf)
	return pdf

def createVelPdf(x, prm):
	u"""パラメータから変位応答を計算します．
	【引数】x: はシミュレーションから取得，prm: 詳細のパラメータ
	【戻り値】pdf: 速度応答分布の値リスト"""
	pdf	= []
	for i in range(0, len(x)):
		tmp_pdf	= 0.
		for ii in range(0, NUM_OF_GAUSS):
			tmp_pdf += prm['a'][ii]*(1./math.sqrt(2.*math.pi)/prm['sigma2'][ii])*math.exp(-1.*(x[i]-prm['mu2'][ii])**2/2./prm['sigma2'][ii]**2)
		pdf.append(tmp_pdf)
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
	simDispX = []
	simDispY = []
	simDispXForPlot	= []
	simDispYForPlot	= []
	lineCount	= 0
	for line_s in open('/usr/local/src/master/dat/l='+str(arg_l)+'/a='+str(arg_a)+'/sim_y1pdf.dat'):
		col_s	= line_s.strip().split(' ')
		simDispX.append(float(col_s[0]))
		simDispY.append(float(col_s[1]))
		if lineCount % 10 == 0:
			simDispXForPlot.append(float(col_s[0]))
			simDispYForPlot.append(float(col_s[1]))
		lineCount	+= 1

	# シミュレーションファイル読み込み
	# simVelX = []
	# simVelY = []
	# for line_s in open('/usr/local/src/master/dat/l='+str(arg_l)+'/a='+str(arg_a)+'/sim_y2pdf.dat'):
	# 	col_s	= line_s.strip().split(' ')
	# 	simVelX.append(float(col_s[0]))
	# 	simVelY.append(float(col_s[1]))

	# プロット用解析解X軸
	anaDispX	= []
	for i in drange(-6, 6, 0.01):
		anaDispX.append(i)
	# anaVelX	= []
	# for i in drange(-12, 12, 0.01):
	# 	anaVelX.append(i)

	# 解析ファイル読み込み
	pop	= []
	for line_a in open('/usr/local/src/master/results/l='+str(arg_l)+'/dat_a='+str(arg_a)+'/ana_gsay1pdf.dat'):
		col_a	= line_a.strip().split(' ')
		ind	= Individual()
		#--- モーメント値 ---
		ind.m	= []
		for i in range(0, NUM_OF_MOMENT):
			ind.m.append(float(col_a[i]))
		#--- 目的関数値 ---
		ind.o	= []
		for i in range(NUM_OF_MOMENT, NUM_OF_MOMENT + NUM_OF_MOMENTEQ):
			ind.o.append(float(col_a[i]))
		#--- パラメータ ---
		simplePrm	= []
		for i in range(NUM_OF_MOMENT + NUM_OF_MOMENTEQ, NUM_OF_MOMENT + NUM_OF_MOMENTEQ + NUM_OF_PARAM):
			simplePrm.append(float(col_a[i]))
		detailPrm	= getDetailParameterFromSimpleNotation(simplePrm)
		#--- 変位 ---
		ind.disp	= anaDispX
		#--- 速度 ---
		# ind.vel	= anaVelX
		#--- 変位応答分布 ---
		dPdfForKL	= createDispPdf(simDispX, detailPrm)
		ind.dPdf	= createDispPdf(anaDispX, detailPrm)
		#--- 速度応答分布 ---
		# vPdfForKL	= createDispPdf(simVelX, detailPrm)
		# ind.dPdf	= createVelPdf(anaVelX, detailPrm)
		#--- 変位応答分布のKLダイバージェンス ---
		ind.dKL	= klDivergence(dPdfForKL, simDispY)
		#--- 速度応答分布のKLダイバージェンス ---
		# ind.vKL	= klDivergence(vPdfForKL, simDispY)
		
		pop.append(ind)
		del ind

	### for 中間発表<<ここから
	### キーを指定してソート
	# pop	= sorted(pop, key=lambda x: x.dKL)
	# # KLダイバージェンスと相関のあるモーメントを見つけよう
	# dKL	= []
	# obj	= np.zeros((NUM_OF_MOMENTEQ, len(pop)))
	# moment	= np.zeros((NUM_OF_MOMENT, len(pop)))
	# min_dKL	= 100
	# min_obj	= [100000]*NUM_OF_MOMENTEQ
	# minObjPop	= [Individual()]*NUM_OF_MOMENTEQ
	# for i in range(0, len(pop)):
	# 	# 散布図プロット用
	# 	dKL.append(pop[i].dKL)
	# 	for ii in range(0, NUM_OF_MOMENT):
	# 		moment[ii, i]	= pop[i].m[ii]
	# 	for ii in range(0, NUM_OF_MOMENTEQ):
	# 		obj[ii, i]	= pop[i].o[ii]
	# 	# KLダイバージェンス最小個体を見つける用
	# 	if pop[i].dKL < min_dKL:
	# 		min_dKL	= pop[i].dKL
	# 		minDKLPop	= pop[i]
	# 	# 各目的関数を最小にしている個体を見つける
	# 	for ii in range(0, NUM_OF_MOMENTEQ):
	# 		if abs(pop[i].o[ii]) < min_obj[ii]:
	# 			min_obj[ii]	= abs(pop[i].o[ii])
	# 			minObjPop[ii]	= pop[i]
	### 各モーメントとの相関ランキング
	# cf	 = []
	# for i in range(0, NUM_OF_MOMENT):
	# 	if np.cov(moment[i]) == 0:
	# 		cf.append(0)
	# 	else:
	# 		cf.append(np.corrcoef(dKL, moment[i])[0,1])
	# 	print(cf[i])
	# for i in range(0, NUM_OF_MOMENTEQ):
	# 	if np.cov(obj[i]) == 0:
	# 		cf.append(0)
	# 	else:
	# 		cf.append(np.corrcoef(dKL, obj[i])[0,1])
	# 	print(cf[i])
	# plt.figure(figsize=(16, 5))
	# plt.plot(range(0, NUM_OF_MOMENT), cf)
	# plt.ylim([-1, 1])
	# plt.grid(which='major',color='black',linestyle='--')
	# # plt.xticks(range(0, NUM_OF_MOMENT), MOMENT_ARRAY)
	# plt.xlabel("Number of moment")
	# plt.ylabel("Correlation coefficient")
	# plt.savefig("/usr/local/src/master/fig/corrcoef.png")
	### KLダイバージェンスが最小の個体
	# print("Min dKL = " + str(min_dKL))
	# for i in range(0, NUM_OF_MOMENTEQ):
	# 	print("Moment No." + str(i) + " = " + str(minDKLPop.o[i]))
	## 散布図をプロット
	# plt.scatter(dKL, moment)
	# plt.xlabel("KL-Divergence")
	# plt.ylabel("Moment")
	# plt.savefig("/usr/local/src/master/fig/dKL-moment_scatter.png")
	## 各目的関数を最小にしている個体情報を表示
	# for i in range(0, len(minObjPop)):
	# 	print("Moment No." + str(i) + ": dKL = " + str(minObjPop[i].dKL) + ", object = " + str(minObjPop[i].o[i]))
	## dKLが上位の個体情報を表示
	# num	= 5
	# for i in range(0, num):
	# 	print("Pop No." + str(i) + ": dKL = " + str(pop[i].dKL))
	# 	plt.plot(range(0, NUM_OF_MOMENTEQ), pop[i].o, label=i)
	# 	plt.xlabel("Number of MomentEQ")
	# 	plt.ylabel("Value of MomentEQ")
	# 	# plt.xticks(np.arange(0, NUM_OF_MOMENTEQ, 1))
	# 	plt.ylim([-140, 20])
	# 	plt.grid(which='major',color='black',linestyle='--')
	# 	plt.legend(loc='lower left')
	# 	plt.savefig("/usr/local/src/master/fig/momentEQ5.png")
	## 指定した目的関数値の和の合計が最小の個体を選ぶ
	# minSumObj	= 100000000.
	# targetInd	= Individual()
	# for i in range(0, len(pop)):
	# 	SumObj	= 0.
	# 	for number in range(0, NUM_OF_MOMENTEQ):
	# 		SumObj	+= pop[i].o[number]**2
	# 	if SumObj < minSumObj:
	# 		targetInd	= pop[i]
	# 		minSumObj	= SumObj
	# print(minSumObj)
	# print("dKL = " + str(targetInd.dKL))
	# # plt.plot(range(0, NUM_OF_MOMENTEQ), targetInd.o)
	# # plt.ylim([-100, 20])
	# plt.plot(simDispXForPlot, simDispYForPlot, 'o')
	# plt.plot(targetInd.disp, targetInd.dPdf)
	# plt.savefig("/usr/local/src/master/fig/pdf.png")
	### for 中間発表>>ここまで


	# プロット
	# KLダイバージェンスの理想は1を下回っているもの
	outputPath	= "/usr/local/src/master/results"
	for i in range(0, len(pop)):
		int_dPdf	= culcIntegralPdf(pop[i].dPdf)
		if (0.95 > int_dPdf) or (int_dPdf > 1.05):
			print("PDF integral violation: " + str(int_dPdf))
			continue
		if pop[i].dKL > 100:
			continue
		plt.xlim([-6, 6])
		plt.ylim([0.00001, 1])
		plt.xlabel("displacement", fontsize=18)
		plt.ylabel("Probability density", fontsize=18)
		plt.grid(which='major',color='black',linestyle='--')
		plt.plot(pop[i].disp, pop[i].dPdf)
		plt.plot(simDispXForPlot, simDispYForPlot, 'o')
		plt.legend(['Analytical solution', 'Simulation'], loc='lower center')
		plt.yscale('log')
		filename	= "KL=" + str(pop[i].dKL) + ", E[y1^2]=" + str(pop[i].m[0]) + ", E[y1^4]=" + str(pop[i].m[3]) + ", E[y1^6]=" + str(pop[i].m[8]) + ".png"
		plt.savefig(outputPath+"/l="+str(arg_l)+"/fig_a="+str(arg_a)+"/"+filename)
		plt.clf()