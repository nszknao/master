# -*- coding: utf-8 -*- 
from scipy import stats

if __name__ == "__main__":
	# シミュレーションファイル読み込み
	simX = []
	simY = []
	for line_s in open('/usr/local/src/master/results/l=0.05/dat_a=0.7/y1_pdf.dat'):
		col_s	= line_s.strip().split(' ')
		simX.append(float(col_s[0]))
		simY.append(float(col_s[1]))

	# 解析ファイルをanaInfoに格納
	# --- anaInfoの形式 ---
	# {
	#  'P': [1,2,...,10],
	#  'X': [1,2,...,10],
	#  'Y': [1,2,...,10]
	# }
	# --- anaInfoの形式 ---
	anaInfo	= {}
	cnt	= -1
	for line_a in open('/usr/local/src/master/results/l=0.05/dat_a=0.7/gsay1pdf.dat'):
		col_a	= line_a.strip().split(' ')
		if col_a[0] == '#':
			cnt	+= 1
			anaInfo[cnt]	= {}
			anaInfo[cnt]['X'] = []
			anaInfo[cnt]['Y'] = []

			col_a.remove('#')
			anaInfo[cnt]['P']	= map(float, col_a)
		else:
			anaInfo[cnt]['X'].append(float(col_a[0]))
			anaInfo[cnt]['Y'].append(float(col_a[1]))

	print "sim:" + str(len(simX)) + ", ana:" + str(len(anaInfo[0]['X']))

	# print(stats.entropy(anaInfo[0]['X'], simX))