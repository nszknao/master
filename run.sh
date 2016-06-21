#!/bin/bash
# ----------
# banalyze.sh : compile bimodal-force.c and brunge.c
# ----------

############ Compile-force #############################################
g++48 -lstdc++ -Wl,--heap,900000000,--stack,536870912 ./src/research.cpp -o research.exe -lm -lgsl -lgslcblas
g++48 -lstdc++ -Wl,--heap,900000000,--stack,536870912 ./src/analysis.cpp -o lsm.exe -lm -lgsl -lgslcblas

############ Reading init_value.txt ####################################
initfile=init_value_disp.txt
for line in `cat ${initfile} | grep -v ^#`
do
 lambda=`echo ${line} | cut -d ',' -f1`
 alpha=`echo ${line} | cut -d ',' -f2`
 lininit=`echo ${line} | cut -d ',' -f3`
 loginit=`echo ${line} | cut -d ',' -f4`

tmp1=`echo "scale=7; 1 / ${lambda}" | bc`
beta2=${tmp1}

echo "lambda = " $lambda "beta2 = " $beta2

##### Handling if "results" exists. #############
results="results"

if [ ! -e ./${results}/ ]
then
mkdir ./${results}/
fi

##### Handling if "results" exists. #############
params="l=${lambda},b2=${beta2}"

X=`echo "${alpha} * 10" | bc`
XX=${X%.*}
echo ${XX}
if [[ $XX -eq 1 ]]
then
if [ -e ./${results}/${params}/ ]
then
rm -Rf ./${results}/${params}/*
else
mkdir ./${results}/${params}/
fi
fi

############ Plot ######################################################
##### Handlings if a directry "dat" exists. #####
dat="dat_a=${alpha}"

if [ ! -e ${dat} ]
then
mkdir ${dat}
fi

##### Handlings if a directry "eps" exists. ######
eps="eps_a=${alpha}"

if [ ! -e ${eps} ] 
then
mkdir ${eps}
mkdir ${eps}/"displacement(log)"
mkdir ${eps}/"displacement(linear)"
mkdir ${eps}/"velocity(log)"
mkdir ${eps}/"velocity(linear)"
fi

##### Handling if "bplot.plt" exists. #############
plot="plot.plt"

if [ -e ${plot} ]
then
rm ${plot}
touch ${plot}
else
touch ${plot}
fi

#### Making dat file only at first ###############
./research.exe ${lambda} ${beta2} ${alpha}

#### Run GSA #####################################
./lsm.exe ${lambda} ${beta2} ${alpha} ${loginit}

##### value of kurtosis ###########################
#exec < "sim_y1_var.dat"
#read kurt_sim
#echo ${kurt_sim}
#
#if [ -e "gsay1pdf.dat" ]
#then
#
#exec < "anl_y1_var.dat"
#read kurt_anl
#echo ${kurt_anl}
#
#var_ratio=`echo "scale=7; ${kurt_sim} / ${kurt_anl}" | bc`
#
#if [ "$(echo "${var_ratio} > 0.9900" | bc)" -eq 1 -a "$(echo "${var_ratio} < 1.010" | bc)" -eq 1 ]
#then

##### A setup of figure ###########################
cat <<EOF >${plot}
set terminal postscript eps enhanced color
set style line 1 lt 1 lw 2
set size 0.45,0.55
unset key
EOF

##入力をプロットする
##cat <<EOF >>${plot}
##set output "${eps}/excitation.eps"
##set xlabel "time{/Symbol-Oblique t}/2{/Symbol-Oblique p}"
##set ylabel "input"
##p [0:50][-50:50] "./t_force.dat" with lines linetype 1 notitle
##EOF
#
##変位応答をプロットする
#cat <<EOF >>${plot}
#set size 0.6,0.4
#set output "${eps}/response.eps"
##set xlabel "time {/Italic-Times t}"
##set ylabel "displacement {/Italic-Times y_1}"
#set ytics 5
#p [0:100][-5:5] "./t_x1.dat" with lines linetype 1 lw 3, 1.754543 with l lt 1 lw 3 lc rgb "blue", -1.754543 with l lt 1 lw 3 lc rgb "blue"
#EOF

##変位応答の確率密度関数(リニア軸)をプロットする
#cat <<EOF >>${plot}
#set size 0.46,0.54 # 図の縦横比
#set output "${eps}/displacement(linear)/pdf_of_y1_mu1=${loginit}.eps"
#set xlabel "displacement {/Italic-Times y_1}"
#set ylabel "probability density {/Italic-Times p_{Y_1}}"
#set xtics 2
#set ytics 0.1
#p [-6:6][0:0.4] "./y1_pdf.dat" with p pt 6 ps 1,\
#		"./gsay1pdf.dat" with lines linetype 1 lw 4 lc rgb "blue" ,\
#		"./Gpdf/y1_Gpdf.dat" with lines linetype 2 lw 4
#EOF

##速度応答の確率密度関数(リニア軸)をプロットする
#cat <<EOF >>${plot}
#set output "${eps}/velocity(linear)/pdf_of_y2_mu1=${loginit}.eps"
#set xlabel "velocity {/Italic-Times y_2}"
#set ylabel "probability density {/Italic-Times p_{Y_2}}"
#set xtics 4
#set ytics 0.1
#p [-12:12][0:0.4]  "./gsay2pdf.dat" with lines linetype 1 lw 2 lc rgb "blue",\
#		"./y2_pdf.dat" every 2 with p  pt 6 ps 1 lc rgb "red",\
#		"./Gpdf/y2_Gpdf.dat" with lines linetype 2 lw 2
#EOF


#変位応答の確率密度関数(対数軸)をプロットする
cat <<EOF >>${plot}
set output "${eps}/displacement(log)/log_pdf_of_y1_mu1=${loginit}.eps"
set xlabel "displacement {/Italic-Times y_1}"
set ylabel "probability density {/Italic-Times p_{Y_1}}"
set xtics 2
set ytics 0.1
set logscale y
set format y "10^{%L}"
p [-6:6][0.00001:1.0] "./y1_pdf.dat" with p pt 6 ps 1 lc rgb "blue",\
		"./gsay1pdf.dat" with lines linetype 1 lw 4 lc rgb "red",\
		"./Gpdf/y1_Gpdf.dat" with lines linetype 2 lw 4
EOF

##速度応答の確率密度関数(対数軸)をプロットする
#cat <<EOF >>${plot}
#set output "${eps}/velocity(log)/log_pdf_of_y2_mu1=${loginit}.eps"
#set xlabel "velocity {/Italic-Times y_2}"
#set ylabel "probability density {/Italic-Times p_{Y_2}}"
#set xtics 4
#set ytics 0.1
#set logscale y
#set format y "10^{%L}"
#p [-12:12][0.00001:1.0]  "./gsay2pdf.dat" with lines linetype 1 lw 2 lc rgb "blue",\
#			"./y2_pdf.dat" every 2 with p pt 6 ps 1 lc rgb "red",\
#			"./Gpdf/y2_Gpdf.dat" with lines linetype 2 lw 2
#EOF

gnuplot ${plot}
echo -e "Plotting is success!!\n"

#rm "gsay1pdf.dat"
#fi
#fi
#### The sign of an end ####################

##### Handling if "eps" exists in "results". #############
if [ ! -e ./${results}/${params}/"eps" ]
then
mkdir ./${results}/${params}/"eps"
fi

##### Handling if "dat" exists in "results". #############
if [ ! -e ./${results}/${params}/"dat" ]
then
mkdir ./${results}/${params}/"dat"
fi

############ Move directory ##############################
mv ./${eps} ./${results}/${params}/"eps"
mv ./*.dat ${dat}/
mv ./${dat} ./${results}/${params}/"dat"


echo -e "a=${alpha}, mu1_width=${loginit} was processed.\n\n"
done
