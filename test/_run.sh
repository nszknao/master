#!/bin/bash

export COMPILER_PATH=/usr/bin

if [ -e gatest.exe ]; then
        rm gatest.exe
fi
############ Compile-force #############################################
# g++ -g -O0 -lstdc++ ../test/testMain.cpp ../ga/src/Ga.cpp ../ga/src/GaCommon.cpp -o gatest.exe -lm -lgsl -lgslcblas -std=c++11
g++48 -Wall -g -O0 -I /opt/shark-2.3.4/usr/local/include ../ga/src/nsga2Sample_Shark.cpp -o gatest.exe -lm -std=c++11 -lshark -lpthread
# g++48 -I /usr/include/paradiseo/moeo -I /usr/include/paradiseo/mo -I /usr/include/paradiseo/eo -Wall -g -O0 ../ga/src/Sch1.cpp -o gatest.exe -std=c++0x -lm -lcma -leo -leoutils -les -lga -lmoeo

if [ -e gatest.exe ]; then
        ./gatest.exe
fi
