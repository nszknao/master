#!/bin/bash

export COMPILER_PATH=/usr/bin

if [ -e gatest.exe ]; then
        rm gatest.exe
fi
############ Compile-force #############################################
#g++ -g -O0 -lstdc++ ../test/testMain.cpp ../ga/src/Ga.cpp ../ga/src/GaCommon.cpp -o gatest.exe -lm -lgsl -lgslcblas -std=c++11
g++48 -Wall -g -O0 ../ga/src/SchwefelEllipsoidCMSA.cpp -o gatest.exe -lm -std=c++11 -lshark

if [ -e gatest.exe ]; then
        ./gatest.exe
fi
