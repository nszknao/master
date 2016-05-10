#!/bin/bash

LD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH

if [ -e gatest.exe ]; then
        rm gatest.exe
fi
############ Compile-force #############################################
g++ -lstdc++ ../ga/GA.cpp main.cpp -o gatest.exe -lm -lgsl -lgslcblas -std=c++11

if [ -e gatest.exe ]; then
        ./gatest.exe
fi