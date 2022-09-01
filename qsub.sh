#!/bin/sh

#a=(0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0 2.0)
#a=(0.001 0.003 0.006 0.01 0.03 0.06 0.1)
#a=(0.001 0.003 0.006 0.01 0.03 0.06 0.1 0.3 0.6 1.0)
#a=(0.09 0.1 0.2 0.4 0.5)
phi=(0.86 0.90 1.2)
P=(1e-6 1e-7 1e-8 1e-9)
#a=(4.0 5.0)
icc -O3 SS_jam_constP.cpp  -o SS_jam_constP.out 

for ((i=0 ; i<3 ; i++))
do
    for ((j=0 ; j<4 ; j++))
    do qsub jam.bat ${phi[i]} ${P[j]}
#do rm ${a[i]}/*
    done
done