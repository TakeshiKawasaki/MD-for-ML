#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M kawasaki@r.phys.nagoya-u.ac.jp
#$ -m ea
#$ -V
#
#$ -q all.q@lemon
#$ -q all.q@mango
#$ -q all.q@orange
#
#$ -N L60Heltz

./SS_jam_constP.out $1 $2