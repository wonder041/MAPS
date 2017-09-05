import fileinput
import itertools

pre='''#!/bin/bash
# Copyright (c) 2017 Yanze Li, Florian Prodinger
# This software is released under the MIT License, see LICENSE.
#
# Name:qsub.sh
# Function:Run pipeline via qsub system
# E-mail:yanzeli@kuicr.kyoto-u.ac.jplace
# Initial version: 2017-05-18

# commits below are for qsub
#
#PBS -q sp
#PBS -j oe
#PBS -o /lustre1/aptmp/ideas2/yanzeli/qlog
#PBS -l mem=100gb
#======================================================
'''


for line_arr in itertools.zip_longest(*(fileinput.input(),)*5):
    print(pre+" & ".join(line for line in line_arr if line))
    break


# line_arr=list(fileinput.input())
# print(line_arr)
# itertools.zip_longest(a)