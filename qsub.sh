#!/bin/bash
# Copyright (c) 2017 Yanze Li, Florian Prodinger
# This software is released under the MIT License, see LICENSE.
#
# Name:qsub.sh
# Function:Run pipeline via qsub system
# E-mail:yanzeli@kuicr.kyoto-u.ac.jplace
# Initial version: 2017-05-18

# commits below are for qsub
#
#PBS -q cdb
#PBS -j oe
#PBS -o /user1/scl1/yanzeli/qsub.out
#PBS -l mem=400gb
#======================================================


#Please change path which in PBS comment part into your path in ABSOLUTE path.
#you can relative path in main text below

# sh ~/aptmp/Paper_pipeline/Scripts/Pipeline.sh -in /aptmp/yanzeli/Paper_pipeline/Sources/ -out /aptmp/yanzeli/Paper_pipeline/Outputs_170901/ -t 28
sh ~/aptmp/Paper_pipeline/Scripts/Pipeline.sh -in /aptmp/yanzeli/Paper_pipeline/Test_Sources/ -out /aptmp/yanzeli/Paper_pipeline/Test_Outputs_0/ -t 28 &
sh ~/aptmp/Paper_pipeline/Scripts/Pipeline.sh -in /aptmp/yanzeli/Paper_pipeline/Test_Sources/ -out /aptmp/yanzeli/Paper_pipeline/Test_Outputs_1/ -t 28