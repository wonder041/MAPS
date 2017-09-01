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

sh /aptmp/yanzeli/Paper_pipeline/Scripts/Pipeline.sh -in /aptmp/yanzeli/Megaviridae/Source_FL/ -out /aptmp/yanzeli/Megaviridae/Outputs_170821/ -t 2
