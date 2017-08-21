#!/bin/bash
# Copyright (c) 2017 Yanze Li
# This software is released under the MIT License, see LICENSE.
#
# Name:BestTrim.sh
# Function:find best trimming method of Megaviridae amplicon
# E-mail:yanzeli@kuicr.kyoto-u.ac.jplace
# Initial version: 2017-05-18

# commits below are for qsub
#
#PBS -q cdb
#PBS -j oe
#PBS -o ~/Mega_pipe_strout
#PBS -l mem=500gb
#==============================================================================


READ_PATTERN[0]="_R1.fastq"
READ_PATTERN[1]="_R2.fastq"

TRIMMOMATIC_PATH="$1 $2 $3"
FLASH_PATH=$4
INPUT_PATH=$5
OUTPUT_DIR=$6
MERGE_ERROR_RATE=$7
MERGE_OVERLAP_MIN=$8

BARCODE=$(echo "${INPUT_PATH}"|rev|cut -d/ -f1|rev)
OUTPUT_PATH=${OUTPUT_DIR}${BARCODE}.fastq
LOG_PATH=${OUTPUT_DIR}${BARCODE}/log

mkdir -p ${OUTPUT_DIR}${BARCODE}/

TRIM_MERGE(){
    mkdir -p ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/
   ${TRIMMOMATIC_PATH} PE -threads 2 -phred33 \
   ${INPUT_PATH}/out.notCombined_1.fastq \
   ${INPUT_PATH}/out.notCombined_2.fastq \
   ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/Trim${READ_PATTERN[0]} \
   /dev/null \
   ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/Trim${READ_PATTERN[1]} \
   /dev/null \
   SLIDINGWINDOW:${SIZE_NOW}:${QUAL_NOW} \
   &> /dev/null
   ${FLASH_PATH} \
    -t 2 \
    -m 20 \
    -M 300 \
    -x ${MERGE_ERROR_RATE} \
    ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/Trim${READ_PATTERN[0]} \
    ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/Trim${READ_PATTERN[1]} \
    -d ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/ > /dev/null

}

SEARCH_BEST_SIZE(){
SIZE_STEP=8
SIZE_END=800
SIZE_ORI=1
NUM_PRE=0
SIZE_NOW=1
SIZE_STEP_SIZE_ORI=`expr $SIZE_STEP \* $SIZE_ORI`

while [[ ture ]] ;do
    if [[ ${SIZE_NOW} -gt 500 ]] ;then
        echo "SIZE over 100 at ${QUAL_NOW} !" >> ${LOG_PATH}
        exit
    fi
        
    if [[ ! -f ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/out.extendedFrags.fastq ]] ;then
        TRIM_MERGE
    fi
    
    NUM_NOW=$(cut -f2 ${OUTPUT_DIR}${BARCODE}/${QUAL_NOW}_${SIZE_NOW}/out.hist|awk '{ sum += $1 }; END { print sum }')
    echo "QUAL $QUAL_NOW SIZE $SIZE_NOW NUM $NUM_NOW" >> ${LOG_PATH}
    
    if [[ $NUM_NOW -lt $NUM_PRE ]] ;then
        if [[ $SIZE_STEP -eq 1 ]] ;then
            SIZE_BEST=`expr $SIZE_NOW - $SIZE_STEP_SIZE_ORI`
            return
        else
            TWO_SIZE_STEP=`expr 2 \* $SIZE_STEP_SIZE_ORI`
            SIZE_END=`expr $SIZE_NOW - $TWO_SIZE_STEP`
            SIZE_STEP=`expr $SIZE_STEP / 2`
            SIZE_ORI=`expr $SIZE_ORI / -1`
            SIZE_STEP_SIZE_ORI=`expr $SIZE_STEP \* $SIZE_ORI`
            
            SIZE_NOW=`expr $SIZE_NOW + $SIZE_STEP_SIZE_ORI`
            NUM_PRE=$NUM_NOW
            continue
        fi
    else
        if [[ $SIZE_NOW -eq 1 ]] ;then
            if [[ $SIZE_ORI -eq -1 ]] ;then
                SIZE_BEST=$SIZE_NOW
                NUM_PRE=$NUM_NOW
                return
            fi
        fi
        SIZE_NOW=`expr $SIZE_NOW + $SIZE_STEP_SIZE_ORI`
        NUM_PRE=$NUM_NOW
    fi
done
}

SEARCH_BEST_QUAL(){
QUAL_STEP=8
QUAL_END=40
QUAL_ORI=1
NUM_QUAL_PRE=0
QUAL_NOW=9
QUAL_STEP_QUAL_ORI=`expr $QUAL_STEP \* $QUAL_ORI`

while [[ ture ]] ;do

    if [[ ${QUAL_NOW} -gt 500 ]] ;then
        echo "QUAL over 100!" >> ${LOG_PATH}
        exit
    fi
        
    SEARCH_BEST_SIZE
    NUM_QUAL_NOW=$NUM_PRE

    echo "BEST SIZE of QUAL $QUAL_NOW is $SIZE_BEST with $NUM_QUAL_NOW" >> ${LOG_PATH}
    
  
    if [[ $NUM_QUAL_NOW -lt $NUM_QUAL_PRE ]] ;then
        if [[ $QUAL_STEP -eq 1 ]] ;then
            QUAL_BEST=`expr $QUAL_NOW - $QUAL_STEP_QUAL_ORI`
            echo "$BARCODE BSET QUAL $QUAL_BEST BEST SIZE $SIZE_BEST_PRE BEST NUMBER $NUM_QUAL_PRE" >> ${LOG_PATH}
            cp ${OUTPUT_DIR}${BARCODE}/${QUAL_BEST}_${SIZE_BEST_PRE}/out.extendedFrags.fastq ${OUTPUT_PATH}
            break
        else
            TWO_QUAL_STEP=`expr 2 \* $QUAL_STEP_QUAL_ORI`
            QUAL_END=`expr $QUAL_NOW - $TWO_QUAL_STEP`
            QUAL_STEP=`expr $QUAL_STEP / 2`
            QUAL_ORI=`expr $QUAL_ORI / -1`
            QUAL_STEP_QUAL_ORI=`expr $QUAL_STEP \* $QUAL_ORI`
            
            QUAL_NOW=`expr $QUAL_NOW + $QUAL_STEP_QUAL_ORI`
            NUM_QUAL_PRE=$NUM_QUAL_NOW
        fi
    else
        if [[ $QUAL_NOW -lt 0 ]] ;then
            QUAL_BEST=`expr $QUAL_NOW - $QUAL_STEP_QUAL_ORI`
            echo "$BARCODE No trim is the Best" >> ${LOG_PATH}
            echo "$BARCODE BSET QUAL $QUAL_BEST BEST SIZE $SIZE_BEST_PRE BEST NUMBER $NUM_QUAL_PRE" >> ${LOG_PATH}
            cp ${OUTPUT_DIR}${BARCODE}/${QUAL_BEST}_${SIZE_BEST_PRE}/out.extendedFrags.fastq ${OUTPUT_PATH}
            break
        fi
        QUAL_NOW=`expr $QUAL_NOW + $QUAL_STEP_QUAL_ORI`
        if [[ $QUAL_NOW -gt 35 ]] ;then
            QUAL_NOW=`expr $QUAL_NOW - $QUAL_STEP_QUAL_ORI`
            TWO_QUAL_STEP=`expr 2 \* $QUAL_STEP_QUAL_ORI`
            QUAL_END=`expr $QUAL_NOW - $TWO_QUAL_STEP`
            QUAL_STEP=`expr $QUAL_STEP / 2`
            QUAL_ORI=`expr $QUAL_ORI / -1`
            QUAL_STEP_QUAL_ORI=`expr $QUAL_STEP \* $QUAL_ORI`
            
            QUAL_NOW=`expr $QUAL_NOW + $QUAL_STEP_QUAL_ORI`
        fi
        NUM_QUAL_PRE=$NUM_QUAL_NOW
    fi
    SIZE_BEST_PRE=$SIZE_BEST
done
}

SEARCH_BEST_QUAL
