#!/bin/bash
# Copyright (c) 2017 Yanze Li, Florian Prodinger
# This software is released under the MIT License, see LICENSE.
#
# Name:Pipeline.sh
# Function:pipeline for megaviridae anplify analysis
# E-mail:yanzeli@kuicr.kyoto-u.ac.jplace
# Initial version: 2017-05-18
# Update 2018-08-03 add chimeric sequence check


# commits below are for qsub
#
#PBS -q cdb
#PBS -j oe
#PBS -o ~/Mega_pipe_strout
#PBS -l mem=400gb
#======================================================

echo "-Start pipeline"
VERSION="2.1"

Usage(){
echo "MegaPipeline V${VERSION}"
echo "Usage under editing"
}

#read arguments
while [ "$1" != "" ]; do
    case $1 in
        -i | -in )              shift
                                SOURCES_DIR=$1
                                if [ "${SOURCES_DIR: -1}" != "/" ];then
                                    SOURCES_DIR="${SOURCES_DIR}/"
                                fi
                                ;;
        -o | -out )             shift
                                OUTPUTS_DIR=$1
                                if [ "${OUTPUTS_DIR: -1}" != "/" ];then
                                    OUTPUTS_DIR="${OUTPUTS_DIR}/"
                                fi
                                ;;
        -t | --threads )        shift
                                THREADS=$1
                                ;;
        --diamond )             DIAMOND=1
                                ;;
        --env_script )          shift
                                ENV_SCRIPT=$1
                                source ${ENV_SCRIPT}
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        -version | --version )  echo "${VERSION}"
                                exit
                                ;;
        * )                     echo "!!Wrong argument!!"
                                usage
                                exit 1
    esac
    shift
done

export PATH=/usr/appli/freeware/bin/:${PATH}

#start time recorder
TIME_P=$(date +%s)

#build environment
{
echo "--Start building environment"

#build Set_default Check file
Set_default_check_file(){
#$1:target variable $2:value
eval $1='$'{$1:-$2}
#if path not exist then exit
eval "if [ ! -f "'$'"$1 ]; then echo !!Envrionment ERROR, $1 = "'$'"$1 not found!!;exit 1;fi"
}

#set main directory
MAIN_DIR=${MAIN_DIR:-"/aptmp/yanzeli/Paper_pipeline/"}

#build sub directorys
SOURCES_DIR=${SOURCES_DIR:-"${MAIN_DIR}Sources/"}
REFERENCES_DIR=${REFERENCES_DIR:-"${MAIN_DIR}References/"}
SCRIPTS_DIR=${SCRIPTS_DIR:-"${MAIN_DIR}Scripts/"}
OUTPUTS_DIR=${OUTPUTS_DIR:-"${MAIN_DIR}Outputs/"}

#build step list
STEP_ARR="A5G40 CUTADAPT_G40 MERGE DEDUPLICATION FAA BLASTP ALIGNMENT"
STEP_COUNTER=1
for STEP_NAME in ${STEP_ARR};do
    OUTPUT_DIR_ARR[${STEP_COUNTER}]=${OUTPUT_DIR_ARR[${STEP_COUNTER}]:-"${OUTPUTS_DIR}${STEP_COUNTER}_${STEP_NAME}/"}
    let "STEP_COUNTER +=1"
done

#build read files pattern
READ_PATTERN[0]=${READ_PATTERN[0]:-"_R1.fastq"}
READ_PATTERN[1]=${READ_PATTERN[1]:-"_R2.fastq"}

#build build-in script path
Set_default_check_file SCR_DECODE_PRIMER ${SCRIPTS_DIR}Decode_primer.py
Set_default_check_file SCR_TRIM_MERGE_O ${SCRIPTS_DIR}Trim_merge_O.py
Set_default_check_file SCR_MERGE_RESURE ${SCRIPTS_DIR}Merge_Rescue.sh
Set_default_check_file SCR_SEQ_CONVERTER ${SCRIPTS_DIR}Seq_convert.py
Set_default_check_file SCR_TRANSLATE_RENAME ${SCRIPTS_DIR}Translate_rename.py
Set_default_check_file SCR_TRIM_COMMON_REGION ${SCRIPTS_DIR}Trim_common_region.py
Set_default_check_file SCR_PPLACER_DECODE ${SCRIPTS_DIR}Pplacer_decode.py
Set_default_check_file SCR_PPLACER_FNA_ID ${SCRIPTS_DIR}Pplacer_fna_id.py

#build build-in reference path
Set_default_check_file REF_PRIMER ${REFERENCES_DIR}primer.fna
Set_default_check_file REF_MIXTURE ${REFERENCES_DIR}mixture.txt
Set_default_check_file REF_POLB_HOMOLOGY_SEARCH ${REFERENCES_DIR}PolB_homology_search.faa
Set_default_check_file REF_DESIGN_PROT_ALN ${REFERENCES_DIR}PolB_design_prot.aln
Set_default_check_file REF_DESIGN_PRIMER_NUCL_ALN ${REFERENCES_DIR}PolB_design_primer_nucl.aln

Set_default_check_file REF_PPLACER_ALN ${REFERENCES_DIR}Pplacer.aln
Set_default_check_file REF_PPLACER_RES ${REFERENCES_DIR}Pplacer.res
Set_default_check_file REF_PPLACER_INFO ${REFERENCES_DIR}Pplacer.info

#set log
PIPELINE_LOG_PATH=${PIPELINE_LOG_PATH:-"${OUTPUTS_DIR}Pipeline.log"}

#build logger
Logger(){
echo $1 |tee -a ${PIPELINE_LOG_PATH}
}

#set threads
THREADS=${THREADS:-"8"}

#make outputs directory
mkdir -p ${OUTPUTS_DIR}

Logger "--Finish building environment"
}




#build software path
{
Logger "--Start building software path"

#build Set_default Check run
Set_default_check_executable(){
#$1:target variable $2:value
eval $1='$'{$1:-$2}
#if path not executable then exit
eval "if [ ! -x "'$'"$1 ]; then echo !!Software path ERROR, $1 = "'$'"$1 not found!!;exit 1;fi"
}

#build software path from user directory
Set_default_check_executable CUTADAPT_PATH /aptmp/yanzeli/miniconda3/envs/cutadapt/bin/cutadapt

#build software path from module
source /usr/share/modules/init/bash
module load "trimmomatic/0.35 pplacer/1.1.alpha19 mafft/7.310 flash/1.2.11 Python/3.6 R/3.3.1 cd-hit/4.6.8 blast+/2.5.0 hmmer/3.1b2 diamond/0.9.9 qiime/1.9.1" > /dev/null

TRIMMOMATIC_PATH=${TRIMMOMATIC_PATH:-"/usr/bin/java -jar /usr/appli/freeware/trimmomatic/0.35/trimmomatic-0.35.jar"}
Set_default_check_executable PPLACER_PATH /usr/appli/freeware/pplacer/1.1.alpha19/pplacer
Set_default_check_executable MAFFT_PATH /usr/appli/freeware/mafft/7.310/bin/mafft
Set_default_check_executable FLASH_PATH /usr/appli/freeware/flash/1.2.11/flash

Set_default_check_executable PYTHON_PATH /usr/appli/freeware/Python/3.6.0/bin/python
Set_default_check_executable R_PATH /usr/appli/freeware/R/3.3.1/bin/R

Set_default_check_executable CDHIT_PATH /usr/appli/freeware/cd-hit/4.6.8/cd-hit
Set_default_check_executable CDHIT_EST_PATH /usr/appli/freeware/cd-hit/4.6.8/cd-hit-est

Set_default_check_executable BLASTP_PATH /usr/appli/freeware/blast/2.5.0+/bin/blastp
Set_default_check_executable BLASTDBCMD_PATH /usr/appli/freeware/blast/2.5.0+/bin/blastdbcmd
Set_default_check_executable MAKEBLASTDB_PATH /usr/appli/freeware/blast/2.5.0+/bin/makeblastdb
Set_default_check_executable DIAMOND_PATH /usr/appli/freeware/diamond/0.9.9/diamond

Set_default_check_executable CHIMERIC_PATH /usr/appli/freeware/qiime/1.9.1/Anaconda2/bin/identify_chimeric_seqs.py

Logger "--Finish building software path"
}

#build function header
Header(){
Logger "--Start ${FUNCNAME[1]}"

STEP_NUMBER=1
for STEP_NAME in ${STEP_ARR};do
    if [ ${STEP_NAME} == ${FUNCNAME[1]} ];then
        break
    fi
    let "STEP_NUMBER += 1"
done

OUTPUT_DIR=${OUTPUT_DIR_ARR[${STEP_NUMBER}]}
let "STEP_NUMBER -= 1"
if [ ${STEP_NUMBER} -gt 0 ];then
    INPUT_DIR=${OUTPUT_DIR_ARR[${STEP_NUMBER}]}
else
    INPUT_DIR=${SOURCES_DIR}
fi

#make output directory
mkdir -p ${OUTPUT_DIR}


Logger "---IN:${INPUT_DIR}"
Logger "---OUT:${OUTPUT_DIR}"

}

#build function footer

Footer(){
Logger "--Finish ${FUNCNAME[1]}"
TIME=$(date +%s)
Logger "--${FUNCNAME[1]} cost `expr \( $TIME - $TIME_P \) / 60`min`expr \( $TIME - $TIME_P \) % 60`s"
TIME_P=${TIME}
}


A5G40(){
Header

BARCODE_ARR=$(ls ${INPUT_DIR}*${READ_PATTERN[0]}|rev |cut -d/ -f1|rev|cut -d_ -f1)

(for BARCODE in ${BARCODE_ARR};do
    R1_INPUT_PATH=${INPUT_DIR}${BARCODE}${READ_PATTERN[0]}
    R2_INPUT_PATH=${INPUT_DIR}${BARCODE}${READ_PATTERN[1]}
    R1_OUTPUT_PATH=${OUTPUT_DIR}${BARCODE}${READ_PATTERN[0]}
    R2_OUTPUT_PATH=${OUTPUT_DIR}${BARCODE}${READ_PATTERN[1]}
    R1_NOPAIR_OUTPUT_PATH=${OUTPUT_DIR}${BARCODE}_np${READ_PATTERN[0]}
    R2_NOPAIR_OUTPUT_PATH=${OUTPUT_DIR}${BARCODE}_np${READ_PATTERN[1]}    
    LOG_PATH=${OUTPUT_DIR}${BARCODE}.log    
    echo "${TRIMMOMATIC_PATH} PE -threads 8 -phred33 -trimlog ${LOG_PATH} \
    ${R1_INPUT_PATH} ${R2_INPUT_PATH} ${R1_OUTPUT_PATH} ${R1_NOPAIR_OUTPUT_PATH} \
    ${R2_OUTPUT_PATH} ${R2_NOPAIR_OUTPUT_PATH} AVGQUAL:5 MINLEN:40 &> /dev/null "
done) > TRIMMOMATIC.com
qsubarraypbs -sync -q sp -l select=1:ncpus=8:mem=80gb TRIMMOMATIC.com

Footer
}

CUTADAPT_G40(){
Header

BARCODE_ARR=$(ls ${INPUT_DIR}*${READ_PATTERN[0]}|grep -v _np|rev |cut -d/ -f1|rev|cut -d_ -f1)

(for BARCODE in $BARCODE_ARR;do
    R1_INPUT_PATH=${INPUT_DIR}${BARCODE}${READ_PATTERN[0]}
    R2_INPUT_PATH=${INPUT_DIR}${BARCODE}${READ_PATTERN[1]}
    R1_OUTPUT_PATH=${OUTPUT_DIR}${BARCODE}${READ_PATTERN[0]}
    R2_OUTPUT_PATH=${OUTPUT_DIR}${BARCODE}${READ_PATTERN[1]}
    LOG_PATH=${OUTPUT_DIR}${BARCODE}.log    
    PRIMER_CMD=$(${PYTHON_PATH} ${SCR_DECODE_PRIMER} ${REF_PRIMER} ${REF_MIXTURE} ${BARCODE})
    echo "${CUTADAPT_PATH} --minimum-length 40 --discard-untrimmed -e 0.1 ${PRIMER_CMD} \
    -o ${R1_OUTPUT_PATH} -p ${R2_OUTPUT_PATH} ${R1_INPUT_PATH} ${R2_INPUT_PATH} &> ${LOG_PATH}"    
done) > CUTADAPT.com
qsubarraypbs -sync -q sp -l select=1:ncpus=1:mem=10gb CUTADAPT.com


#delete empty files
find ${OUTPUT_DIR} -name "*" -type f -size 0c | xargs -n 1 rm &> /dev/null

Footer
#Explanation:
   #BARCODE - merged primer set
   #BARCODE_ARR - a list (array) with all the primers
   #PSS - single primer was used

   #DIR_ARR[3]=3_Cutadapted/
   #PP = primer pair, the loop generates PP01 - PPn and the sequence  
   #PL = primer left, acceses $MAIN_DIR"References/" and obtains the primer pair information
   #PR_RC = primer right RC-> i reverse code?
#CUTADAPT
   #first the "left" primer is cut, hence the first file is called _L
   #then the second primer is cut
   # the input file is just put into the command line without anything written before it
   # -o for output; there are two output files: removal of 3' and 5' primers.
   # -a for the adapter that will be removed 
   # -g ^  for the anchored 5' adapter
   # -m minimum length
   # --discard-untrimmed reads without adapter are 'thrown away'
   # &> (bash command) redirects output to a file
   # && (bash command) if first CUTADAPt call was succesfull also do the second (check: || , && , ;) 
#Generated Files
   #the command generates 4 files
   # ***fastq is the output
   # ***log contains the the unused sequences, they had no primer
}


#this function uses FLASH
MERGE(){
Header

MERGE_ERROR_RATE=0.1
MERGE_OVERLAP_MIN=100

MERGE_I_DIR=${OUTPUT_DIR}Merge_I/
MERGE_O_DIR=${OUTPUT_DIR}Merge_O/
MERGE_RESCUE_DIR=${OUTPUT_DIR}Merge_Rescue/
CHIMERA_DIR=${OUTPUT_DIR}Chimera/

mkdir -p ${MERGE_I_DIR} ${MERGE_O_DIR} ${MERGE_RESCUE_DIR} ${CHIMERA_DIR}

#merge i
BARCODE_ARR=$(ls ${INPUT_DIR}*${READ_PATTERN[0]}|grep -v _np|rev |cut -d/ -f1|rev|cut -d_ -f1)
(for BARCODE in $BARCODE_ARR;do
    
    R1_INPUT_PATH=${INPUT_DIR}${BARCODE}${READ_PATTERN[0]}
    R2_INPUT_PATH=${INPUT_DIR}${BARCODE}${READ_PATTERN[1]}
    
    mkdir -p ${MERGE_I_DIR}${BARCODE}
    LOG_PATH=${MERGE_I_DIR}${BARCODE}/log
    echo "$FLASH_PATH -t 8 -m ${MERGE_OVERLAP_MIN} -M 300 -x ${MERGE_ERROR_RATE} \
    ${R1_INPUT_PATH} ${R2_INPUT_PATH} -d ${MERGE_I_DIR}${BARCODE}/ &> ${LOG_PATH}"
done) > FLASH.com
qsubarraypbs -sync -q sp -l select=1:ncpus=8:mem=80gb FLASH.com

#Explanation: flash
   #takes two input files and merges them.
   #path to software
   # -t number of working threads
   # -m min overlap
   # -M max overlap
   # -x max mismatch density, e.g. error rate
   # -O allows outies (see manual for innies and outies)
   # now the files are passed
   # -d directory for outputfiles

#Merge O
BARCODE_ARR=$(ls -p ${MERGE_I_DIR}|rev |cut -d/ -f2|rev)
(for BARCODE in $BARCODE_ARR;do
    R1_INPUT_PATH=${MERGE_I_DIR}${BARCODE}/out.notCombined_1.fastq
    R2_INPUT_PATH=${MERGE_I_DIR}${BARCODE}/out.notCombined_2.fastq
    
    mkdir -p ${MERGE_O_DIR}${BARCODE}
    LOG_PATH=${MERGE_O_DIR}${BARCODE}/log
    echo "${FLASH_PATH} -t 8 -m ${MERGE_OVERLAP_MIN} -M 300 -x ${MERGE_ERROR_RATE} \
    -O ${R1_INPUT_PATH} ${R2_INPUT_PATH} -d ${MERGE_O_DIR}${BARCODE}/ &> ${LOG_PATH} \
    && ${PYTHON_PATH} ${SCR_TRIM_MERGE_O} ${MERGE_O_DIR}${BARCODE}/out.extendedFrags.fastq \
    ${R1_INPUT_PATH} ${R2_INPUT_PATH} ${MERGE_O_DIR}${BARCODE}.fastq"
        
done) > FLASH2.com
qsubarraypbs -sync -q sp -l select=1:ncpus=8:mem=80gb FLASH2.com

#Merge Rescue
BARCODE_ARR=$(ls -p ${MERGE_O_DIR}|rev |grep ^/|cut -d/ -f2|rev)

(for BARCODE in $BARCODE_ARR;do
    echo "sh ${SCR_MERGE_RESURE} ${TRIMMOMATIC_PATH} ${FLASH_PATH} ${MERGE_O_DIR}${BARCODE} ${MERGE_RESCUE_DIR} ${MERGE_ERROR_RATE} ${MERGE_OVERLAP_MIN}"
    
    # TRIMMOMATIC_PATH="$1 $2 $3"
    # FLASH_PATH=$4
    # INPUT_PATH=$5
    # OUTPUT_DIR=$6
    # MERGE_ERROR_RATE=$7
    # MERGE_OVERLAP_MIN=$8

    
done) > SCR_MERGE_RESURE.com
qsubarraypbs -sync -q sp -l select=1:ncpus=4:mem=40gb SCR_MERGE_RESURE.com



#chimera check, merge results and convert into fna
BARCODE_ARR=$(ls ${INPUT_DIR}*${READ_PATTERN[0]}|grep -v _np|rev |cut -d/ -f1|rev|cut -d_ -f1)
(for BARCODE in $BARCODE_ARR;do
    echo "cat ${MERGE_I_DIR}${BARCODE}/out.extendedFrags.fastq ${MERGE_O_DIR}${BARCODE}.fastq \
    ${MERGE_RESCUE_DIR}${BARCODE}.fastq > ${CHIMERA_DIR}${BARCODE}.fastq \
    && ${PYTHON_PATH} ${SCR_SEQ_CONVERTER} ${CHIMERA_DIR}${BARCODE}.fastq ${CHIMERA_DIR}${BARCODE}.fna \
    && ${MAKEBLASTDB_PATH} -in ${CHIMERA_DIR}${BARCODE}.fna -dbtype nucl -hash_index -parse_seqids \
    && ${CHIMERIC_PATH} -m usearch61 --suppress_usearch61_ref -i ${CHIMERA_DIR}${BARCODE}.fna -o ${CHIMERA_DIR}${BARCODE}/ \
    && ${BLASTDBCMD_PATH} -db ${CHIMERA_DIR}${BARCODE}.fna -entry_batch ${CHIMERA_DIR}${BARCODE}/non_chimeras.txt -out ${OUTPUT_DIR}${BARCODE}.fna"
done) > CAT.com
qsubarraypbs -sync -q sp -l select=1:ncpus=1:mem=10gb CAT.com

Footer
}



DEDUPLICATION(){
Header

BARCODE_ARR=$(ls ${INPUT_DIR}*.fna|rev |cut -d/ -f1|rev|cut -d. -f1)
(for BARCODE in ${BARCODE_ARR};do
    INPUT_FILE=${INPUT_DIR}${BARCODE}.fna
    OUTPUT_FILE=${OUTPUT_DIR}${BARCODE}.fna
    LOG_PATH=${OUTPUT_DIR}${BARCODE}.log
    
    echo "${CDHIT_EST_PATH} -i ${INPUT_FILE} -o ${OUTPUT_FILE} \
	-G 0 -aS 1 -c 1 -n 11 -d 0 -M 8000 -T 8 > ${LOG_PATH}"
done) > CDHIT_EST.com
qsubarraypbs -sync -q sp -l select=1:ncpus=8:mem=80gb CDHIT_EST.com

Footer 

#Explanation:
  # -i input file    -o output file
  # -G 0  use global sequence identity. makes the next command necessary
  # -aS alignment coverage for the shorter sequence must be 100%, hence the function name
  # -c 1 100% identity
  # -n word length
  # -d take the whole description of the fasta file
  # -M how muc memory will be used
  # -T how many threads(cpus) 
  # > save output to log
}


FAA(){
Header

BARCODE_ARR=$(ls ${INPUT_DIR}*.fna|rev |cut -d/ -f1|rev|cut -d. -f1)

(for BARCODE in $BARCODE_ARR;do
	echo "$PYTHON_PATH ${SCRIPTS_DIR}Translate_rename.py ${INPUT_DIR}$BARCODE.fna \
	${OUTPUT_DIR}$BARCODE.fna ${OUTPUT_DIR}$BARCODE.faa \
	&& $MAKEBLASTDB_PATH -dbtype prot -in ${OUTPUT_DIR}$BARCODE.faa \
	-parse_seqids -hash_index &> ${OUTPUT_DIR}$BARCODE.faa.log \
    && $MAKEBLASTDB_PATH -dbtype nucl -in ${OUTPUT_DIR}$BARCODE.fna \
	-parse_seqids -hash_index &> ${OUTPUT_DIR}$BARCODE.fna.log"

done) > TRANSLATE_RENAME.com
qsubarraypbs -sync -q sp -l select=1:ncpus=1:mem=10gb TRANSLATE_RENAME.com

Footer

#Explanation: FAA calls a Python3 script in Lis folder. This script has two functions:
  #rename takes two args and is later called to rename an array, creates seqr
  #translate the created seqr in redard of stop codons and frame shifts
 
#First command:	#Lis script takes three args, the input (array 5).fna, the fna output
		#and of course the translated .faa file.
#Sec. command:	#calls makeblastdb. this is explained well in the documentation
		#makes a prot db from .faa file.
		#-parse_seqids -hash_index are needed to call the db from the consol
		# this creates .phd, .phi, .phr ... 8 binary files
#thrd commansd:	#the > stores all the echoed info in a log file.
}


BLASTP(){
Header

# make REF_POLB_HOMOLOGY_SEARCH
if [ "$DIAMOND" = "1" ];then
    if [[ ! -s "${REF_POLB_HOMOLOGY_SEARCH}.dmnd" ]] ;then
        ${DIAMOND_PATH} makedb --in ${REF_POLB_HOMOLOGY_SEARCH} --db ${REF_POLB_HOMOLOGY_SEARCH} > ${REF_POLB_HOMOLOGY_SEARCH}.dmnd.log
    fi
else
    if [[ ! -s "${REF_POLB_HOMOLOGY_SEARCH}.phd" ]] ;then
        ${MAKEBLASTDB_PATH} -dbtype prot -in ${REF_POLB_HOMOLOGY_SEARCH} -parse_seqids -hash_index > ${REF_POLB_HOMOLOGY_SEARCH}.log
    fi
fi


BARCODE_ARR=$(ls ${INPUT_DIR}*.faa|rev |cut -d/ -f1|rev|cut -d. -f1)
(for BARCODE in $BARCODE_ARR;do
    if [ "$DIAMOND" = "1" ];then
        echo -n "${DIAMOND_PATH} blastp -p 4 -e 1e-5 -k 1 -d ${REF_POLB_HOMOLOGY_SEARCH} \
        -q ${INPUT_DIR}${BARCODE}.faa -o ${OUTPUT_DIR}${BARCODE}.res > ${OUTPUT_DIR}${BARCODE}.log"
    else
        echo -n "${BLASTP_PATH} -query ${INPUT_DIR}${BARCODE}.faa -db $REF_POLB_HOMOLOGY_SEARCH \
        -out ${OUTPUT_DIR}${BARCODE}.res -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 2"
    fi
	echo "&& cat ${OUTPUT_DIR}$BARCODE.res|grep MEGA|cut -f1|sort -u > ${OUTPUT_DIR}$BARCODE.tit \
    && test -s ${OUTPUT_DIR}$BARCODE.tit \
    && $BLASTDBCMD_PATH -db ${INPUT_DIR}$BARCODE.faa \
    -entry_batch ${OUTPUT_DIR}$BARCODE.tit -out ${OUTPUT_DIR}$BARCODE.faa \
	&& cat ${OUTPUT_DIR}$BARCODE.tit |cut -dF -f1|sort -u > ${OUTPUT_DIR}${BARCODE}_fna.tit \
	&& $BLASTDBCMD_PATH -db ${INPUT_DIR}$BARCODE.fna \
	-entry_batch ${OUTPUT_DIR}${BARCODE}_fna.tit -out ${OUTPUT_DIR}$BARCODE.fna"
done) > DIAMOND_BLASTP.com
qsubarraypbs -sync -q sp -l select=1:ncpus=4:mem=40gb DIAMOND_BLASTP.com

Footer

#diamond blastp --more-sensitive -p 8 -e 1e-5 -k 1 -d ../../References/ref_for_blastp.diamond.dmnd -q S0-PP1.faa -o S0-PP1.faa.m8

#Explanation per line
  #calls blastp
  #-query "input file"
  #-db database for search. In this case polymerase B from different organisms
  #-out output file
  #-evalue limit for evalue, -outfmt 6 = tabular
  #-max_target_seq = maximum number of seq to keep default 500, here only 1 is kept
  # && links commands and executes them when the first was successful
  #only writes the first entry of every line in a new tit file
  #&&
  #calls blast again
  # -db is the created aa sequence
  # -entry_batch is a list of files with the titel only
  # -> blast makes a file which has title and the sequence 
}


ALIGNMENT(){
Header

BARCODE_ARR=$(ls ${INPUT_DIR}*.faa|rev |cut -d/ -f1|rev|cut -d. -f1)

(for BARCODE in $BARCODE_ARR;do
    echo "${MAFFT_PATH} --quiet --thread 20 --6merpair --addfragments ${INPUT_DIR}$BARCODE.faa \
    ${REF_PPLACER_ALN} > ${OUTPUT_DIR}${BARCODE}.combo.fasta \
    && ${PPLACER_PATH} -j 20 --verbosity 0 -o ${OUTPUT_DIR}${BARCODE}.combo.jplace \
    -t ${REF_PPLACER_RES} -s ${REF_PPLACER_INFO} ${OUTPUT_DIR}${BARCODE}.combo.fasta > /dev/null \
    && ${PYTHON_PATH} ${SCR_PPLACER_DECODE} ${OUTPUT_DIR}${BARCODE}.combo.jplace > ${OUTPUT_DIR}${BARCODE}_Pplacer.tit \
    && cut -dF -f1 ${OUTPUT_DIR}${BARCODE}_Pplacer.tit > ${OUTPUT_DIR}${BARCODE}_Pplacer_fna.tit \
    && ${MAKEBLASTDB_PATH} -dbtype prot -in ${INPUT_DIR}$BARCODE.faa -parse_seqids -hash_index > ${INPUT_DIR}$BARCODE.faa.log \
    && ${MAKEBLASTDB_PATH} -dbtype nucl -in ${INPUT_DIR}$BARCODE.fna -parse_seqids -hash_index > ${INPUT_DIR}$BARCODE.fna.log \
    && ${BLASTDBCMD_PATH} -db ${INPUT_DIR}$BARCODE.faa \
    -entry_batch ${OUTPUT_DIR}${BARCODE}_Pplacer.tit -out ${OUTPUT_DIR}${BARCODE}_Pplacer.faa \
    && ${BLASTDBCMD_PATH} -db ${INPUT_DIR}$BARCODE.fna \
    -entry_batch ${OUTPUT_DIR}${BARCODE}_Pplacer_fna.tit -out ${OUTPUT_DIR}${BARCODE}_Pplacer.fna \
    && ${MAFFT_PATH} --quiet --thread 20 --6merpair --addfragments ${OUTPUT_DIR}${BARCODE}_Pplacer.faa \
    ${REF_DESIGN_PROT_ALN} > ${OUTPUT_DIR}$BARCODE.aln \
    && ${PYTHON_PATH} ${SCR_TRIM_COMMON_REGION} ${REF_DESIGN_PRIMER_NUCL_ALN} ${OUTPUT_DIR}$BARCODE.aln \
    ${OUTPUT_DIR}${BARCODE}_Pplacer.fna ${OUTPUT_DIR}${BARCODE}_Trimmed.faa ${OUTPUT_DIR}${BARCODE}_Trimmed.fna"
done) > MAFFT.com
qsubarraypbs -sync -q sp -l select=1:ncpus=20:mem=200gb MAFFT.com

Footer
}

A5G40
CUTADAPT_G40
MERGE
DEDUPLICATION
FAA
BLASTP
ALIGNMENT


Logger "-Finish Pipeline"