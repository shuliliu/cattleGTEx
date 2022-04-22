#!/bin/bash


###########################################################################################
## Step 1. covert bam to junction files
###########################################################################################
## custom time command
timecmd=/usr/bin/time
## General setup, error behaviour, paths, utilities and script inputs.
source ~/common_general_setup.sh
source ~/common_logging_setup.sh

logfile=$srs_name.Splicing.log
DIR_log=$DIR_work/log

DIR_bam=$DIR_work/SRA/${srs_name}
DIR_splicing=$DIR_work/Splicing
DIR_splicing_ds=$DIR_work/Splicing_ds

mkdir -p $DIR_splicing
mkdir -p $DIR_splicing_ds

[ -z $srs_name ] && srs_name=$1
##########################################
log_newstage "Splicing 1: Splicing."
export PATH=$PATH:$LeafCutterDIR/scripts
$timecmd sh bam2junc.sh $DIR_bam/${srs_name}-STARAligned.sortedByCoord.out.bam $DIR_splicing/${srs_name}_leafCutter.junc \
        && log_info "bam2junc.sh ok" || log_error "bam2junc.sh failed."
[ -e $DIR_splicing/${srs_name}_leafCutter.junc ] && log_success "Splicing 1: Splicing ok"  # echo "Splicing 1: Splicing ok" >> $logfile   

###########################################################################################
###########################################################################################



###########################################################################################
## Step 2. Intron clustering
###########################################################################################
## custom time command
timecmd=/usr/bin/time
## General setup, error behaviour, paths, utilities and script inputs.

DIR_splicing=$DIR_work/Splicing
mkdir -p $DIR_splicing

# OUT_clustering=Buffalo_RNASeq_tissues
[ -z $OUT_clustering ] && OUT_clustering=$1

cd $DIR_splicing
ls | grep "junc$" > juncfiles.txt
$timecmd python $LeafCutterDIR/clustering/leafcutter_cluster.py -j $DIR_splicing/juncfiles.txt -m 50 -o $OUT_clustering -l 500000
###########################################################################################
###########################################################################################

