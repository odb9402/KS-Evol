#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
file_dir=$1
src_dir=$DIR"/Monte-Carlo_Based_KS-Test/"
final_output=$2

###### 01.Make KS-test input ###############################
#
# To generate the column of the output file correctly, the
# file names should be formatted as like "XXXXX_t0.depth".
# Where t0 is the specific time point.
#
############################################################
echo -e "$(date)> Step 1: Make KS_test Input . . . \n"
python $src_dir"prog1_0_ksInput.py" $1

###### 02.Binomial sampling  ###############################
#
# It generates KS distances between p0 based on binomial
# samples.
#
############################################################
echo -e "$(date)>: Step 2: Binomial sampling for bootstrapping to make null hypothesis. . .\n"
python $src_dir"prog2_binSampling.py" 0_result.tsv 1_binSampling.tsv

###### 03.Make histogram of binomial samples  ##############
#
############################################################
echo -e "$(date)> Step 3: Make histograms of binomial samples to calculate p-values of KS distances. . .\n"
python $src_dir"prog2_ksTest.py" 0_result.tsv 2_ksTest.tsv
python $src_dir"prog2_ksStar_Dist.py" 1_binSampling.tsv 2_p0_histogram.json

###### 04.Calculate p-values  ##############################
#
# It calculates p-values based on the location of KS-distances
# of allele frequencies among the histogram of binomial samples
#
############################################################
echo -e "$(date)> Step 4: Make p-values of each KS distances for all time-points. . . \n"
python $src_dir"prog3_ksPercentile.py" 2_p0_histogram.json 2_ksTest.tsv 3_ksPvals.tsv

###### 05.Make KS-test input ###############################
#
# It takes the maximum p-values for
#
############################################################
echo -e "$(date)> Step 5:Convert allele depths data to summerize KS distances for all time-points. . . \n"
python $src_dir"prog4_Result.py" 0_result.tsv 4_ksResult.tsv

###### 06.Make KS-test input ###############################
#
############################################################
echo -e "$(date)> Step 6: Merge and take the maximum distance from allele depths data for all time-points. . . \n"
python $src_dir"prog5_mergeKS.py" 4_ksResult.tsv 3_ksPvals.tsv 5_merged_ks.tsv


awk -F "," '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' 5_merged_ks.tsv | sort -rk9 | awk '{print NR,$0}' > $final_output

