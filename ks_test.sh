#!/bin/bash

function description
{
        echo -e "\t<USAGE of KS_test>"
        echo -e "\t:\n"
        echo -e "\t[ks_test.sh -n <number of samplings> -t <timepoints> <input> <output>]\n"
        echo -e "\t:[-n] The number of sampling (Default : 50000). \n"
        echo -e "\t:[-t] Timepoints of input depth file.\n"
}


while [ -n "$1" ]
do
    case $1 in

    -h)
        description
        exit 1
        ;;

    -n)
        sampleNum=$2
        shift 2
        ;;
    -t)
        timePoint=$2
        shift 2
        ;;

    -*)
        echo "you used wrong options $1"
        description
        exit 1
        ;;
     *)
        break
        ;;
esac
done

if [ $? != 0 ]; then description; fi
if [ -z "$sampleNum" ]; then sampleNum=50000; fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
file_dir=$1
src_dir=$DIR"/Monte-Carlo_Based_KS-Test/"
final_output=$2

if [ -z "$1" ]; then description;
echo -e "The input file name is missed."; exit 1; fi
if [ -z "$2" ]; then description;
echo -e "The output file name is missed."; exit 1; fi
if [ -z "$timePoint" ]; then description;
echo -e "The number of time points is missed."; exit 1; fi

echo -e " The number of sampling for each allele is : ${sampleNum}"
echo -e " The number of timepoints is : ${timePoint}"
echo -e " The input directory name is : $1"
echo -e " The output file name is : $2"

###### 01.Make KS-test input ###############################
#
# To generate the column of the output file correctly, the
# file names should be formatted as like "XXXXX_t0.depth".
# Where t0 is the specific time point.
#
############################################################
#echo -e "$(date)> Step 1: Make KS_test Input . . . \n"
#python $src_dir"prog1_0_ksInput.py" -i $1

###### 02.Binomial sampling  ###############################
#
# It generates KS distances between p0 based on binomial
# samples.
#
############################################################
echo -e "$(date)>: Step 2: Binomial sampling for bootstrapping to make null hypothesis. . .\n"

n_proc=22
n_proc_10=$(expr 10 + $n_proc)

line=$(wc -l < 0_result.tsv)
line_tp=$(expr $line '/' $timePoint)

line_tp_each=$(expr $line_tp '/' $n_proc)
remain_lines=$(expr $line_tp '%' $n_proc)

split_lines=$(expr $line_tp_each '*' $timePoint)

split -l $split_lines -d --numeric-suffixes=10 0_result.tsv

seq 10 $n_proc_10 | parallel $src_dir"prog2_binSampling" -n $sampleNum -t $timePoint -i x{} -o 1_binSampling_{}.tsv

cat 1_binSampling_*.tsv > 1_binSampling.tsv
#SET=$(seq 10 $n_proc_10)
#for i in $SET
#do
#    $src_dir"prog2_binSampling" -n $sampleNum -t $timePoint -i x$i -o 1_binSampling_$i.tsv &
#done
#wait

###### 03.Make histogram of binomial samples  ##############
#
############################################################
echo -e "$(date)> Step 3: Make histograms of binomial samples to calculate p-values of KS distances. . .\n"
python $src_dir"prog2_ksTest.py" -t $timePoint -i 0_result.tsv -o 2_ksTest.tsv
python $src_dir"prog2_ksStar_Dist.py" -i 1_binSampling.tsv -o 2_p0_histogram.json

###### 04.Calculate p-values  ##############################
#
# It calculates p-values based on the location of KS-distances
# of allele frequencies among the histogram of binomial samples
#
############################################################
echo -e "$(date)> Step 4: Make p-values of each KS distances for all time-points. . . \n"
python $src_dir"prog3_ksPercentile.py" -n $sampleNum -hist 2_p0_histogram.json -k 2_ksTest.tsv -o 3_ksPvals.tsv

###### 05.Make KS-test input ###############################
#
# It takes the maximum p-values for
#
############################################################
echo -e "$(date)> Step 5:Convert allele depths data to summerize KS distances for all time-points. . . \n"
python $src_dir"prog4_Result.py" -i 0_result.tsv -o 4_ksResult.tsv

###### 06.Make KS-test input ###############################
#
############################################################
echo -e "$(date)> Step 6: Merge and take the maximum distance from allele depths data for all time-points. . . \n"
python $src_dir"prog5_mergeKS.py" -t $timePoint -i 4_ksResult.tsv -p 3_ksPvals.tsv -o 5_merged_ks.tsv

awk -F "," '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' 5_merged_ks.tsv | sort -gk9 | awk '{print NR,$0}' > $final_output

rm x*
rm 1_binSampling_*.tsv