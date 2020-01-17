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
echo -e "$(date)>: Make KS test depths from $1 . . . \n"
if [[ $1 == *.sync ]]
then
    echo -e "$(date)>: $1 is a mpileup data ( .sync ) \n"
    python $src_dir"/../sync2AD.py" $1
    echo -e "$(date)>: Make KS_test Input . . . \n"
    python $src_dir"makeKSInput.py" -i $(dirname $1)

elif [[ $1 == *.vcf ]]
then
    echo -e "$(date)>: $1 is a vcf data ( .vcf ) \n"
    python vcf2AD.py $1
else
    echo -e "$(date)>: $1 is not both of .vcf or .sync files."
    #exit
fi

###### 02.Binomial sampling and get p-values ###############
#
# It generates KS distances between p0 based on binomial
# samples.
#
############################################################
echo -e "$(date)>: Get p-values for each allele using bootstrapping. . .\n"
$src_dir"getPvals" -n $sampleNum -t $timePoint -i ks_input.tsv -o ks_pvals.tsv

#n_proc_raw=$(nproc)
#n_proc=$(expr $n_proc_raw '-' 1)
#if [[ "$n_proc_raw" > 22 ]]
#then
#    n_proc=22
#fi

#n_proc_10=$(expr 10 + $n_proc)

#line=$(wc -l < 0_result.tsv)
#line_tp=$(expr $line '/' $timePoint)

#line_tp_each=$(expr $line_tp '/' $n_proc)
#remain_lines=$(expr $line_tp '%' $n_proc)

#split_lines=$(expr $line_tp_each '*' $timePoint)

#split -l $split_lines -d --numeric-suffixes=10 ks_input.tsv

#seq 10 $n_proc_10 | parallel -k $src_dir"prog2_binSampling" -n $sampleNum -t $timePoint -i x{} -o 1_ks_pvals_{}.tsv

#for ((i=10;i<=$n_proc_10;i++))
#do
#    $src_dir"getPvals" -n $sampleNum -t $timePoint -i x$i -o ks_pvals_$i.tsv &
#done
#wait

#for ((i=10;i<=$n_proc_10;i++))
#do
#    cat ks_pvals_$i.tsv >> ks_pvals.tsv
#    rm ks_pvals_$i.tsv
#    rm x$i
#done
###### 03.Make KS-test input ###############################
#
# It takes the maximum p-values for
#
############################################################
echo -e "$(date)> Step 3:Convert allele depths data to summerize KS distances for all time-points. . . \n"
python $src_dir"summarizeKS.py" -i ks_input.tsv -o ks_input_summarized.tsv

###### 04.Make KS-test input ###############################
#
############################################################
echo -e "$(date)> Step 4: Merge and take the maximum distance from allele depths data for all time-points. . . \n"
python $src_dir"makeKSOutput.py" -t $timePoint -i ks_input_summarized.tsv -p ks_pvals.tsv -o ks_output.tsv

awk '{FS="\t"; OFS="\t"; print $1,$2,$4,$5,$6,$7,$9}' ks_output.tsv | sort -gk7 | awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,NR-1}' > $final_output

awk '{FS="\t"; OFS="\t"; print $1,$2,$7}' $final_output > $final_output".pval"
