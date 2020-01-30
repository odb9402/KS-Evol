#!/bin/bash

function description
{
        echo -e "\t<USAGE of preprocessing for KS-test>"
        echo -e "\t:\n"
        echo -e "\t[allele_depth_pipe.sh -r <reference sequence> -d <working directory> -o <output prefix of sync file>]\n"
        echo -e "\t:[-d] The working directory which includes the raw sequencing data. \n"
        echo -e "\t:Note that the raw sequencing data have to be formatted as 'XXXX_timepoint.fastq'\n"
}


while [ -n "$1" ]
do
    case $1 in

    -h)
        description
        exit 1
        ;;

    -r)
        reference=$2
        shift 2
        ;;
    -d)
        workingDir=$2
        shift 2
        ;;
    -o)
        outputPrefix=$2
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
if [ -z "$reference" ]; then description;
echo -e "The reference sequence -r is missed."; exit 1; fi
if [ -z "$workingDir" ]; then description;
echo -e "The working director -d is missed."; exit 1; fi

#####01.Sequence allignment###############################
# 1. Sequence alignment with bwa.
# 2. Convert file formats of sam to bam.
#
# Parameter : [1](Reference fasta file)
#             [2](Workind directory)
#
##########################################################
echo "Reference file $reference."
echo "Working directory with fastq files $workingDir."

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cutoff_freq="0.05"

for fastq in $workingDir/*.fastq
do
    echo "Bwa for $fastq"
    prefix="${fastq%.fastq}"
    bwa mem -t 28 $reference $fastq > "${prefix}.sam"
done

for sam in $workingDir/*.sam
do
    echo "Samtools for $sam"
    prefix="${sam%.sam}"
    samtools view -b $sam > "${prefix}.bam"
done

for bam in $workingDir/*.bam
do
    echo "Bamtools sort for $bam"
    prefix="${bam%.bam}"
    bamtools sort -in $bam -out "${prefix}_sorted.bam"
done

mkdir $workingDir/tmp
mv $workingDir/*_sorted.bam $workingDir/tmp
rm $workingDir/*.bam
mv $workingDir/./tmp/*_sorted.bam .
rm -r $workingDor/tmp

for bam in $workingDir/*.bam
do
    removed="${bam%_sorted.bam}"
    mv $bam "${removed}.bam"
done

echo "Merge alignment files . . ."
python $DIR/align2sync.py -i $workingDir -r $reference > "{$workingDir}/{$outputPrefix}.sync"

echo "Convert the merged alignment file to .sync format. . ."
java -ea Xmx8g -jar $DIR/mpileup2sync.jar --input "{$workingDir}/{$outputPrefix}.sync" --output "{$workingDir}/{$outputPrefix}.sync" --threads 16

echo "Filtering alleles . . ."
python2 $DIR/filtAF.py -i "{$workingDir}/{$outputPrefix}.sync" -o "{$workingDir}/{$outputPrefix}_filt1.sync"

python2 $DIR/filtZeroAF.py -i "{$workingDir}/{$outputPrefix}_filt1.sync" -o "{$workingDir}/{$outputPrefix}_filt2.sync"

#####02.SNP calling from alignment #########################
# Freebayes to call varient for every bam file in the
# current directory.
#
# Vcftools to filter out indels and take only SNPs.
#
# Parameter : [1](Reference fasta file)
#
############################################################

#echo "Reference file $1."

#for bam in *.bam
#do
#    echo "freebayes for $bam"
#    prefix="${bam%.bam}"
#    freebayes -f $1 $bam >"${prefix}.vcf"
#done

#for vcf in *.vcf
#do
#    echo "Filtering indels for $vcf"
#    prefix="${vcf%.vcf}"
#    vcftools --vcf $vcf --remove-indels --recode --recode-INFO-all --out "${prefix}"
#done

#mkdir tmp
#mv *.recode.vcf tmp
#rm *.vcf

#mv ./tmp/*.vcf .
#rm -r tmp

#vcf-merge *.vcf > merged.vcf


#####02.SNP calling from alignment #########################
# Make allele depth file using vcf files.
############################################################

#for vcf in *.vcf
#do
#    echo "Make allele depth for $vcf"
#    prefix="${vcf%.vcf}"
#    python vcf2AD.py $vcf > "${prefix}.depth"
#done


