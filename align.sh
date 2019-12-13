#bwa mem -t 28 $1 $2 > $3
#bwa mem -t 28 LYZE01.fa SRR940358.fastq > LTE_BYB1_H10_825.bam
#bwa mem -t 28 LYZE01.fa SRR940360.fastq > LTE_BYB1_H10_665.sam
#bwa mem -t 28 LYZE01.fa SRR940378.fastq > LTE_BYB1_H10_240.sam
#bwa mem -t 28 LYZE01.fa SRR940424.fastq > LTE_BYB1_H10_505.sam
#bwa mem -t 28 LYZE01.fa SRR940455.fastq > LTE_BYB1_H10_745.sam
#bwa mem -t 28 LYZE01.fa SRR940482.fastq > LTE_BYB1_H10_140.sam
#bwa mem -t 28 LYZE01.fa SRR940520.fastq > LTE_BYB1_H10_415.sam
#bwa mem -t 28 LYZE01.fa SRR940653.fastq > LTE_BYB1_H10_585.sam
#bwa mem -t 28 LYZE01.fa SRR940654.fastq > LTE_BYB1_H10_0.sam
#bwa mem -t 28 LYZE01.fa SRR940672.fastq > LTE_BYB1_H10_335.sam
#bwa mem -t 28 LYZE01.fa SRR940928.fastq > LTE_BYB1_H10_1000.sam
#bwa mem -t 28 LYZE01.fa SRR940958.fastq > LTE_BYB1_H10_910.sam


##########################################################
# 1. Sequence alignment with bwa.
# 2. Convert file formats of sam to bam.
#
# Parameter : [1](Reference fasta file)
#
##########################################################
echo "Reference file $1."

for fastq in *.fastq
do
    echo "Bwa for $fastq"
    prefix="${fastq%.fastq}"
    bwa mem -t 28 $1 $fastq > "${prefix}.sam"
done

for sam in *.sam
do
    echo "Samtools for $sam"
    prefix="${sam%.sam}"
    samtools view -b $sam > "${prefix}.bam"
done

for bam in *.bam
do
    echo "Bamtools sort for $bam"
    prefix="${bam%.bam}"
    bamtools sort -in $bam -out "${prefix}_sorted.bam"
done

mkdir tmp
mv *_sorted.bam tmp
