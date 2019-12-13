#####01.Sequence allignment###############################
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
rm *.bam
mv ./tmp/*_sorted.bam .
rm -r tmp

#####02.SNP calling from alignment #########################
# Freebayes to call varient for every bam file in the
# current directory.
#
# Vcftools to filter out indels and take only SNPs.
#
# Parameter : [1](Reference fasta file)
#
############################################################

echo "Reference file $1."

for bam in *.bam
do
    echo "freebayes for $bam"
    prefix="${bam%.bam}"
    freebayes -f $1 $bam >"${prefix}.vcf"
done

for vcf in *.vcf
do
    echo "Filtering indels for $vcf"
    prefix="${vcf%.vcf}"
    vcftools --vcf $vcf --remove-indels --recode --recode-INFO-all --out "${prefix}"
done

mkdir tmp
mv *.recode.vcf tmp
rm *.vcf

mv ./tmp/*.vcf .
rm -r tmp

vcf-merge *.vcf > merged.vcf


#####02.SNP calling from alignment #########################
# Make allele depth file using vcf files.
############################################################

for vcf in *.vcf
do
    echo "Make allele depth for $vcf"
    prefix="${vcf%.vcf}"
    python vcf2AD.py $vcf > "${prefix}.depth"
done


