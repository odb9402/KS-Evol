##########################################################
# Freebayes to call varient for every bam file in the
# current directory.
#
# Vcftools to filter out indels and take only SNPs.
#
##########################################################

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
