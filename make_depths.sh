for vcf in *.vcf
do
    echo "Make allele depth for $vcf"
    prefix="${vcf%.vcf}"
    python vcf2AD.py $vcf > "${prefix}.depth"
done
