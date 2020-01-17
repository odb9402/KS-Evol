DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
for vcf in *.vcf
do
    echo "Make allele depth for $vcf"
    prefix="${vcf%.vcf}"
    python $DIR"/vcf2AD.py" $vcf > "${prefix}.depth"
done
