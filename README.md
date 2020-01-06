# Allele Frequency Hypothesis Test Using Monte-Carlo based Kolmogorov-Smirnov test

The genomic varient calling pipeline using the monte-carlo based Kolmogorov-Smirnov(KS) test. 



## Quick start

Since it uses **GSL**(GNU Scientific Library) to optimize bootstrapping processes, it can be executed **only for Linux** environments.

### 0. Install

```shell
./install.sh
```



### 1. For a mpileup (.sync) input

```
./ks_test.sh -n <the number of samples> -t <timepoints> <input> <output>

e.g >
./ks_test.sh -n 100000 -t 9 sample_data.sync sample_output.tsv
```



### 2. For fastq, fasta inputs

If you have several **fastq sequence files for each time point** and their reference file, it gives the pipeline to handle the whole processes including **sequence alignments, variant calling and its statistical test** with KS distance.

To use the pipeline, all the fastq files should be included at the working directory so that the pipeline automatically detects the fastq files.

```
./allele_depth_pipe.sh <reference>
```




## Data description

### 1. Overview

The universal meaning of the allele frequency is just a biological derivation rate. However, for the SNP, it could be the changing rate of each letter from a sequence.

E. g) AAAATCC -> AAATTCA  ( SNP for the 4'th letter and 8'th letter) 

	allele frequency for the sequence -> 2/8 = 0.25.

Q ) How about multiple sequences?

AATCCGG AATCCGA

ATTCGGC ATTCGGC

ATTGCGG GGTGCGG

3/21 ~ 0.142 ??

However, what we want to do is **giving statistical score ( confidences, p-value, ... etc ) that can represent the significance of the allele frequency** based on some statistical measurement as like p-value of Poisson distribution.

After giving scores for each SNP, **we can rank them** with the score and decide the threshold value to determine which SNPs are important or not based on the threshold.

The evaluation process of the statistical test is just comparing true significant SNP with the significant SNP from the test based on our hypothesis.

 

*Observed allele counts over time t1~tn at a bi-allelic marker with alleles major and minor for the m-th SNP*

| Allele( m`th ) |    t1    | ...  | ...  | ...  |    tn    |
| :------------: | :------: | ---- | ---- | ---- | :------: |
|     Major      | $a_{m1}$ |      |      |      | $a_{mn}$ |
|     Minor      | $b_{m1}$ |      |      |      | $b_{mn}$ |
|     Total      | $n_{m1}$ |      |      |      | $n_{mn}$ |

Note that the number of all of these allele is M, which means there are M number of tables like the above.
