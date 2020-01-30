import argparse
import time

parser = argparse.ArgumentParser(description='This script filter out the SNPs if the maximum allele frequency of the SNPs less than specific threshold.')
parser.add_argument('-i','--input', dest='input', type=str, required=True, help='Input file containing allele counts in a synchronized file.')

parser.add_argument('-o','--output', dest='output', type=str, required=True, help='Output file containing the filtered data in the same format.')
parser.add_argument('-t','--treshold', type=float, default=0.05, help='Allele frequency cut-off for the each variant')
args = parser.parse_args()

infile=open(args.input,'r')
outfile=open(args.output,'w')

start = time.time()

for l in infile:
        a=l.split()
        ref_allele = a[2]
        all_depths = 0.
        
        counts_sum = {'A':0,'T':0,'C':0,'G':0,'N':0,'D':0}
        for item in a[3:]:
            counts = [int(x) for x in item.split(':')]
            for i in range(6):
                counts_sum[list(counts_sum.keys())[i]] += counts[i]
                all_depths += counts[i]
            
        sorted_counts = sorted(counts_sum.items(), key= lambda x : x[1], reverse=True)
        
        if sorted_counts[0][0] == ref_allele:
            freq = sorted_counts[1][1] / all_depths
        else:
            freq = sorted_counts[0][1] / all_depths
        
        #print 'Freq : ', freq, 'all_depths : ', all_depths , counts_sum, ref_allele, sorted_counts
        
        if freq > args.threshold and (1-freq) > args.threshold:
            print >> outfile, "\t".join(a)

infile.close()
outfile.close()
