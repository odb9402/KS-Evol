import argparse
import time

parser = argparse.ArgumentParser(description='This script filter out the rows from the syncronized file which contain counts with only zeros at a certain time-replicate position (column).')
parser.add_argument('-i','--input', dest='input', type=str, required=True, help='Input file containing allele counts in a synchronized file.')

'''
input format:

3L      8724429 G       63:0:0:39:0:0   279:0:0:179:0:0 147:0:0:85:0:0  84:0:0:107:0:0  ...
3L      8724539 A       21:0:137:0:0:0  59:0:337:0:0:0  38:0:193:1:0:0  65:0:158:0:0:0  ...
3L      8736299 A       129:0:0:124:0:0 225:0:0:221:0:0 94:0:0:110:0:0  0:0:0:0:0:0     ...

Note: The 3rd row will be filtered out since there are only zeros in the last column.

'''

parser.add_argument('-o','--output', dest='output', type=str, required=True, help='Output file containing the filtered data in the same format.')
args = parser.parse_args()

infile=open(args.input,'r')
outfile=open(args.output,'w')

start = time.time()

for l in infile:
        a=l.split()
        l=[]
        counts_sum = [0,0,0,0,0,0]
        for item in a[3:]:
            counts = [int(x) for x in item.split(':')]
            for i in range(6):
                counts_sum[i] += counts[i]
                
        num_zero = 0
        for c_sum in counts_sum:
            if c_sum == 0:
                num_zero += 1
        
        if num_zero != 5:
            print >> outfile, "\t".join(a)

infile.close()
outfile.close()
