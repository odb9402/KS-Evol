import sys
#input='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.res.2.csv'
#pInput='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.pVal.csv'
#output='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.res.ks.csv'
import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-i", "--input", help="KS-distances input file name")
arg_parser.add_argument("-p", "--pvalueInput", help="P-values file name")
arg_parser.add_argument("-o", "--output", help="An output file name")
arg_parser.add_argument("-t", "--timepoints", help="Number of timepoint")

args = arg_parser.parse_args()

input= args.input
pInput= args.pvalueInput
output= args.output

timePoint = int(args.timepoints) 

iFile = open(input,'r')
pFile = open(pInput, 'r')
oFile = open(output, 'w')

pVals = pFile.readlines()

oFile.write('chr\tloc\ttimepoint\tmajor\tminor\ttotal_depth\tminor_depth\tp0\tpValue,')
for i in range(timePoint - 1):
    oFile.write(str(i))
    oFile.write('\t')
oFile.write(str(timePoint-1))
oFile.write('\n')
p=0
iFile.readline()

while (True):
    line = []

    time = iFile.readline()
    tmp = time.strip().split('\t')

    if(time == ''): break

    line = tmp[0:5]
    line.append(tmp[8])
    line.append(tmp[9])
    line.append(tmp[11])
    try:
        line.append(pVals[p].rstrip())
    except IndexError:
        continue
    line.append(tmp[10])

    p = p + 1

    for i in range(0, timePoint-1):
        time = iFile.readline()
        if(time == ''): break
        tmp = time.strip().split('\t')

        line.append(tmp[10])

    oFile.write('\t'.join(line))
    oFile.write('\n')


