#inputFile='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.merge.csv'
#outputFile='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.merge.2.csv'

import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-i", "--input", help="KS-distances input file name")
arg_parser.add_argument("-p", "--pvalueInput", help="P-values file name")
arg_parser.add_argument("-o", "--output", help="An output file name")
arg_parser.add_argument("-t", "--timepoints", help="Number of timepoint")

args = arg_parser.parse_args()

TIMEPOINT=args.timepoints
TIMELIST=["t{}".format(i) for i in range(TIMEPOINT)]
file = open(args.input,'r')
oFile = open(outputFile,'w')
time=0

header = file.readline()

while(True):

    chunk=[]

    for i in TIMELIST:
        line = file.readline()
        if(line == ''):
            exit()
        tmp = line.strip().split('\t')

        time=int(tmp[2])
        if(i==time):
            chunk.append(line)
        else:
            break

    if(len(chunk) == 12):
        for j in chunk:
            oFile.write(j)

oFile.close()


