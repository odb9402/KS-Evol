import numpy.random as rd
import sys
import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-t", "--timepoints", help="Number of timepoint")
arg_parser.add_argument("-i", "--input", help="An input file name")
arg_parser.add_argument("-o", "--output", help="An output file name")

args = arg_parser.parse_args()

timePoint = int(args.timepoints)

inFileName = args.input 
inFile = open(inFileName, 'r')

outFileName = args.output
outFile = open(outFileName, 'w')

while(True):
    ref, alt = [], []

    line = inFile.readline()
    if( line == ''):    break
    line = line.strip().split('\t')

    ref.append(int(line[5]))
    alt.append(int(line[6]))

    i = 0
    while (i < timePoint-1):
        line = inFile.readline()
        if(line==''):   break
        line = line.strip().split('\t')

        ref.append(int(line[5]))
        alt.append(int(line[6]))
        i = i +1


    p_0 = sum(alt) / ( sum(ref) + sum(alt) )
    p_mt = []

    for a_mt, r_mt in zip(alt,ref):
        mt = a_mt + r_mt
        if( mt == 0):
            p_mt.append(0)
        else:
            p_mt.append( (a_mt / (a_mt + r_mt)) )

    ks = [abs(p_0 - p) for p in p_mt ]

    maxKS = max(ks)
    val = '%.4f'%maxKS
    outFile.write(str(val) + '\n')

outFile.close()
inFile.close()
