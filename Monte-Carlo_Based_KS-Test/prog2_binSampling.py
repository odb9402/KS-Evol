import numpy.random as rd
import sys
import progressbar
import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-t", "--timepoints", help="Number of timepoint")
arg_parser.add_argument("-i", "--input", help="An input file name")
arg_parser.add_argument("-o", "--output", help="An output file name")
arg_parser.add_argument("-n", "--sampleNum", default=50000, help="A sampling size")

args = arg_parser.parse_args()

SSIZE = int(args.sampleNum)

inFileName = args.input 
print("Count number of lines. . .")
num_lines = sum(1 for line in open(inFileName))
print("Number of lines of {} : {}".format(inFileName, num_lines))
inFile = open(inFileName, 'r')

timePoint = int(args.timepoints)
outFileName = args.output
outFile = open(outFileName, 'w')

prog = 0
bar = progressbar.ProgressBar(max_value=num_lines)

write_strings = ""

while(True):
    bar.update(prog)
    ref, alt = [], []
    i = 1

    line = inFile.readline()
    if( line == ''):
        break

    line = line.strip().split('\t')
    ref.append(int(line[5]))
    alt.append(int(line[6]))


    while (i < timePoint):
        line = inFile.readline()
        line = line.strip().split('\t')
        ref.append(int(line[5]))
        alt.append(int(line[6]))
        i = i +1

    p_0 = sum(alt) / ( sum(ref) + sum(alt) )
    p_mt = [ ( a_mt / (a_mt + r_mt) ) for a_mt, r_mt in zip(alt, ref) ]

    ks = [abs(p_0 - p_tmp) for p_tmp in p_mt ]
    maxKS = max(ks)

    binSampling = []

    for time in range(timePoint):
        totalN = ref[time] + alt[time]
        sampling = rd.binomial(n=totalN, p=p_0, size=SSIZE)
        toProb = [(x / totalN) for x in sampling]
        binSampling.append(toProb)

    ks_new = []

    for i in range(SSIZE):
        tmp = []
        for j in range(timePoint):
            tmp.append(abs(binSampling[j][i] - p_0 ))

        cut = format(max(tmp), '.4f')
        ks_new.append(cut)

    for ksStar in ks_new:
        write_strings = write_strings + ksStar + '\t'
        #outFile.write(ksStar+'\t')
    #outFile.write('\n')
    write_strings = write_strings + '\n'
    
    if prog % 10000 == 0:
        outFile.write(write_strings)
        write_strings = ""
    prog = prog + 1

inFile.close()
