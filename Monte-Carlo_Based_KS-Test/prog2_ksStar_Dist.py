import json
import os
import collections
import sys
import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-i", "--input", help="An input file name")
arg_parser.add_argument("-o", "--output", help="An output file name")

args = arg_parser.parse_args()

inputDir = args.input 
outputDir = args.output

def list_to_Dict(li):
    dct = {}
    for item in li:
        if (item in dct):
            dct[item] = dct[item] + 1
        else:
            dct[item] = 1


    return  collections.OrderedDict(sorted(dct.items()))

input = open(inputDir, 'r')
output = open(outputDir, 'w')

for line in input:
    line = line.strip().split()

    newLine = []

    for n in line:
        newLine.append(float(n))

    histogram = list_to_Dict(newLine)
    output.write(json.dumps(histogram))
    output.write('\n')

output.close()
