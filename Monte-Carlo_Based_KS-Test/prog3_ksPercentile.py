import json
import numpy as np
import sys
import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-h", "--histogram", help="A histogram input file name")
arg_parser.add_argument("-k", "--ksFile", help="A ks-test input file name")
arg_parser.add_argument("-o", "--output", help="An output file name")
arg_parser.add_argument("-n", "--sampleNum", default=50000, help="A sampling size")

args = arg_parser.parse_args()

SAMPLING_SIZE=int(args.sampleNum)

def jsonToDict(json_str):
    json_data = json.loads(json_str)
    return json_data

def ks_loci(ks, ndd_keys):
    loci = 0

    ks = float(ks)
    for item in ndd_keys:
        tmp = float(item)
        if( ks <= tmp):
            loci = item
            break

    return loci

ndFile = args.histogram
ksFile = args.ksFile

ndLines = open(ndFile, 'r', encoding="UTF-8")
ksLines = open(ksFile, 'r')
oFileName = args.output 
oFile = open(oFileName, 'w')

j = 0

for ks in ksLines:
    nd = ndLines.readline()
    if nd == "":
        continue
    ndd = jsonToDict(nd)
    ndd_keys = list(ndd.keys())

    loci = ks_loci(ks, ndd_keys)
    if ( loci == 0 ):
        res = 0
    else:
        start =ndd_keys.index(loci)
        val = 0
        for i in range(start, len(ndd_keys)):
            val += ndd[ ndd_keys[i] ]

        res = val / SAMPLING_SIZE


    oFile.write( str(res) + '\n')

oFile.close()



