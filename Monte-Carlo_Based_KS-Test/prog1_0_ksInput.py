import os
import pandas as pd
import glob
import sys
import argparse

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-i", "--input", help="An input directory name")

args = arg_parser.parse_args()

dir_name = args.input

fileNames = glob.glob(dir_name + "/*.depth")

baseCount_head = ['CHROM', 'POS', 'REF', 'ALT', 'R_CNT', 'A_CNT']
timeList = []

for i in range(0,len(fileNames)):
    print("Make a KS-Test input for {}".format(fileNames[i]))
    #tmp = fileNames[i].split('_')[4]
    tmp =fileNames[i].rsplit('/',1)[1].split('_')[-1].split('.')[0]

    time = pd.read_csv(fileNames[i], sep='\t', header=None, names=baseCount_head)
    time['TIME'] = tmp
    time = time[['CHROM', 'POS', 'TIME', 'REF', 'ALT', 'R_CNT', 'A_CNT']]

    timeList.append(time)

print("{} files total : It might imply the dataset has {} timepoints.".format(
    len(fileNames),len(fileNames)))

refAlt = pd.concat(timeList, ignore_index=True)
refAlt = refAlt.sort_values(by = ['CHROM', 'POS', 'TIME'])

outName = '0_result.tsv'
print(" Merged KS-test input is writing. . .")
refAlt.to_csv(outName, sep='\t', header=None, index=None)
