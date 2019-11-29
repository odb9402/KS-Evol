import os
import pandas as pd
import glob

fileNames = glob.glob("*.depth")

baseCount_head = ['CHROM', 'POS', 'REF', 'ALT', 'R_CNT', 'A_CNT']
timeList = []

for i in range(0,len(fileNames)):
    tmp = fileNames[i].split('_')[3]

    time = pd.read_csv(fileNames[i], sep='\t', header=None, names=baseCount_head)
    time['TIME'] = tmp
    time = time[['CHROM', 'POS', 'TIME', 'REF', 'ALT', 'R_CNT', 'A_CNT']]

    timeList.append(time)


refAlt = pd.concat(timeList, ignore_index=True)
refAlt = refAlt.sort_values(by = ['CHROM', 'POS', 'TIME'])

outName = '../result.csv'
refAlt.to_csv(outName, sep='\t', header=False, index=None)

