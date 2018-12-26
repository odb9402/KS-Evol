import os
import pandas as pd

inputDir='E:/evolutinoarypattern/RMB2_B10.snpCalling'

fileNames = os.listdir(inputDir)

baseCount_head = ['CHROM', 'POS', 'REF', 'ALT', 'R_CNT', 'A_CNT']
timeList = []

for i in range(0,len(fileNames)):
	name = inputDir + '/' + fileNames[i]
	tmp = fileNames[i].split('-')[2]
	tmp = tmp.split('_')[0]
	
	time = pd.read_csv(name, sep='\t', header=None, names=baseCount_head)
	time['TIME'] = tmp
	time = time[['CHROM', 'POS', 'TIME', 'REF', 'ALT', 'R_CNT', 'A_CNT']]
	
	timeList.append(time)
	

refAlt = pd.concat(timeList, ignore_index=True)
refAlt = refAlt.sort_values(by = ['CHROM', 'POS', 'TIME'])

outName = 'E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.merge.csv'
refAlt.to_csv(outName, sep='\t', header=True, index=None)

