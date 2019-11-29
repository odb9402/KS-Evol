import pandas as pd
import numpy as np
import sys

#input = "E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.merge.2.csv"
#output = "E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.res.2.csv"
input = sys.argv[1] 
output = sys.argv[2] 
data = pd.read_csv( input,
					sep='\t',
					names=['CHROM', 'POS', 'TIMEPOINT', 'REF', 'ALT', 'R_CNT', 'A_CNT'])
					
data['T_CNT'] = data['R_CNT'] + data['A_CNT']
data = data.join(data.groupby('POS')['T_CNT'].sum(), on='POS', rsuffix='_ALL')
data = data.join(data.groupby('POS')['A_CNT'].sum(), on='POS', rsuffix='_ALL')
data['PMT'] = data['A_CNT'] / data['T_CNT']
data['P0'] = data['A_CNT_ALL'] / data['T_CNT_ALL']
	
data.to_csv(output, sep=',', index=False)




