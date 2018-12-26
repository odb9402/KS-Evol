input='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.res.2.csv'
pInput='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.pVal.csv'
output='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.res.ks.csv'

timePoint = 12

iFile = open(input,'r')
pFile = open(pInput, 'r')
oFile = open(output, 'w')

pVals = pFile.readlines()


oFile.write('CHROM,POS,TIMEPOINT,REF,ALT,TOTAL,ALT_ALL,P_ALL,pValue,P0,P140,P240,P335,P415,P505,P585,P665,P745,P825,P910,P1000')
oFile.write('\n')
p=0
iFile.readline()
while (True):
	line = []
	
	time = iFile.readline()
	tmp = time.strip().split(',')
	
	if(time == ''):	break
	
	line = tmp[0:5]
	line.append(tmp[8])
	line.append(tmp[9])
	line.append(tmp[11])
	line.append(pVals[p].rstrip())
	line.append(tmp[10])

	p = p + 1
	
	for i in range(0, timePoint-1):
		time = iFile.readline()
		if(time == ''):	break
		tmp = time.strip().split(',')

		line.append(tmp[10])
		
	oFile.write(','.join(line))
	oFile.write('\n')


