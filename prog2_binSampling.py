import numpy.random as rd

SSIZE = 100000


inFIleName = "E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.merge.2.csv"
inFile = open(inFIleName, 'r')

timePoint = 12
outFileName = "E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.binDist.csv"
outFile = open(outFileName, 'w')

while(True):

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
		outFile.write(ksStar+'\t')
	outFile.write('\n')
		
inFile.close()
