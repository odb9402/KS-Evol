import json
import os
import collections


inputDir = 'E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.binDist.csv'
outputDir = 'E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.KSHist.csv'

def list_to_Dict(li):
	dct = {}
	for item in li:
		if (item in dct):
			dct[item] = dct[item] + 1
		else:
			dct[item] = 1
	

	return	collections.OrderedDict(sorted(dct.items()))

input = open(inputDir, 'r')
output = open(outputDir, 'w')

for line in input:
	line = line.strip().split()

	newLine = []
	
	for n in line:
		newLine.append(float(n) )
	
	histogram = list_to_Dict(newLine)
	output.write(json.dumps(histogram))
	output.write('\n')
	
output.close()
