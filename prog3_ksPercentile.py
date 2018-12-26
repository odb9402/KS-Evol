import json
import numpy as np

SAMPLING_SIZE=100000

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

ndFile = 'E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.KSHist.csv'
ksFile = 'E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.ks.csv'
	
ndLines = open(ndFile, 'r')
ksLines = open(ksFile, 'r')
oFileName = 'E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.pVal.csv'
oFile = open(oFileName, 'w')

for ks in ksLines:
	nd = ndLines.readline()	
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

		
		
