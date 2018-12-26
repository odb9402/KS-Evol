CHROM_LIST= {'chrI':'LYZE01000001.1',
             'chrII':'LYZE01000002.1',
             'chrIII':'LYZE01000003.1',
             'chrIV':'LYZE01000004.1',
             'chrV':'LYZE01000005.1',
             'chrVI':'LYZE01000006.1',
             'chrVII':'LYZE01000007.1',
             'chrVIII':'LYZE01000008.1',
             'chrIX':'LYZE01000009.1',
             'chrX':'LYZE01000010.1',
             'chrXI':'LYZE01000011.1',
             'chrXII':'LYZE01000012.1',
             'chrXIII':'LYZE01000013.1',
             'chrXIV':'LYZE01000014.1',
             'chrXV':'LYZE01000015.1',
             'chrXVI':'LYZE01000016.1'}



gffPATH='E:/evolutinoarypattern/w303.gff.sorted/w303.gene.gff'
posPATH='E:/evolutinoarypattern/w303.gff.sorted/merging_total_pos.txt'
outPATH='E:/evolutinoarypattern/w303.gff.sorted/BYS2_E03_total.gene.txt'

mine = open(posPATH, 'r')
gff = open(gffPATH, 'r')
out = open(outPATH, 'w')
geneLine = gff.readlines()

for line in mine:
	
	tmp = line.strip().split('\t')
	pos = int(tmp[1])
	chrom = CHROM_LIST[tmp[0].strip()]

	
	find =0
	for i in range(0, len(geneLine)):
		gene = geneLine[i]
		tmp =	gene.strip().split(' ')
		
		gChrom = tmp[0]
		
		if(gChrom == chrom ):
			start = int(tmp[3])
			end = int(tmp[4])

			if( pos >= start and pos <= end):
				if(find > 0):
					res = '\t' + tmp[8] 
				else:
					res = chrom + '\t' + str(pos) + '\t' + tmp[8]
				find += 1
				out.write(res)
		
	if(find == 0):
		res = chrom + '\t' + str(pos) + '\n'
		out.write(res)
	else:
		out.write('\n')
out.close()


