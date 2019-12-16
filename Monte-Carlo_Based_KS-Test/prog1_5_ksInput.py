inputFile='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.merge.csv'
outputFile='E:/evolutinoarypattern/RMB2_B10.freebayes_snps/RMB2_B10.merge.2.csv'

TIMEPOINT=12
TIMELIST=[0,140,240,335,415,505,585,665,745,825,910,1000]
file = open(inputFile)


time=0

header = file.readline()

oFile = open(outputFile,'w')

while(True):

    chunk=[]

    for i in TIMELIST:
        line = file.readline()
        if(line == ''):
            exit()
        tmp = line.strip().split('\t')

        time=int(tmp[2])
        if(i==time):
            chunk.append(line)
        else:
            break

    if(len(chunk) == 12):
        for j in chunk:
            oFile.write(j)

oFile.close()


