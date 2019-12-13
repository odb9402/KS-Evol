import sys

output_format=sys.argv[1].rsplit(".")[0] + "_{}.depth"

allele_index={'A':0, 'T':1, 'C':2, 'G':3, 'N':4, 'D':5}

input=open(sys.argv[1], 'r')
header=input.readline()
header=header.strip().split('\t')
counts_all=header[3:]
time_points=len(counts_all)
input.close()

input=open(sys.argv[1], 'r')
output_files = [open(output_format.format("t"+str(x)), 'w') for x in range(time_points)]

### REf:: Major , ALT:: minor
for line in input:
    if (line[0]=='#'):
            pass;
    else:
        line = line.strip().split('\t')
        chrom=line[0]
        if chrom == 'chr?':
            continue
        pos=line[1]
        major_allele=line[2]
        counts_all=line[3:]

        for i in range(time_points):
            counts = list(map(lambda x: int(x), counts_all[i].split(':')))
            major_depth = counts[allele_index[major_allele]]

            max_counts=0
            max_index=-1
            for j in range(len(counts)):
                if j == allele_index[major_allele]:
                    continue
                if counts[j] > max_counts:
                    max_counts=counts[j]
                    max_index=j
            minor_depth = max_counts
            minor_allele = list(allele_index.keys())[max_index]

            output_files[i].write(chrom+"\t"+pos+"\t"+major_allele+"\t"+minor_allele+"\t"+str(major_depth)+"\t"+str(minor_depth)+"\n")

