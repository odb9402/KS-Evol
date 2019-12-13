import sys

CHROM_LIST= {'LYZE01000001.1':'chr1',
             'LYZE01000002.1':'chr2',
             'LYZE01000003.1':'chr3',
             'LYZE01000004.1':'chr4',
             'LYZE01000005.1':'chr5',
             'LYZE01000006.1':'chr6',
             'LYZE01000007.1':'chr7',
             'LYZE01000008.1':'chr8',
             'LYZE01000009.1':'chr9',
             'LYZE01000010.1':'chr10',
             'LYZE01000011.1':'chr11',
             'LYZE01000012.1':'chr12',
             'LYZE01000013.1':'chr13',
             'LYZE01000014.1':'chr14',
             'LYZE01000015.1':'chr15',
             'LYZE01000016.1':'chr16',
             'LYZE01000017.1':'chr?',
             'LYZE01000018.1':'chr?',
             'LYZE01000019.1':'chr?',
             'LYZE01000020.1':'chr?',
             'LYZE01000021.1':'chr?',
             }

input=open(sys.argv[1], 'r')


### REf:: Major , ALT:: minor

for line in input:
        if (line[0]=='#'):
                pass;
        else:
                line = line.strip().split('\t')
                chrom= CHROM_LIST[ line[0] ]
                if chrom == 'chr?':
                    continue
                pos=line[1]
                ref=line[3]
                alt=line[4]
                if(len(ref) > 1 or len(alt) > 1):
                        pass;
                else:
                        info = line[9].split(':')
                        depth = info[1]
                        allele = info[2].split(',')
                        rDepth = allele[0]
                        aDepth = allele[1]
                        print(chrom,'\t',pos,'\t',ref,'\t',alt,'\t',rDepth,'\t',aDepth)



