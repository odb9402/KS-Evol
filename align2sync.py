import glob
import argparse
import subprocess as sp

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-i", "--input", help="An input directory name")
arg_parser.add_argument("-r", "--reference", help="A reference fastq name")

args = arg_parser.parse_args()

bam_files = glob.glob("./{}/*.bam".format(args.input)).sort(key= lambda x: int(x.rsplit('.')[0].split('_')[-1]))
bam_strs = ''
for bam in bam_files:
    bam_strs = bam+" "

command = "samtools mpileup -f {} {}".format(args.reference, bam_strs)
sp.call([command], shell=True)

