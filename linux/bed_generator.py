#!/usr/bin/python

from argparse import ArgumentParser
import sys
import os
import glob

from Bio import SeqIO

def abspath(path):
	return os.path.abspath(path)
	
	
cur_dir = abspath("./")
parser = ArgumentParser(description="Automatically builds bed file")
parser.add_argument("--ipath", required=True, help="Path to your reference fasta set" , type=str)
parser.add_argument("--opath", required=True, help="Path to your reference fasta set default=%s"%(cur_dir) , type=str)
parser.add_argument("--prefix" , required=True , help= "Prefix to bed file", type=str)
args = parser.parse_args()

ipath = abspath(args.ipath)
opath = abspath(args.opath)
prefix = args.prefix

falist = glob.glob(ipath + "/*.fa")


bed_file_name = opath + "/" + prefix + ".bed"
bed_file = open(bed_file_name , "w")

for fname in falist:
	with open(fname) as reffile:
		for record in SeqIO.parse(reffile,"fasta"):
			seq = str(record.seq)
			id = record.id
			bed_file.write(id + "\t0\t" + str(len(seq)) + "\n")
			
bed_file.close()
