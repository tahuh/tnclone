#!/usr/bin/python

import sys
import os
import argparse

from Bio import SeqIO

def absp(path) :
	return os.path.dirname(os.path.abspath(path)) + "/" + path.split("/")[-1]

usags = "python buildConfigFile.py [option] <outfile>"

parser = argparse.ArgumentParser(description="Make config file for de novo assembly")

parser.add_argument('-k', "--ksize" , type=int, help="Kmer size")
parser.add_argument('-b' ,"--bed" , type=str , help="Bed file")
parser.add_argument('-r' , "--ref" , type=str , help="Reference")
parser.add_argument('-c',"--cov", type=int,help = "Kmer coverage frequency")
parser.add_argument('-n', "--samples",type=int,help="Number of samples")
parser.add_argument('-o' , "--out",type=str,help="configuration file")

args = parser.parse_args()

ksize = args.ksize
bedname = args.bed
ref = args.ref
cov = args.cov
N = args.samples
out = args.out
outfile = open(absp(out), "aw")
# print "Output to %s"%(absp(out))
# print "Parsing reference"
REFS = {}
for record in SeqIO.parse(open(ref) , "fasta") :
	refid = record.id
	seq =str(record.seq)
	REFS[refid] = seq
	
BEDS = {}
for line in open(bedname) :
	data=  line.rstrip().split("\t")
	#print line
	I = data[0]
	st = int(data[1])
	ed = int(data[2])
	
	try:
		BEDS[I].append((st,ed))
	except KeyError :
		BEDS[I] = [(st,ed)]

# sys.stderr.write("[ConfigBuilder] Build start...\n")
for n in range(1,N+1) :
#	outfile.write(refid + "_" + str(n) + "\t")
	for k in REFS.keys() :
		refid = k
		seq = REFS[refid]
		try:
			region = BEDS[refid]
			for r in region :
				st = r[0]
				ed = r[1]
				outfile.write(refid + "_" + str(n) + "\t")
				outfile.write(str(ksize) + "\t")
				outfile.write(seq[st:st+ksize] + "\t")
				outfile.write(seq[ed - ksize:ed] + "\t")
				outfile.write(str(cov) + "\n")
		except KeyError :
			continue

outfile.close()
# sys.stderr.write("[ConfigBuilder] Done building...\n")
	
