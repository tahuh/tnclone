#!/usr/bin/python

usage = "python handle_multiple.py [final.fa] <1.f1> <2.fq> <bwa out path>"

import sys
import os
import subprocess
import re

from Bio import SeqIO

if len(sys.argv) < 2 :
	print usage
	exit(-1)

try:
	final = sys.argv[1]
	dest = final + "2"
	os.rename(final,dest)
	f = open(dest)
	if os.path.getsize(dest) == 0 :
		fsize0 = True
	else:
		fsize0 = False
except IOError :
	print "Missing file %s...."%(final)
	exit(-1)

def bwa_index(ref) :
	return "./bwa index %s"%(ref)

def bwa_mem(ref, fq1, fq2) :
	gene = fq1.split("/")[-1].split("_")[0] + fq1.split("/")[-1].split("_")[1]
	opath = sys.argv[-1]
	return "./bwa mem -M -R @RG\t@ID:foo\tSM:bar %s %s %s > %s"%(ref,fq1,fq2,opath + "/" + gene  + ".bwa.sam") , opath + "/" + gene + ".bwa.sam"


if fsize0 == False:
	bwa_index_s = bwa_index(dest)
	bwa_mem_script , osamname = bwa_mem(dest,sys.argv[2],sys.argv[3])

	subprocess.call(bwa_index_s , shell=True)
	subprocess.call(bwa_mem_script, shell=True)

else:
	exit(-1)

sam = open(osamname)

def get_cigar_alpha(c) :
	return re.split("[0-9]+",c)[1:]

d = {}

for line in sam :
	if line.startswith("@") : continue

	data = line.rstrip().split("\t")

	rname = data[2]
	pos = int(data[3])

	cigar = data[5]
	alphas = get_cigr_alpha(cigar)
	
	if rnmae == '*' : continue

	if not d.has_key(rname) :
		d[rname] = 0
	if 'M' in alphas :
		d[rname] += 1


### Record Sequence

SEQ = {}

for record in SeqIO.parse(f, "fasta") :
	seq = str(record.seq)
	id = record.id
	SEQ[id] = seq

### Compare results

tot = 0

for k ,v in d.items() :
	tot += v

for k in d.keys() :
	d[k] = 1.0 * d[k] / tot 

freq = 0.9

keys = d.keys()
for k in keys :
	if d[k] < freq :
		d.pop(k)

o = open(final , "w")

for k in d.keys() :
	o.write(">" + k + "\n" + SEQ[k] + "\n")
o.close()

