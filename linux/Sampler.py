#!/usr/bin/python

"""
Sampler

Sampling code was originally written by B.H

Class implementation was carried by S.H.
"""

import sys
import random
import itertools
import os
import gzip
import HTSeq

class Sampler(object):
	def __init__(self, fq1, fq2,sample_ratio=0.1):
		P = "/".join(os.path.abspath(fq1).split("/")[:-1])
		#print P
		self.fq1 = fq1

		self.fq2 = fq2

		#print self.fq2
		self.o1 = P + "/" + fq1.split("/")[-1]+ ".dwn.fastq"
		self.o2 = P + "/" + fq2.split("/")[-1]+ ".dwn.fastq"
		#print self.o1
		self.rate = sample_ratio

	def sample(self):
		fraction = float(self.rate)
		in1 = iter(HTSeq.FastqReader(self.fq1))
		in2 = iter(HTSeq.FastqReader(self.fq2))
		o1 = open(self.o1, "w")
		o2 = open(self.o2, "w")

		for read1 , read2 in itertools.izip(in1, in2):
			if random.random() < fraction :
				read1.write_to_fastq_file(o1)
				read2.write_to_fastq_file(o2)

		o1.close()
		o2.close()
