#!/usr/bin/python

"""
Author : Sunghoon Heo
"""
import os
from Bio import SeqIO

class SequencingReadParser:
	def __init__(self, fwdpath, revpath) :
		self.fwdread = []
		self.revread = []
		self.fwdpath = os.path.dirname(os.path.abspath(fwdpath)) + "/" + fwdpath.split("/")[-1]
		self.revpath = os.path.dirname(os.path.abspath(revpath)) + "/" + revpath.split("/")[-1]
	
	def readFwd(self) :
		F = open(self.fwdpath)
		for record in SeqIO.parse(F, "fastq") :
			seq = str(record.seq)
			self.fwdread.append(seq)
		F.close()

	def readRev(self) :
		R = open(self.revpath)
		for record in SeqIO.parse(R, "fastq") :
			seq = str(record.seq)
			self.revread.append(seq)
		R.close()
