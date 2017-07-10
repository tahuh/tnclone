#!/usr/bin/python

"""
Author : Sunghoon Heo
"""
from DNAStrManipulate import convert2barr, barr2seq

class KFS :
	def __init__(self,oname,ksize = 63) :
		self.kfs = {}
		self.ksize = ksize
		self.oname = oname
	def computeFreq(self, read) :
		for i in range(len(read) - self.ksize+1) :
			kmer = read[i:i+self.ksize]
			# converted = convert2barr(kmer)
			# if converted == None : continue
			converted = kmer
			try:
				self.kfs[converted] += 1
			except KeyError :
				self.kfs[converted] = 1

	def writeRecord(self) :
		ofile = open(self.oname , "w")
		for k , v in self.kfs.items() :
			ofile.write(k + "\t" + str(v) + "\n")
		ofile.close()
