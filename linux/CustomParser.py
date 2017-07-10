#!/usr/bin/python

class FastaParser:
	def __init__(self ,fname) :
		self.fname = fname
		self.F = None
		self.id = None
		self.seq = None

	def open(self) :
		self.F = open(self.fname)
	def close(self) :
		self.F.close()

	def read(self) :
		self.id = self.F.readline().rstrip()
		self.seq = self.F.readline().rstrip()

	
