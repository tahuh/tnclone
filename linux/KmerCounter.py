#!/usr/bin/python

from collections import Counter
from Seq import reverse_complement , generate_kmers
class KmerCounter:
	"""
	pobj1 = SequenceParser object1 (must be ready to parse)
	pobj2 = SequenceParser object2 (must be ready to parse)
	k = kmer size
	"""
	def __init__(self, pobj1, pobj2, k):
		self.p1 = pobj1
		self.p2 = pobj2
		self.k = k
		self.kmers = []
	def Count(self):
		for id, seq, qual in self.p1.parse():
			rc = reverse_complement(seq)
			for km1 , km2 in zip(list(generate_kmers(seq, self.k)) , list(generate_kmers(rc,self.k))):
				self.kmers.append(km1)
				self.kmers.append(km2)

		for id, seq, qual in self.p2.parse():
			rc = reverse_complement(seq)
			for km1, km2 in zip(list(generate_kmers(seq, self.k)), list(generate_kmers(rc, self.k))):
				self.kmers.append(km1)
				self.kmers.append(km2)

		self.count = dict(Counter(self.kmers))
		del self.kmers	
	def GetCountMap(self):
		return self.count

	def GetCountOfKmer(self, kmer):
		return self.count[kmer]
