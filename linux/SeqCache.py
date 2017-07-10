#!/usr/bin/python

"""
SeqCache.py

This not really a 'cache' but stores some assembly sequence information
for TnClone's analysis


Author : Sunghoon Heo
"""

class SeqRecorder(object):
	def __init__(self,id=None,seq=None,qual=None):
		self.seq = seq
		self.id = id
		self.qual = qual

	def get_seq(self):
		return self.seq
	
	def get_id(self):
		return self.id
	
	def get_qual(self):
		return self.qual

	def set_seq(self, seq):
		self.seq = seq

	def set_qual(self, qual):
		self.qual = qual

	def set_id(self, id):
		self.id = id
 
class SeqCache(object):
	def __init__(self, seq_record):
		if type(seq_record) == list or type(seq_record) == tuple:
			self.record = seq_record
			self.is_list = True
		else:
			self.record = seq_record
			self.is_list = False

	def seq(self):
		if self.is_list:
			ret = map(lambda x : x.get_seq() , self.record)
		else:
			ret = self.record.get_seq()

		return ret

	def id(self):
		if self.is_list:
			ret = map(lambda x : x.get_id(), self.record)
		else:
			ret = self.record.get_id()

		return ret

	def get_qual(self):
		if self.is_list:
			ret = map(lambda x : x.get_qual(), self.record)
		else:
			ret = self.record.get_qual()

		return ret
