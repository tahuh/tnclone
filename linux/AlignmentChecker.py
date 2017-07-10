#!/usr/bin/python

"""
AlignmentChecker.py

Checks INDEL presence in aligned sequence before testing contig's stats

Author : Sunghoon Heo

"""
import re
from SeqCache import SeqCache , SeqRecorder

class sub_recorder(object):
	def __init__(self , data):
		self.query = data[0]
		self.flag = data[1]
		self.rname = data[2]
		self.pos = int(data[3])
		self.mapq = int(data[4])
		self.cigar = data[5]
		self.rnext = data[6]
		self.pnext = data[7]
		self.slen = int(data[8])
		self.seq = data[9]
		try:
			self.qual  = data[10]
		except IndexError:
			self.qual = None
		try:
			self.misc = data[11]
		except IndexError:
			self.misc = None
	def reform(self):
		reform_str = []
		reform_str.append(self.query)
		reform_str.append(self.flag)
		reform_str.append(self.rname)
		reform_str.append(self.pos)
		reform_str.append(self.mapq)
		reform_str.append(self.cigar)
		reform_str.append(self.rnext)
		reform_str.append(self.pnext)
		reform_str.append(slen)
		reform_str.append(seq)
		if self.qual != None:
			reform_str.append(self.qual)
		if self.misc != None:
			reform_str.append(self.misc)
			
		return "\t".join(reform_str)
class SamReader(object):

	def __init__(self, samfile) :
		self.sam = samfile
		self.headers = []
		self.header_read = False
	def open(self):
		self.sam_file = open(self.sam)

	def read_header(self):
		for line in self.sam_file:
			if line.startswith("@"):
				self.headers.append(line.rstrip())
				continue
			else:
				self.sam_file.seek(0)
				break
		self.header_read = True

	def read_record(self):
		for line in self.sam_file:
			
			if line.startswith("@") : 
				if self.header_read == False:
					self.headers.append(line.rstrip())
					continue
				else:
					continue

			data_field = {}
			data = line.rstrip().split("\t")

			rec = sub_recorder(data)

			yield rec

	def show_header(Self):
		if self.headers == [] :
			raise AttributeError("No header read. Please use read_header method")
		else:
			return self.headers

	def rewind(self):
		self.sam_reader.seek(0)

	def close(self):
		self.sam_file.close()

	def cigar_alpha(self,cigar):
		alphas = re.split("[MIDNSHPX=]+" , cigar)[:-1]
		return alphas

	def cigar_chars(self, cigar):
		chars = re.split("[0-9]+" , cigar)[1:]
		return chars

class IndelMarker(object):
	def __init__(self, sam_reader):
		self.reader = sam_reader
		self.cache = None
	def mark(self):
		seq_recs = []
		for sam_rec in self.reader.read_record():
			query = sam_rec.query
			seq = sam_rec.seq
			cigar = sam_rec.cigar
			qual = sam_rec.qual
			cigar_chars = self.reader.cigar_chars(cigar)
			if ('I' in cigar_chars) or ('D' in cigar_chars):
				query = query + "|INDEL"
			
			seq_rec = SeqRecorder(id=query, seq=seq,qual=qual)
			seq_recs.append(seq_rec)

		
		self.cache = SeqCache(seq_recs)

