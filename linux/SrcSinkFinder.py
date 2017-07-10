#!/usr/bin/python

import sys
import os

class SrcSinkFinder:
	@staticmethod
	def config(path):
		return os.path.dirname(os.path.abspath(path)) + "/" + path.split("/")[-1]
	def __init__(self, assemResult, Ofile,src,snk,editor) :
		self.In = assemResult
		self.Out = Ofile
		self.src = src
		self.snk = snk
		self.contigs = {}
		self.editor = editor
	def loadResult(self) :
		#self.editor.append("[TnClone::SrcSnkFinder] Src&Snk In file : %s\n"%(self.config(self.In)))
		#sys.stderr.write("[SH ASSEMBLER MACHINERY] Src&Snk In file : %s\n"%(self.config(self.In)))
		with open(self.config(self.In)) as I :
			header = ''
			seq = ''
			for line in I :
				if line.startswith(">") :
					header = line.rstrip()
				else:
					seq = line.rstrip()
					self.contigs[header] = seq
	def judge(self, seq) :
		return (seq[:len(self.src)] == self.src) and (seq[-len(self.snk):] == self.snk)
	def run(self) :
		self.loadResult()
		o = open(self.config(self.Out) , "w")
		#self.editor.append("[TnClone::SrcSnkFinder] Src&Snk Out file : %s\n"%(self.config(self.Out)))
		#sys.stderr.write("[SH ASSEMBLER MACHINERY] Src&Snk Out file : %s\n"%(self.config(self.Out)))
		npass = 0
		for ID , SEQ in self.contigs.iteritems() :
			reslt = self.judge(SEQ)
			if(reslt) :
				if ID.startswith(">contig"):
					npass += 1
					o.write(ID + "\n" + SEQ + "\n")
		o.close()
		#self.editor.append("[TnClone::SrcSnkFinder] Src&Snk %d / %d passed\n"%(npass,len(self.contigs.keys())))
		#sys.stderr.write("[SH ASSEMBLER MACHINERY] Src&Snk %d / %d passed\n"%(npass,len(self.contigs.keys())))
