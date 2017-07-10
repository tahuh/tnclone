#!/usr/bin/python

from DeBruijn import DeBruijnFromSequencingRead as DBFSR
from Traverse2 import Walker2
from FastQParser import SequencingReadParser as SRP
from GraphModifier import GraphModifier as GM
from DumpIntermediate import Dump

class Assembler(object) :
	def __init__(self, kfs, src, snk, outname, kmer_size, ffq, rfq, tmp) :
		self.kfsname = kfs
		self.src = src
		self.snk = snk
		self.outname = outname
		self.kmer_size = kmer_size
		self.ffq = ffq
		self.rfq = rfq
		self.tmp = tmp
		self.graph = None
	def buildGraph(self) :
		graph = DBFSR(k=self.kmer_size)
		
		graph.loadKFS(self.kfsname)
		srp = SRP(self.ffq , self.rfq)
		srp.readFwd()
		srp.readRev()
		graph.loadArchive(srp.fwdread, rc = 'f')
		graph.loadArchive(srp.revread, rc = 'r')
		del srp

		graph.build()
		self.graph = graph

	def dump2Tmp(self, tmpfile) :
		dumper = Dump(tmpfile)
		dumper.open()
		dumper.write(self.graph.bdg)
		dumper.close()

	def modifyGraph(self, tmpfile) :
		gm = GM()
		reshaped_graph = gm.reshapeGraphWithDepthRatio(tmpfile)
		gm.detect_dead_end(reshaped_graph)
		self.graph.dbg = reshaped_graph

	def assemble(self):
		
		walker = Walker2(self.graph.bdg , self.graph.kfsMap, self.src, self.snk)
		return walker

	def clearAssembler(self) :
		self.graph.clearArchive()
		self.graph.clearKFSMap()
