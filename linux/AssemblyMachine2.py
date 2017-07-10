#!/usr/bin/python

from DeBruijn import DeBruijnFromSequencingRead as DBFSR
from Traverse2 import Walker2
from FastQParser import SequencingReadParser as SRP
#from GraphModifier import GraphModifier as GM
from DumpIntermediate import Dump
#from Modifier import Modifier
from Modifier import Modifier2 as Modifier

class Assembler(object) :
	def __init__(self, kfs, src, snk, outname, kmer_size, ffq, rfq, tmp, mink_occ) :
		self.kfsname = kfs
		self.src = src
		self.snk = snk
		self.outname = outname
		self.kmer_size = kmer_size
		self.ffq = ffq
		self.rfq = rfq
		self.tmp = tmp
		self.graph = None
		self.occ = mink_occ
	def set_src(self, s) :
		self.src = s
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
		modifier = Modifier(self.graph.bdg, self.occ, 0.05)
		try:
			modifier.enforce_modify()
		except ValueError:
			raise ValueError("TnClone raised ValueError in AssemblyMachinery due to no graph presence.")
		self.graph.bdg = modifier.graph

		"""
		modifier = Modifier(self.graph.bdg, 3, 0.05)
		modifier.remove_dead_ends()
		modifier.remove_nodes_depthcut()
		modifier.remove_nodes_edge_freq()
		self.graph.bdg = modifier.graph
		
		gm = GM()
		reshaped_graph = gm.reshapeGraphWithDepthRatio(tmpfile)
		gm.detect_dead_end(reshaped_graph)
		self.graph.dbg = reshaped_graph
		"""
	def assemble(self):
		
		walker = Walker2(self.graph.bdg , self.graph.kfsMap, self.src, self.snk)
		return walker

	def clearAssembler(self) :
		self.graph.clearArchive()
		self.graph.clearKFSMap()
