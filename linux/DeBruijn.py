#!/usr/bin/python

"""
Author: Sunghoon Heo

Modification note
05-Apr-2016 Depth Calculation Code block upadted
"""

from PrimitiveDataType import Node, Edge
from PrimitiveDataType import BiDiGraph as BDG
from DNAStrManipulate import convert2barr, barr2seq
### General constants
A=0x0
T=0x3
G=0x1
C=0x2
MASK=0x3

class DeBruijnFromSequencingRead(BDG):
	def __init__(self,k=63) :
		BDG.__init__(self)
		self.reads = []
		self.ksize = k
	def loadKFS(self , sai) :
		self.kfsMap = {}
		with open(sai) as saifile :
			for record in saifile :
				data = record.rstrip().split("\t")
				kmer , freq =  data[0] , int(data[1])
				# self.kfsMap[convert2barr(kmer)] = freq
				self.kfsMap[kmer] = freq
	@staticmethod			
	def RC(string) :
		switcher = dict(A='T',T='A',G='C',C='G',N='N')
		complement = ''
		for base in string :
			complement = switcher[base] + complement
		return complement
		
	def loadArchive(self, archive , rc = 'f') :
		if rc == 'f' :
			self.reads.extend(archive)
		elif rc == 'r' :
			for seq in archive :
				self.reads.append(self.RC(seq))
	@staticmethod
	def chop(string,ksize) :
		cutsize = ksize + 1
		for i in range(len(string) - cutsize+1) :
			yield string[i:i+cutsize]
		
	

	def build(self) :
		""" Construct directed de Bruijn graph """
		
		for read in self.reads :
			for k1mer in self.chop(read,self.ksize) :
				# leftK = convert2barr(k1mer[:-1])
				# rightK = convert2barr(k1mer[1:])
				leftK = k1mer[:-1]
				rightK = k1mer[1:]
				if leftK == None : continue
				if rightK == None : continue
				# Search if start node has enougph coverage
				
				if not self.kfsMap.has_key(leftK) : continue
				if not self.kfsMap.has_key(rightK) : continue
				
				# Check if left kmer is already inserted in graph
				
				if not self.bdg.has_key(leftK)  :
					# Not in graph
					self.bdg[leftK] = {"in" : 0 , "out" : 0 , "depth" : self.kfsMap[leftK] ,"edges" : [] }
					# Judge if right K-mer passed coverage constraint
					if not self.kfsMap.has_key(rightK) : continue
					self.bdg[leftK]["edges"].append(rightK)
					self.bdg[leftK]["out"] += 1
					if self.bdg.has_key(rightK) :
						# if current graph already has right k-mer
						self.bdg[rightK]["in"] += 1
					#	self.bdg[rightK]["depth"] += 1
					else:
						self.bdg[rightK] = {"in" : 1, "out" : 0 , "depth" : self.kfsMap[rightK] , "edges" : [] }
				else:
					# Present in graph
					# Judge if right K-mer passed coverage constraint
					# self.bdg[leftK]["depth"] += 1
					if not self.kfsMap.has_key(rightK) : continue
					
					if rightK not in self.bdg[leftK]["edges"] :
						# Reduces duplicate insertion for path find
						self.bdg[leftK]["edges"].append(rightK)
						self.bdg[leftK]["out"] += 1

					if rightK in self.bdg :
					#	self.bdg[rightK]["depth"] += 1
						self.bdg[rightK]["in"] += 1
					else:
						self.bdg[rightK] = {"in":1,"out":0,"depth": self.kfsMap[rightK], "edges" : []}
						
	def clearArchive(self) :
		del self.reads[:]
	
	def clearKFSMap(self) :
		self.kfsMap.clear()
		
class NaiveDeBruijn(BDG) :
	def __init__(self) :
		BDG.__init__(self)
	
	def load(self,kindexfname) :
		SAIresult = open(kindexfname)
		for line in SAIresult :
			data = line.rstrip().split("\t")
			idx , kmer , freq = int(data[0]) , data[1] , int(data[2])
			lkmer = kmer[:-1]
			rkmer = kmer[1:]
			try:
				self.bdg[lkmer].add(rkmer)
			except KeyError :
				self.bdg[lkmer] = set([rkmer])
			
			try :
				self.bdg[rkmer].add(lkmer)
			except KeyError :
				self.bdg[rkmer] = set([lkmer])
			"""
			if self.bdg.has_key("TCCCATCGACCCACGATTTTTCATACGCTTGCATGAATCAGCGCTTTGGCGGGACGACGTTG") :
				print "Exist"
			"""
		SAIresult.close()
class DeBruijnGraph(BDG):
	def __init__(self) :
		BDG.__init__(self)
		
	def load(self, kindexfname) :
		# line specification
		# index kmer freq
		# Input file is SAI result file
		SAIresult = open(kindexfname)
		
		for line in SAIresult :
			data = line.rstrip().split("\t")
			idx , kmer , freq = int(data[0]) , data[1] , int(data[2])
			lkmer = kmer[:-1]
			rkmer = kmer[1:]
			"""
			if not self.bdg.hash_key(lkmer):
				# If there is no current node as start node in graph
				# Add this node to graph first
				lkmerNode = Node(lkmer , freq)
				
				
			"""
			if self.bdg.has_key(lkmer): # check if already this k-mer is loaded
				
				if rkmer in self.bdg[lkmer] :
					continue
				else :
					# if right k-mer is not present in graph
					node = Node(rkmer,freq)
					lkmerNode = self.nodes[self.nodeIdx-1]
					node.sortIdx = idx
					node.freq = freq
					# insert node to graph
					self.pushNode(node)
					
					# insert edge to graph
					edge = Edge(lkmerNode , node)
					self.pushEdge(edge)
					
					# Connect two nodes
					self.bdg[lkmer].append(rkmer)
					self.boolTable[rkmer] = False
					lkmerNode.outdeg += 1
					node.indeg += 1
					self.nvertex += 1
					self.nedges += 1
					
					
			else:
				# Newly introduced k-1 mer node
				lkmerNode = Node(lkmer,freq)
				self.pushNode(lkmerNode)
				self.bdg[lkmer] = [rkmer]
				rkmerNode = Node(rkmer,freq)
				self.pushNode(rkmerNode)
				edge = Edge(lkmerNode, rkmerNode)
				self.pushEdge(edge)
				lkmerNode.outdeg += 1
				rkmerNode.indeg += 1
				lkmerNode.edges = [rkmerNode]
				
				lkmerNode.sortIdx = idx
				rkmerNode.sortIdx = idx
				lkmerNode.freq = freq
				rkmerNode.freq = freq
				self.nvertex += 2
				self.nedges += 1
				self.boolTable[lkmer] = False
				self.boolTable[rkmer] = False
			
		SAIresult.close()
	
	
		
