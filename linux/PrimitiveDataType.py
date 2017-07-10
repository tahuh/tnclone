#!/usr/bin/python
"""
Author : Sunghoon Heo
"""
class Node:
	def __init__(self,kmer,freq) :
		self.kmer = kmer
		self.freq = freq
		self.edges = None # designate end node
		self.indeg = 0
		self.outdeg = 0
		self.idx = 0      # index of node
		self.sortIdx = 0
		visit = False
	def __str__(self) :
		return self.kmer
	def __hash__(self) :
		return hash(self.kmer)

class Edge:
	def __init__(self, node1, node2) :
		self.aparture = node1
		self.departure = node2
		self.idx = 0
		self.lidx = node1.idx # from node index
		self.ridx = node2.idx # to node index

class BiDiGraph :
	def __init__(self) :
		self.bdg = {}
		self.nodes = []  # something like vector<Node> in C++
		self.edges = []  # something like vector<Edge> in C++
		self.boolTable = {}
		self.nvertex = 0
		self.nedges = 0
		self.nodeIdx = 0 # Current node index if given node puts into the graph
		self.edgeIdx = 0 # Current edge index if given edge is assigned to the graph
		
	def pushNode(self,node) :
		node.idx = self.nodeIdx
		self.nodes.append(node)
		self.nodeIdx += 1
		
	def pushEdge(self, edge) :
		edge.idx = self.edgeIdx;
		self.edges.append(edge)