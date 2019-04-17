#!/usr/bin/python

from Seq import reverse_complement, generate_kmers

class Graph:
	def __init__(self, p1obj, p2obj, k):
		self.k = k
		self.graph = {}
		self.p1 = p1obj
		self.p2 = p2obj

	def __getitem__(self, item):
		return self.graph[item]
	def pop(self,item):
		self.graph.pop(item)
	def has_key(self, item):
		return True if item in self.graph else False
	def items(self):
		for k , v in self.graph.items():
			yield (k,v)

	def keys(self):
		return self.graph.keys()

	def values(self):
		return self.graph.values()
	def nodes(self):
		return self.graph.keys()
		
	def add_vertex(self, v):
		try:
			self.graph[v]["depth"] += 1
		except KeyError:
			self.graph[v] = { "depth" : 1, "edges" : set() }

	def add_edge(self, v1, v2):
		self.graph[v1]["edges"].add(v2)
	def nodes_incoming_to_current_node(self, node):
		nodes = set()
		for n in [ x + node[:-1] for x in "ATGCN" ]:
			if n not in self.graph:
				continue
			if node in self.graph[n]["edges"]:
				nodes.add(n)
		return nodes
	def nodes_outgoing_from_current_node(self, node):
		return self.graph[node]["edges"]
	def build_naive(self):
		for _ , seq , _ in self.p1.parse():
			if len(seq) < self.k +1 :
				continue
			rc = reverse_complement(seq)
			for km1 , km2 in zip(list(generate_kmers(seq , self.k + 1)) , list(generate_kmers(seq, self.k + 1))):
				lkmer1 = km1[:-1]
				rkmer1 = km1[1:]
				lkmer2 = km2[:-1]
				rkmer2 = km2[1:]
				self.add_vertex(lkmer1); self.add_vertex(rkmer1)
				self.add_vertex(lkmer2); self.add_vertex(rkmer2)
				self.add_edge(lkmer1, rkmer1)
				self.add_edge(lkmer2, rkmer2)
		for _ , seq , _ in self.p2.parse():
			if len(seq) < self.k +1 :
				continue
			rc = reverse_complement(seq)
			for km1 , km2 in zip(list(generate_kmers(seq , self.k + 1)) , list(generate_kmers(seq, self.k + 1))):
				lkmer1 = km1[:-1]
				rkmer1 = km1[1:]
				lkmer2 = km2[:-1]
				rkmer2 = km2[1:]
				self.add_vertex(lkmer1); self.add_vertex(rkmer1)
				self.add_vertex(lkmer2); self.add_vertex(rkmer2)
				self.add_edge(lkmer1, rkmer1)
				self.add_edge(lkmer2, rkmer2)
	def num_nodes(self):
		return len(self.graph.keys())
	def is_island(self, n):
		if n not in self.graph:
			return False
		if len(self.graph[n]["edges"]) == 0 and len(self.nodes_incoming_to_current_node(n)) == 0 :
			return True
		else:
			return False
	def is_outer_tip(self, n):
		if len(self.graph[n]["edges"]) == 0 and len(self.nodes_incoming_to_current_node(n)) != 0 :
			return True
		else:
			return False
	def is_incoming_tip(self, n):
		if len(self.graph[n]["edges"]) != 0 and len(self.nodes_incoming_to_current_node(n)) == 0 :
			return True
		else:
			return False
	def is_tip(self, n):
		if n not in self.graph:
			return False
		return self.is_outer_tip(n) or self.is_incoming_tip(n)
	def RemoveNodeSet(self, nodeset):
		for n in nodeset:
			try:
				incomings = self.nodes_incoming_to_current_node(n)
				outgoings = self.nodes_outgoing_from_current_node(n)
			except KeyError:
				continue
			for i in incomings :
				try:
					self.graph[i]["edges"].remove(n)
					#self.graph[i]["depth"] -= self.graph[n]["depth"]
				except KeyError:
					continue
			try:
				self.graph.pop(n)
			except KeyError:
				continue