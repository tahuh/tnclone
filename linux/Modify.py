#!/usr/bin/python

class MinimumDepthKmerCollector:
	"""
		count_obj : KmerCounter object
		cut       : cutoff value (minimum k-mer occurence)
	"""
	def __init__(self, count_obj, cut):
		self.obj = count_obj
		self.cutoff = cut
		self.kmers_to_reduce = {}

	def Action(self):
		Counts = self.obj.GetCountMap()
		for km , v in Counts.items():
			if v <= self.cutoff:
				self.kmers_to_reduce[km] = False

	def IsThisKmerTargetOfRemoval(self, km):
		try:
			self.kmers_to_reduce[km]
			return True
		except KeyError:
			return False

	def IsThisPreviouslyRemoved(self, km):
		try:
			return self.kmers_to_reduce[km]
		except KeyError:
			return None
				
				
class GraphModificationEngine:
	"""
		graph - Graph object
		obj   - MinimumDepthKmerCollector object
		nton  - Node to node edge freq
		multi - multiple edge
	"""
	def __init__(self, graph, obj, nton, multi):
		nodeset = set(obj.kmers_to_reduce.keys())
		graph.RemoveNodeSet(nodeset)
		self.g = graph
		self.nton = nton
		self.multi = multi
	def Action(self):
		MAXITERATION = 1000
		it = 0
		while 1:
			tip_nodes = set()
			island_nodes = set()
			for n in self.g.nodes():
				if self.g.is_tip(n):
					tip_nodes.add(n)
				if self.g.is_island(n):
					island_nodes.add(n)
			self.g.RemoveNodeSet(island_nodes)
			self.g.RemoveNodeSet(tip_nodes)
			
			edge_freq_removed = self.edge_freq_removal()
			
			if len(tip_nodes) ==0 and len(island_nodes) == 0 and edge_freq_removed == False:
				break
			it += 1
			if len(self.g.nodes()) == 0 :
				break
			if it > MAXITERATION:
				break
	def ShowGraph(self):
		return self.g

	def edge_freq_removal(self):
		nton_list = set()
		multi_list = set()
		total_nodes = self.g.nodes()
		multi_perf = False
		nton_perf = False
		for n in total_nodes:
			### Iterate over all nodes
			outgoing_nodes = list(self.g.nodes_outgoing_from_current_node(n))
			
			### First remove multiple nodes
			multiple_edges_depths = []
			multiple_total = 0
			for o in outgoing_nodes:
				multiple_edges_depths.append(self.g[o]["depth"])
			multiple_total = sum(multiple_edges_depths)
			if multiple_total != 0:
				for i, x in enumerate(multiple_edges_depths):
					multiple_edges_depths[i] = multiple_edges_depths[i] / float(multiple_total)
					if multiple_edges_depths[i] < self.multi:
						multi_list.add(outgoing_nodes[i])
				if len(multi_list) != 0 :
					self.g.RemoveNodeSet(multi_list)
					multi_perf = True
			multi_list = set()
			
			### From alived nodes perform node to node strategy
			outgoing_nodes = list(self.g.nodes_outgoing_from_current_node(n))
			multiple_edges_depths = []
			for o in outgoing_nodes:
				multiple_edges_depths.append(self.g[o]["depth"])
				
			my_depth = self.g[n]["depth"]
			if my_depth != 0:
				for i, x in enumerate(multiple_edges_depths):
					multiple_edges_depths[i] = multiple_edges_depths[i] / float(my_depth)
					if multiple_edges_depths[i] < self.nton:
						nton_list.add(outgoing_nodes[i])
				if len(nton_list) != 0 :
					self.g.RemoveNodeSet(nton_list)
					nton_perf = True
			nton_list = set()
		if nton_perf == False and multi_perf == False:
			return False
		else:
			return True