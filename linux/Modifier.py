#!/usr/bin/python

"""
class Modifier
Author : sunghoon heo
"""

import numpy as np

class Modifier2(object):
	def __init__(self, graph, node_depth, edge_freq=0.05):
		self.graph = graph
		self.node_depth = node_depth
		self.allele_freq = edge_freq


	def enforce_modify(self):
		island_pass = False
		dead_end_pass = False
		low_depth_pass = False
		low_allele_pass = False

		while 1:
			if len(self.graph.keys()) == 0 :
				print "Too much graph reducing resulted in graph elimination. Do DEBUG!!!"
				print "Possible option one can use is excute TnClone until Trimming process."
				print "Then plot diagnosis plot by pressing diagnosis button next to run button"
				raise ValueError("TnClone raised ValueError while reducing graph.")

			island_pass = self.remove_islands()
			dead_end_pass = self.remove_dead_end()
			low_depth_pass = self.remove_below_depth_cut_node()
			low_allele_pass = self.remove_below_allele_freq_edges()

			if island_pass and dead_end_pass and low_depth_pass and low_allele_pass :
				break

			


	def remove_islands(self):
		island_nodes = []
		### Search step
		for node in self.graph:
			edges = self.graph[node]["edges"]
			if len(edges) != 0 :
				### To accelerate algorithm
				### If we compute out degree first and it is not 0 then it is not island
				continue
			incoming_nodes = self._gen_target_incomings(node)
			
			if len(incoming_nodes) == 0 and len(edges) == 0:
				island_nodes.append(node)

		### Removing step
		if len(island_nodes) == 0: return 1 ### To very there is no islands

		for island in island_nodes :
			self._remove_vertex(island)
		return -1


	def remove_dead_end(self):
		### Search step
		dead_end_nodes = []
		for node in self.graph:
			if self._is_dead_end(self.graph, node):
				### Add to removal storage
				dead_end_nodes.append(node)
				### Search for incoming nodes
				incoming_nodes = self._gen_target_incomings(node)
				
				### Now it is time to remove connection information from incoming nodes
				for incoming in incoming_nodes:
					self._remove_edge_nodes_connected_to_incoming_nodes(incoming, node)

		### Remove dead end vertices from graph

		if len(dead_end_nodes) == 0 : return 1 ## Verify there is no dead ends

		for dead_end in dead_end_nodes:
			self._remove_vertex(dead_end)

		return -1

	def remove_below_depth_cut_node(self):
		### Search step
		below_depth_cut_nodes = []
		for node in self.graph:
			if self._pass_depth_cut(self.graph, node) : continue
			
			### Add to storage
			below_depth_cut_nodes.append(node)
			incoming_nodes = self._gen_target_incomings(node)
			for incoming in incoming_nodes:
				self._remove_edge_nodes_connected_to_incoming_nodes(incoming, node)

		if len(below_depth_cut_nodes) == 0 : return 1 ### Verify that there is no low depth node

		for low_depth_node in below_depth_cut_nodes:
			self._remove_loop(low_depth_node)
			self._remove_vertex(low_depth_node)
		return -1 ### Verify that analysis has done

	def remove_below_allele_freq_edges(self):
		### Search step
		selected_edges = []
		for node in self.graph:
			edges = self.graph[node]["edges"]

			edge_depth_set = []
			total_depth = 0
			for edge in edges :
				edge_depth = self.graph[edge]["depth"]
				total_depth += edge_depth
				edge_depth_set.append(edge_depth)

			### convert to numpy array for vectorized computation
			numpy_edge_depth_set = np.array(edge_depth_set ,dtype=np.float32)
			numpy_edge_depth_set /= total_depth

			for z in zip(edges , numpy_edge_depth_set):
				## z = ( edge , ratio ) , tuple
				ratio = z[1]
				if self._pass_edge_ratio_cut(ratio) : continue

				### Remove connection information

				self.graph[node]["edges"].remove(z[0])
				#if len(selected_edges) == 0 :
				if z[0] not in selected_edges:
					selected_edges.append(z[0])

		if len(selected_edges) == 0 : return 1 ### verify there is no edges below allele freq cut

		for node in selected_edges :
			for incoming in self._gen_target_incomings(node):
				self._remove_edge_nodes_connected_to_incoming_nodes(incoming, node)
			self._remove_loop(node)
			self._remove_vertex(node)	
		
		return -1 ### analysis completed

	### Accesory functions

	def _remove_loop(self, vertex) :
		for e in self.graph[vertex]["edges"]:
			### remove connection to prevent loop structure
			if e == vertex :
				self.graph[vertex]["edges"].remove(e)

	def _remove_vertex(self, vertex):
		self.graph.pop(vertex)
		return
	def _remove_edge_nodes_connected_to_incoming_nodes(self, incoming, node):
		if node in self.graph[incoming]["edges"]:
			### To ensure that there is connection between incmoning and target node
			self.graph[incoming]["edges"].remove(node)
		return

	def _pass_depth_cut(self, graph, node):
		return graph[node]["depth"] > self.node_depth
				
	def _pass_edge_ratio_cut(self, ratio):
		return ratio > self.allele_freq

	def _gen_target_incomings(self,node) :
		incoming_nodes = []
		for incoming in self._gen_incomings(node):
			try:
				assert self.graph.has_key(incoming)
				incoming_nodes.append(incoming)
			except AssertionError:
				continue
		return incoming_nodes

	def _gen_incomings(self , node):
		ret = []
		frag = node[:-1]
		for b in "ATGC" :
			yield b + frag

	def _is_dead_end(self, graph, node):
		### dead end : No out going nodes present for given node
		return len(graph[node]["edges"]) == 0
	
""" Below class will be deprecated soon """
class Modifier(object):
	def __init__(self, graph, node_depth, edge_freq=0.05):
		self.graph = graph
		self.node_depth = node_depth
		self.allele_freq = edge_freq

	def remove_dead_ends(self):
		nodes_to_remove = []
		for node in self.graph:
			edges = self.graph[node]["edges"]
			if len(edges) == 0 :
				is_island = False
				cnt = 0
				tmp_nodes = []
				for b in "ATGC":
					node_from = b + node[:-1]
					if not self.graph.has_key(node_from):
						cnt += 1
					else:
						tmp_nodes.append(node_from)

				if cnt == 4 :
					is_island= True

				if is_island :
					nodes_to_remove.append(node)

				else:
					for x_node in tmp_nodes:
						try:
							self.graph[x_node]["edges"].remove(node)
							nodes_to_remove.append(node)
						except ValueError:
							continue

		already_removed = []
		for target in nodes_to_remove:
			try:
				self.graph.pop(target)
				already_removed.append(target)
			except KeyError:
				if target in already_removed :
					continue
				else:
					print "ARKKKKKKKKKKKKKKKKKKK"
					

	def remove_nodes_depthcut(self):
		targets = []
		for node in self.graph:
			depth = self.graph[node]["depth"]
			if depth < self.node_depth:
				targets.append(node)

		for target in targets :
			cnt = 0
			tmp_nodes = []
			for b in "ATGC" :
				node_from = b + target[:-1]
				if self.graph.has_key(node_from):
					tmp_nodes.append(node_from)

			for tmp in tmp_nodes:
				try:
					self.graph[tmp]["edges"].remove(target)
				except ValueError:
					continue

			self.graph.pop(target)

	def remove_nodes_edge_freq(self):
		for node in self.graph:
			edges = self.graph[node]["edges"]
			e_depths = []
			for e in edges :
				e_depth = self.graph[e]["depth"]
				e_depths.append(e_depth)

			numerator = sum(e_depths)
			for i in xrange(len(e_depths)):
				e_depths[i] = e_depths[i] / float(numerator)

			rm_edge_target = []
			for j , e_freq in enumerate(e_depths):
				if e_freq < self.allele_freq:
					rm_edge_target.append(edges[j])

			for target in rm_edge_target:
				self.graph[node]["edges"].remove(target)
				other_connections = []
				for b in "ATGC":
					new_node = b + node[:-1]
					if self.graph.has_key(new_node):
						other_connections.append(new_node)
				for other in other_connections:
					try:
						self.graph[other]["edges"].remove(target)
					except ValueError:
						continue

