#!/usr/bin/python

"""
Author : Sunghoon Heo
"""

import os
import random

import numpy as np

class GraphModifier(object):
	def __init__(self) :
		pass
	
	def detect_dead_end(self, G) :
		# G is a de Bruijn Graph
		dead_ends = []
		before = ''
		current = ''
		nxt = ''
		# Selects a node randomly
		# Dead end criterion
		# if there is a node which has no outgoing edge then discard
		iter = 0
		marked_dead = []
		for node in G :
			edges = G[node]
			if len(edges) == 0 :
				marked_dead.append(node)
				
		for marked in marked_dead :
			try:
				G.pop(marked)
			except KeyError :
				continue
		# while iter < 1000000 :
			# before = current = nxt = random.choice(G.keys())
			# if len(G[current]["edges"]) == 0 :
				# G.pop(current)
				# iter += 1
				# continue
			# else:
				# for edge in G[current]["edges"] :
					# nxt = edge
					# iter += 1
					# if len(G[nxt]["edges"]) == 0 :
						# G.pop(nxt)
				# #self.detect_dead_end(G)
				# continue
				
	def reshapeGraphWithDepthRatio(self, tmp_file , ratio = 0.30) :
		# For memory issue, previous graph must be cleared up
		# Less than 5 % is discarded due to the Sanger  sequencing signal error
		reshaped_graph = {}
		# Read tmp file
		with open(tmp_file) as f :
			for record in f :
				L = record.rstrip().split("\t")
				if len(L) < 4 :
					# In the case when no edges connected
					continue
				node = L[0]
				node_depth = L[1] 
				edges = L[2]
				edge_depths = L[3]
				if not reshaped_graph.has_key(node) :
					reshaped_graph[node] = {"edges" : [] , "depth" : int(node_depth) , "in" : 0 , "out" : 0}
				# computing ratio for each edge_depths
				edge_set = edges.split(",")
				
				edge_depth_set = [int(x) for x in edge_depths.split(",")]
				
				numpy_edge_depth_set = np.array(edge_depth_set , dtype=np.int32)
				total_sum = np.sum(edge_depth_set,dtype=np.float32)
				numpy_ratio = numpy_edge_depth_set / total_sum
				
				for i , r in enumerate(numpy_ratio) :
					if r < ratio :
						# discard this edge since invalid
						continue
					else:
						edge_selected = edge_set[i]
						edge_depth_selected = edge_depth_set[i]
						edge_full_seq = node[1:] + edge_selected
						
						reshaped_graph[node]["edges"].append(edge_full_seq) # Add to current node
						reshaped_graph[node]["out"] += 1
						try :
							reshaped_graph[edge_selected]
							reshaped_graph[edge_selected]["depth"] += edge_depth_selected
							reshaped_graph[edge_selected]["in"] += 1
						except KeyError :
							reshaped_graph[edge_selected] = {"depth" : edge_depth_selected , "in" : 1 , "out" : 0 , "edges" : [] }
							
		return reshaped_graph
