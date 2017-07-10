#!/usr/bin/python
"""
Author : Sunghoon Heo
"""

# Modification report
# 2016-02-25
# Since graph is now using merged graph, in DFS, using graph directly
# 2016-04-05
# Average path depth calculation code applied

import sys
import numpy as np

class Walker2:
	def __init__(self,graphType,covMap,start=None,end=None, mode='dfs',max_iter=1000000) :
		self.mode = mode
		self.graphType = graphType
		self.start = start
		self.end = end
		self.covMap = covMap
		self.max_iter = max_iter
	def walk(self) :
		if self.mode == 'dfs' :
			paths = []
			for path in self.DFS(self.graphType,self.start,self.end,self.covMap,self.max_iter) :
				paths.append(path)
			return paths
		elif self.mode == 'dfsm' :
			paths = []
			for path in self.DFSMergeResult(self.graphType,self.start,self.end,self.covMap) :
				paths.append(path)
			return paths
		elif self.mode == 'dfss' :
			paths = []
			for path in self.DFSSingleSrc(self.graphType, self.start) :
				paths.append(path)
			return paths
		else :
			return self.walkWithSeedAndEnd(self.graphType,self.start,self.end)
			
	def walkWithSeedAndEnd(self, graphInstance, start ,end,path=[]) :
		# Walk through graph using paths
		
		nodes = graphInstance.bdg[start]
		path = path + [start]
		paths = []
		paths.append(path)
		#sys.stderr.write("[INDEL CALLER] RECURSIVE PATH TRACKER : %s\n"%(start))
		for node in nodes :
			if node == end :
				path.append(node)
				paths.append(path)
				return paths
			if node not in path :
				newpaths = self.walkWithSeedAndEnd(graphInstance,node ,end , path)
				for newpath in newpaths :
					paths.append(newpath)
					
		return paths
	
	def DFS(self, graphInstance, start , end, covMap,max_iter) :
		#http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
		# cov = covMap[start]
		# stack = [(start,[start],cov)]
		try:
			stack = [(graphInstance[start]["depth"],start,[start],1, [graphInstance[start]["depth"]], {start:1.0})]  # total depth , node , path , node used , depth list, ratio list
		except KeyError:
			stack = []
			yield ("NOSEED")
		iter = 0
		while stack :
			# (vertex , path, cov) = stack.pop()
			(depth, vertex, path, node_used, d_list, r_list) = stack.pop()
			# print Statement below added 2016-02-20 for debuggung
			#sys.stderr.write("[INDEL CALLER] DFS PATH TRACKER : %s\n"%(vertex))
			
			#2016-02-25
			if not graphInstance.has_key(vertex) : 
			# if not graphInstance.has_key(vertex) :
				# yield path , cov
				#yield depth, path , node_used,d_list, r_list , "BROKEN" , iter
				pass
				#print path
				#sys.stderr.write("#"*100 + "\n")
				#break
			#2016-02-25
			# if graphInstance[vertex] == [] :
				# yield path
			# for next in set(graphInstance[vertex]) - set(path) :
			if iter > max_iter :
				yield depth, path , node_used ,d_list , r_list , "BROKEN" , iter
				#print path
				break
			branch = set(graphInstance[vertex]["edges"])
			#print branch
			
			has_branch = False
			
			
			edge_sum = 0.0
			ed_set = {}
			for e in graphInstance[vertex]["edges"] :
				ed = graphInstance[e]["depth"]
				ed_set[e] = ed
				edge_sum += ed

			for e in ed_set :
				ed_set[e] = ed_set[e] / edge_sum

			if len(graphInstance[vertex]["edges"]) == 0:
				#yield depth, path , node_used,d_list, r_list , "BROKEN" , iter
				pass
				#print path
				#break

			for next in set(graphInstance[vertex]["edges"]) - set(path) :
				iter += 1
				

				if next == end :
				# if end in next:
					#sys.stderr.write("[DFS] HIT : %s\n" %(next))
					# cov += covMap[next]
					# yield path + [next] , cov
					depth += graphInstance[next]["depth"] ## add up all depth
					node_used += 1
					d_list.append(graphInstance[next]["depth"])
					r_list[next] = ed_set[next]
					yield depth,path + [next],node_used , d_list, r_list
				else:
					#sys.stderr.write("[DFSTRACKER] %s\n"%(next))
					# cov += covMap[next]
					# stack.append((next,path + [next],cov))
					node_used += 1
					
					try:
						depth += graphInstance[next]["depth"]
						d_list.append(graphInstance[next]["depth"])
						r_list[next] = ed_set[next]
						"""
						if has_branch :
							r_list.append(ed_set[next])
							print next , ed_set[next]
						"""
						stack.append((depth , next , path + [next], node_used , d_list, r_list))
					except KeyError :
						continue
	def DFSMergeResult(self,graphInstance, start , end, covMap) :
		#http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
		# cov = covMap[start]
		# stack = [(start,[start],cov)]
		stack = [(start,[start])]
		while stack :
			# (vertex , path, cov) = stack.pop()
			(vertex,path) = stack.pop()
			# print Statement below added 2016-02-20 for debuggung
			#sys.stderr.write("[INDEL CALLER] DFS PATH TRACKER : %s\n"%(vertex))
			
			#2016-02-25
			if not graphInstance.bdg.has_key(vertex) : 
			# if not graphInstance.has_key(vertex) :
				# yield path , cov
				yield path
				#sys.stderr.write("#"*100 + "\n")
				continue
			#2016-02-25
			if graphInstance.bdg[vertex] == [] :
				yield path
			# for next in set(graphInstance[vertex]) - set(path) :
			for next in set(graphInstance.bdg[vertex]) - set(path) :
				
				# if next == end :
				if end in next:
					#sys.stderr.write("[DFS] HIT : %s\n" %(next))
					# cov += covMap[next]
					# yield path + [next] , cov
					yield path + [next]
				else:
					#sys.stderr.write("[DFSTRACKER] %s\n"%(next))
					# cov += covMap[next]
					# stack.append((next,path + [next],cov))
					stack.append((next , path + [next]))
	def DFSSingleSrc(self, graphInstance, start) :
		visited , stack = set(), [ start ]
		while stack :
			vertex = stack.pop()
			if vertex not in visited :
				visited.add(vertex)
				if len(graphInstance.bdg[vertex]) == 0 :
					yield visited
				else:
					stack.extend(set(graphInstance.bdg[vertex]) - visited)
	#	return list(visited)
			
