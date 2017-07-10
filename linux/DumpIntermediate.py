#!/usr/bin/python

class Dump:
	def __init__(self,file_name) :
		self.fname = file_name
		self.F = None
	def open(self) :
		self.F = open(self.fname, "w")
	def close(self) :
		self.F.close()
	def write(self, graph):
		for cur_node in graph:
			
			depth = graph[cur_node]["depth"]
			edges = graph[cur_node]["edges"]
			try:
				edge_depths = [str(graph[e]["depth"]) for e in edges ]
				self.F.write(cur_node + "\t" + str(depth))
			
				edge_bases = [x[-1] for x in edges]
				#print cur_node, str(depth) , ",".join(edge_bases), ",".join(edge_depths)
				self.F.write("\t" + ",".join(edge_bases))
				self.F.write("\t" + ",".join(edge_depths))
				self.F.write("\n")
			except KeyError:
				self.F.write(cur_node + "\t" + str(0) + "\n")
			
				
