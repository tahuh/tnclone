#!/usr/bni/python

import subprocess

class ErrorSeedCorrector(object):
	def __init__(self, graph, seed, path, log, max_mis=3,use_depth_profile=False):
		self.graph = graph
		self.seed = seed
		self.max_mis = max_mis
		self.corrected_seed = None
		self.path = path
		#self.nvidia = nvidia_enable
		self.seed_file = self.path + "/" + "seed.fa"
		self.use_depth_profile = use_depth_profile
		self.log = log

	def select(self):
		if self.use_depth_profile :
			self._correct_wrt_depth()
		else:
			self._correct_wrt_max_mis()
	def _correct_wrt_depth(self) :
		#self.log.append("[LOGGER] Using depth criteria for seed selection. This is so heuristic. Please be aware....\n")
		keyset = self.graph.keys()
		cnt_set = {}
		self.seed_file_out = open(self.seed_file , "w")
		for k in key_set :
			if cnt_set.has_key(k) : continue
			else: cnt_set[k] = self.graph[k]["depth"]
		max_node_depth = max(cnt_set.values())
		factor = 0.7
		min_node_depth = int(factor * max_node_depth)
		### select node

		for k in cnt_set :
			if cnt_set[k] >= min_node_depth and cnt_set[k] <= max_node_depth:
				self.seed_file_out.write(k + "\n")
		self.seed_file_out.close()
		self.log.append("[LOGGER] seed file can be found at %s\n"%(self.seed_file))
	def _correct_wrt_max_mis(self):
		#self.log.append("[LOGGER] Using mis-match criteria for seed selection.\n")
		keyset = self.graph.keys()
		### Very greedy
		seed_file = open(self.seed_file , "w")
		#if not self.nvidia:
		#self.log.append("[LOGGER] Searching all possible seeds. This might take a while\n")
		for node in keyset :
			cnt = self._string_cmp(self.seed, node)
			if cnt <= self.max_mis :
				seed_file.write(node + "\n")
		#else:
		#	self.log.append("[LOGGER] NVidia graphics card is available. Your device information is listed below.\n")
		#	device_query = "./deviceQuery > deviceInfo.txt"
		#	subprocess.call(device_query)
		#	for line in open("./deviceInfo.txt") :
		#		self.log.append(line)

			### Dump all nodes
		#	dmp_file = open("./cuda_dump.txt","w")
		#	for k in keyset :
		#		dmp_file.write(k + "\n")
			

		#	self._fast_string_compare("./cuda_dmp.txt" , self.seed_file)
		#self.log.append("[LOGGER] seed file can be found at %s\n"%(self.seed_file))

	def _string_cmp(self, s1, s2) :
		if s1 == s2 :
			return 0
		else:
			cnt = 0
			for z in zip(s1,s2):
				if z[0] != z[1] : cnt += 1

			return cnt

	#def _fast_string_compare(self, in1, out):
		### At here we will use awesome thing
		### The data paralization!!!
		### The GPU thing!!!!
	#	script = "./kmer_compare -k %s -s %s -o %s"%(self.seed, in1, out)
	#	subprocess.call(script , shell = True)
		
