#!/usr/bin/python

"""
Stats.py

for stats, see manuscript

Assumes that TnClone already performed internal alignment between sequences


Author : Sunghoon Heo
"""

from collections import Counter
import copy
import numpy as np
from scipy.stats import chi2_contingency as chi2_test

class ContigStats2(object):
	def __init__(self, seq_set, ksize, graph):
		### seq_set must be in the form of (id , seq) which id is marked with indel or not
		self.ksize = ksize
		self.seqs = seq_set
		self.graph = graph
	def kmerize(self) :
		for i in xrange(len(self.seqs) - self.ksize - 1):
			yield self.seqs[i:i+self.ksize+1]


	def path_score(self):
		path_score = 0.0
		n_junction = 0
		e_freq_ratios = []
		for k1mer in self.kmerize():
			lkmer = k1mer[:-1]
			rkmer = k1mer[1:]

			connected_nodes = self.graph[lkmer]["edges"]
			if len(connected_nodes) >= 2:
				n_junction += 1
				e_freqs = {}
				for e in connected_nodes:
					e_freqs[e] = self.graph[e]["depth"]
				summation = sum(e_freqs.values())
				for e in connected_nodes:
					
					e_freqs[e] = e_freqs[e] / float(summation)
				e_freq_ratios.append(e_freqs[rkmer])

		if n_junction == 0 :
			path_score = 1.0
		else:
			path_score = sum(e_freq_ratios) / n_junction

		return path_score
	def cv(self):
		cv = 0.0
		depths = []
		
		for i in xrange(len(self.seqs) - self.ksize):
			kmer = self.seqs[i:i+self.ksize]
			depths.append(int(self.graph[kmer]["depth"]))

		avg = np.mean(depths)

		std = np.std(depths)

		cv = std / avg

		return avg , cv	
class ContigStats(object):
	def __init__(self, seq_set, ksize, graph):
		### seq_set must be in the form of (id , seq) which id is marked with indel or not
		self.ksize = ksize
		self.seqs = seq_set
		self.graph = graph
		self.p_val = None
		self.path_scores = []
		self.indel_seqs = self._filter_out_indel_seqs()
		self.rows , self.cols, self.junction_index_storage = self._set_dimension()
		self.major_minor_allele_dic = self._get_major_minor_allele()
	#### 2017-03-16 renewd version
	
	def _get_major_minor_allele(self):
		#### junction required
		if len(self.seqs) == 0 :
			return {}
			
		seq_target =self.seqs[0][1][1]
		ret = {}
		for idx in self.junction_index_storage:
			transformed_idx = idx - 1
			my_kmer = seq_target[transformed_idx:transformed_idx + self.ksize]
			next_kmers = self.graph[my_kmer]["edges"]
			ds = []
			for e in next_kmers:
				ds.append(self.graph[e]["depth"])
			max_val = max(ds)
			max_val_idx = ds.index(max_val)
			min_val = min(ds)
			min_val_idx = ds.index(min_val)
			ret[idx] = [max_val , min_val]
		return ret
	def _filter_out_indel_seqs(self):
		### Search for non-indel marked sequences
		ids = map(lambda x : x[0] , self.seqs)
		indel_index = []
		for i, id in enumerate(ids):
			if self._is_indel_marked(id):
				indel_index.append(id)
		indel_seqs = []	
		if len(indel_index) != 0 :
			for idx, id in enumerate(indel_index) :
				indel_seqs.append((ids[idx] , self.seqs[idx]))
				
		### refactoring seqs
		retain = []
		for i in range(len(self.seqs)):
			if self.seqs[i][0] in indel_index:
				continue
			else:
				retain.append((ids[i],self.seqs[i]))
		self.seqs = retain
		
		return indel_seqs
		
	def _is_indel_marked(self, id):
		if id[-5:] == "INDEL" :
			return True
		else:
			return False
			
	def _build_data_array(self):
		### Build data matrix
		
		data_array = np.zeros(shape=(self.rows,self.cols))
		for col_idx, index in enumerate(self.junction_index_storage):
			for i in range(self.rows):
				kmer_ = self.seqs[i][1][1][index - self.ksize + 1 : index + 1]
				#print i , kmer_
				depth_ = self.graph[kmer_]["depth"]
				data_array[i,col_idx] = np.float32(depth_)
		return data_array
		
	def _set_dimension(self):
#		print self.seqs
		ids_retain = map(lambda x : x[0] , self.seqs)
		seqs_retain = map(lambda x : x[1][1] , self.seqs)
		if len(seqs_retain) == 0 :
			self.p_val = 'Nan'
			return 0 , 0, []
		
			##Z = zip(seqs_retain)
		junction_index_storage = []
		for base_idx in xrange(len(seqs_retain[0])):
			tmp = []
			for seq in seqs_retain:
				tmp.append(seq[base_idx])
			c = Counter(tmp)
			if len(c.items()) == 1:
				continue
			else:
				junction_index_storage.append(base_idx)
		"""
		for i,z in enumerate(Z):
			c = Counter(z)
			
			nums = len(c.items())
			
			if nums == 1:
				continue
			else:
				junction_index_storage.append(i)
		"""		
		rows = len(seqs_retain)
		cols = len(junction_index_storage)
		
		return rows, cols , junction_index_storage
		
	def compute_path_scores(self) :
		### TODO Must be fixed
		if len(self.seqs) == 0 :
			self.path_scores = [-1]
			return
		if len(self.seqs) == 1 :
			self.path_scores = [1.0]
			return
		
		if self.junction_index_storage == []:
			return [-1]
		self.path_scores = []
		for seq in self.seqs:
			### We are computing path score here
			ratios = []
			ksize = self.ksize
			for junction in self.junction_index_storage:
				
				kmer = seq[1][1][junction-ksize:junction]
				my_depth = self.graph[kmer]["depth"]
				kmer_before = seq[1][1][junction-1-ksize:junction-1]
				edges = self.graph[kmer_before]["edges"]
				depths = 0
				for e in edges :
					e_depth = self.graph[e]["depth"]
					depths += e_depth
				ratio = (1.0 * my_depth) / depths ### Ratio for current node
				ratios.append(ratio)
			n = len(ratios)
			S = sum(ratios)
			ps = S / n
				
			self.path_scores.append(ps)
				
		# data_array = self._build_data_array()
		
		# columnwise_sum = np.sum(data_array, axis=0,dtype=np.float32)
		# path_scores = []
		# for i in range(self.rows):
			# row_data = data_array[i]
			# for col_idx in range(len(row_data)):
				# column_sum = columnwise_sum[col_idx]
				# row_data[col_idx] = row_data[col_idx] / column_sum
			# path_scores.append(np.mean(row_data))
			
		# self.path_scores = path_scores
		
		
	def chi_test(self):
		#data_array = self._build_data_array()
		
		if len(self.junction_index_storage) == 0 :
			return "Nan"
		data_array = np.zeros(shape=(2,len(self.junction_index_storage)))
		pos_set = self.major_minor_allele_dic.keys()
		if len(pos_set) == 0 :
			return "Nan"
		pos_set.sort()
		for i , k in enumerate(pos_set):
			numeric = self.major_minor_allele_dic[k]
			data_array[0,i] = numeric[0]
			data_array[1,i] = numeric[1]
		p_val = chi2_test(data_array)
		self.p_val = p_val[1]
		return p_val[1]
		
	def cv(self,seq):
		cv = 0.0
		depths = []
		
		for i in xrange(len(seq) - self.ksize):
			kmer = seq[i:i+self.ksize]
			depths.append(int(self.graph[kmer]["depth"]))

		avg = np.mean(depths)

		std = np.std(depths)

		cv = std / avg

		return avg , cv	
	#### Functions below do not use
	# def kmerize(self) :
		# for i in xrange(len(self.seq) - self.ksize - 1):
			# yield self.seq[i:i+self.ksize+1]


	# def path_score(self):
		# path_score = 0.0
		# n_junction = 0
		# e_freq_ratios = []
		# for k1mer in self.kmerize():
			# lkmer = k1mer[:-1]
			# rkmer = k1mer[1:]

			# connected_nodes = self.graph[lkmer]["edges"]
			# if len(connected_nodes) >= 2:
				# n_junction += 1
				# e_freqs = {}
				# for e in connected_nodes:
					# e_freqs[e] = self.graph[e]["depth"]
				# summation = sum(e_freqs.values())
				# for e in connected_nodes:
					
					# e_freqs[e] = e_freqs[e] / float(summation)
				# e_freq_ratios.append(e_freqs[rkmer])

		# if n_junction == 0 :
			# path_score = 1.0
		# else:
			# path_score = sum(e_freq_ratios) / n_junction

		# return path_score
