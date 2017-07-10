#!/usr/bin/python

class PvalReporter(object):
	def __init__(self, pval_info_as_dict, ohandle):
		### ohandle must be opened
		self.d = pval_info_as_dict
		self.F = ohandle

	def flush(self):
		gene_set = self.d.keys()
		gene_set.sort()
		self.F.write("[P-Value Stat]\n")
		for g in gene_set:
			number_set = map(lambda x : int(x) , self.d[g].keys())
			number_set.sort()
			self.F.write(g)
			for num in number_set:
				s_num = str(num)
				self.F.write("\t" + str(self.d[g][s_num][0]))
			self.F.write("\n")
