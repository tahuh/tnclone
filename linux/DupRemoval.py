#!/usr/bin/python
import sys
import os
import re
from Bio import SeqIO
class DupRemovalTool:
	def __init__(self,In,Out,editor,mode=True):
		self.In = os.path.dirname(os.path.abspath(In)) + "/" + In.split("/")[-1]
		self.out = os.path.dirname(os.path.abspath(Out)) + "/" + Out.split("/")[-1]
		self.contigs = {}
		self.origin = 0
		self.editor = editor
		self.mode = mode
	def load(self) :
		#self.editor.append("[TnClone::DupRmv] Infile : %s\n"%(self.In))
		#sys.stderr.write("[SH DOWNSTREAM TOOLBOX::DupRmv] Infile : %s\n"%(self.In))
		with open(self.In) as In :
			for record in SeqIO.parse(In,"fasta") :
				self.origin += 1
				ID = record.id
				SEQ = str(record.seq)
				
				self.contigs[SEQ] = ID


	def run(self) :
		self.load()
		o = open(self.out , "w")
		#self.editor.append("[TnClone::DupRmv] Outfile : %s\n"%(self.out))
		#sys.stderr.write("[SH DOWNSTREAM TOOLBOX::DupRmv] Outfile : %s\n"%(self.out))
		
		#### Select max and mins sample
		key_set = self.contigs.keys()
		if len(key_set) == 1 :
			o.write(">" + self.contigs[key_set[0]] + "\n" + key_set[0] + "\n")
		elif len(key_set) == 0 :
			pass
		else:
			"""
			### DEPRECATED CODE after inspection of DATA
			### select max and min
			entries = list(self.contigs.iteritems()) ## [ ( seq , id ) , ... ]
			ids = []
			for e in entries :
				id = e[1]; ids.append(id)
			path_scores = []
			for id in ids :
				#print id
				if self.mode:
					if "INDEL" in id:
						if "path_score" in id:
							contig_number, sample_name, indel , depth, ps, cv = re.split("[:|]+" , id)[1::2]
							path_scores.append(float(ps))
						else:
							if "INDEL" in id:
								contig_number, sample_name = re.split("[:|]+" , id)[1::2]
							else:
								contig_number, sample_name = re.split("[:|]+" , id)[1::2]
					else:
						if "path_score" in id:
							contig_number, sample_name, depth, ps, cv = re.split("[:|]+" , id)[1::2]
							path_scores.append(float(ps))
						else:
							contig_number, sample_name = re.split("[:|]+" , id)[1::2]
					#contig_number, depth, path_score, cv , sample_name = re.split("[:|]+" , id)[1::2]
				else:
					break
			if len(path_scores) != 0 :
				min_path = min(path_scores); max_path = max(path_scores)
				min_index = path_scores.index(min_path)
				max_index = path_scores.index(max_path)
				selected_min = entries[min_index]
				selected_max = entries[max_index]
				### write MAX
				o.write(">" + selected_max[1] + "\n" + selected_max[0] + "\n")
				o.write(">" + selected_min[1] + "\n" + selected_min[0] + "\n")
			else:
				for seq , id in self.contigs.items():
					o.write(">"  + id + "\n" + str(seq) + "\n")
			
			"""		
		
			for k , v in self.contigs.iteritems() :
	#sys.stderr.write( k)
				o.write(">" + v + "\n" + k + "\n")
		
		#sys.stderr.write("[SH DOWNSTREAM TOOLBOX::DupRmv] %d / %d passed\n"%(self.origin,len(self.contigs.keys())))
		
		o.close()
		# if len(key_set) > 1 :
			# self.editor.append("[TnClone::DupRmv] %d / %d passed.\n"%(self.origin,len(self.contigs.keys())))
		# else:
			# self.editor.append("[TnClone::DupRmv] Single clone have detected. 1 contig passed.\n")
