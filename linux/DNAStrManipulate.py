#!/usr/bin/python
import sys
def convert2barr(kmer) :
	A = 0x0
	T = 0x3
	G = 0x1
	C = 0x2
	byte = len(kmer) / 4 + (len(kmer) % 4 != 0)
	insert = [0] * byte
	suc = True
	for i in range(len(kmer)) :
		base = kmer[i]
		shift = 6 - 2 * ( i % 4 )
		if base == 'A' :
			insert[i/4] |= A << shift
		elif base == 'T' :
			insert[i/4] |= T << shift
		elif base == 'C' :
			insert[i/4] |= C << shift
		elif base == 'G':
			insert[i/4] |= G << shift
		else:
			#sys.stderr.write("Invalid base to convert : %c [ NONE of A,T,G,C ]\n Perhaps base 'N,R,Y?'"%(base))
			suc = False
			break
		if shift == 0 :
			shift = 6
		else:
			shift = shift - 2
	if suc :
		return tuple(insert)
	else:
		return None

def barr2seq(barr,length):
	A = 0x0
	T = 0x3
	G = 0x1
	C = 0x2
	seq = []
	for i in range(length) :
		shift = 6 - 2 * (i%4)
		mask = MASK << shift
		base = (barr[i/4]&mask) >> shift
		if base ==A :
			seq.append('A')
		elif base == 'T':
			seq.append('T')
		elif base == 'G':
			seq.append('G')
		elif base == 'C':
			seq.append('C')

	return ''.join(seq)
# def barr2seq(barr, length) :
	# seq = []
	# MASK = 0x3
	# for i in range(length) :
		# shift = 6 - 2 * (i%4)
		# mask = MASK << shift
		
		# base = (barr[i/4] & mask) >> shift
		
		# if base == A:
			# seq.append('A')
			
		# elif base == T:
			# seq.append('T')
			
		# elif base == G:
			# seq.append('G')
			
		# elif base == C:
			# seq.append('C')
			
	# return "".join(seq)


def gen_rnd_dna(len) :
	bases = ['A','T','G','C']
	seq = []
	for i in range(len) :
		seq.append(random.choice(bases))
	return ''.join(seq)