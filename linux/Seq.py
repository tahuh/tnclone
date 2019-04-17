#!/usr/bin/python

BASE_CONV = { x : y for x , y in zip('ATGCN' , 'TACGN')}
def reverse_complement(seq):
	ret = []
	for b in seq :
		ret.append(BASE_CONV[b])
	ret = ret[::-1]
	return ''.join(ret)

def generate_kmers(seq, k):
	ret = []
	for i in range(len(seq) - k + 1):
		ret.append(seq[i:i+k])
	return ret


