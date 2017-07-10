#!/usr/bin/python

faname = "./empty2.fa"

from SequenceParser import FAParser

parser = FAParser(faname)

parser.open()

cnt = 0

for id , desc , seq in parser.parse():
	cnt += 1

print cnt
