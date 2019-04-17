#!/usr/bin/python

class AssemblyEngine:
	def __init__(self, seed, terminal, graph):
		self.g = graph
		self.seed = seed
		self.terminal = terminal
	
	def Action(self):
		""" Performs DFS """
		if self.seed not in self.g.nodes():
			yield 'NO SEED'
		stack = [ ( self.seed , [ self.seed ] ) ]
		MAXITER = 100000
		it = 0
		while stack :
			vertex , path = stack.pop()
			try:
				for n in self.g[vertex]["edges"] - set(path):
					if n == self.terminal :
						yield path + [n]
					else:
						stack.append( ( n , path + [n] ) )
				it += 1
				if it > MAXITER:
					break
			except KeyError:
				yield 'BROKEN'
				break