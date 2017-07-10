#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

infile = open("/Data2/SHN/BE/Plasmids/161025_Nextseq_2/LbNBE_GGS1/1.pileup")
posits = []
depths = []

for line in infile:
	tmp = line.strip().split("\t")
	posits.append(int(tmp[1]))
	depths.append(int(tmp[3]))
	
xmin = min(posits)
xmax = max(posits)
	
fig = plt.figure(figsize=(10,8.5))

subplt = fig.add_subplot(111)
subplt.set_ylabel('Coverage')
subplt.set_xlabel('Coordinates')

axes = plt.gca()
axes.set_xlim([xmin,xmax])

subplt.plot(posits, depths,lw=2)
subplt.axhline(y=np.mean(depths), ls = 'dashed', color = 'red', lw = 3)
subplt.fill_between(posits, depths, facecolor='blue', alpha=0.5)


for label in subplt.get_yticklabels():
    label.set_visible=(False)


fig.suptitle('Diagnostic plot for evaluating coverage of given reference\n(Dashed line indicates average coverage)')	
plt.savefig("Coverage_diagnostic.jpeg")
