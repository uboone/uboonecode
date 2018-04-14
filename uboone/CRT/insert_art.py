#!/usr/bin/env python
import sys, getopt
import os

infile = sys.argv[1]
with open(infile,"r+") as f:
	lines = f.readlines()
	for i in range(len(lines)):
		if 'CRTRawtoCRTHit' in lines[i]:
			lines.insert(i+1,'"family": "art",')
	print "\n".join(lines)
