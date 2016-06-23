#!/usr/bin/env python
# This script takes as its single argument a path to a directory (archetypally Rosetta/demos) and checks every .md file in that directory to see if it has keywords, defined as "a line starting with KEYWORDS: and containing at least two more things." Designed to fail noisily if not, as per a test.
# author Frank David Teets
# Written during XRW2016
import os
import sys
import shutil
demos = []
for root, directories, filenames in os.walk(sys.argv[1]):
	for filename in filenames:
		if filename.endswith(".md") and not(filename.endswith("Home.md")):
			demos.append(os.path.join(root, filename))
for demoname in demos:
	has_keywords = False
	current_demo = open(demoname, 'r')
	for line in current_demo:
		if ("KEYWORDS:" in line):
			tokens = line.split()
			if(len(tokens)>2):
				has_keywords = True
	if not(has_keywords):
		print str(demoname)  +" has no keywords!"
		os._exit(1)
