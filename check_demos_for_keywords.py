#!/usr/bin/env python
# This script takes as its single argument a path to a directory (archetypally Rosetta/demos) and checks every .md file in that directory to see if it has keywords, defined as "a line starting with KEYWORDS: and containing at least two more things." Designed to fail noisily if not, as per a test.
# author Frank David Teets
# Written during XRW2016
import os
import sys
import shutil
demos = []
approved_keywords = []
approved_keywords.append("KEYWORDS:")
keywords = open(sys.argv[2], 'r')
for line in keywords:
	if not(line.startswith("#")):
		tokens = line.split()
		if(len(tokens) > 0):
			approved_keywords.append(tokens[0])
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
				for keyword in tokens:
					if keyword not in approved_keywords:
						print str(demoname) + " has unapproved keyword "+ keyword
						os._exit(1)
	if not(has_keywords):
		print str(demoname)  +" has no keywords!"
		os._exit(1)
