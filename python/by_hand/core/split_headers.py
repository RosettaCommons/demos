# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
import sys,os,re
t = open("/Users/sheffler/svn/branches/mini/demo/python/by_hand/core/core.hh").xreadlines()
r = re.compile("#include .core/(.*?)/(.*\.hh)")
h = dict()
all_h = list()
for l in t:
  m = r.match(l)
  if m:
    ns = m.group(1)
    if not ns in h: h[ns] = [ '#ifndef core_%s\n#define core_%s\n'%(ns,ns) ]
    h[ns].append( '#include "core/%s/%s"'%(ns,m.group(2) ) )
  elif l.strip() and not l.count("coreall") and not l=="#endif" and not l.startswith('//'): all_h.append(l)

for k,v in h.items():
  f = open("/Users/sheffler/svn/branches/mini/demo/python/by_hand/core/core_"+k+'.hh','w')
  print >>f, v[0]
  for l in all_h: print >>f, l
  for l in v[1:]: print >>f, l
  print >>f, '\n#endif'
  f.close()
