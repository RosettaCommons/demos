###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-


import argparse
import textwrap
import library

class ExampleArgumentParser(argparse.ArgumentParser):
	def __init__(self, **kwargs):
		if 'examples' in kwargs:
			self.examples=kwargs['examples']
			del kwargs['examples']
		else:
			self.examples=None
#		print kwargs
		if 'aliases' in kwargs:
			self.aliases=kwargs['aliases']
			del kwargs['aliases']
		else:
			self.aliases=None

		argparse.ArgumentParser.__init__(self, **kwargs)

	def format_help(self):
#		return argparse.ArgumentParser.format_help(self)
		help_str=argparse.ArgumentParser.format_help(self)
		if self.examples:
			help_str=help_str+"\nexamples:"
			for ex in self.examples:
				if library.obj_is_list(ex):
					help_str=help_str+"\n\n  "+ex[0]
					help_str=help_str+"\n        "
					help_str=help_str+"\n        ".join(textwrap.wrap(ex[1],100))
				else:
					help_str=help_str+"\n  "+ex+""
			help_str=help_str+"\n\n"

		if self.aliases:
			self.aliases.remove(self.prog.replace('.py',''))
		if self.aliases and len(self.aliases)>0:
			if len(self.aliases)>1:
				s='s'
				es='es'
			else:
				s=''
				es=''
			help_str=help_str+("\nalias%s:\nThe same application can be called using the name%s: ")%(es,s)+" ".join(self.aliases)+"\n"

		return help_str%dict(prog=self.prog)



