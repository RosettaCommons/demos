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
# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#                 Biopython License Agreement
#
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.
#
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE.

"""Turn an mmCIF file into a dictionary."""

import os.path
import warnings
import PDB.mmCIF.MMCIFlex
from UserDict import UserDict


class MMCIF2Dict(UserDict):
    # The token identifiers
    NAME=1
    LOOP=2
    DATA=3
    SEMICOLONS=4    
    DOUBLEQUOTED=5
    QUOTED=6
    SIMPLE=7

    def __init__(self, filename):
        # this dict will contain the name/data pairs 
        self.data={}
        # entry for garbage
        self.data[None]=[]
        if not os.path.isfile(filename):
            raise IOError("File not found.")
        PDB.mmCIF.MMCIFlex.open_file(filename)
        self._make_mmcif_dict()
        PDB.mmCIF.MMCIFlex.close_file()

    def _make_mmcif_dict(self): 
        # local copies
        NAME=self.NAME
        LOOP=self.LOOP
        DATA=self.DATA
        SEMICOLONS=self.SEMICOLONS
        DOUBLEQUOTED=self.DOUBLEQUOTED
        QUOTED=self.QUOTED
        SIMPLE=self.SIMPLE
        get_token=PDB.mmCIF.MMCIFlex.get_token
        # are we looping?
        loop_flag=0
        # list of names in loop
        temp_list=[]
        # last encountered name
        current_name=None
        # get first token/value pair
        token, value=get_token()
        # print token, value
        mmcif_dict=self.data
        # loop until EOF (token==0)
        while token:
            if token==NAME:
                if loop_flag:
                    # Make lists for all the names in the loop
                    while token==NAME:
                        # create  a list for each name encountered in loop
                        new_list=mmcif_dict[value]=[]
                        temp_list.append(new_list)
                        token, value=get_token()  
                        # print token, value
                    loop_flag=0         
                    # nr of data items parsed
                    data_counter=0
                    # corresponding data name
                    pos=0
                    nr_fields=len(temp_list)
                    # Now fill all lists with the data
                    while token>3:
                        pos=data_counter%nr_fields
                        data_counter=data_counter+1
                        temp_list[pos].append(value)
                        token, value=get_token()  
                        # print token, value
                    if pos!=nr_fields-1:
                        warnings.warn("ERROR: broken name-data pair "
                                      "(data missing)!", RuntimeWarning)
                    # The last token was not used, so
                    # don't set token to None! (this means the 
                    # last parsed token goes through the loop again)
                else:   
                    # simple name-data pair (no loop)
                    # so next token should be the data
                    next_token, data=get_token()  
                    # print token, value
                    mmcif_dict[value]=data
                    if next_token<4:
                        warnings.warn("ERROR: broken name-data pair "
                                      "(name-non data pair)!", RuntimeWarning)
                        # print token, value
                    else:   
                        # get next token
                        token=None
            elif token==LOOP:
                loop_flag=1
                temp_list=[]
                # get next token
                token=None
            elif token==DATA:
                mmcif_dict[value[0:5]]=value[5:]
                token=None
            else:
                # we found some complete garbage
                warnings.warn("ERROR: broken name-data pair "
                              "(missing name)!\n%s %s" % (token, value),
                              RuntimeWarning)
                mmcif_dict[None].append(value)
                # get next token
                token=None
            if token==None:
                token, value=get_token()
                # print token, value

    def __getitem__(self, key):
        return self.data[key]


if __name__=="__main__":

    import sys

    if len(sys.argv)!=2:
        print "Usage: python MMCIF2Dict filename."

    filename=sys.argv[1]    

    mmcif_dict=MMCIF2Dict(filename)

    entry = ""
    print "Now type a key ('q' to end, 'k' for a list of all keys):"
    while(entry != "q"):
        entry = raw_input("MMCIF dictionary key ==> ")    
        if entry == "q":
            sys.exit()
        if entry == "k":
            for key in mmcif_dict:
                print key
            continue
        try:
            value=mmcif_dict[entry]
            if type(value)==type([]):
                for item in value:
                    print item
            else:
                print value
        except KeyError:
            print "No such key found."

