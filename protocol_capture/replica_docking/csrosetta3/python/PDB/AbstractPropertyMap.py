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

"""Class that maps (chain_id, residue_id) to a residue property."""


class AbstractPropertyMap:
    def __init__(self, property_dict, property_keys, property_list):
        self.property_dict=property_dict
        self.property_keys=property_keys
        self.property_list=property_list

    def _translate_id(self, entity_id):
        return entity_id

    def __contains__(self, id):
        """True if the mapping has a property for this residue.

        Example:
            >>> if (chain_id, res_id) in apmap:
            ...     res, prop = apmap[(chain_id, res_id)]

        @param chain_id: chain id
        @type chain_id: char 

        @param res_id: residue id
        @type res_id: char 
        """
        translated_id = self._translate_id(id)
        return (translated_id in self.property_dict)

    def __getitem__(self, key):
        """
        Return property for a residue.

        @param chain_id: chain id
        @type chain_id: char 

        @param res_id: residue id
        @type res_id: int or (char, int, char) 

        @return: some residue property 
        @rtype: anything (can be a tuple)
        """
        translated_id=self._translate_id(key)
        return self.property_dict[translated_id]

    def __len__(self):
        """
        Return number of residues for which the property is available.

        @return: number of residues
        @rtype: int
        """
        return len(self.property_dict)

    def has_key(self, id):
        """True if the mapping has a property for this residue.

        (Obsolete; use "id in mapping" instead.)

        Example:

            >>> if apmap.has_key((chain_id, res_id)):
            ...     res, prop = apmap[(chain_id, res_id)]

        Is equivalent to:

            >>> if (chain_id, res_id) in apmap:
            ...     res, prop = apmap[(chain_id, res_id)]

        @param chain_id: chain id
        @type chain_id: char 

        @param res_id: residue id
        @type res_id: char 
        """
        import warnings
        warnings.warn("This function is obsolete; use 'id in mapping' instead", PendingDeprecationWarning)
        return (id in self)

    def keys(self):
        """
        Return the list of residues.

        @return: list of residues for which the property was calculated
        @rtype: [(chain_id, res_id), (chain_id, res_id),...] 
        """
        return self.property_keys

    def __iter__(self):
        """
        Iterate over the (entity, property) list. Handy alternative to 
        the dictionary-like access.

        Example:
            >>> for (res, property) in iter(map):
            ...     print res, property

        @return: iterator
        """
        for i in range(0, len(self.property_list)):
            yield self.property_list[i]


class AbstractResiduePropertyMap(AbstractPropertyMap):
    def __init__(self, property_dict, property_keys, property_list):
        AbstractPropertyMap.__init__(self, property_dict, property_keys, 
                property_list)

    def _translate_id(self, ent_id):
        chain_id, res_id=ent_id
        if isinstance(res_id, int):
            ent_id=(chain_id, (' ', res_id, ' '))
        return ent_id


class AbstractAtomPropertyMap(AbstractPropertyMap):
    def __init__(self, property_dict, property_keys, property_list):
        AbstractPropertyMap.__init__(self, property_dict, property_keys, 
                property_list)

    def _translate_id(self, ent_id):
        if len(ent_id)==4:
            chain_id, res_id, atom_name, icode=ent_id
        else:
            chain_id, res_id, atom_name=ent_id
            icode=None
        if isinstance(res_id, int):
            ent_id=(chain_id, (' ', res_id, ' '), atom_name, icode)
        return ent_id

