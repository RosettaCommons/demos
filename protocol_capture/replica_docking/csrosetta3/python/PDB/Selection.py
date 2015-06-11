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

"""Selection of atoms, residues, etc."""

from PDB.Entity import Entity
from PDB.PDBExceptions import PDBException


entity_levels=["A", "R", "C", "M", "S"]


def uniqueify(items):
    """Return a list of the unique items in the given iterable.

    Order is NOT preserved.
    """
    return list(set(items))


def get_unique_parents(entity_list):
    """Translate a list of entities to a list of their (unique) parents.""" 
    parents = [entity.get_parent() for entity in entity_list]
    return uniqueify(parents)


def unfold_entities(entity_list, target_level):
    """
    Unfold a list of entities to a list of entities of another 
    level.  E.g.:

    list of atoms -> list of residues
    list of modules -> list of atoms
    list of residues -> list of chains

    o entity_list - list of entities or a single entity
    o target_level - char (A, R, C, M, S)
    """
    if not target_level in entity_levels:
        raise PDBException("%s: Not an entity level." % target_level)
    if isinstance(entity_list, Entity):
        # single entity
        entity_list=[entity_list]
    # level of entity list
    level=entity_list[0].get_level()
    for entity in entity_list:
        if not (entity.get_level()==level):
            raise PDBException("Entity list is not homogeneous.")
    target_index=entity_levels.index(target_level)
    level_index=entity_levels.index(level)
    if level_index==target_index:
        # already right level
        return entity_list
    if level_index>target_index:
        # we're going down, e.g. S->A
        for i in range(target_index, level_index):
            new_entity_list=[]
            for entity in entity_list:
                new_entity_list=new_entity_list+entity.get_list()
            entity_list=new_entity_list
    else:
        # we're going up, e.g. A->S
        for i in range(level_index, target_index):
            new_entity_list=[]  
            for entity in entity_list:
                parent=entity.get_parent()
                new_entity_list.append(parent)
            # find unique parents
            entity_list=uniqueify(new_entity_list)
    return entity_list

