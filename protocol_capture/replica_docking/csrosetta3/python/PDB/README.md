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
KEYWORDS: GENERAL GENERAL
This module is copied from biopython-1.57 and has been cleaned up in a hack and slash fashion to reduce external dependencies. This means that probably many parts of the API 
are now broken. We only use this for PDB import/export and could probably reduce the number of files here drastically or use a different PDB reader/write entirely.

If you want to use any functionality of biophython we suggest to install the full package and import via Bio.PDB. 


Here is the Biophython License Statement

                Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

Here is a list of Contributors to the origian BioPython package: 
CONTRIBUTORS
============

This is a list of people who have made contributions to Biopython.
This is certainly not comprehensive, and if you've been overlooked
(sorry!), please mention it on the development mailing list.
People are listed alphabetically by surname.

Cecilia Alsmark <Cecilia.Alsmark at domain ebc.uu.se>
Tiago Antao <tiagoantao at gmail.com>
Sebastian Bassi <sbassi at domain asalup.org>
Bill Barnard <bill at domain barnard-engineering.com>
Yves Bastide <ybastide at domain irisa.fr>
Paul T. Bathen
Yair Benita <Y.Benita at domain pharm.uu.nl>
Peter Bienstman <Peter.Bienstman at domain rug.ac.be>
Jose Blanca
Bob Bussell <rgb2003 at domain med.cornell.edu>
Diego Brouard <diego at domain conysis.com>
James Casbon <j.a.casbon at domain qmul.ac.uk>
Hye-Shik Chang <perky at domain fallin.lv>
Jeffrey Chang <jchang at domain smi.stanford.edu>
Brad Chapman <chapmanb at domain arches.uga.edu>
Peter Cock <p.j.a.cock at googlemail dot com>
Marc Colosimo <mcolosimo at domain mitre.org>
Andres Colubri <andres dot colubri at gmail dot com>
Cymon J Cox <cymon at domain duke.edu>
Gavin E Crooks <gec at domain compbio.berkeley.edu>
Andrew Dalke <dalke at domain acm.org>
Michiel de Hoon <mdehoon at domain c2b2.columbia.edu>
Bart de Koning <bratdaking gmail>
Sjoerd de Vries <sjoerd at domain nmr.chem.uu.nl>
Nathan J. Edwards <nje5 at edu domain georgetown>
Kyle Ellrott
Jeffrey Finkelstein <jeffrey.finkelstein at domain gmail.com>
Iddo Friedberg <idoerg at domain burnham.org>
Bertrand Frottier <bertrand.frottier at domain free.fr>
Phillip Garland <pgarland at gmail>
Walter Gillett < https://github.com/wgillett >
Frederik Gwinner
Jason A. Hackney <jhackney at domain stanford.edu>
Thomas Hamelryck <thamelry at domain binf.ku.dk>
Michael Hoffman <hoffman+biopython at domain ebi.ac.uk>
Thomas Holder
Yu Huang <krocea at domain yahoo.com.cn>

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
Kevin Jacobs <jacobs at bioinformed dot com>
Diana Jaunzeikare
Joanna & Dominik Kasprzak
Frank Kauff <fkauff at domain duke.edu>
Siong Kong
Andreas Kuntzagk <andreas.kuntzagk at domain mdc-berlin.de>
Michal Kurowski <michal at domain genesilico.pl>
Uri Laserson <laserson at Massachusetts Institute of Technology dot edu>
Chris Lasher <chris.lasher at gmail.com>
Gaetan Lehman <gaetan.lehmann at domain jouy.inra.fr>
Katharine Lindner <katel at domain worldpath.net>
Erick Matsen <surname at fhcrc dot org>
Tarjei Mikkelsen <tarjei at domain genome.wi.mit.edu>
Konstantin Okonechnikov <k.okonechnikov at domain gmail.com>
Cheng Soon Ong <chengsoon.ong at tuebingen.mpg.de>
Anne Pajon <ap one two at sanger ac uk>
Claude Paroz <claude at two (as digit) xlibre dot net>
Andrea Pierleoni <andrea at the Italian domain biocomp dot unibo>
Mike Poidinger <Michael.Poidinger at domain eBioinformatics.com>
Leighton Pritchard <lpritc at domain scri.sari.ac.uk>
Wolfgang Schueler <wolfgang at domain proceryon.at>
Peter Slickers <piet at domain clondiag.com>
Thomas Sicheritz-Ponten <thomas at domain cbs.dtu.dk>
Frederic Sohm <fsms at domain users.sourceforge.net>
Joao Rodrigues <anaryin at the domain gmail dot com>
Thomas Rosleff Soerensen <rosleff at domain mpiz-koeln.mpg.de>
Eric Talevich <eric.talevich at domain gmail.com>
Bartosz Telenczuk <bartosz.telenczuk at domain gmail.com>
Carlos Rios Vera <crosvera at domain gmail.com>
Johann Visagie <wjv at domain cityip.co.za>
Dan Vogel <dmv at domain andrew.cmu.edu>
David Weisman <david.weisman at domain acm.org>
Bartek Wilczynski <bartek at domain rezolwenta.eu.org>
David Winter <david dot winter at gmail dot com> 
Hongbo Zhu
Christian Zmasek
Harry Zuzan <iliketobicycle at domain yahoo.ca>


