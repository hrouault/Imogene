#    
# Copyright (C) 2006-2011 Hervé Rouault <rouault@lps.ens.fr>
# Copyright (C) 2009-2011 Marc Santolini <santolin@lps.ens.fr>
# 
# This file is part of Imogene.
# 
# Imogene is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Imogene is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Imogene; see the file COPYING  If not,
# see <http://www.gnu.org/licenses/>.

BUILT_SOURCES = getalign.py mapping.py extractbackground.py weblogo-display.py\
				emf2fasta.pyx 

nodist_python_PYTHON = getalign.py mapping.py extractbackground.py\
					   weblogo-display.py emf2fasta.pyx 
CLEANFILES = $(python_PYTHON)
EXTRA_DIST = getalign.py.in mapping.py.in extractbackground.py.in\
			 weblogo-display.py.in emf2fasta.pyx.in 

do_subst = sed -e 's,[@]PYTHON[@],$(PYTHON),g'\
		       -e 's,[@]pythondir[@],$(pythondir),g'\
               -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

getalign.py: getalign.py.in
	$(do_subst) < $< > $@
	chmod +x $@

mapping.py: mapping.py.in
	$(do_subst) < $< > $@
	chmod +x $@

extractbackground.py: extractbackground.py.in
	$(do_subst) < $< > $@
	chmod +x $@

weblogo-display.py: weblogo-display.py.in
	$(do_subst) < $< > $@
	chmod +x $@

emf2fasta.pyx: emf2fasta.pyx.in
	$(do_subst) < $< > $@
	chmod +x $@
