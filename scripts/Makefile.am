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

bin_SCRIPTS = getalign extractbackground weblogo-display
python_PYTHON = mapping.py

CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = getalign.py.in mapping.py.in extractbackground.py.in\
			 weblogo-display.py.in emf2fasta.pyx.in emf2fasta.pyx

do_subst = $(SED) -e 's,[@]PYTHON[@],$(PYTHON),g'\
		       -e 's,[@]pythondir[@],$(pythondir),g'\
               -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

getalign: getalign.py.in
	$(do_subst) < $< > $@
	chmod +x $@

mapping.py: mapping.py.in
	$(do_subst) < $< > $@

extractbackground: extractbackground.py.in
	$(do_subst) < $< > $@
	chmod +x $@

weblogo-display: weblogo-display.py.in
	$(do_subst) < $< > $@
	chmod +x $@

emf2fasta.pyx: emf2fasta.pyx.in
	$(do_subst) < $< > $@

emf2fasta.c: emf2fasta.pyx
	$(PYTHON) $(srcdir)/setup.py build_ext --inplace

pyexec_LTLIBRARIES = emf2fasta.la
emf2fasta_la_SOURCES = emf2fasta.c
emf2fasta_la_LDFLAGS = -avoid-version -module
emf2fasta_la_CPPFLAGS = -I$(PYTHON_INCLUDE_DIR)
