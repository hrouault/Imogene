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
# along with Imogene; see the file COPYING  If not, see <http://www.gnu.org/licenses/>.

nodist_python_PYTHON = getalign-droso.py getalign-eutherian.py mapping.py emf2fasta.pyx extractbackground-droso.py extractbackground-eutherian.py weblogo-display.py
CLEANFILES = $(python_PYTHON)
EXTRA_DIST = getalign-droso.py.in getalign-eutherian.py.in mapping.py.in emf2fasta.pyx.in extractbackground-droso.py.in extractbackground-eutherian.py.in weblogo-display.py.in

do_subst = sed -e 's,[@]PYTHON[@],$(PYTHON),g' -e 's,[@]pythondir[@],$(pythondir),g' -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

getalign-droso.py: getalign-droso.py.in
	$(do_subst) < $(srcdir)/getalign-droso.py.in > getalign-droso.py
	chmod +x getalign-droso.py

getalign-eutherian.py: getalign-eutherian.py.in
	$(do_subst) < $(srcdir)/getalign-eutherian.py.in > getalign-eutherian.py
	chmod +x getalign-eutherian.py

emf2fasta.pyx: emf2fasta.pyx.in
	$(do_subst) < $(srcdir)/emf2fasta.pyx.in > emf2fasta.pyx

do_substmap = sed -e 's,[@]pkgdatadir[@],$(pkgdatadir),g' 

mapping.py: mapping.py.in
	$(do_substmap) < $(srcdir)/mapping.py.in > mapping.py
	chmod +x mapping.py

weblogo-display.py: weblogo-display.py.in
	$(do_substmap) < $(srcdir)/weblogo-display.py.in > weblogo-display.py
	chmod +x weblogo-display.py

extractbackground-droso.py: extractbackground-droso.py.in
	$(do_subst) < $(srcdir)/extractbackground-droso.py.in > extractbackground-droso.py
	chmod +x extractbackground-droso.py

extractbackground-eutherian.py: extractbackground-eutherian.py.in
	$(do_subst) < $(srcdir)/extractbackground-eutherian.py.in > extractbackground-eutherian.py
	chmod +x extractbackground-eutherian.py
