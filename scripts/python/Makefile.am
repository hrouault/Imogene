
python_PYTHON = getalign-droso.py getalign-eutherian.py mapping.py emf2fasta.pyx
CLEANFILES = $(python_PYTHON)
EXTRA_DIST = getalign-droso.py.in getalign-eutherian.py.in mapping.py.in emf2fasta.pyx.in

do_subst = sed -e 's,[@]pythondir[@],$(pythondir),g' -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

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
