

python_PYTHON = maskalignedgenomes.py mapping.py

do_subst = sed -e 's,[@]pythondir[@],$(pythondir),g' -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

maskalignedgenomes.py: maskalignedgenomes.py.in
	$(do_subst) < $(srcdir)/maskalignedgenomes.py.in > maskalignedgenomes.py
	chmod +x maskalignedgenomes.py

do_substmap = sed -e 's,[@]pkgdatadir[@],$(pkgdatadir),g' 

mapping.py: mapping.py.in
	$(do_substmap) < $(srcdir)/mapping.py.in > mapping.py
	chmod +x mapping.py
