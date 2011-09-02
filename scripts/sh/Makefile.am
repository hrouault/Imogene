
bin_SCRIPTS = generate-droso-alignlist.sh generate-eutherian-alignlist.sh
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = generate-droso-alignlist.sh.in generate-eutherian-alignlist.sh.in

do_subst = sed -e 's,[@]pythondir[@],$(pythondir),g' -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

generate-droso-alignlist.sh: generate-droso-alignlist.sh.in Makefile
	$(do_subst) < $(srcdir)/generate-droso-alignlist.sh.in > generate-droso-alignlist.sh
	chmod +x generate-droso-alignlist.sh

generate-eutherian-alignlist.sh: generate-eutherian-alignlist.sh.in Makefile
	$(do_subst) < $(srcdir)/generate-eutherian-alignlist.sh.in > generate-eutherian-alignlist.sh
	chmod +x generate-eutherian-alignlist.sh
