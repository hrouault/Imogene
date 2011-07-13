
bin_SCRIPTS = generate-droso-alignlist.sh generate-eutherian-alignlist.sh
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = generate-droso-alignlist.sh.in

do_subst = sed -e 's,[@]pythondir[@],$(pythondir),g' 

generate-droso-alignlist.sh: generate-droso-alignlist.sh.in Makefile
	$(do_subst) < $(srcdir)/generate-droso-alignlist.sh.in > generate-droso-alignlist.sh
	chmod +x generate-droso-alignlist.sh
