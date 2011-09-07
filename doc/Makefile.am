
BUILT_SOURCES = extract.help distinfo.help genmot.help

extract.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene extract --help > $@

distinfo.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene distinfo --help > $@

genmot.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene genmot --help > $@

info_TEXINFOS = imogene.texi
imogene_TEXINFOS = fdl.texi extract.texi extract.help distinfo.texi distinfo.help genmot.texi genmot.help
