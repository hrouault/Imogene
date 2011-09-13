
BUILT_SOURCES = extract.help distinfo.help genmot.help display.help

extract.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene extract --help > $@

distinfo.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene distinfo --help > $@

genmot.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene genmot --help > $@

display.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene display --help > $@

info_TEXINFOS = imogene.texi
imogene_TEXINFOS = fdl.texi extract.texi extract.help distinfo.texi distinfo.help genmot.texi genmot.help install.texi display.texi display.help test.texi
