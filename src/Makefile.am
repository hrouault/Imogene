## Process this file with automake to produce Makefile.in

bin_PROGRAMS = imogene

BUILT_SOURCES = extract_cmdline.c distinfo_cmdline.c genmot_cmdline.c imogene_cmdline.c display_cmdline.c

imogene_SOURCES = imogene.hpp imogene.cpp extract.cpp extract.hpp genmot.cpp genmot.hpp help.cpp help.hpp scangen.cpp scangen.hpp sequence.cpp sequence.hpp motif.cpp motif.hpp tree.cpp tree.hpp random.cpp random.hpp vectortypes.cpp const.cpp const.hpp vectortypes.hpp distinfo.cpp distinfo.hpp display.cpp display.hpp imogene_cmdline.c distinfo_cmdline.c extract_cmdline.c display_cmdline.c genmot_cmdline.c
imogene_CPPFLAGS = -DDATA_PATH='"${pkgdatadir}"'

extract_cmdline.c: extract.ggo
	gengetopt --input=$(top_srcdir)/src/extract.ggo --func-name=extract_cmdline_parser \
   	--arg-struct-name=extract_args_info --file-name=extract_cmdline

distinfo_cmdline.c: distinfo.ggo
	gengetopt --input=$(top_srcdir)/src/distinfo.ggo --func-name=distinfo_cmdline_parser \
   	--arg-struct-name=distinfo_args_info --file-name=distinfo_cmdline

display_cmdline.c: display.ggo
	gengetopt --input=$(top_srcdir)/src/display.ggo --func-name=display_cmdline_parser \
   	--arg-struct-name=display_args_info --file-name=display_cmdline

imogene_cmdline.c: imogene.ggo
	gengetopt --input=$(top_srcdir)/src/imogene.ggo --func-name=imogene_cmdline_parser \
   	--arg-struct-name=imogene_args_info --file-name=imogene_cmdline

genmot_cmdline.c: genmot.ggo
	gengetopt --input=$(top_srcdir)/src/genmot.ggo --func-name=genmot_cmdline_parser \
   	--arg-struct-name=genmot_args_info --file-name=genmot_cmdline

