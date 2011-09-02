## Process this file with automake to produce Makefile.in

#noinst_LIBRARIES = libdistinfo.a
#libdistinfo_a_SOURCES = distinfo.cpp distinfo.hpp distinfo_cmdline.c distinfo_cmdline.h

bin_PROGRAMS = imogene ##distinfo
##bin_PROGRAMS = scangen distinfo 
## 
imogene_SOURCES = imogene.hpp imogene.cpp extract.cpp extract.hpp genmot.cpp genmot.hpp help.cpp help.hpp scangen.cpp scangen.hpp sequence.cpp sequence.hpp motif.cpp motif.hpp tree.cpp tree.hpp random.cpp random.hpp vectortypes.cpp const.cpp const.hpp vectortypes.hpp extract_cmdline.c extract_cmdline.h genmot_cmdline.c genmot_cmdline.h distinfo.cpp distinfo.hpp distinfo_cmdline.c distinfo_cmdline.h
imogene_CPPFLAGS = -DDATA_PATH=${pkgdatadir}
#imogene_LDADD = libdistinfo.a

## distinfo_SOURCES =  distinfo.cpp distinfo.hpp motif.cpp motif.hpp sequence.cpp sequence.hpp tree.cpp tree.hpp const.cpp const.hpp vectortypes.cpp vectortypes.hpp

## scangen_SOURCES = scangen.cpp scangen.hpp motif.cpp motif.hpp sequence.cpp sequence.hpp tree.cpp tree.hpp random.cpp random.hpp vectortypes.cpp vectortypes.hpp const.hpp const.cpp cmdline.c cmdline.h
## ##scangen_CPPFLAGS = -O3
## 
#distinfo_SOURCES = distinfo.cpp distinfo.hpp motif.cpp motif.hpp vectortypes.cpp vectortypes.hpp const.hpp const.cpp distinfo_cmdline.c distinfo_cmdline.h
#distinfo_CPPFLAGS = -DDATA_PATH=${pkgdatadir}
#distinfo_LDADD = libdistinfo.a
## ##distinfo_CPPFLAGS = -O3 
