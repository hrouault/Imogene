## Process this file with automake to produce Makefile.in

noinst_LIBRARIES = libdistinfo.a
libdistinfo_a_SOURCES = distinfo.cpp distinfo.hpp

bin_PROGRAMS = imogene_genmot ##distinfo
##bin_PROGRAMS = scangen distinfo 
## 
imogene_genmot_SOURCES = imogene-genmot.hpp imogene-genmot.cpp sequence.cpp sequence.hpp motif.cpp motif.hpp tree.cpp tree.hpp random.cpp random.hpp vectortypes.cpp const.cpp const.hpp vectortypes.hpp imogene-genmot_cmdline.c imogene-genmot_cmdline.h
imogene_genmot_CPPFLAGS = -DDATA_PATH=${pkgdatadir}
imogene_genmot_LDADD = libdistinfo.a

## distinfo_SOURCES =  distinfo.cpp distinfo.hpp motif.cpp motif.hpp sequence.cpp sequence.hpp tree.cpp tree.hpp const.cpp const.hpp vectortypes.cpp vectortypes.hpp

## scangen_SOURCES = scangen.cpp scangen.hpp motif.cpp motif.hpp sequence.cpp sequence.hpp tree.cpp tree.hpp random.cpp random.hpp vectortypes.cpp vectortypes.hpp const.hpp const.cpp cmdline.c cmdline.h
## ##scangen_CPPFLAGS = -O3
## 
## distinfo_SOURCES = distinfo.cpp distinfo.hpp motif.cpp motif.hpp sequence.cpp sequence.hpp tree.cpp tree.hpp random.cpp random.hpp vectortypes.cpp vectortypes.hpp const.hpp const.cpp distinfo_cmdline.c distinfo_cmdline.h
## ##distinfo_CPPFLAGS = -O3 
