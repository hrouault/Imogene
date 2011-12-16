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
# along with Imogene; see the file COPYING  If not, see
# <http://www.gnu.org/licenses/>.
#
## Process this file with automake to produce Makefile.in

dist_noinst_DATA = display.ggo distinfo.ggo extract.ggo genmot.ggo scangen.ggo

display_cmdline.c: display.ggo
	  gengetopt --input=$< --func-name=display_cmdline_parser\
   	--arg-struct-name=display_args_info --file-name=display_cmdline

distinfo_cmdline.c: distinfo.ggo
	  gengetopt --input=$< --func-name=distinfo_cmdline_parser\
   	--arg-struct-name=distinfo_args_info --file-name=distinfo_cmdline

extract_cmdline.c: extract.ggo
	  gengetopt --input=$< --func-name=extract_cmdline_parser\
   	--arg-struct-name=extract_args_info --file-name=extract_cmdline

genmot_cmdline.c: genmot.ggo
	  gengetopt --input=$< --func-name=genmot_cmdline_parser\
   	--arg-struct-name=genmot_args_info --file-name=genmot_cmdline

scangen_cmdline.c: scangen.ggo
	  gengetopt --input=$< --func-name=scangen_cmdline_parser\
   	--arg-struct-name=scangen_args_info --file-name=scangen_cmdline

bin_PROGRAMS = imogene

imogene_SOURCES = distinfo_cmdline.c distinfo_cmdline.h\
				  extract_cmdline.c extract_cmdline.h\
				  display_cmdline.c display_cmdline.h\
				  genmot_cmdline.c genmot_cmdline.h\
				  scangen_cmdline.c scangen_cmdline.h\
				  imogene.hpp imogene.cpp\
				  extract.cpp extract.hpp\
				  genmot.cpp genmot.hpp\
				  help.cpp help.hpp\
				  scangen.cpp scangen.hpp\
				  sequence.cpp sequence.hpp\
				  motif.cpp motif.hpp\
				  tree.cpp tree.hpp\
				  random.cpp random.hpp\
				  vectortypes.cpp vectortypes.hpp\
				  const.cpp const.hpp\
				  distinfo.cpp distinfo.hpp\
				  display.cpp display.hpp
imogene_CPPFLAGS = -DDATA_PATH='"${pkgdatadir}"' -DPYTHON_PATH='"${pythondir}"'
