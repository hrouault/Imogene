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
# along with Imogene; see the file COPYING  If not, see <http://www.gnu.org/licenses/>.
#
## Process this file with automake to produce Makefile.in

bin_PROGRAMS = imogene

BUILT_SOURCES = extract_cmdline.c distinfo_cmdline.c genmot_cmdline.c imogene_cmdline.c display_cmdline.c scangen_cmdline.c
   
imogene_SOURCES = imogene.hpp imogene.cpp extract.cpp extract.hpp genmot.cpp genmot.hpp help.cpp help.hpp scangen.cpp scangen.hpp sequence.cpp sequence.hpp motif.cpp motif.hpp tree.cpp tree.hpp random.cpp random.hpp vectortypes.cpp const.cpp const.hpp vectortypes.hpp distinfo.cpp distinfo.hpp display.cpp display.hpp imogene_cmdline.c distinfo_cmdline.c extract_cmdline.c display_cmdline.c genmot_cmdline.c scangen_cmdline.c
imogene_CPPFLAGS = -DDATA_PATH='"${pkgdatadir}"' -DPYTHON_PATH='"${pythondir}"'
	   

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

scangen_cmdline.c: scangen.ggo
	gengetopt --input=$(top_srcdir)/src/scangen.ggo --func-name=scangen_cmdline_parser \
   	--arg-struct-name=scangen_args_info --file-name=scangen_cmdline

