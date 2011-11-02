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

BUILT_SOURCES = extract.help distinfo.help genmot.help display.help

extract.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene extract --help > $@

distinfo.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene distinfo --help > $@

genmot.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene genmot --help > $@

display.help: $(top_builddir)/src/imogene
	$(top_builddir)/src/imogene display --help > $@

#AM_MAKEINFOHTMLFLAGS = --document-language="en"

info_TEXINFOS = imogene.texi
imogene_TEXINFOS = fdl.texi extract.texi extract.help distinfo.texi distinfo.help genmot.texi genmot.help install.texi display.texi display.help test.texi
