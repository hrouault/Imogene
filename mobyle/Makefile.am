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

xmldir = $(pkgdatadir)/xml

nodist_xml_DATA = imogene.xml imogene_env.xml

do_subst = sed -e 's,[@]PYTHON[@],$(PYTHON),g' -e 's,[@]pythondir[@],$(pythondir),g' -e 's,[@]bindir[@],$(bindir),g' -e 's,[@]PACKAGE_VERSION[@],$(PACKAGE_VERSION),g' -e 's,[@]PACKAGE_NAME[@],$(PACKAGE_NAME),g' -e 's,[@]PACKAGE[@],$(PACKAGE),g'

imogene.xml: imogene.xml.in
	$(do_subst) < $(srcdir)/imogene.xml.in > imogene.xml

imogene_env.xml: imogene_env.xml.in
	$(do_subst) < $(srcdir)/imogene_env.xml.in > imogene_env.xml

htmldir = $(pkgdatadir)/help

html_DATA = mobyle-help.html

EXTRA_DIST =  help.md

