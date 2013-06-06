#!@PYTHON@
# -*- coding: utf-8 -*-
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

import argparse
import submodule

# Command line parser
parser = argparse.ArgumentParser(description='Genmot wrapper for galaxy.')

parser.add_argument('-w', '--width', help='Width of the motifs')

parser.add_argument('-t',
                    '--threshold',
                    help='Threshold used for motif scanning')

parser.add_argument('-u',
                    '--threshold_disp',
                    help='Threshold used for motif display')

parser.add_argument('-x',
                    '--neighbext',
                    help='Extent of the motif search within an alignment')

parser.add_argument('-e',
                    '--evolutionary-model',
                    help='Evolutionary model used for motif generation' \
                         '(1=felsen, 2=halpern)')

parser.add_argument('-s', '--species', help='Species selected') 


args = parser.parse_args()

submodule.call



