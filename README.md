Overview
========

**Imogene** works in two different modes. In the *genmot* mode, motifs are predicted on the basis of a set of functionally related enhancers. In the *scangen* mode, enhancers are predicted from a set of motifs.

It is programmed in C++, accompanied by a bunch of python help script. In addition, an xml description page is provided for the mobyle system.


----------

Installation
============

You can refer to the INSTALL file for detailled installation instructions. 

Quick install
-------------

"""bash
mkdir build
cd build
../configure
make
make install
"""

Once the package is installed, the first thing to do is to download the required genomes. To help in that task, scripts have been added in the script folder.
execute `$prefix/bin/generate-droso-align.py` (for drosophila genomes)

Files will be automatically installed in `$(prefix)/share/imogene`.


Third-party software
====================

This package includes third party software in the folder mobyle under the GNU LGPLv2.1 license. Please refer to the README file in that folder for more details.