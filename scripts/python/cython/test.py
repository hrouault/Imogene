#!/usr/bin/env python
#try:
#   import Cython
#   print "Using Cython for emf2fasta conversion..."
#   print
##   import os
##   os.system("/usr/bin/env python setup.py build_ext -inplace")
#   import pyximport
#   pyximport.install()
#   import emf2fasta
#except:
#   print "Using Python for emf2fasta conversion..."
#   print
#   execfile("emf2fasta.pyx")



import time
import os

start=time.time()

print "Using C++ for emf2fasta conversion..."
print

os.system("./emf2fasta")

print
print "it took", time.time() - start, "seconds."
print

start=time.time()

import Cython
print "Using Cython for emf2fasta conversion..."
print
import pyximport
pyximport.install()
import emf2fasta

print
print "it took", time.time() - start, "seconds."
print

start=time.time()

print "Using Python for emf2fasta conversion..."
print
execfile("emf2fasta.pyx")

print
print "it took", time.time() - start, "seconds."
print
