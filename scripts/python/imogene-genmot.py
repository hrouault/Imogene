#!/usr/bin/python

import os
import argparse # unfortunately this is not supported prior to version 2.7

# create the top-level parser
parser = argparse.ArgumentParser(description='Genome analysis for the inference of gene enhancers.',prog="@progname@")

subparsers = parser.add_subparsers(help='sub-command help')

# create the parser for the "genmot" command
parser_genmot = subparsers.add_parser('genmot', help='genmot help')
parser_genmot.add_argument('motif_width', type=int, help='motif width help')

# create the parser for the "scangen" command
parser_scangen = subparsers.add_parser('scangen', help='scangen help')
parser_scangen.add_argument('--baz', choices='XYZ', help='baz help')

## Version information
parser.add_argument('--version', action='version', version='%(prog)s @version@')

## Width of the generated motifs
parser.add_argument('--width', type=int,help='Width of the motifs')

args=parser.parse_args()



## Hard coded variables
extraction_cutoff=2*args.motif_width

conca_droso=0.3 # A concentration
conca_eutherian=0.263 # A concentration

## Arguments pour genmot

## Print version information in a dedicated file
try:
   finfo=os.open("parameters_run.txt",os.O_CREAT | os.O_EXCL | os.O_WRONLY)
except:
   print "The file parameters_run.txt already exists"
   print "Cannot record parameters. Exiting."
   exit()

os.write(finfo,'Extraction_cuttoff '+str(extraction_cutoff))

os.close(finfo)


## Print parameters for the execution of imogene-genmot in the same file

## Execute imogene-motgen

## Genereate pictures of the motifs

def coordtoseq(coordinates):
   for i in coordinates:

