#!/usr/bin/python

import argparse # unfortunately this is not supported prior to version 2.7

# create the top-level parser
parser = argparse.ArgumentParser(description='Generate a list of motif from an initial series of alignments',prog="@progname@")

subparsers = parser.add_subparsers(help='sub-command help')

# create the parser for the "genmot" command
parser_genmot = subparsers.add_parser('genmot', help='genmot help')
parser_genmot.add_argument('motif_width', type=int, help='motif width help')

# create the parser for the "scangen" command
parser_scangen = subparsers.add_parser('scangen', help='scangen help')
parser_scangen.add_argument('--baz', choices='XYZ', help='baz help')

## Version information
parser.add_argument('--version', action='version', version='%(prog)s @version@')

parser.parse_args()

## Test parser
##parser.parse_args(['--help'])


## Arguments pour genmot

## Print version information in a dedicated file
try:
   finfo=open("parameters_run.txt",'w')
   finfo.write("Essai")
except:
   print "The file parameters_run.txt already exists"
   print "Cannot record parameters. Exiting."
   exit()


## Print parameters for the execution of imogene-genmot in the same file

## Execute imogene-motgen

## Genereate pictures of the motifs

