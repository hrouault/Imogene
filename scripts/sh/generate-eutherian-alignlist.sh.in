#!/bin/sh

find $PWD -name "*.fa" > foo

sed "s/.*chr//;s/\// /;s/-/ /;s/\.fa//" <foo > bar
paste bar foo| tr " " "\t" > files.dat
rm foo
rm bar
