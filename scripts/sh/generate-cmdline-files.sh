#!/bin/sh


gengetopt --input=genmot.ggo --func-name=genmot_cmdline_parser \
   --arg-struct-name=genmot_args_info --file-name=genmot_cmdline

gengetopt --input=extract.ggo --func-name=extract_cmdline_parser \
   --arg-struct-name=extract_args_info --file-name=extract_cmdline\

gengetopt --input=distinfo.ggo --func-name=distinfo_cmdline_parser \
   --arg-struct-name=distinfo_args_info --file-name=distinfo_cmdline\
