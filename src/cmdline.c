/*
  File autogenerated by gengetopt version 2.22.4
  generated with the following command:
  gengetopt --input=extract.ggo --func-name=extract_cmdline_parser 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "Alignment extraction.";

const char *gengetopt_args_info_usage = "Usage: " EXTRACT_CMDLINE_PARSER_PACKAGE " [OPTIONS]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help         Print help and exit",
  "  -V, --version      Print version and exit",
  "  -i, --input=str    Coordinate file",
  "  -o, --output=str   Output folder with fasta formatted alignments  \n                       (default=`align')",
  "  -s, --species=str  Species selected  (possible values=\"droso\", \n                       \"eutherian\")",
    0
};

typedef enum {ARG_NO
  , ARG_STRING
} extract_cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
extract_cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct extract_cmdline_parser_params *params, const char *additional_error);

static int
extract_cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);

const char *extract_cmdline_parser_species_values[] = {"droso", "eutherian", 0}; /*< Possible values for species. */

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->input_given = 0 ;
  args_info->output_given = 0 ;
  args_info->species_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->input_arg = NULL;
  args_info->input_orig = NULL;
  args_info->output_arg = gengetopt_strdup ("align");
  args_info->output_orig = NULL;
  args_info->species_arg = NULL;
  args_info->species_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->input_help = gengetopt_args_info_help[2] ;
  args_info->output_help = gengetopt_args_info_help[3] ;
  args_info->species_help = gengetopt_args_info_help[4] ;
  
}

void
extract_cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(EXTRACT_CMDLINE_PARSER_PACKAGE_NAME) ? EXTRACT_CMDLINE_PARSER_PACKAGE_NAME : EXTRACT_CMDLINE_PARSER_PACKAGE),
     EXTRACT_CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  extract_cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
}

void
extract_cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
extract_cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);
}

void
extract_cmdline_parser_params_init(struct extract_cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct extract_cmdline_parser_params *
extract_cmdline_parser_params_create(void)
{
  struct extract_cmdline_parser_params *params = 
    (struct extract_cmdline_parser_params *)malloc(sizeof(struct extract_cmdline_parser_params));
  extract_cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
extract_cmdline_parser_release (struct gengetopt_args_info *args_info)
{

  free_string_field (&(args_info->input_arg));
  free_string_field (&(args_info->input_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->species_arg));
  free_string_field (&(args_info->species_orig));
  
  

  clear_given (args_info);
}

/**
 * @param val the value to check
 * @param values the possible values
 * @return the index of the matched value:
 * -1 if no value matched,
 * -2 if more than one value has matched
 */
static int
check_possible_values(const char *val, const char *values[])
{
  int i, found, last;
  size_t len;

  if (!val)   /* otherwise strlen() crashes below */
    return -1; /* -1 means no argument for the option */

  found = last = 0;

  for (i = 0, len = strlen(val); values[i]; ++i)
    {
      if (strncmp(val, values[i], len) == 0)
        {
          ++found;
          last = i;
          if (strlen(values[i]) == len)
            return i; /* exact macth no need to check more */
        }
    }

  if (found == 1) /* one match: OK */
    return last;

  return (found ? -2 : -1); /* return many values or none matched */
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  int found = -1;
  if (arg) {
    if (values) {
      found = check_possible_values(arg, values);      
    }
    if (found >= 0)
      fprintf(outfile, "%s=\"%s\" # %s\n", opt, arg, values[found]);
    else
      fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
extract_cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", EXTRACT_CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->input_given)
    write_into_file(outfile, "input", args_info->input_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->species_given)
    write_into_file(outfile, "species", args_info->species_orig, extract_cmdline_parser_species_values);
  

  i = EXIT_SUCCESS;
  return i;
}

int
extract_cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", EXTRACT_CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = extract_cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
extract_cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  extract_cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
extract_cmdline_parser (int argc, char **argv, struct gengetopt_args_info *args_info)
{
  return extract_cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
extract_cmdline_parser_ext (int argc, char **argv, struct gengetopt_args_info *args_info,
                   struct extract_cmdline_parser_params *params)
{
  int result;
  result = extract_cmdline_parser_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      extract_cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
extract_cmdline_parser2 (int argc, char **argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct extract_cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = extract_cmdline_parser_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      extract_cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
extract_cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  int result = EXIT_SUCCESS;

  if (extract_cmdline_parser_required2(args_info, prog_name, 0) > 0)
    result = EXIT_FAILURE;

  if (result == EXIT_FAILURE)
    {
      extract_cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
extract_cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error)
{
  int error = 0;
  FIX_UNUSED (additional_error);

  /* checks for required options */
  if (! args_info->input_given)
    {
      fprintf (stderr, "%s: '--input' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  if (! args_info->species_given)
    {
      fprintf (stderr, "%s: '--species' ('-s') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error = 1;
    }
  
  
  /* checks for dependences among options */

  return error;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see extract_cmdline_parser_params.check_ambiguity
 * @param override @see extract_cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               extract_cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  if (possible_values && (found = check_possible_values((value ? value : default_value), possible_values)) < 0)
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: %s argument, \"%s\", for option `--%s' (`-%c')%s\n", 
          package_name, (found == -2) ? "ambiguous" : "invalid", value, long_opt, short_opt,
          (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: %s argument, \"%s\", for option `--%s'%s\n", 
          package_name, (found == -2) ? "ambiguous" : "invalid", value, long_opt,
          (additional_error ? additional_error : ""));
      return 1; /* failure */
    }
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };


  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
extract_cmdline_parser_internal (
  int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct extract_cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    extract_cmdline_parser_init (args_info);

  extract_cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "input",	1, NULL, 'i' },
        { "output",	1, NULL, 'o' },
        { "species",	1, NULL, 's' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVi:o:s:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          extract_cmdline_parser_print_help ();
          extract_cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          extract_cmdline_parser_print_version ();
          extract_cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'i':	/* Coordinate file.  */
        
        
          if (update_arg( (void *)&(args_info->input_arg), 
               &(args_info->input_orig), &(args_info->input_given),
              &(local_args_info.input_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input", 'i',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output folder with fasta formatted alignments.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, "align", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 's':	/* Species selected.  */
        
        
          if (update_arg( (void *)&(args_info->species_arg), 
               &(args_info->species_orig), &(args_info->species_given),
              &(local_args_info.species_given), optarg, extract_cmdline_parser_species_values, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "species", 's',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", EXTRACT_CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */



  if (check_required)
    {
      error += extract_cmdline_parser_required2 (args_info, argv[0], additional_error);
    }

  extract_cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  return 0;

failure:
  
  extract_cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
