/*
 * Copyright (C) 2006-2011 Hervé Rouault <rouault@lps.ens.fr>
 * Copyright (C) 2009-2011 Marc Santolini <santolin@lps.ens.fr>
 *
 * This file is part of Imogene.
 *
 * Imogene is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Imogene is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Imogene; see the file COPYING  If not, see <http://www.gnu.org/licenses/>.
 */


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <cstring>

using namespace std;

#include "config.h"

#include "const.hpp"
#include "sequence.hpp"
#include "motif.hpp"
#include "vectortypes.hpp"
#include "sequence.hpp"
#include "imogene.hpp"
#include "distinfo.hpp"
#include "display.hpp"
#include "extract.hpp"
#include "genmot.hpp"
#include "test.hpp"
#include "help.hpp"
#include "scangen.hpp"


const char usage_string[] =
    "imogene [--version] [--help]\n"
    "           <command> [<args>]";

const char more_info_string[] =
    "See 'imogene help <command>' for more information on a specific command.";


static int handle_options(char ** *argv, int * argc, int * envchanged)
{
    char ** orig_argv = *argv;
    while (*argc > 0) {
        char * cmd = (*argv)[0];
        if (cmd[0] != '-')
            break;
        if (!strcmp(cmd, "--help") || !strcmp(cmd, "--version")) {
            break;
        } else {
            fprintf(stderr, "Unknown option: %s\n", cmd);
            fprintf(stderr, "Usage : %s\n", usage_string);
        }
        (*argv)++;
        (*argc)--;
    }
    return (*argv) - orig_argv;
}

struct cmd_struct {
    const char * cmd;
    int (*fn)(int, char **);
    const char * help;
};

static int run_builtin(struct cmd_struct * p, int argc, char ** argv)
{
    int status, help;
    //	struct stat st;
    const char * prefix;
    prefix = NULL;
    help = argc == 2 && !strcmp(argv[1], "-h");
    status = p->fn(argc, argv);
    if (status)
        return status;
    return 0;
}

static struct cmd_struct commands[] = {
    { "distinfo", cmd_distinfo, "Computes the distance between two motifs." },
    { "display", cmd_display, "Displays motifs on sequences." },
    { "extract", cmd_extract, "Extracts alignments from coordinates." },
    { "genmot", cmd_genmot, "generate motifs" },
    { "test", cmd_test, "Run tests" },
    { "help", cmd_help, "Help message" },
    { "scangen", cmd_scangen, "infere CRMs" },
    { "version", cmd_version, "Print Imogene version" }
};


static void
handle_command(int argc, char ** argv)
{
    char * cmd = argv[0];
    for (unsigned int i = 0; i < ARRAY_SIZE(commands); i++) {
        struct cmd_struct * p = commands + i;
        if (strcmp(p->cmd, cmd))
            continue;
        exit(run_builtin(p, argc, argv));
    }
    cout <<  "Usage : " << usage_string << endl;
    list_common_cmds_help();
    cout << "\n" << more_info_string << endl;
}

static inline void mput_char(char c, unsigned int num)
{
    while (num--)
        putchar(c);
}


void list_cmds_help(void)
{
    unsigned int longest = 0;
    for (unsigned int i = 0; i < ARRAY_SIZE(commands); i++) {
        if (longest < strlen(commands[i].cmd))
            longest = strlen(commands[i].cmd);
    }
    puts("The available Imogene commands are:");
    for (unsigned int i = 0; i < ARRAY_SIZE(commands); i++) {
        printf("   %s   ", commands[i].cmd);
        mput_char(' ', longest - strlen(commands[i].cmd));
        puts(commands[i].help);
    }
}


void
print_reportbugs()
{
    cout << endl;
    cout << "Maintained by Hervé Rouault rouault@lps.ens.fr and" << endl;
    cout << "Marc Santolini santolin@lps.ens.fr" << endl;
    cout << "Report bugs to one of us." << endl;
}

void
print_copyright()
{
    cout << " Copyright (C) 2003-2011 Hervé Rouault" << endl;
    cout << endl;
    cout << " Imogene is free software: you can redistribute it and/or modify" << endl;
    cout << " it under the terms of the GNU General Public License as published by" << endl;
    cout << " the Free Software Foundation, either version 3 of the License, or" << endl;
    cout << " (at your option) any later version." << endl;
    cout << endl;
    cout << " Imogene is distributed in the hope that it will be useful," << endl;
    cout << " but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
    cout << " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
    cout << " GNU General Public License for more details." << endl;
    cout << endl;
    cout << " You should have received a copy of the GNU General Public License" << endl;
    cout << " along with Imogene.  If not, see <http://www.gnu.org/licenses/>." << endl;
    cout << endl;
    cout << " Written by Hervé Rouault and Marc Santolini." << endl;
}

#define is_dir_sep(c) ((c) == '/' || (c) == '\\')

char *
extract_argv0_cmd(char * argv0)
{
    char * slash;
    if (!argv0 || !*argv0)
        return NULL;
    slash = argv0 + strlen(argv0);
    while (argv0 <= slash && !is_dir_sep(*slash))
        slash--;
    if (slash >= argv0) {
        return slash + 1;
    }
    return argv0;
}

int prefixcmp(char * str, const char * prefix)
{
    for (; ; str++, prefix++)
        if (!*prefix)
            return 0;
        else if (*str != *prefix)
            return (unsigned char) * prefix - (unsigned char) * str;
}


/**
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  imogene-genmot main file, responsible for the de novo inference of
 *  motifs
 * =====================================================================================
 */
int
main(int argc, char ** argv)
{
    char * cmd;
    cmd = extract_argv0_cmd(argv[0]);
    if (!cmd)
        strncpy(cmd, "imogene-help", 15);
    if (!prefixcmp(cmd, "imogene-")) {
        cmd += 8;
        argv[0] = cmd;
        handle_command(argc, argv);
        return EXIT_FAILURE;
    }
    argv++;
    argc--;
    handle_options(&argv, &argc, NULL);
    if (argc > 0) {
        if (!strcmp(argv[0], "--help") || !strcmp(argv[0], "--version")) {
            argv[0] += 2;
        }
    } else {
        cout <<  "Usage : " << usage_string << endl;
        list_common_cmds_help();
        cout << "\n" << more_info_string << endl;
        return EXIT_FAILURE;
    }
    cmd = argv[0];
    if (!prefixcmp(cmd, "imogene-")) {
        cmd += 8;
        argv[0] = cmd;
        handle_command(argc, argv);
        return EXIT_FAILURE;
    }
    cmd = argv[0];
    handle_command(argc, argv);
    fprintf(stderr, "Failed to run command '%s': %s\n",
            cmd, strerror(errno));
    return EXIT_FAILURE;
}

