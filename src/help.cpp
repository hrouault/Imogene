/*
 * =====================================================================================
 *
 *       Filename:  help.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06.08.2011 14:18:23
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <iostream>

#include "config.h"

using namespace std;

struct cmdname_help {
   char name[16];
   char help [80];
};

static struct cmdname_help common_cmds[] = {
   {"extract", "Extract an alignment from a coordinate file"},
   {"distinfo", "Distance between PWMs"},
   {"genmot", "Generating motifs de novo"}
};

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

static inline void mput_char(char c, unsigned int num)
{
	while(num--)
		putchar(c);
}

void
list_common_cmds_help(void)
{
   int i, longest = 0;

	for (i = 0; i < ARRAY_SIZE(common_cmds); i++) {
		if (longest < strlen(common_cmds[i].name))
			longest = strlen(common_cmds[i].name);
	}

   cout << "\n" << "The most commonly used imogene commands are:" << endl;
	for (i = 0; i < ARRAY_SIZE(common_cmds); i++) {
		printf("   %s   ", common_cmds[i].name);
		mput_char(' ', longest - strlen(common_cmds[i].name));
		puts(common_cmds[i].help);
	}
}
//void
//list_common_cmds_help(void)
//{
//	int i, longest = 0;
//
//	for (i = 0; i < ARRAY_SIZE(common_cmds); i++) {
//		if (longest < strlen(common_cmds[i].name))
//			longest = strlen(common_cmds[i].name);
//	}
//
//	puts("The most commonly used git commands are:");
//	for (i = 0; i < ARRAY_SIZE(common_cmds); i++) {
//		printf("   %s   ", common_cmds[i].name);
//		mput_char(' ', longest - strlen(common_cmds[i].name));
//		puts(common_cmds[i].help);
//	}
//}


int
cmd_help(int argc, char **argv)
{
   int nongit;
	const char *alias;
	enum help_format parsed_help_format;
	load_command_list("git-", &main_cmds, &other_cmds);

	argc = parse_options(argc, argv, prefix, builtin_help_options,
			builtin_help_usage, 0);
	parsed_help_format = help_format;

	if (show_all) {
		printf("usage: %s\n\n", git_usage_string);
		list_commands("git commands", &main_cmds, &other_cmds);
		printf("%s\n", git_more_info_string);
		return 0;
	}

	if (!argv[0]) {
		printf("usage: %s\n\n", git_usage_string);
		list_common_cmds_help();
		printf("\n%s\n", git_more_info_string);
		return 0;
	}

	setup_git_directory_gently(&nongit);
	git_config(git_help_config, NULL);

	if (parsed_help_format != HELP_FORMAT_NONE)
		help_format = parsed_help_format;

	alias = alias_lookup(argv[0]);
	if (alias && !is_git_command(argv[0])) {
		printf("`git %s' is aliased to `%s'\n", argv[0], alias);
		return 0;
	}

	switch (help_format) {
	case HELP_FORMAT_NONE:
	case HELP_FORMAT_MAN:
		show_man_page(argv[0]);
		break;
	case HELP_FORMAT_INFO:
		show_info_page(argv[0]);
		break;
	case HELP_FORMAT_WEB:
		show_html_page(argv[0]);
		break;
	}

	return 0;
}

int
cmd_version(int argc, char **argv)
{
   cout << PACKAGE_NAME << " version " << PACKAGE_VERSION << endl;
	return 0;
}
