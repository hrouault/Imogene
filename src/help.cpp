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

	puts("The most commonly used imogene commands are:");
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
   cout << "help function\n";
}
