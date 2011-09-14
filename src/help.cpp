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
 * along with Imogene.  If not, see <http://www.gnu.org/licenses/>.
 *
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
#include <cstring>

#include "config.h"
#include "imogene.hpp"

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
		cout << c;
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
		cout << "   " << common_cmds[i].name << "   ";
		mput_char(' ', longest - strlen(common_cmds[i].name));
		cout << common_cmds[i].help << endl;
	}
}

int
cmd_help(int argc, char **argv)
{
   cout <<  "Usage : " << usage_string << endl;
   list_common_cmds_help();
   cout << "\n" << more_info_string << endl;

	return 0;
}

int
cmd_version(int argc, char **argv)
{
   cout << PACKAGE_NAME << " version " << PACKAGE_VERSION << endl;
	return 0;
}
