/*
 * =====================================================================================
 *
 *       Filename:  display.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06.08.2011 13:03:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef Display_H
#define Display_H

class svg
{
public:
	int xsize;
	int ysize;
	int xoffset;
	int yoffset;
	int pos;
	
	svg();
};

int cmd_extract(int argc, char **argv);


#endif // Display_Hpp
