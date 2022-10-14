/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : Knox_Stddef.h
//  Purpose   : Set of routines and old habit data types
//
//  Developer : David Knox (david.knox@colorado.edu) Jan 2011
//  Copyright : Copyright (C) 2007-2011 David Knox
//
//  Web site  : http://parsinsert.sourceforge.net/
//
/////////////////////////////////////////////////////////////////////////////////////////////////// 
//	This file is part of ParsInsert.
//
//    ParsInsert is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ParsInsert is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ParsInsert.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////////////////////////
#if !defined(__KNOX_STDDEF_H__)
#define __KNOX_STDDEF_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <string>
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////

typedef vector<string> CStringList;

#define PVOID void*
#define LPSTR char *
#define LPCSTR const char *
#define BYTE unsigned char
#define WORD unsigned short
#define DWORD long unsigned int
#define BOOL int
#define TRUE 1
#define FALSE 0

#ifndef INT_MAX
#define INT_MAX 20000000
#endif

#define DEBUG TRUE
///////////////////////////////////////////////////////////////////////////////////////////////////

//
// Function that will only print if in DEBUG mode
#define TRACE if (DEBUG) printf

//
// Application Log: Manages a log file.  Display functions will print to 
//					multiple locations (Stdout, Application log)
//					File is always current (Flush after each write)
//
extern FILE		*appLog;
extern void 	OpenAppLog(LPCSTR filename, LPCSTR mode="a+");
extern void 	CloseAppLog();

extern void 	Display(LPSTR format, ...);
extern void 	DisplayT(LPSTR format, ...);
extern void 	DisplayL(LPSTR format, ...);

// Helper Functions

//
//  Read Next Line from file into buffer
//
extern BOOL		ReadNextLine(LPSTR line, int len, FILE *f);

// Trim the string, removing unwanted characters from the ends.
extern void		Trim(string& s, LPCSTR unwanted);	
extern void		TrimLeft(string& s, LPCSTR unwanted);	
extern void		TrimRight(string& s, LPCSTR unwanted);	

// convert the string to all same case
extern void		ToLower(string& s);
extern void		ToUpper(string& s);

///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
