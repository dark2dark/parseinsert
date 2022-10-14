/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : Knox_Stddef.cpp
//  Purpose   : Set of routines and old habit data types
//
//  Developer : David Knox (david.knox@colorado.edu) Jan 2011
//  Copyright : Copyright (C) 2007-2011 David Knox
//
//  Web site  : http://parsinsert.sourceforge.net/
//
/////////////////////////////////////////////////////////////////////////////////////////////////// 
//	  This file is part of ParsInsert.
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
#include "Knox_Stddef.h"
#include <ctype.h>

// File open for the application log
// File is always current (Flushed after each write)
FILE				*appLog = NULL;

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Application Log: Manages a log file.  Display functions will print to 
//					multiple locations (Stdout, Application log)
//

void OpenAppLog(LPCSTR filename, LPCSTR mode)
	{
	if (appLog != NULL)
		CloseAppLog();
		
	appLog = fopen(filename, mode);
	}

void CloseAppLog()
	{
	if (appLog != NULL)
		fclose(appLog);
		
	appLog = NULL;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Display functions will write to multiple locations (Stdout, Applog, ...)

void DisplayT(LPSTR format, ...)
	{
	char		msg[64*1024];
	va_list args;
	va_start(args, format);     /* Initialize variable arguments. */
	vsprintf(msg, format, args);
	printf("%s", msg);
	fflush(stdout);
	//TRACE("%s", msg);
	va_end(args);              /* Reset variable arguments.      */
	}

void Display(LPSTR format, ...)
	{
	char		msg[64*1024];
	va_list args;
	va_start(args, format);     /* Initialize variable arguments. */
	vsprintf(msg, format, args);
	printf("%s", msg);
	va_end(args);              /* Reset variable arguments.      */
	}

void DisplayL(LPSTR format, ...)
	{
	char		msg[64*1024];
	va_list args;
	va_start(args, format);     /* Initialize variable arguments. */
	vsprintf(msg, format, args);
	printf("%s", msg);
	if (appLog != NULL)
		{
		int		last = strlen(msg);
		if (last > 0) --last;

		fprintf(appLog, "%s%s", msg, (msg[last]!='\n'?"\n":""));
		fflush(appLog);
		}

	va_end(args);              /* Reset variable arguments.      */
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Read the next line from a file.  Skips empty lines or lines begining 
//		with "//" or "#".  Removes the trailing CR/LF/space.
//
BOOL ReadNextLine(LPSTR line, int len, FILE *f)
	{
	while ((f != NULL) && !feof(f) && fgets(line, len-1, f))
		{
		// skip empty lines, commented lines, lines without enough data
		if (strlen(line) == 0)	continue;
		if ((line[0] == '/') && (line[1] == '/')) continue;
		if (line[0] == '#') continue;

		// remove ending CR/LF and spaces
		for (int i=strlen(line) ; i > 0 ; i--)
			{
			if (   (line[i-1] == '\r') 
				|| (line[i-1] == '\n') 
				|| (isspace(line[i-1])))
				line[i-1] = 0;
			else 
				break;
			}
		return TRUE;
		}

	return FALSE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Trim unwanted characters from ends of strings
//
void Trim(string& s, LPCSTR unwanted)
	{
	TrimLeft(s, unwanted);
	TrimRight(s, unwanted);
	}
	
void TrimLeft(string& s, LPCSTR unwanted)
	{
	int start = s.find_first_not_of(unwanted);
	s = s.substr(start, s.length());
	}
	
void TrimRight(string& s, LPCSTR unwanted)
	{
	int		end = s.length();
	while ((end > 0) &&  (strchr(unwanted, s[end-1]) != NULL))
		--end;
	s = s.substr(0, end);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convert string to all one case	
//
void ToLower(string& s)
	{
	for (int i=0 ; i < s.length() ; ++i)
		s[i] = tolower(s[i]);
	}
	
void ToUpper(string& s)
	{
	for (int i=0 ; i < s.length() ; ++i)
		s[i] = toupper(s[i]);
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
