/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : AttrList.h
//  Purpose   : Class to keep set of key,value pairs
//				Access methods to get/set using different data types.
//				All information is stored as strings and converted as required.  This will
//				lead to slow execution for variables being acessed often.  This class is 
//				intended for long term storage of attributes.
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
#if !defined(__ATTRLIST_H__)
#define __ATTRLIST_H__

#include <string>
#include <map>

#include "Knox_Stddef.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
typedef pair<string,string>			CAttrEntry;
typedef map<string,string>			CAttrDict;
typedef CAttrDict::const_iterator	CAttrDictIter;

class CAttrList
{
	public:
		CAttrDict	m;

		LPCSTR		operator[](LPCSTR key) const;

		int			GetInt(LPCSTR key, int _default=INT_MAX);
		DWORD		GetHex(LPCSTR key, DWORD _default=INT_MAX);
		BOOL		GetBOOL(LPCSTR key, int _default=-1);
		double		GetDouble(LPCSTR key, double _default=INT_MAX);
		string		GetString(LPCSTR key, LPCSTR _default="");

		BOOL		Add(LPCSTR key, int value);
		BOOL		Add(LPCSTR key, double value);
		BOOL		Add(LPCSTR key, LPCSTR value);
		BOOL		Add(LPCSTR key, string& value);
		BOOL		AddHex(LPCSTR key, DWORD value);

		BOOL		Add(const CAttrList& list);
		
		BOOL		Clear();
};


///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
