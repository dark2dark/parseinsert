/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : AttrList.cpp
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

#include "AttrList.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// AttrList.cpp: implementation of the CAttrList class.
//
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////////////////////////////////

LPCSTR CAttrList :: operator[](LPCSTR key) const
	{
	CAttrDictIter iter;
	iter = m.find(key);
	if (iter == m.end())
		return NULL;

	return (*iter).second.c_str();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CAttrList :: Clear()
	{
	m.clear();

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CAttrList :: GetInt(LPCSTR key, int _default)
	{
	LPCSTR	value = m[key].c_str();

	if ((value == NULL) || (value == ""))
		return _default;

	return atoi(value);
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////

DWORD CAttrList :: GetHex(LPCSTR key, DWORD _default)
	{
	LPCSTR	value = m[key].c_str();

	if ((value == NULL) || (value == ""))
		return _default;

	return strtol(value, NULL, 0);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CAttrList :: GetBOOL(LPCSTR key, int _default)
	{
	LPCSTR	value = m[key].c_str();

	if ((value == NULL) || (value == ""))
		return _default;

	return (atoi(value) != 0);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

double CAttrList :: GetDouble(LPCSTR key, double _default)
	{
	LPCSTR	value = m[key].c_str();

	if ((value == NULL) || (value == ""))
		return _default;

	return atof(value);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

string CAttrList :: GetString(LPCSTR key, LPCSTR _default)
	{
    CAttrDictIter iter;
    iter = m.find(key);
    
	string	value;
	if (iter != m.end())
		value = iter->second;
    else if (_default != NULL)
        value = _default;
        
	return value;
	}


///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CAttrList :: Add(LPCSTR key, int value)
	{
	char buf[64];
	sprintf(buf,"%d",value);		// itoa(value, buf, 10);

	m[key] = buf;
	
	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CAttrList :: Add(LPCSTR key, double value)
	{
	char buf[64];
	sprintf(buf, "%f", value);

	m[key] = buf;

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CAttrList :: Add(LPCSTR key, LPCSTR value)
	{
	m[key] = value;

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CAttrList :: Add(LPCSTR key, string& value)
	{
	m[key] = value;

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CAttrList :: AddHex(LPCSTR key, DWORD value)
	{
	char buf[64];
	sprintf(buf, "0x%08X", value);
	m[key] = buf;

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CAttrList :: Add(const CAttrList &attrs)
	{
	CAttrDictIter iter;
	
	for (iter=attrs.m.begin() ; iter != attrs.m.end() ; ++iter)
		{
		m[(*iter).first] = (*iter).second;
		}

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

