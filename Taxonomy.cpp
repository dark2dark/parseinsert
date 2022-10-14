/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : Taxonomy.cpp
//  Purpose   : Handles the taxonomy functions for ParsInsert
//
//              CTaxEntry: Class to manage heirarchy of taxonomic classifications
//				           Count of number of occurances of given name at given rank
//				           Keeps list of all entries that were defined to lower ranks
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
#include "Taxonomy.h"

///////////////////////////////////////////////////////////////////////////////////////////////////

CTaxEntry :: CTaxEntry(LPCSTR _name)
    {
	name = _name;
	count = 0;
	entries.clear();
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

CTaxEntry :: ~CTaxEntry()
    {
	// release all the subentries
    for (CTaxEntryList::iterator iter=entries.begin() ; iter != entries.end() ; ++iter)
        {
        delete (*iter).second;
        }
    entries.clear();
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CTaxEntry :: Add(LPCSTR tax, int votes)
    {
	// Add the first name in the taxonomy to this entry
	count += votes;
	if (strlen(tax) == 0)
		return TRUE;
    
	// strip next rank off string
	string		subentry = tax;  // assume whole string is this name
	string		rest = "";
    int         len = strcspn(tax, ";/");
    if (len < subentry.length())
        {
        rest = subentry.substr(len+1);
        subentry = subentry.substr(0, len);
        }
        
	Trim(subentry, " \r\n\t");
    
	// lookup in table to get entry
	CTaxEntry *ptr = entries[subentry];
	if (ptr == NULL)
		{
        ptr = new CTaxEntry(subentry.c_str());
        entries[subentry] = ptr;
		}
    
	// add rest to that entry
	ptr->Add(rest.c_str(), votes);
    
	return TRUE;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CTaxEntry :: Display(LPCSTR leader, int levels)
    {
    if (levels <= 0) return FALSE;
    --levels;
    
	if (leader != NULL)
		DisplayL("%s[%s] %d\n", leader, name.c_str(), count);
    
	string		indent = leader;
	indent += "   ";
    
	// display all the subentries
    for (CTaxEntryList::iterator iter=entries.begin() ; iter != entries.end() ; ++iter)
		{
        string      key  = (*iter).first;
        CTaxEntry	*ptr = (*iter).second;
        
        ptr->Display(indent.c_str(), levels);
		}
    
	return TRUE;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

int CTaxEntry :: FindBest(int threshPercent, CStringList& taxList)
	{
	// Find the best taxonomy with at least thresh count
	taxList.clear();
    
	int			thresh = count * threshPercent / 100;
    
    for (CTaxEntryList::iterator iter=entries.begin() ; iter != entries.end() ; ++iter)
		{
        string          key  = (*iter).first;
        CTaxEntry       *ptr = (*iter).second;
            
        ptr->FindBestSubentry(thresh, "", taxList);
		}
    
	return taxList.size();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CTaxEntry :: FindBestSubentry(int thresh, LPCSTR lineage, CStringList& taxList)
	{
	// Find the best taxonomy with at least thresh count
	int			ntax = 0;
	string		myLineage = lineage;
    
	if (!myLineage.empty())
        myLineage += ";";
    myLineage += name;
    
    for (CTaxEntryList::iterator iter=entries.begin() ; iter != entries.end() ; ++iter)
		{
            string      key  = (*iter).first;
            CTaxEntry	*ptr = (*iter).second;
            
            if (ptr->count < thresh)
                continue;
            
            ++ntax;
            ptr->FindBestSubentry(thresh, myLineage.c_str(), taxList);
		}
    
	if (ntax == 0)
		taxList.push_back(myLineage); // if no other lineage added, add my lineage
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
