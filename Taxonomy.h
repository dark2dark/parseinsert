/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : Taxonomy.h
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
#if !defined(__TAXONOMY_H__)
#define __TAXONOMY_H__

#include "Knox_Stddef.h"
#include <map>

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////////////

class CTaxEntry;
typedef map<string, CTaxEntry*> CTaxEntryList;

class CTaxEntry		// Class to manage heirarchy of taxonomic classifications
	{
	public:
		string				name;		// name of this item
		int					count;		// number of times this item seen
		CTaxEntryList		entries;	// list of entries for each sub-type
        
	public:
		CTaxEntry(LPCSTR _name);
		~CTaxEntry();
        
		BOOL		Add(LPCSTR tax, int votes=1);		
       				// Add the first name in the taxonomy to this entry
        
		int			FindBest(int threshPrecent, vector<string>& taxList);
			        // Find the best taxonomy with at least threshold count
        
		BOOL		Display(LPCSTR leader, int levels=100);	
			        // Show entries in hierarical list
        
	protected:
        
		void		FindBestSubentry(int thresh, LPCSTR lineage, CStringList& taxList);
					// Find the best taxonomy with at least thresh count
					// Recursive internal use
					// thresh is a raw count
	};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#endif
