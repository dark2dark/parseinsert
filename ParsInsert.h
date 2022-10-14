/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : ParsInsert.h
//  Purpose   : Parsimonious Insertion of unclassified seqeunces into Phylogenetic Tree
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
#ifndef PARS_INSERT_H

#include "SeqList.h"
#include "ParsimonySet.h"

#define RC_ERROR_OPTIONS                -1000   // Error in the command line options
#define RC_ERROR_CMDLINE                -1001   // Error not enough parameters supplied
#define RC_ERROR_MASKFILE               -1002   // Error accessing mask file
#define RC_ERROR_TAXFILE                -1003   // Error accessing taxonomy file
#define RC_ERROR_TREEFILE               -1004   // Error accessing newick tree file
#define RC_ERROR_TREESEQ                -1005   // Error accessing tree seq for tree
#define RC_ERROR_TREESEQ_MISSING        -1006   // Error could not find tree leaf in sequences
#define RC_ERROR_NO_INSERT_SEQ          -1007   // Error no insertion sequences to process
#define RC_ERROR_
#define RC_ERROR_
#define RC_ERROR_

///////////////////////////////////////////////////////////////////////////////////////////////////

class CInsertPos 	// class stores a possible seqeunce insertion location into a tree
	{
	public:
		CSequenceItem		*seq;		// sequence being inserted
		int					nSites;		// number of active sites being scored
		CParsimonySet		pars;		// parsimony sequence for insertion
		CBestLocation		best;		// location with best match
        
	public:
		CInsertPos()
        {
			seq    = NULL;
			nSites = 0;
        }
        
		~CInsertPos()
        {
			seq = NULL;
			best.list.clear();
        }
	};
typedef vector<CInsertPos*> CInsertPosArray;

///////////////////////////////////////////////////////////////////////////////////////////////////

#endif //PARS_INSERT_H
