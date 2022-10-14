/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : SeqList.h
//  Purpose   : Handles the sequence file functions for ParsInsert
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
#if !defined(__SEQLIST_H__)
#define __SEQLIST_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Knox_Stddef.h"
#include <vector>
#include <string>
#include <map>

using namespace std;

#define MAX_SEQ_SIZE (64*1024)

class CPNode;

///////////////////////////////////////////////////////////////////////////////////////////////////

class CBestLocationEntry
	{
	public:

		CPNode		*node;			// pointer to node with best match
		int			score;			// score of match at this position
		string		tax;			// taxonomy assignd to this position
		int			levels;			// number of ranks in taxonomy

	public:
		CBestLocationEntry(int s, CPNode *n, LPCSTR t)
			{
			score  = s;
			node   = n;
			tax    = t;
			levels = 0;

			if (!tax.empty())
				{
				for (int i=0; i < tax.length() ; ++i)
					if (tax[i] == ';')
						++levels;
				}
			}
	};

typedef vector<CBestLocationEntry> CBestList;

///////////////////////////////////////////////////////////////////////////////////////////////////

class CBestLocation
	{
	public:
		static int			default_len;

		CBestList			list;		// list of best matches
		int					len;		// number of matches to keep
	public:
		CBestLocation(int _len=-1)
			{
			if (_len > 0)
				len = _len;
			else
				len = default_len;
			}

		void Add(int score, CPNode *node, LPCSTR t);

		int  WorstScore();
	};

///////////////////////////////////////////////////////////////////////////////////////////////////

class CSequenceItem	
	{
	public:
		long		offset;		// offset in sequence file where sequence entry begins
		string		name;		// id from the fasta file
		int			len;		// length of the sequence data (sequence only)
		LPSTR		data;		// sequence string
		LPSTR		hdr;		// header line from fasta file

	public:
		CSequenceItem();
		CSequenceItem(LPCSTR _name, long _offset, int _len, LPCSTR _data=NULL);
		~CSequenceItem();

		LPCSTR		GetSeqData();

		BOOL		AllocateSeqData();
		BOOL		ReleaseSeqData();

		BOOL		AllocateSeqHeader(int size);
		BOOL		ReleaseSeqHeader();

		LPCSTR		ReadSeqHeader(FILE *f);
		LPCSTR		ReadSeqData(FILE *f);
		BOOL		WriteSeqData(FILE *f);

	};

typedef map<string,CSequenceItem*> 	CSequenceList;
typedef CSequenceList::iterator		CSequenceListIter;

///////////////////////////////////////////////////////////////////////////////////////////////////

class CSequenceFile
	{
	public:
		FILE				*f;					// file to read data
		string				fname;				// name of file opened

		int					seqCount;			// number of sequences used in array
		int					seqN;				// number of items allocated in array
		
		CSequenceItem*		*seqArray;			// array of sequence numbers
		CSequenceList		seqlist;			// map of name to sequence

		CSequenceListIter	seqIter;			// iterator to step thru map
		int					posArray;			// iterator position in array

        static int          progressCount;		// number of sequences between progress reports
        static int			verbose;
        
	public:
		CSequenceFile(LPCSTR _fname=NULL, int size=MAX_SEQ_SIZE, int hashsize=10*1024);
		~CSequenceFile();

		CSequenceItem * GetSequence(LPCSTR name);
		CSequenceItem * GetSequenceHeader(LPCSTR name);

		BOOL			Open(LPCSTR _fname, LPCSTR mode="r");
		void			Close();

		int				GetCount()
							{
							return seqCount;
							}

		int				ReadSequenceIndexFile(LPCSTR filename);
		int				ReadSequenceFile(int type);
		CSequenceItem	*GetSeq(LPCSTR name);
		BOOL			WriteSequence(CSequenceItem *seq);

		void			ResetSeqIterator();
		CSequenceItem	*GetNextSeq();

		int				ReadTaxonomyFile();

	};

#endif // __SEQLIST_H__
