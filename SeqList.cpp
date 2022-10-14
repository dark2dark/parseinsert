/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : SeqList.cpp
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
#include "SeqList.h"
#include "PNode.h"
#include <vector>
#include <algorithm>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////

bool operator<(const CBestLocationEntry& x, const CBestLocationEntry& y)
	{
	// sort by score first
    if (x.score != y.score)
		return x.score < y.score;

	// check type of object
	if (x.node->IsBranchNode() && !y.node->IsBranchNode())
		return false;
	if (!x.node->IsBranchNode() && y.node->IsBranchNode())
		return true;

	// whoever's taxonomy is longer is first in list
	if (x.tax.empty() && !y.tax.empty())
		return false;
	if (!x.tax.empty() && y.tax.empty())
		return true;

	return x.levels > y.levels;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int			CBestLocation::default_len = 10;

///////////////////////////////////////////////////////////////////////////////////////////////////

int  CBestLocation :: WorstScore()
	{
	// return the current worst score
	if (list.size() < len)
		return INT_MAX;
	else
		return (list.end()-1)->score;	
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CBestLocation :: Add(int score, CPNode *node, LPCSTR taxonomy)
	{
	list.push_back(CBestLocationEntry(score, node, taxonomy));
	sort(list.begin(), list.end());
	if (list.size() > len)
		list.pop_back();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceItem :: CSequenceItem()
	{
	name.clear();
	offset = -1;
	len    = -1;
	data   = NULL;
	hdr    = NULL;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceItem :: CSequenceItem(LPCSTR _name, long _offset, int _len, LPCSTR _data)
	{
	name   = _name;
	offset = _offset;
	len    = _len;
	data   = NULL;
	hdr    = NULL;

	if (_data != NULL)
		{
		if (strlen(_data) > len)
			len = strlen(_data);
		if (len > 0)
			{
			AllocateSeqData();
			strcpy(data, _data);
			}
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceItem :: ~CSequenceItem()
	{
	if (data != NULL)
		ReleaseSeqData();

	if (hdr != NULL)
		ReleaseSeqHeader();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

LPCSTR	CSequenceItem :: GetSeqData()
	{
	return data;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CSequenceItem :: AllocateSeqData()
	{
	if (data != NULL)
		ReleaseSeqData();

	if (len > 0)
		{
		data = new char[len+1];
		memset(data, 0, len+1);
		}

	return (data != NULL);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL	CSequenceItem :: ReleaseSeqData()
	{
	if (data != NULL)
		{
		delete [] data;
		data = NULL;
		}

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CSequenceItem :: AllocateSeqHeader(int size)
	{
	if (hdr != NULL)
		ReleaseSeqHeader();

	hdr = new char[size];
	hdr[0] = 0;

	return (hdr != NULL);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL	CSequenceItem :: ReleaseSeqHeader()
	{
	if (hdr != NULL)
		{
		delete [] hdr;
		hdr = NULL;
		}

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

LPCSTR	CSequenceItem :: ReadSeqHeader(FILE *f)
		{
		if (f == NULL) return NULL;

		fseek(f, offset, SEEK_SET);
		if (ftell(f) != offset) return NULL;

		// assuming FASTA file format
		// skip the first line, seq description
		char		buffer[96*1024];

		fgets(buffer, sizeof(buffer)-1, f);
//		if (buffer[0] != '>') return NULL;

		if (strlen(buffer) > 0)
			{
			AllocateSeqHeader(strlen(buffer)+1);
			if (hdr != NULL)
				strcpy(hdr, buffer);
			}

		//TRACE("READ %d chars",strlen(data));
		return hdr;
		}

///////////////////////////////////////////////////////////////////////////////////////////////////

LPCSTR	CSequenceItem :: ReadSeqData(FILE *f)
		{
		if (f == NULL) return NULL;

		fseek(f, offset, SEEK_SET);
		if (ftell(f) != offset) return NULL;

		// assuming FASTA file format
		// skip the first line, seq description
		char		buffer[MAX_SEQ_SIZE];

		fgets(buffer, sizeof(buffer)-1, f);
//		if (buffer[0] != '>') return NULL;

		if (strlen(buffer) > 0)
			{
			AllocateSeqHeader(strlen(buffer)+1);
			if (hdr != NULL)
				strcpy(hdr, buffer);
			}

		if (len > 0)
			{
			AllocateSeqData();

			// while (we have not found the next entry)
			//		add new data to sequence
			while ((fgets(buffer, sizeof(buffer)-1, f) != NULL) && (buffer[0] != '>'))
				{
				//TRACE("===%s", buffer);
				// add alignment data to sequence string
				int			i;
				int			j;
				for (i=0 ; (i < sizeof(buffer)) && isspace(buffer[i]) && (buffer[i] != 0) ; ++i)
					; /* find first non-whitespace */
				for (j=i ; (j < sizeof(buffer)) && !isspace(buffer[j]) && (buffer[j] != 0) ; ++j)
					; /* find next whitespace */
				buffer[j] = 0;

				if (strlen(data)+strlen(&buffer[i]) <= len)	
					strcat(data, &buffer[i]);
				}
			}

		//TRACE("READ %d chars",strlen(data));
		return data;
		}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CSequenceItem :: WriteSeqData(FILE *f)
	{
	if (f == NULL)	return FALSE;
	if (data == NULL) return FALSE;

	offset = ftell(f);
	fprintf(f, ">%s" , name.c_str());
	for (int i=0 ; i < len ; i+=50)
		{
		if (i%50 == 0)
			fprintf(f, "\n          ");
		int count = min(len-i, 50);
		fwrite(&data[i], 1, count, f);
		}

	fprintf(f, "\n");
	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int CSequenceFile::progressCount 	= 20000;
int CSequenceFile::verbose 			= 0;

///////////////////////////////////////////////////////////////////////////////////////////////////
    
CSequenceFile :: CSequenceFile(LPCSTR _fname, int size, int hashsize)
	{
	seqCount = 0;
	f        = NULL;

	seqArray = NULL;
	seqN     = size;
	if (seqN > 0)
		{
		int			nbytes = seqN * sizeof(CSequenceItem*);
		seqArray = (CSequenceItem**)malloc(nbytes);
		memset(seqArray, 0, nbytes);
		}

	//seqlist.InitHashTable(hashsize);

	if (_fname != NULL)
		{
		Open(_fname);
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceFile :: ~CSequenceFile()
	{
	if (f != NULL)
		fclose(f);
	f = NULL;

	// delete all the seqs
	for (int i=0 ; i < seqN ; ++i)
		if (seqArray[i] != NULL)
			delete seqArray[i];
	
	if (seqArray != NULL)
		{
		delete seqArray;
		seqArray = NULL;
		}
////	memset(seqArray, 0, sizeof(seqArray));

	for (CSequenceListIter iter=seqlist.begin() ; iter != seqlist.end() ; ++iter)
		{
		delete ((*iter).second);
		}

	seqlist.clear();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceItem * CSequenceFile :: GetSequenceHeader(LPCSTR name)
	{
	CSequenceItem		*seq = this->GetSeq(name);

	if ((seq != NULL) && (seq->hdr == NULL))
		seq->ReadSeqHeader(f);
	
	return seq;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceItem * CSequenceFile :: GetSequence(LPCSTR name)
	{
	CSequenceItem		*seq = this->GetSeq(name);

	if ((seq != NULL) && (seq->data == NULL))
		seq->ReadSeqData(f);
	
	return seq;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CSequenceFile :: Open(LPCSTR _fname, LPCSTR mode)
	{
	if (f != NULL)
		Close();
	fname = _fname;
	f = fopen(fname.c_str(), mode);

	return (f != NULL);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CSequenceFile :: Close()
	{
	if (f != NULL)
		fclose(f);
	f = NULL;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CSequenceFile :: ResetSeqIterator()
	{
	seqIter  = seqlist.begin();
	posArray = 0;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceItem*	CSequenceFile :: GetNextSeq()
	{
	if (seqIter != seqlist.end())
		{
		CSequenceItem	*item = (*seqIter).second;
		++seqIter;
		return item;
		}

	while (posArray < seqN)
		{
		CSequenceItem	*seq = seqArray[posArray];
		++posArray;
		if (seq != NULL)
			return seq;
		}

	return NULL;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CSequenceFile :: ReadSequenceIndexFile(LPCSTR filename)
	{
	FILE		*f = fopen(filename, "r");
	if (f == NULL) return 0;

	// Index file format
	//     offset  len  name

	char				buffer[MAX_SEQ_SIZE];

	//int					seqlen = 0;
	CSequenceItem		*seq = NULL;
	int					count = 0;

	while (fgets(buffer, sizeof(buffer)-1, f) != NULL)
		{
		char *sep = strchr(buffer, '\t');
		if (sep == NULL) continue;

		*sep = 0;
		++sep;

		long offset = atol(buffer);

		int len = atoi(sep);
		sep = strchr(sep, '\t');
		if (sep == NULL) continue;

		*sep = 0;
		++sep;

		char	*name = sep;

		sep = strchr(sep, '\n');
		if (sep == NULL) continue;
		*sep = 0;
            
        sep = strchr(name, '\r');
        if (sep != NULL)
            *sep = 0;
            
		seq = new CSequenceItem(name, offset, len);
		if (isdigit(buffer[0]))
			{
			int id = atoi(name);
			if ((id > 0) && (id < sizeof(seqArray)/sizeof(CSequenceItem*)))
				{
				if (seqArray[id] != NULL)
					{
					delete seqArray[id];
					--seqCount;
				//	TRACE(" *** Replaced %d with new location\n", id);
					}
				seqArray[id] = seq;
				++seqCount;
				}
			else
				seqlist[seq->name] = seq;
			}
		else
			{
			seqlist[seq->name] = seq;
			}

		++count;
		if (count%progressCount == 0)
			{
			if (verbose)
				DisplayT("Loading sequence %d [%s]\n", count, name);
			else
				DisplayT("Loading sequence %d\n", count);
			}
		}

	DisplayT("Loaded %d sequences from [%s]\n", count, filename);
	fclose(f);

	return count;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CSequenceFile :: ReadTaxonomyFile()
	{
	if (f == NULL) return 0;

	fseek(f, 0, SEEK_SET); // rewind
	
	char				buffer[MAX_SEQ_SIZE];

	int					offset = ftell(f);
	CSequenceItem		*seq = NULL;
	
	seqCount = 0;

	while (fgets(buffer, sizeof(buffer)-1, f) != NULL)
		{
		int 		j;
		for (j=1 ; (j < sizeof(buffer)) && !isspace(buffer[j]) && (buffer[j] != 0) ; ++j)
			; /* find next whitespace */

		buffer[j] = 0;

		seq = new CSequenceItem(buffer, offset, 0);
		seqlist[seq->name] = seq;
		if (++seqCount%progressCount == 0)
			if (verbose)
				DisplayT("Loading sequence %d [%s]\n", seqCount, seq->name.c_str());
			else
				DisplayT("Loading sequence %d\n", seqCount);

		offset = ftell(f);
		}

	return seqCount;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CSequenceFile :: ReadSequenceFile(int type)
	{
	if (f == NULL) return 0;

	fseek(f, 0, SEEK_SET); // rewind
	
	char				buffer[MAX_SEQ_SIZE];

	int					seqlen = 0;
	int					offset = ftell(f);
	CSequenceItem		*seq = NULL;
	
	//seqCount = 0;

	while (fgets(buffer, sizeof(buffer)-1, f) != NULL)
		{
		if (buffer[0] == '>')
			{
			// found a new entry
			// update last entry
			if (seq != NULL)
				{
				seq->len = seqlen;
				if (seqCount%progressCount == 0)
					DisplayT("Loading sequence %d [%s]\n", seqCount, seq->name.c_str());
				}

			// extract taxa name
			// find first whitespace
			int			i;
			int			j;
			for (i=1 ; (i < sizeof(buffer)) && (isspace(buffer[i])||((buffer[i] == '>'))) && (buffer[i] != 0) ; ++i)
				; /* find first non-whitespace */
			for (j=i ; (j < sizeof(buffer)) && !isspace(buffer[j]) && (buffer[j] != 0) ; ++j)
				; /* find next whitespace */
			buffer[j] = 0;

			seq = new CSequenceItem(&buffer[i], offset, 0);
			if (seqlist.find(seq->name) != seqlist.end())
				DisplayL("Duplicate sequence found [%s]\n", seq->name.c_str());
			else
				++seqCount;
			seqlist[seq->name] = seq;

			seqlen = 0;
			}
		else
			{
			// update the sequence length
			int			i;
			int			j;
			for (i=0 ; (i < sizeof(buffer)) && isspace(buffer[i]) && (buffer[i] != 0) ; ++i)
				; /* find first non-whitespace */
			for (j=i ; (j < sizeof(buffer)) && !isspace(buffer[j]) && (buffer[j] != 0) ; ++j)
				; /* find next whitespace */					
			seqlen += (j - i + 1);
			}

		offset = ftell(f);
		}

	if (seq != NULL)
		{
		seq->len = seqlen;
		DisplayT("Loading sequence %d [%s]\n", seqCount, seq->name.c_str());
		}

	return seqCount;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceItem* CSequenceFile :: GetSeq(LPCSTR name)
	{
	int					id   = atoi(name);
	CSequenceItem		*seq = NULL;
	if ((id > 0) && (id < sizeof(seqArray)/sizeof(CSequenceItem*)))
		seq = seqArray[id];

	if (seq == NULL)
		seq = (CSequenceItem*)seqlist[name];

	return seq;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CSequenceFile :: WriteSequence(CSequenceItem *seq)
	{
	if (f == NULL) return FALSE;
	if (seq == NULL) return FALSE;

	fseek(f, 0, SEEK_END);
	seq->WriteSeqData(f);

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
