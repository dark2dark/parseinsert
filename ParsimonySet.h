/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : ParsimonySet.h
//  Purpose   : Handles the parsimony functions for ParsInsert
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
#if !defined(__PARS_SET_H__)
#define __PARS_SET_H__

#include "Knox_Stddef.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////// 

class CParsimonySet
	{
	public:
		BYTE		*data;
		int			len;
		int			start;
		int			end;

        int         nSeg;
        BYTE        *segCounts;

	public:
		CParsimonySet()
			{
			data  = NULL;
			len   = 0;
			start = -1;
			end   = -1;
            
            nSeg = 0;
            segCounts = NULL;
			}

		CParsimonySet(const CParsimonySet *other)
			{
			data  = NULL;
			len   = other->len;
			start = other->start;
			end   = other->end;

            nSeg = 0;
            segCounts = NULL;
            
			Allocate(len);
			memcpy(data, other->data, len);
			}

		~CParsimonySet()
			{
			Release();
			len = 0;
			}

		void Release()
			{
			if (data != NULL)
				delete [] data;
			data = NULL;
            if (segCounts != NULL)
                delete [] segCounts;
            segCounts = NULL;
			}

		BOOL Allocate(int _len)
			{
			if (data != NULL)
				Release();

			len = _len;
			data = new BYTE[len];

			memset(data, 0, sizeof(BYTE)*len);

			return TRUE;
			}

		BOOL Convert(LPCSTR seq)
			{
			if (seq == NULL) return FALSE;

			Allocate(strlen(seq));

			// check to see if this is bits or nucleotides
			int			span = 0;
			for (LPCSTR s=seq ; *s != 0 ; ++s)
				if (isxdigit(*s))
					++span;
			BOOL		bits = (span > len/2);

			for (int i=0 ; *seq != 0 ; ++seq,++i)
				{
				BYTE			val = 0;

				if (bits)
					{
					switch (toupper(*seq))
						{
						case '1': val = 1; break;
						case '2': val = 2; break;
						case '3': val = 3; break;
						case '4': val = 4; break;
						case '5': val = 5; break;
						case '6': val = 6; break;
						case '7': val = 7; break;
						case '8': val = 8; break;
						case '9': val = 9; break;
						case 'A': val = 10; break;
						case 'B': val = 11; break;
						case 'C': val = 12; break;
						case 'D': val = 13; break;
						case 'E': val = 14; break;
						case 'F': val = 15; break;
						}
					}
				else
					{
					/*
						A 	Adenosine
						C 	Cytosine
						G 	Guanine
						T 	Thymidine
						U 	Uracil
						R 	G A (puRine)
						Y 	T C (pYrimidine)
						K 	G T (Ketone)
						M 	A C (aMino group)
						S 	G C (Strong interaction)
						W 	A T (Weak interaction)
						B 	G T C (not A) (B comes after A)
						D 	G A T (not C) (D comes after C)
						H 	A C T (not G) (H comes after G)
						V 	G C A (not T, not U) (V comes after U)
						N 	A G C T (aNy)
						X 	masked
						- 	gap of indeterminate length
					*/
#define rA 0x01
#define rC 0x02
#define rG 0x04
#define rT 0x08
					switch (toupper(*seq))
						{
						case 'A': val = rA;					break;
						case 'C': val = rC;					break;
						case 'G': val = rG;					break;
						case 'T': val = rT;					break;
						case 'U': val = rT;					break;
						case 'R': val = rA | rG;			break;
						case 'Y': val = rC | rT;			break;
						case 'K': val = rG | rT;			break;
						case 'M': val = rA | rC;			break;
						case 'S': val = rG | rC;			break;
						case 'W': val = rA | rT;			break;
						case 'B': val = rG | rC | rT;		break;
						case 'D': val = rA | rG | rT;		break;
						case 'H': val = rA | rC | rT;		break;
						case 'V': val = rA | rC | rG;		break;
						case 'N': val = rA | rC | rG | rT;	break;
						case 'X': val = 0;					break;
						case '.': break;
						case '-': break;
						default:
							if (DEBUG)	printf("%c", toupper(*seq));
						}
					}

				data[i] = val;
				}

            SetEnds();
            
			return TRUE;
			}

		BOOL Union(CParsimonySet *other)
			{
			if (other == NULL)   return FALSE;
			if (other->len == 0) return FALSE;

			if (data == NULL)
				Allocate(other->len);

			for (int i=0 ; (i < len) && (i < other->len) ; ++i)
				{
				data[i] |= other->data[i];
				}

			return TRUE;
			}

		BOOL Intersect(CParsimonySet *other)
			{
			if (other == NULL) return FALSE;

			if (data == NULL)
				{
				Allocate(other->len);
				memcpy(data, other->data, other->len);
				}
			else
				for (int i=0 ; (i < len) && (i < other->len) ; ++i)
					{
					data[i] &= other->data[i];
					}

			return TRUE;
			}

		int Set(CParsimonySet *unionSet, CParsimonySet *intersectSet)
			{
			int				cost = 0;

			if (unionSet == NULL) return cost;
			if (intersectSet == NULL) return cost;

			if (data == NULL)
				Allocate(unionSet->len);

			for (int i=0 ; (i < len) && (i < unionSet->len) && (i < intersectSet->len) ; ++i)
				{
				data[i] = intersectSet->data[i];
				if (data[i] == 0)
					{
					if (unionSet->data[i] != 0) // check for both being gaps
						{
						cost += 1;
						data[i] = unionSet->data[i];
						}
					}
				}
                
            SetEnds();

			return cost;
			}
            
        int Force(CParsimonySet *parentSet)
            {
			if (parentSet == NULL)  return 0;            
			if (data == NULL)       return 0;
            
       //     LPSTR        str;
       //     parentSet->BuildString(str, -1);
       //     DisplayL("Parent: [%200.200s]\n", str);
       //     BuildString(str, -1);
       //     DisplayL("Child:  [%200.200s]\n", str);

            int count = 0;
			for (int i=0 ; (i < len) && (i < parentSet->len) && (i < len) ; ++i)
				{
                BYTE    orig = data[i];
               // BYTE    p    = parentSet->data[i];
               // if ((p != orig) && p)
               //    ++count;
                BYTE    d    = parentSet->data[i] & data[i];
                if (d != 0)
                    data[i] = d;
                if (data[i] != orig)
                    ++count;
				}
                
       //     DisplayL("Updated:[%200.200s]\n\n", str);
       //     free(str);
           return count;
            }
        
        void SetEnds()
            {
            start = -1;
            end   = len;
            
            if (data == NULL) return;
            
            for (int i=0 ; i < len ; ++i)
                {
                if (data[i] != 0)
                    {
                    if (start < 0)
                        start = i;
                    end = i;
                    }
                }
            }

		LPSTR BuildString(LPSTR &str, int str_len)
			{
			if ((str == NULL) || (str_len < 1))
				{
				str_len = len + 1;
				str = (LPSTR)calloc(str_len,1);
				}
				
			str[0] = 0;
			
			for (int i=0 ; i < len ; ++i)
				{
				char x[10];
				sprintf(x,"%X", data[i]);
				strcat(str, x);
				}
				
			return str;
			}

		void BuildSets(LPSTR str, int str_len)
			{
			char		x[10];
			
			str[0] = 0;
			strcat(str, "[");

			for (int i=0 ; i < len ; ++i)
				{
				x[0] = 0;
				if (data[i]&1)
					strcat(x, "A");
				if (data[i]&2)
					strcat(x, "C");
				if (data[i]&4)
					strcat(x, "G");
				if (data[i]&8)
					strcat(x, "T");
				if (data[i] == 0)
					strcat(x, ".");

				if (strlen(x) > 1)
					{
					strcat(str, "{");
					strcat(str, x);
					strcat(str, "}");
					}
				else
					strcat(str, x);
				}
			
			strcat(str, "]");
			}

// SEGMENT_SIZE must be 256 or less to use BYTE counts
#define SEGMENT_SIZE 16
		int BuildSegCounts(char *mask=NULL)
            {
			if (segCounts != NULL)
				{
                delete [] segCounts;
                segCounts = NULL;
                nSeg = 0;
                }
            
            if (data == NULL)
                return -1;
                
			nSeg = (len + SEGMENT_SIZE-1) / SEGMENT_SIZE;
			segCounts = new BYTE[nSeg];
            
			memset(segCounts, 0, sizeof(BYTE)*nSeg);
            
            int         seg = 0;
            for (int j=0 ; (j < len) ; j+=SEGMENT_SIZE)
                {
                // count items in segment
                int     n = 0;
                for (int k=0 ; (k < SEGMENT_SIZE) && (k+j < len) ; ++k)
                    {
                    if ((mask == NULL) || (mask[j+k] != 0))
                        if (data[j+k] != 0)
                            ++n;
                    }
                if (n > 255)
                    n = 255;
                segCounts[seg++] = n;
                }  
            
			return nSeg;
            }

        int CompareSegments(CParsimonySet *other)
            {
            // compare the counts, looking for known error
            int         s   = (start + SEGMENT_SIZE-1) / SEGMENT_SIZE;
            int         e   = (end + 1)   / SEGMENT_SIZE;
        //    BYTE		*p1   = &segCounts[s];
        //    BYTE		*p2   = &other->segCounts[s];
            int         indel = 0;
            
            for (int k=s ; (k < nSeg) && (k < e) ; ++k)
                {
                indel += abs(segCounts[k] - other->segCounts[k]);
        //        indel += abs(*p1-*p2);
        //        ++p1;
        //        ++p2;
                }
            
            return indel;
            }
            
        void TraceSegments(LPCSTR label)
            {
            if (segCounts == NULL) return;
            
            DisplayL("%s [", label);
            for (int k=0 ; k < nSeg ; ++k)
                DisplayL(" %3d", segCounts[k]);
            DisplayL("]\n");
            }
            
		void Trace(LPCSTR label, LPCSTR str)
			{
			if (DEBUG)
				{
				printf("%s", label);
				while (strlen(str) > 500)
					{
					printf("%.500s", str);
					str += 500;
					}
				printf("%s\n", str);
				}
			}
	};

/////////////////////////////////////////////////////////////////////////////////////////////////// 

typedef map<string,CParsimonySet*> 	CParsimonyList;
typedef CParsimonyList::iterator	CParsimonyListIter;

/////////////////////////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////// 
#endif
