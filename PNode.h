/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : PNode.h
//  Purpose   : Phylogenetic Tree Node
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
#if !defined(__PNODE_H__)
#define __PNODE_H__

#include "Knox_Stddef.h"
#include "AttrList.h"
#include "ParsimonySet.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

extern LPCSTR globalAttrs[];
extern LPCSTR treeAttrs[];
extern LPCSTR branchAttrs[];
extern LPCSTR branchAttrs2[];
extern LPCSTR applyAttrs[];

/////////////////////////////////////////////////////////////////////////////////////////////////// 

extern BOOL ValidRegEx(LPCSTR expr);
extern BOOL MatchRegEx(LPCSTR expr, LPCSTR text);

/////////////////////////////////////////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////// 

class CPNode;
class CPTree;

typedef vector<CPNode*> CPNodeList;
typedef vector<CPNode*>::iterator CPNodeListIter;

class CPNode 
	{
	public:

		string				title;			// name of the node
		int					id;				// internal id
		CPTree				*tree;			// tree the node belongs to
		CPNode				*parent;		// parent node of this node
		
		CPNodeList			edges;			// children
		CAttrList			attrs;			// list of attributes assigned to this node

	public:

		CPNode(CPTree *ptree, LPCSTR v="");
		~CPNode();

		void				AddEdge(CPNode *other);
		CPNode*				FindEdgeNode(int e);
		void				InvertEdges(BOOL recurse=FALSE);	// invert the edge list
		inline BOOL			IsBranchNode() {return (edges.size() > 0); }

		BOOL				GetComment(LPCSTR  attr, string& str, LPCSTR _default=NULL);
		void				SetAttr(LPCSTR k, string& v);
		void				SetAttr(LPCSTR k, LPCSTR v);
		void				AddAttrs(LPCSTR v);
		BOOL				ReadAttrs(LPCSTR filename);
		BOOL				WriteAttrs(FILE *f);
		BOOL				Display(LPCSTR key);
		BOOL				WriteNode(FILE *f, int indent=0, BOOL includeAttrs=FALSE, BOOL withComment=FALSE);
		string&				CleanString(string& s, BOOL searchForBootstrap=FALSE);

		void				ResetParent(CPNode *p);
								// traverse the tree setting parent to correct items
								// used by tree inversion

		int					GetAncestors(CPNodeList &list, BOOL fromBottom=TRUE);
								// return the ancestors in an ordered list
                                // fromBottom - if TRUE, inserts each parent at end of list 

		int					GetDecendants(CPNodeList &list, BOOL leafOnly=TRUE);
								// return the nodes the decendants
								// leafOnly - if TRUE, only add the leaf id's to list
	};

/////////////////////////////////////////////////////////////////////////////// 

typedef BOOL (* ProgressFunction)(LPCSTR msg, int curInProcess, int totalToProcess, PVOID data);

class CPTree
	{
	public:	
		string					filename;			// file tree read from or wriiten to
		CPNode					*root;				// root node of the tree
		CAttrList				attrs; 				// list of attributes for tree

		int						internal;			// internal id number for last created node
		CPNodeList				nodeList;			// list of nodes in the tree
///		map<string,CPNode*>		nameTable;			// lookup table for node names

		//
		// Variables for parsing
		//
		LPCSTR					buffer;				// points to current location in the buffer to get chars
		FILE					*bufFile;			// file to read from 
		int						lineCount;			// current position in lines of the file
		int						colCount;			// Current column within the line
		char					errorReason[2048];

		// Status reporting 
		ProgressFunction		Progress;
		int						totalProcess;
		int						nProcess;
		LPCSTR					endParseStr;		// str position when parsing complete

		//
		// Class variables
		//
		static int				PRECISION;  		// number of decimal points to use when printing
	public:
		CPTree();
		~CPTree();

		void				ShowError(int codeline, LPCSTR msg);

		BOOL				WriteNewickTree(LPCSTR filename, CPNode *base, BOOL withComment=FALSE);

		CPNode *			Decode(LPCSTR &str);		// Decode encoded node description
		void				Parse(LPCSTR str);
		CPNode*				ParseTree(LPCSTR &str);
		void				ReadComment(LPCSTR &v, string &s);
		static string&		QuoteString(string& comment, char quote='\"');
		void				GetQuotedString(LPCSTR &str, string &v);
		void				ReadToken(LPCSTR &str, string &v);
		BOOL				ReadTreeAttrs();
		BOOL				ReadNodeAttrs();
		BOOL				WriteAttrs();
		void				SetAttr(LPCSTR k, LPCSTR v);
		BOOL				RemoveAttribute(LPCSTR attr, BOOL all);
								// Remove the given attribute.
								// all - if TRUE, remove from all nodes, 
								//		 otherwise from only from global TREE attrs						
		void				ProcessAttrList(LPCSTR str, CPNode *node);
		
		CPNode*				Find(LPCSTR v);
		CPNode*				FindCommonAncestor(CPNode *n1, CPNode *n2);
								// returns the nearest common ancestor to given nodes

		// Progress Functions Access
		BOOL				SetProgressFunction(ProgressFunction callBack);
		BOOL				ClearProgressFunction();
		BOOL				ShowProgress(LPCSTR msg, int n, int total, PVOID data);
	};


#endif
