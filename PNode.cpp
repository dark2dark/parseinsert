/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : PNode.cpp
//  Purpose   : Phylogenetic Tree Node
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

#include "PNode.h"
#include "Attrs.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class variables
//
int				CPTree::PRECISION = 5;  		// number of decimal points to use when printing

///////////////////////////////////////////////////////////////////////////////////////////////////

string& MakeCSVQuoted(string &q_str, LPCSTR others=NULL);
string& MakeQuoted(string &q_str, LPCSTR others=NULL);

///////////////////////////////////////////////////////////////////////////////////////////////////

CPTree :: CPTree()
	{
	root = NULL;

	internal = 0;
	
	nodeList.clear();

	Progress = NULL;
	totalProcess = 0;
	nProcess = 0;

	bufFile = NULL;
	}

CPTree :: ~CPTree()
	{
	CPNodeListIter iter = nodeList.begin();
	while (iter != nodeList.end())
		{
		if ((*iter) != NULL)
			{
			(*iter)->tree = NULL;
			delete (*iter);
			}
		++iter;
		}

	nodeList.clear();
	root = 0;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL	CPTree :: SetProgressFunction(ProgressFunction callBack)
	{
	Progress = callBack;
	return TRUE;
	}

BOOL	CPTree :: ClearProgressFunction()
	{
	Progress = NULL;
	return TRUE;
	}

BOOL	CPTree :: ShowProgress(LPCSTR msg, int n, int total, PVOID data)
	{
	if (Progress != NULL)
		return (*Progress)(msg, n, total, data);

	printf("Progress: (%d of %d) %s\n", n, total, msg);

	return TRUE;
	}
	
void CPTree::ShowError(int codeline, LPCSTR msg)
	{
	printf("ERROR [CodeLine:%d] %s [Input line:%d, col: %d]\n", codeline, msg, lineCount, colCount);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CPTree :: WriteNewickTree(LPCSTR filename, CPNode *base, BOOL withComment)
	{
	if (filename == NULL) 
		return FALSE;

	if (base == NULL)
		return FALSE;

	FILE		*f = fopen(filename, "wb+");
	if (f == NULL) 
		return FALSE;

	base->WriteNode(f, 0, FALSE, withComment);

	fprintf(f, ";\n");

	fclose(f);

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL IsAllDigits(LPCSTR s)
	{
	while ((*s != 0) && isdigit(*s))
		++s;

	return (*s == 0);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPTree :: GetQuotedString(LPCSTR &str, string &v)
	{
	// search for ending quote
	char quote = *str;

	v += *str;	// add beginning quote

	// read until next matching quote
	do 	{
		++str;
		++nProcess;
		v += *str;
		if ((*str) == '\\')
			{
			++str;
			v += *str;
			}
		if ( ((*str) == '\'') && (str[1] == '\''))
			{
			++str;
			++nProcess;
			v += *str; // skip 2nd quote
				
			++str;
			++nProcess;
			if ((*str) != '\'')
				{ // get the next char after it
				v += *str;
				}
			}
		} while ( ((*str) != quote) && (*str != 0) ); 
	
	if (*str == quote)
		{
		++str;	// skip ending quote
		++nProcess;
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPTree :: ReadComment(LPCSTR &v, string &comment)
	{
	// skip blanks
	while (*v == ' ')
		{
		++v;
		++nProcess;
		}

	if (*v == '[')
		{
		++v; // skip opening
		++nProcess;
		while ((*v != 0) && (*v != ']'))
			{
			if (*v == '\"') 
				GetQuotedString(v, comment);
			else
				{
				comment += *v;
				++v;
				++nProcess;
				}
			}
		if (*v != 0)
			{
			++v; // skip ending
			++nProcess;
			}
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPTree :: ReadToken(LPCSTR &str, string &v)
	{
	LPCSTR		delim = " ,:[]()";

	while ((*str != 0) && (strchr(delim, *str) == NULL) )
		{
		if ((*str == '\'') || (*str == '\"'))
			{
			GetQuotedString(str, v);
			}
		else
			{
			v += *str;
			++str;
			++nProcess;
			}
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

string& CPTree :: QuoteString(string& comment, char quote)
	{
///	comment.Replace("\'", "\'\'");
///	comment.Replace("\"", "\\\"");

	if (!comment.empty())
		{
		comment.insert(comment.begin(), quote);
		comment += quote;
		}

	return comment;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPTree::Parse(LPCSTR str)
	{
    BOOL        done = FALSE;
    // strip initial arb comment
    while (!done && (*str != 0))
        {
        switch (*str)
            {
            default: done = TRUE;
            case ' ': 
                {
                ++str;
                }
                break;
                
            case '[':
                {
                while ((*str != ']') && (str != 0))
                    ++str;
                }
                break;
			}
        }
        
	totalProcess = strlen(str);
	nProcess     = 0;
	endParseStr = str + totalProcess;

	root = ParseTree(str);
	if ((nProcess < totalProcess) && (*str == ',')) // more to process
		{
		CPNode		*child = root;
		root = new CPNode(this);
		
		int			prev_processed = 0;
		while ((nProcess < totalProcess) && (child != NULL))
			{
			if (*str == ',') ++str;
			
			root->AddEdge(child);
			//DisplayL("Processed %d chars of %d, [%-100.100s]\n", nProcess, totalProcess, str);
			prev_processed = nProcess;
			child = ParseTree(str);
			if (nProcess == prev_processed)
				{
				DisplayL("Could not process all the chars in NEWICK file\n");
				break;
				}
			}
		}
	DisplayL("Processed %d chars of %d\n", nProcess, totalProcess);
	}


///////////////////////////////////////////////////////////////////////////////////////////////////

CPNode* CPTree::ParseTree(LPCSTR &str)
	{
	while (*str != 0)
		{
		//TRACE("[%d] [%s]\n", __LINE__, str);
		
		nProcess = totalProcess - (endParseStr - str); // chars processed is total chars - chars left to process
		switch (*str)
			{
			default:	{
						//   read next token
						string		v;
						v.clear();
						while ((*str != 0) && (*str != ',') && (*str != ')') && (*str != ';'))
							{
							if ((*str == '\'') || (*str == '\"'))
								{
								GetQuotedString(str, v);
								}
							else
								{
								v += *str;
								++str;
	//							++nProcess;
								}
							}
				//		TRACE("[%d] Create new node [%08X][%s]\n", __LINE__, this, v.c_str());
						return new CPNode(this, v.c_str());
						}

			case ' ':	++str;
	//					++nProcess;
						break;

			case ':':   {
				//		TRACE("[%d] [:]\n", __LINE__);
						// find attrribute data upto ',' or ')'
							string		v;
							v.clear();
							while ((*str != 0) && (*str != ',') && (*str != ')') && (*str != ';'))
								{
								if ((*str == '\'') || (*str == '\"'))
									GetQuotedString(str, v);
								else
									{
									v += *str;
									++str;
	//								++nProcess;
									}
								}
						}
						break;
						
			case ';':   { // end of the data
				//		TRACE("[%d] [;]\n", __LINE__);
						nProcess = totalProcess;
						str = endParseStr;
						}
						break;

			case '[':   
						{
				//		TRACE("[%d] [\\[]\n", __LINE__);
						string		comment;
						ReadComment(str, comment);
						if (comment[0] == '{')
							{
							//add attrs to tree
							ProcessAttrList(comment.c_str(), NULL);
							}
						else
							{
							// strip quotes around whole comment
							if ((comment[0] == '\"') && (comment[comment.length()-1] == '\"'))
								comment = comment.substr(1,comment.length()-2); 

							SetAttr(ATTR_COMMENT, comment.c_str());
							}
						}
						break;

			case '(':	{
				//		TRACE("[%d] [%s]\n", __LINE__, str);
						BOOL	complete = FALSE;
						CPNode	*R = new CPNode(this);
						++str;
	//					++nProcess;

						while (!complete)
							{
							CPNode		*S = ParseTree(str);

							R->AddEdge(S);

							if ((*str == ')') || (*str == 0))
								complete = TRUE; 

							++str;
	//						++nProcess;
							}

						// find attrribute data upto ',' or ')'
						string		v;
						v.clear();
						while ((*str != 0) && (*str != ',') && (*str != ')') && (*str != ';'))
							{
							if ((*str == '\'') || (*str == '\"'))
								GetQuotedString(str, v);
							else
								{
								v += *str;
								++str;
	//							++nProcess;
								}
							}

				//		TRACE("[%d] ATTRS[%s]\n", __LINE__, v.c_str());
						R->AddAttrs(v.c_str());

						return R;
						}
						
			case ')':	{
						// This is for the root
				//		TRACE("[%d] Close[%s]\n", __LINE__, str);
						++str;

						// find attrribute data upto ',' or ')'
						string		v;
						v.clear();
						while ((*str != 0) && (*str != ',') && (*str != ')') && (*str != ';'))
							{
							if ((*str == '\'') || (*str == '\"'))
								GetQuotedString(str, v);
							else
								{
								v += *str;
								++str;
	//							++nProcess;
								}
							}

				//		TRACE("[%d] CLOSE ATTRS[%s]\n", __LINE__, v.c_str());
						root->AddAttrs(v.c_str());
						}
			}
		}
	return NULL;
	}


///////////////////////////////////////////////////////////////////////////////////////////////////

void CPTree :: SetAttr(LPCSTR k, LPCSTR v)
	{
	attrs.m[k] = v;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPTree :: ProcessAttrList(LPCSTR str, CPNode *node)
	{
	int		pos;

	while (*str != '\0')
		{
		switch (*str)
			{
			case '{':	{
						// add new attribute
						pos = strcspn(++str, "=}");
				//		++nProcess;

						string		key(str, pos);
						str += pos; // skip to next char
				//		nProcess += pos;

						pos = strcspn(++str, "}");
				//		++nProcess;

						string		value(str, pos);
						str += pos; // skip to next char
				//		nProcess += pos;

						++str; // skip closing brace
				//		++nProcess;

						//TRACE("[%s] %s = %s\n", title, key, value);
						if (node == NULL)
							SetAttr(key.c_str(), value.c_str());
						else
							node->SetAttr(key.c_str(), value.c_str());
						}
						break;

			default:
						++str;
				//		++nProcess;
						break;
			}
		}
	}

//////////////////////////////////////////////////////////////////////

CPNode* CPTree :: Find(LPCSTR v)
    {
	CPNodeListIter iter = nodeList.begin();
	while (iter != nodeList.end())
		{
        //CPNode		*n = *iter;
        if ((*iter)->title.compare(v) == 0)
            return *iter;
        ++iter;
		}
	return NULL;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

CPNode :: CPNode(CPTree *ptree, LPCSTR v)
	{
	char 		str[2048];
	tree = ptree;

	if (tree->internal%100 == 0)
		{
		string		partial = v;
		partial = partial.substr(0, 100);

///		char		msg[2048];
///		sprintf(msg,"Adding node %d {%s}\n", tree->internal, partial);
///		tree->ShowProgress(msg, tree->nProcess, tree->totalProcess, NULL);
		}

	if ((v == NULL) || (strlen(v) == 0))
		{
		id = ++(tree->internal);
		sprintf(str, "%s_%04d", "I", id);
		title = str;
		}
	else
		{
		id = ++(tree->internal);
		// copy upto the ':' or '['
		string		name;
		name.clear();
		while ((*v != 0) && (*v != ':') && (*v != '['))
			{
			if ((*v == '\'') || (*v == '\"'))
				tree->GetQuotedString(v, name);
			else
				{
				name += *v;
				++v;
				}
			}
		AddAttrs(v);

		if (!name.empty() && (name[0] == '\''))
			{
			// this is from ARB.  Strip the quotes, use chars upto first ',' as name
			string			realName;
			string			comment;

			int first = name.find_first_not_of("\' ");
			int last  = name.find_last_not_of("\' ");
			name = name.substr(first, last-first+1);

			int n = name.find(",");
			realName = name.substr(0,n);
			comment  = name.substr(n+1);
			name = realName;

			if (!comment.empty())
				{
				// strip quotes around whole comment
				if ((comment[0] == '\"') && (comment[comment.length()-1] == '\"'))
					comment = comment.substr(1,comment.length()-2); 
				attrs.Add(ATTR_COMMENT, comment);
				}
			}

		int first = name.find_first_not_of(" ");
		int last  = name.find_last_not_of(" ");
		name = name.substr(first, last-first+1);
		attrs.Add(ATTR_NAME, name);
		if (name.empty())
			{
			sprintf(str, "%s_%04d", "I", ++(tree->internal));
			title = str;
			}
		else 
			title = name;
		}

	parent = NULL;
	
	tree->nodeList.push_back(this);
	}

CPNode :: ~CPNode()
	{
	// find position in tree's node list and remove it
	if (tree != NULL)
		{
		CPNodeList::iterator iter = tree->nodeList.begin();
		
		while ((iter != tree->nodeList.end()) && (*iter != this))
			++iter;
		if (iter != tree->nodeList.end())
			tree->nodeList.erase(iter);
		}
		
	id = -1;
	parent = NULL;
	attrs.Clear();
	edges.clear();
	title = "DELETED";
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPNode :: AddEdge(CPNode *other)
	{
	edges.push_back(other);
	other->parent = this;
//	TRACE("Edge from [%s] to [%s]\n", title, other->title);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPNode :: InvertEdges(BOOL recurse)
	{
	// invert the edge list
	CPNodeList		orig;
	orig.insert(orig.begin(), edges.begin(), edges.end());
	
	edges.clear();
	
	for (CPNodeListIter iter=orig.end() ; iter != orig.begin() ; --iter)
		edges.push_back(*iter);
		
	if (recurse)
		{
		for (CPNodeListIter iter=edges.begin() ; iter != edges.end() ; ++iter)
				(*iter)->InvertEdges(TRUE);
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPNode :: ResetParent(CPNode *p)
	{
	// traverse the tree setting parent to correct items
	// used by tree inversion
	parent = p;

	for (CPNodeListIter iter=edges.begin() ; iter != edges.end() ; ++iter)
			(*iter)->ResetParent(this);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

string& CPNode :: CleanString(string& s, BOOL searchForBootstrap)
	{
	if (s.empty()) return s;

	// strip quotes, convert escaped chars
	int			pos;
	
	do {
		pos = s.find("\\\"");
		if (pos != -1)
			s.replace(pos, 2, "\""); // backslash double quote to double quote
		} while (pos >= 0);

	do {
		pos = s.find("\\\\");
		if (pos != -1)
			s.replace(pos, 2, "\\"); // double backslash to single backslash
		} while (pos >= 0);
		
	do {
		pos = s.find("\'\'");
		if (pos != -1)
			s.replace(pos, 2, ""); // two single quotes to nothing
		} while (pos >= 0);
		
	do {
		pos = s.find("\\\'");
		if (pos != -1)
			s.replace(pos, 2, "");	// backslash single quote to nothing
		} while (pos >= 0);
		
	// Remove beginning and ending double quotes
	if (s[0] == '\"')
		s = s.substr(1);
	if (s[s.length()-1] == '\"')
		s = s.substr(0, s.length()-1);

	// Remove beginning and ending single quotes
	if (s[0] == '\'')
		s = s.substr(1);
	if (s[s.length()-1] == '\'')
		s = s.substr(0, s.length()-1);

	if (searchForBootstrap)
		{
		int			pos = s.find_first_not_of("0123456789:", 0);
		string 		boot = s.substr(0,pos);
		int			n = boot.find(":");
		
		if (s.length() == boot.length())
			n = s.length(); // only digits, no colon

		if (n > 1) 
			{
			// bootstrap data on leading edge
			boot = boot.substr(0,n);
			attrs.Add(ATTR_BOOTSTRAP_DATA, boot);
			if (s.length() == n)
				s.clear();
			else
				s = s.substr(n+1);
			}
		}

	return s;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////

void CPNode :: SetAttr(LPCSTR k, LPCSTR v)
    {
	attrs.Add(k,v);
    }

void CPNode :: SetAttr(LPCSTR k, string& v)
    {
	attrs.Add(k,v);
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

void CPNode :: AddAttrs(LPCSTR v)
	{
	while (*v != 0)
		{
		if (*v == ':')
			{
			++v; // consume the ':'
			int			i = strcspn(v, "[");
			string		dist(v, i);
			v += i; 

			Trim(dist, " \n\r\t");
			attrs.Add(ATTR_DISTANCE, dist);
			}
		else if (*v == '[')
			{
			++v; // skip '['
			int			i = strcspn(v, "]");
			string		comment(v,i);
			v += i; 
			if (*v == ']')
				++v;

			if (comment.length() > 0)
				{
				// strip quotes around whole comment
				if ((comment[0] == '\"') && (comment[comment.length()-1] == '\"'))
					comment = comment.substr(1,comment.length()-2); 

				if ((comment[0] == '{') && (comment[comment.length()-1] == '}'))
					tree->ProcessAttrList(comment.c_str(), this);
				else if (IsBranchNode())
					{
					attrs.Add(ATTR_TAXLABEL, CleanString(comment, TRUE));
					}
				else
					attrs.Add(ATTR_COMMENT, CleanString(comment));
				}
			}
		else if ((*v == '\'') || (*v == '\"'))
			{
			string		comment;
			tree->GetQuotedString(v, comment);

			// strip quotes around whole comment
			if ((comment[0] == '\"') && (comment[comment.length()-1] == '\"'))
				comment = comment.substr(1,comment.length()-2); 

			if (comment.length() > 0)
				{
				// strip quotes around whole comment
				if ((comment[0] == '\"') && (comment[comment.length()-1] == '\"'))
					comment = comment.substr(1,comment.length()-2); 

				if ((comment[0] == '{') && (comment[comment.length()-1] == '}'))
					tree->ProcessAttrList(comment.c_str(), this);
				else if (IsBranchNode())
					attrs.Add(ATTR_TAXLABEL, CleanString(comment, TRUE));
				else
					attrs.Add(ATTR_COMMENT, CleanString(comment));
				}
			}
		else if (IsBranchNode())
			{
			// read title of this item
			while (isspace(*v)) ++v;		// skip white space
			int			i = strcspn(v, "[;:),\'\"");
			string		iname(v, i);
			v += i; 

			if ((iname.length() > 1) && (iname[0] == 'I') && (iname[1] == '_'))
				title = iname;
			else if (iname.length() > 0)
				attrs.Add(ATTR_TAXLABEL, CleanString(iname, TRUE));
			}
		else
			++v; // consume the chars up to [:
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CPNode :: WriteAttrs(FILE *f)
	{
	// for each attr
	//		add key,value to file

	// for each edge
	//		edge->ReadAttrs(filename);

	fprintf(f, "[%s]\n", title.c_str());

	CAttrDictIter iter;
	
	for (iter=attrs.m.begin() ; iter != attrs.m.end() ; ++iter)
		{
		if (((*iter).first.length() > 2) && ((*iter).first[0] == '-') && ((*iter).first[1] == '-')) continue;

		if ((*iter).second.length() > 0)
			fprintf(f, "%s=\"%s\"\n", (*iter).first.c_str(), (*iter).second.c_str());
		}

	for (CPNodeListIter iter=edges.begin() ; iter != edges.end() ; ++iter)
		(*iter)->WriteAttrs(f);

	return TRUE;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CPNode :: WriteNode(FILE *f, int indent, BOOL includeAttrs, BOOL withComment)
	{
	char		str[2048];

	if (f == NULL)
		return FALSE;

	if (IsBranchNode())
		{
		fprintf(f, "(");
		indent += 1;
		fprintf(f, "\n%*s", indent, "");
		for (CPNodeListIter iter=edges.begin() ; iter != edges.end() ; ++iter)
			{
			if (iter != edges.begin())
				{
				fprintf(f, ",");
				fprintf(f, "\n%*s", indent, "");
				}
			(*iter)->WriteNode(f, indent, includeAttrs, withComment);
			}
		indent -= 1;
		fprintf(f, "\n%*s", indent, "");
		fprintf(f, ")");
		}

	string		name = attrs.GetString(ATTR_NAME,"");

	if (includeAttrs)
		{
		string		boot = attrs.GetString(ATTR_BOOTSTRAP_DATA);

		if (IsBranchNode())
			{
			name = attrs.GetString(ATTR_TAXLABEL);
			if (!boot.empty())
				{
				// prepend the boot to the name, place in quotes
				if (name.empty())
					name = boot;
				else
					{
					sprintf(str, "\'%s:%s\'", boot.c_str(), name.c_str());
					name = str;
					}
				}

			fprintf(f, "%s", name.c_str()); 	
			}
		else
			{
			fprintf(f, "%s", name.c_str()); 	

			if (!boot.empty())
				fprintf(f, "[%s]", boot.c_str()); 	
			}

		string		comment;
		GetComment(ATTR_COMMENT, comment);
		tree->QuoteString(comment);

		if (!comment.empty())
			fprintf(f, " [%s]", comment.c_str()); 	

		string		key;
		string		value;
		string		attrStr;

		CAttrDictIter iter;
		for (iter=attrs.m.begin() ; iter != attrs.m.end() ; ++iter)
			{
			key = (*iter).first;
			value = (*iter).second;
			if ((key.length() > 2) && (key[0] == '-') && (key[1] == '-')) continue;
			//ignore attributes printed elsewhere
			if (key.compare(ATTR_NAME) == 0)			continue;
			if (key.compare(ATTR_DISTANCE) == 0)		continue;
			if (key.compare(ATTR_TAXLABEL) == 0)		continue;
			if (key.compare(ATTR_COMMENT) == 0)			continue;
			if (key.compare(ATTR_BOOTSTRAP_DATA) == 0)	continue;

			if (attrStr.empty())
				attrStr += '[';

			sprintf(str, "{%s=%s}", key.c_str(), value.c_str());
			attrStr += str;
			}

		if (!attrStr.empty())
			{
			attrStr += ']';
			}
		}
	else if (withComment)
		{
		string		comment;
		GetComment(ATTR_COMMENT, comment);
		if (!comment.empty())
			{
			name += ", ";
			name += comment;
			tree->QuoteString(name, '\'');
			}
		fprintf(f, "%s", name.c_str()); 	
		}
	else
		fprintf(f, "%s", name.c_str()); 	


	double		dist = 	attrs.GetDouble(ATTR_DISTANCE, 0);
	fprintf(f, ":%.*f", CPTree::PRECISION, dist); 	

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CPNode :: GetComment(LPCSTR  attr, string &str, LPCSTR _default)
	{
	str = attrs.GetString(attr, _default);
	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CPNode :: GetAncestors(CPNodeList &list, BOOL fromBottom)
	{
	// return the names of the ancestors
	list.clear();
	int		n = 0;

	CPNode		*p = parent;
	while (p != NULL)
		{
        if (fromBottom)
            list.push_back(p);
        else
            list.insert(list.begin(), p);
            
		p = p->parent;
		++n;
		}

	return n;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CPNode :: GetDecendants(CPNodeList &list, BOOL leafOnly)
	{
	// return the names of the decendants

	// for each edge
	//	add edge name
	//	add all decendants of edge
	
	for (CPNodeListIter iter=edges.begin() ; iter != edges.end() ; ++iter)
		{
		CPNode *node = (*iter);
		
		if (node->IsBranchNode())
			{
			if (!leafOnly)
				list.push_back(node);
			node->GetDecendants(list, leafOnly);
			}
		else
			list.push_back(node);
		}

	return list.size();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

string& MakeCSVQuoted(string &q_str, LPCSTR others)
	{
	// Escape the following chars
	//		" double quote

	string		local = q_str;
	LPCSTR		s = local.data();

	q_str = "\"";

	while (*s != 0)
		{
		if (*s == '\"')
			q_str += "\"\"";
		else 
			{
//			if (others && (strchr(others, *s) != NULL))
//				q_str += "\\"; // escape the char
				
			q_str += *s;
			}
		++s;
		}

	q_str += "\"";

	return q_str;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

string& MakeQuoted(string &q_str, LPCSTR others)
	{
	// Escape the following chars
	//		/ backslash
	//		" double quote
	//		< open 
	//		> close

	string		local = q_str;
	LPCSTR		s = local.data();

	q_str = "\"";

	while (*s != 0)
		{
		if (*s == '\\')
			q_str += "\\\\";
		else if (*s == '\"')
			q_str += "\\\"";
		else if (*s == '<')
			q_str += "\\<";
		else if (*s == '>')
			q_str += "\\>";
		else 
			{
			if (others && (strchr(others, *s) != NULL))
				q_str += "\\"; // escape the char
				
			q_str += *s;
			}
		++s;
		}

	q_str += "\"";

	return q_str;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

