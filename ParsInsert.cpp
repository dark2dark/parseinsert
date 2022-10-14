/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : ParsInsert.cpp
//  Purpose   : Parsimonious Insertion of unclassified seqeunces into Phylogenetic Tree
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
//
// Sample command for parsimonious insertion:
//		-x rdp.taxonomy -t core_set.tree -s core_set_aligned.fasta unclassified_seq.fasta  
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "Knox_Stddef.h"
#include "ParsInsert.h"
#include "PNode.h"
#include "SeqList.h"
#include "Taxonomy.h"
#include "Attrs.h"

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////

#define APPNAME						"ParsInsert"
#define VERSION						APPNAME " Version 1.04 " __DATE__ " " __TIME__
#define COPYRIGHT \
"ParsInsert: Parsimonious Insertion of unclassified sequences into phylogenetic tree\n"\
"            Copyright (C) 2007-2011 David Knox\n"\
"            This program comes with ABSOLUTELY NO WARRANTY\n"\
"            This is free software, and you are welcome to redistribute under certain conditions - see license."

#define USAGE	"Application:%s\n\n%s [options]  <insert_sequences>      \n"    \
				"      Parsimonious Insertion of Sequences into Given Tree                  \n\n"  \
				"     -m <mask file>           - read mask from this file                   \n"    \
                "     -s <tree sequences>      - read core tree sequences from this file    \n"    \
                "                                (default: PI_Tree.fasta)                   \n"    \
                "     -t <tree file>           - read core tree from this file              \n"    \
                "                                (default: PI_Tree.tree)                    \n"    \
                "     -x <tree taxonomy>       - read core tree taxomony from this file     \n"    \
                "     -o <output file>         - output taxonomy for each insert sequence to this file\n"    \
                "                                (default: PI_Results.log)                  \n"    \
				"     -l[-|<log file>]         - create log file (default is ParsInsert.log)\n"    \
				"     -n#                      - number of best matches to display          \n"    \
				"     -c#                      - percent threshold cutoff                   \n"    \
				"     -p                       - print node comments in newick file         \n"    \
				"     -D#                      - print branch lengths using # decimal places\n"    \
				"\n"  
//              "     -p <parsimony file>      - read parsimony values from this file       \n"    

///////////////////////////////////////////////////////////////////////////////////////////////////

FILE				*newseqs      = NULL;			// sequence file of unclassified insertions

BOOL				createLog     = TRUE;
int					fast		  = 0;

int                 nCalcCost     = 0;				// counters for cost calc stats
int                 nCalcCostFull = 0;				// counters for cost calc stats
BOOL                useCalcCost   = 4;				// calculation method to be used

char                *statsName      = "PI_Results.log";
char                *logFilename    = "ParsInsert.log";
char                *treeFilename   = "PI_Tree.tree";
char                *seqFilename    = "PI_Tree.fasta";
char                *fullTreeName   = NULL;

char                *files[64]      = {64*0};
int                 nfiles          = 0;

BOOL                precision       = FALSE;
int					scoreThresh 	= 80; 			// default to 80% or better to add item	

///////////////////////////////////////////////////////////////////////////////////////////////////

FILE				*stats          = NULL;
time_t              start           = 0;

long long			allocatedMemory = 0;

CSequenceFile		*treeSeqfile   = NULL;
CSequenceFile		*insertSeqfile = NULL;
CPTree              *fullTree      = NULL;

// Variables used during the parsimony calculation
//
volatile int		nSeqRemaining   = 0;    // used to show progress during parsimony calculation
CParsimonyList      parsimonyList;          // parsimonious sequence for common ancestors
CParsimonyList      unionList;              // union of all possible nucleotides at each position
CSequenceFile		*taxonomy       = NULL;
CPTree              *tree           = NULL;
char                *taxName        = NULL;
char                *maskName       = NULL;
char                *parsName       = NULL;

map<string,CPNode*>	coreNames;

CInsertPosArray		inserts;

char				mask[MAX_SEQ_SIZE];
bool				useMask         = FALSE;

int					testing         = 1;
BOOL                fitch           = TRUE;
BOOL                jukes_cantor    = TRUE;
BOOL				withComments	= FALSE;

///////////////////////////////////////////////////////////////////////////////////////////////////

enum				{RANK_DOMAIN, RANK_PHYLUM, RANK_CLASS, RANK_ORDER, RANK_FAMILY, RANK_GENUS, RANK_SPECIES, RANK_N};
int					rankCounts[RANK_N][4];
LPCSTR				rankNames[] = {	"Domain", 
									"Phylum", 
									"Class", 
                                    "Order",
									"Family", 
									"Genus", 
									"Species", 
									NULL};

class CRankCounts
    {
    public:
        int					rankCounts[RANK_N][4];
    };
    
typedef map<string,CRankCounts*>	CMapRankCounts;
CMapRankCounts      				phylumRankCounts;

///////////////////////////////////////////////////////////////////////////////////////////////////

typedef	vector<CInsertPos*>				CTreeInsertList;
typedef map<CPNode*,CTreeInsertList*>	CTreeInserts;

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceFile*	ReadSequenceFile(LPCSTR filename);
CSequenceFile*	ReadTaxonomyFile(LPCSTR filename);
LPCSTR			GetTaxonomy(LPCSTR name, CPNode *node=NULL);

CPTree*			ReadNewickTree(LPCSTR filename);
int				LoadParsimonyIndex();			// loads a parsimony index
int				LoadSequenceIndex();			// loads a sequence index
int				BuildSequenceIndex();			// build seq index from sequence file

BOOL			RemoveUnwantedNodes(CPNode *node, CSequenceFile *seqfile);
BOOL			SetInternalTaxonomy(CPTree *tree);
BOOL            SetMyTaxonomy(CPNode *node);
BOOL            CheckLeafTaxonomy(CPTree *tree);

int             CheckTreeSequences(CPTree *tree, CSequenceFile *treeSeqfile);
int             ParsimonyInsertion();

int             CreateParsimonySets(LPCSTR parsName, CPTree *tree);
int				CalcAncestorParsimony(CPTree *tree);	// Build parsimony values (write to file)
int				CalcNodeParsimony(CPNode *node);
int             ForceNodeParsimony(CPNode *node);
int             WriteParsimonySet(CPNode *node, CSequenceFile *parsimony);

void			FindInsertLocations(CSequenceFile *insertSeqfile);
void			FindInsertLocations(CSequenceFile *treeSeqfile, CPTree *tree, CSequenceFile *insertSeqfile, LPCSTR seqname);
void			FindInsertLocations(CSequenceFile *treeSeqfile, CPTree *tree, CSequenceFile *insertSeqfile, CInsertPosArray &inserts);
void			FindInsertLocations_Serial(CSequenceFile *treeSeqfile, CPTree *tree, CSequenceFile *insertSeqfile, CInsertPosArray &inserts);

int				CalcCost(CParsimonySet &next_pars, CParsimonySet &node_pars, int &diffs, int &partials, int &indel, int thresh=INT_MAX);
int				CalcCost_NEW4(CParsimonySet &next_pars, CParsimonySet &node_pars, int &diffs, int &partials, int &indel, int thresh=INT_MAX);
int				CalcCost_ORIG(CParsimonySet &next_pars, CParsimonySet &node_pars, int &diffs, int &partials, int &indel, int thresh=INT_MAX);
void            BuildScoreMatrix();

void			DisplayInsertLocations(CInsertPosArray &inserts, BOOL all);
void            DisplayBestInsertLocations(CInsertPos *item, BOOL all);
void			OutputInsertLocations(CInsertPosArray &inserts, BOOL all);
void 			DisplayRankCounts(LPCSTR header, int counts[][4]);
void			ReleaseInserts(CInsertPosArray &inserts);
CParsimonySet* 	GetParsimonySet(CPNode *node, BOOL create=FALSE);
int 			GetTaxonomyMatch(string &taxA, string &taxB, string &marker, int counts[][4]);
void            GetRankFromTax(int rank, string &tax, string &rankStr, BOOL only=TRUE);

CTreeInserts*	BuildTreeInsertLists(CInsertPosArray &inserts);
BOOL			BuildInsertTree(CTreeInserts* insertLists, CPTree *tree);
CPNode* 		AddChild(CPTree *tree, CPNode *node, CInsertPos *in);
double 			CalcDist(CPNode *parent, CInsertPos* in, int &diffs, int &parts, int &indel);

BOOL			BuildInsertTree(CInsertPosArray &inserts, CPTree *tree);
BOOL			AddTaxonomyToTree(CPTree *tree, CSequenceFile *taxonomy);

BOOL			ReadMask(LPCSTR filename);

void			CompareTaxonomy(LPCSTR tax1, LPCSTR tax2);

void            TestInsertPositions();
BOOL            FullTree_Initialization();
double          ComparePosition(LPCSTR name, CPNode *parent, CPTree *fullTree, CPTree *coreTree);

///////////////////////////////////////////////////////////////////////////////////////////////////

inline LPCSTR GetNodeName(CPNode *node)
    {
    string		name = node->attrs.GetString(ATTR_NAME,"");
    
    if (name.empty()) 
        name = node->title;
    
    return name.c_str();
    }


///////////////////////////////////////////////////////////////////////////////////////////////////
// The one and only application object

int main(int argc, char* argv[])
	{
	int nRetCode = 0;

	memset(rankCounts, 0, sizeof(rankCounts));

	DisplayL("%s\n", COPYRIGHT);

	int i = 1;
	while (i < argc)
		{
        if (argv[i][0] == '-')
            {
            // Process argument
            switch (argv[i][1])
                {
                case 'c': // cutoff threshold percent
                    {
                    int v;
                    if (argv[i][2] != 0)
                        v = atoi(&argv[i][2]);
                    else if (i < argc)
                        v = atoi(argv[++i]);
                    scoreThresh = v;
                    }
                    break;
                        
                case 'd': // distance calc method
                    {
                    int v;
                    if (argv[i][2] != 0)
                        v = atoi(&argv[i][2]);
                    else if (i < argc)
                        v = atoi(argv[++i]);
                    useCalcCost = v;
                    }
                    break;
                        
                case 'D': // distance calc method
                    {
                    int v;
                    if (argv[i][2] != 0)
                        v = atoi(&argv[i][2]);
                    else if (i < argc)
                        v = atoi(argv[++i]);
                    if (v >= 0)
                    	{
                    	DisplayL("Setting Tree Floating Point Precision to %d\n", v);
                    	CPTree::PRECISION = v;
                    	}
                    }
                    break;
                        
                case 'f':
                    {
                    int v = atoi(&argv[i][2]);
                    fast = v;
                    }
                    break;
                    
                case 'F':
                    {
                    fitch = (argv[i][2] != '-');
                    }
                    break;

                case 'j':
                    {
                    jukes_cantor = (argv[i][2] != '-');
                    }
                    break;
                        
                case 'l':
                    createLog = (argv[i][2] != '-');
                    if (argv[i][2] != 0)
                        logFilename = &argv[i][2];
                    else if ((i < argc) && (argv[i][2] != '-'))
                        logFilename = argv[++i];
                    break;

                case 'n': // number of matches to display
                    {
                    int v = atoi(&argv[i][2]);
                    if (v > 0)
                        CBestLocation::default_len = v;
                    }
                    break;
                    
                case 'm': // mask input file 
                    {
                    if (argv[i][2] != 0)
                        maskName = &argv[i][2];
                    else if (i < argc)
                        maskName = argv[++i];
                    }
                    break;
                    
                case 'o': // output file 
                    {
                    if (argv[i][2] != 0)
                       statsName = &argv[i][2];
                    else if (i < argc)
                       statsName = argv[++i];
                    }
                    break;
                        
                        
                case 'p':
                    {
                    withComments = (argv[i][2] != '-');
                    }
                    break;

                case 'P':
                    {
                    precision = (argv[i][2] != '-');
                    }
                    break;

                case 's': // tree sequence file 
                    {
                    if (argv[i][2] != 0)
                       seqFilename = &argv[i][2];
                    else if (i < argc)
                       seqFilename = argv[++i];
                    }
                    break;
                        
                case 't': // tree file 
                    {
                    if (argv[i][2] != 0)
                       treeFilename = &argv[i][2];
                    else if (i < argc)
                       treeFilename = argv[++i];
                    }
                    break;
                                                
                case 'x': // taxonomy input file 
                    {
                    if (argv[i][2] != 0)
                        taxName = &argv[i][2];
                    else if (i < argc)
                        taxName = argv[++i];
                    }
                    break;
                       
                default:
                    DisplayL("invalid command option '%c'\n", argv[i][1]);
                    nRetCode = RC_ERROR_OPTIONS;
                    break;
                }
            }
        else
            {
            files[nfiles++] = argv[i];
            }
        i++;
        }

	if (createLog)
		OpenAppLog(logFilename, "a+");

	// Write out our version and the command line
	string			cmd;
	for (int k=0 ; k < argc ; ++k)
		{
		if (!cmd.empty())
			cmd += " ";
		cmd += argv[k];
		}
	DisplayL("%s\nCommand Line:[%s]\n", VERSION, cmd.c_str());
    
    if (nfiles < 1)
        {
        DisplayL("ERROR - missing files to process\n");
        DisplayT(USAGE, VERSION, APPNAME);
        nRetCode = RC_ERROR_CMDLINE;
        }

	start = time(NULL);

	if (treeFilename != NULL)
		{
		DisplayL(" tree file: %s\n", treeFilename);
		}
		
	if (seqFilename != NULL)
		{
		DisplayL(" tree sequence file: %s\n", seqFilename);
		}
		
	if (maskName != NULL)
		{
		DisplayL(" mask file: %s\n", maskName);
		
		if (!ReadMask(maskName))
			{
			DisplayT("Cannot open mask file [%s]", maskName);
			nRetCode = RC_ERROR_MASKFILE;
			}

		useMask = TRUE;
		}

	if (taxName != NULL)
		{
		DisplayL(" taxomony index file: %s\n", taxName);
		taxonomy = ReadTaxonomyFile(taxName);
		if (taxonomy->f == NULL)
			{
			DisplayT("Cannot open taxomony index FASTA file [%s]", taxName);
			nRetCode = RC_ERROR_TAXFILE;
			}
        else
            DisplayL(" taxomony index file: %d taxa\n", taxonomy->GetCount());
		}

    BuildScoreMatrix();
        
    if (nRetCode == 0)
        {
        nRetCode = ParsimonyInsertion();
        }
        
	// release the data
        
    if (treeSeqfile != NULL)
        delete treeSeqfile;
    
    if (insertSeqfile != NULL)
        delete insertSeqfile;
    
    if (tree != NULL)
        delete tree;
    
    if (stats != NULL)
        fclose(stats);
        
	if (taxonomy != NULL)
		delete taxonomy;

///	ReleaseInserts(inserts);

	DisplayL("Process Completed: %d sec\n", time(NULL)-start);
	CloseAppLog();

	return nRetCode;
	}


///////////////////////////////////////////////////////////////////////////////////////////////////

int ParsimonyInsertion()
    {
    int             ret = 0;
    time_t			stage;

    // Algorithm
    /*
     1.	If sequence index available
     load sequence index
     else
     create sequence index from sequence FASTA file
     
     2.	Read Tree file
     
     3.	If parsimony index available
     load parsimony index
     else if parsimony FASTA file availble
     create parsimony index from parsimony FASTA file
     else
     Create parsimony FASTA file and index file
     
     4.	If new sequence index available
     Load new sequence index
     else 
     Create new sequence index from new sequence FASTA file
     
     5.	for each group of new sequences that fit into memory
     for each node in tree
     for each seq in group
     calc score to convert node into seq
     Add score,node pair to best for seq
     Keep top K scores
     for each seq in group
     write top scoring positions and Taxonomy
     */

    stats  = fopen(statsName, "w");

    // Read Tree File
    DisplayL("Reading Newick Tree: [%s]\n", treeFilename);
    stage = time(NULL);
    tree = ReadNewickTree(treeFilename);
    if (tree == NULL)
        {
        DisplayL("Error accessing Newick Tree file [%s]\n", treeFilename);
        return RC_ERROR_TREEFILE;
        }
    else
        DisplayL("Newick Tree Read Completed: %d taxa, %d sec\n", (tree->nodeList.size()+1)/2, time(NULL)-stage);

	// Build list of core leaf names
	CPNodeList		&coreNodes = tree->nodeList;
	for (CPNodeListIter iter=coreNodes.begin() ; iter != coreNodes.end() ; ++iter)
		{
		CPNode		*node = *iter;
		string 		name = GetNodeName(node);
		if (!node->IsBranchNode())
			{
			if (coreNames.find(name) != coreNames.end())
				DisplayL("Found duplicate name: [%s]/n", name.c_str());
			coreNames[name] = node;
			}
		}
		
    // Read Sequence File
    DisplayL("Tree Sequence Reading: [%s]\n", seqFilename);
    stage = time(NULL);
    treeSeqfile = ReadSequenceFile(seqFilename);
    ret = CheckTreeSequences(tree, treeSeqfile);
    if (ret != 0)
        return ret;
    DisplayL("Tree Sequence Read Completed: %d sequences, %d sec\n", treeSeqfile->seqlist.size(), time(NULL)-stage);
    
    // read sequences to be inserted
    DisplayL("Reading sequences for insertion: [%s]\n", files[0]);
    stage = time(NULL);
    insertSeqfile = new CSequenceFile(files[0], 0);
    insertSeqfile->ReadSequenceFile(1);
    if (insertSeqfile->GetCount() == 0)
        {
        DisplayL("No insertion sequences found\n");
        return RC_ERROR_NO_INSERT_SEQ;
        }
    DisplayL("Insert Sequecne Read Completed: %d sequences, %d sec\n", insertSeqfile->GetCount(), time(NULL)-stage);
        
    ///			Only calc the ones actually needed
    DisplayL("Checking Taxa Taxonomy\n");
    CheckLeafTaxonomy(tree);
    DisplayL("Setting Internal Taxonomy\n");
    stage = time(NULL);
    SetInternalTaxonomy(tree);
    DisplayL("Internal Taxonomy Completed: %d sec\n", time(NULL)-stage);
    
        // Set parsimony either from file or building from sequences
    DisplayL("Setting internal parsimony sequences ...\n");
    stage = time(NULL);
    CreateParsimonySets(parsName, tree);
    DisplayL("Parsimony Calculation Completed: %d sec\n", time(NULL)-stage);

    if (fullTreeName != NULL)
        FullTree_Initialization();

    DisplayL("[%d sec elapsed] Searching for insertion points:\n", time(NULL)-start);
    //DisplayL("[%d fast]\n", fast);
    if (useCalcCost == 4)
        {
        CPNodeList		&nodes = tree->nodeList;
        for (CPNodeListIter iter=nodes.begin() ; iter != nodes.end() ; ++iter)
            {
            CPNode          *node = *iter;
            CParsimonySet	*node_pars  = GetParsimonySet(node);
            
            if (node_pars->segCounts == NULL)
                node_pars->BuildSegCounts((useMask ? mask : NULL));
            }
        }
    stage = time(NULL);
    switch (fast)
        {
        case 0:
            {
            FindInsertLocations_Serial(treeSeqfile, tree, insertSeqfile, inserts);
            }
            break;
        default:
            {
            FindInsertLocations(treeSeqfile, tree, insertSeqfile, inserts);
            }
            break;
        }                   
    DisplayL("[%d sec elapsed] Insert Locations Completed: %d sec\n", time(NULL)-start, time(NULL)-stage);
    time_t  insertTime = time(NULL)-stage;

    ///    DisplayInsertLocations(inserts, TRUE);
    memset(rankCounts, 0, sizeof(rankCounts));
    OutputInsertLocations(inserts, FALSE);

    // Place inserts into tree
    AddTaxonomyToTree(tree, taxonomy);
    //tree->WriteNewickTree("AfterAddTax.tree", tree->root, TRUE);

    CTreeInserts		*insertmap = BuildTreeInsertLists(inserts);
    BuildInsertTree(insertmap, tree);

    string			treename = insertSeqfile->fname;
    int				n = treename.find_last_of('.');
    if (n > 0)
        treename = treename.substr(0, n) + ".tree";
    DisplayL("Writing Newick tree file: [%s]\n", treename.c_str());
    tree->WriteNewickTree(treename.c_str(), tree->root, withComments);

    DisplayRankCounts("Rank Matches:", rankCounts);
    
    if (testing > 2)
	    DisplayL("CalcCost calls[%d] = %d / %d (%d)\n", useCalcCost, nCalcCostFull, nCalcCost,nCalcCost-nCalcCostFull);
	if (insertTime > 0)
	    DisplayL("Insert Time = %d, %d per hour\n", insertTime, 3600 * inserts.size() / insertTime );

	if (testing > 2)
		{
		// display the rank counts for each phylum
		CMapRankCounts::iterator   iter;
		for (iter=phylumRankCounts.begin() ; iter != phylumRankCounts.end() ; ++iter)
			{
			string      phy  = iter->first;
			DisplayRankCounts(phy.c_str(), (iter->second)->rankCounts);
			}
		}
		
    return ret;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

int CheckTreeSequences(CPTree *tree, CSequenceFile *treeSeqfile)
    {
    // Check to see that we have a sequence for each leaf of the tree
    int         ret = 0;
    
    if (treeSeqfile->seqlist.size() == 0)
        {
        DisplayL("Error accessing Tree Sequencefile [%s]\n", treeSeqfile->fname.c_str());
        ret = RC_ERROR_TREESEQ;
        }
        
    CPNodeList		&nodes = tree->nodeList;
    for (CPNodeListIter iter=nodes.begin() ; iter != nodes.end() ; ++iter)
        {
        CPNode		*node = *iter;
                
        if (node == NULL) continue;
        if (node->IsBranchNode()) continue;
            
        // look up the sequence name in the sequence file
        string		name = GetNodeName(node);
        
        // get sequence from the sequence file
        CSequenceItem *seq = treeSeqfile->GetSeq(name.c_str());
        if ((seq == NULL) || (seq->len == 0))
            {
            DisplayL("ERROR - tree member sequence is not available [%s]\n", name.c_str());
            ret = RC_ERROR_TREESEQ_MISSING;
            }
        }
           
    return ret;
    }
    
///////////////////////////////////////////////////////////////////////////////////////////////////

int CreateParsimonySets(LPCSTR parsName, CPTree *tree)
    {
    // Create parsimony sets from FASTA file
    //
    if (parsName != NULL)
        {
        CSequenceFile		*parsimony    = NULL;
        
        DisplayL(" parsimony index file: %s\n", parsName);
        parsimony = ReadSequenceFile(parsName);
        if (parsimony->f == NULL)
            {
            DisplayT("Cannot open parsimony index FASTA file [%s]", parsName);
            return -10;
            }
        DisplayL(" parsimony file: %d internal nodes\n", parsimony->GetCount());
        
        // loop creating the parsimony sets for each node (leaf and internal)
        parsimony->ResetSeqIterator();
        CSequenceItem	*seq = parsimony->GetNextSeq();
        while (seq != NULL)
            {                            
            seq->ReadSeqData(parsimony->f);
            
            CParsimonySet		*node_pars = NULL;
            
            node_pars = (CParsimonySet*)parsimonyList[seq->name];
            if (node_pars == NULL)
                {
                // does not exist, create new entry
                //	DisplayL("Create new CParsimonySet for %s\n", name.c_str());
                
                node_pars = new CParsimonySet();
                node_pars->Convert(seq->data);
                char    str[100*1024];
                node_pars->BuildSets(str,sizeof(str));
                node_pars->Trace(seq->name.c_str(), str);
                parsimonyList[seq->name] = node_pars;
                }
            else
                node_pars->Convert(seq->data);
            
            seq->ReleaseSeqData();
            seq = parsimony->GetNextSeq();
            }
        
        if (parsimony != NULL)
            delete parsimony;
        }
    else
        {
        // create parsimony for tree from leaf sequences
        DisplayL("Calculating internal parsimony sequences ...\n");
        CalcAncestorParsimony(tree);
        }
        
    return 0;
    }
            
///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceFile*	ReadTaxonomyFile(LPCSTR filename)
	{
	CSequenceFile	*seqfile = new CSequenceFile();

	string			fname = filename;
	if (fname.substr(fname.length()-4) == ".idx")
		{
		// read this index file
		seqfile->ReadSequenceIndexFile(filename);

		fname = fname.substr(0, fname.length()-4);
		seqfile->Open(fname.c_str(), "r");
		}
	else
		{
		// read the file and create index
		seqfile->Open(filename);
		seqfile->ReadTaxonomyFile();	
		}

	return seqfile;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CSequenceFile*	ReadSequenceFile(LPCSTR filename)
	{
	CSequenceFile	*seqfile = new CSequenceFile();

	string			fname = filename;
	if (fname.substr(fname.length()-4) == ".idx")
		{
		// read this index file
		seqfile->ReadSequenceIndexFile(filename);

		fname = fname.substr(0, fname.length()-4);
		seqfile->Open(fname.c_str(), "r");
		}
	else
		{
		// read the file and create index
		seqfile->Open(filename);
		seqfile->ReadSequenceFile(1);	
		}

	return seqfile;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL ReadMask(LPCSTR filename)
	{
	if (filename == NULL) return FALSE;

	FILE			*f = fopen(filename, "r");
	if (f == NULL) return FALSE;

	memset(mask, 0, sizeof(mask));
	int				m = 0;

	char			buffer[MAX_SEQ_SIZE];
	while (ReadNextLine(buffer, sizeof(buffer), f))
		{
		for (int i=0 ; (i < strlen(buffer)) && (m < sizeof(mask)) ; ++i)
			{
            switch (buffer[i])
                {
                case '0':
                case '1':
                    mask[m] = buffer[i] - '0';
                    break;
                
                case '.':
                case '-':
                    mask[m] = 0;
                    break;
            
                default:
                    mask[m] = 1;
                    break;
                }
            ++m;
			}
		}

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CPTree*	ReadNewickTree(LPCSTR filename)
	{
	TRACE("opening source file: [%s]\n", filename);
	
	FILE				*f = fopen(filename, "r");
	if (f == NULL)	
		{
		Display("Could not open source file: [%s]\n", filename);
		return NULL;
		}

	static		char		newick[50*1024*1024];	// should get the file size for the buffer
	int			len=0;

	CPTree *newickTree = new CPTree();

	DisplayT("Opening file: %s\t...", filename);

	if (f != NULL)
		{
		char buffer[32*1024];
		newickTree->filename = filename;
		//newickTree->ReadTreeAttrs(); // get defaults

		memset(newick, 0, sizeof(newick));
		int			maxlinelen = 0;
		int			linecount  = 0;
		while (ReadNextLine(buffer, sizeof(buffer), f) > 0)
			{
			if (buffer[0] == '#') continue;
			int blen = strlen(buffer);
			// strip the leading whitespace
			int			skip;
			for (skip=0 ; (skip < blen) && isspace(buffer[skip]) ; ++skip)
				; // skip whitespace
			LPCSTR buf = buffer + skip;
			memcpy(&newick[len], buf, blen-skip);
			len += blen-skip;
			//ASSERT(len < 50*1024*1024);
			maxlinelen = max(maxlinelen, blen);
			}

		if (++linecount%100 == 0)
			Display("Opening file: %s\tLines:%d", filename, linecount);

		fclose(f);
		}

	DisplayT("Parsing Newick file: %s\t...", filename);
	newickTree->Parse(newick);

	return newickTree;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL CheckLeafTaxonomy(CPTree *tree)
    {
	// calc the taxonomy of the tree
	if (taxonomy == NULL) return FALSE;
	
	CPNodeList		missing;
    int             n = 0;
	CPNodeList		&nodes = tree->nodeList;
	for (CPNodeListIter iter=nodes.begin() ; iter != nodes.end() ; ++iter)
		{
        CPNode		*node = *iter;
        string      tax;
        
        if (node == NULL) continue;
        if (node->IsBranchNode()) continue;
        
        // check for taxonomy
        CSequenceItem *seq = taxonomy->GetSequence(node->title.c_str());
        if (seq != NULL)
            {
            tax = seq->hdr;
            Trim(tax, "\n\r\t ");
            // find first \t rest of line is taxonomy
            int pos = tax.find_first_of('\t');
            if (pos >= 0)
                tax = tax.substr(pos+1,tax.length());
            seq->ReleaseSeqData();
            }
                  
        if (tax.empty())
            {
            //DisplayL("\tNO TAXONOMY for [%s]\n", node->title.c_str());
            missing.push_back(node);
            ++n;
            }
        }
        
    if (n != 0)
        {
        DisplayL("%d of %d nodes are missing taxonomy\n\t", n, nodes.size());
  
		int			count = 0;
	    for (CPNodeListIter iter=missing.begin() ; iter != missing.end() ; ++iter)
			{
   			CPNode		*node = *iter;
   	   		DisplayL(" %-12.12s%s", node->title.c_str(), (++count%8 == 0 ? "\n\t" :""));
   		    }
        DisplayL("\n");
		}
		
    return (n == 0);
    }
            
///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL SetInternalTaxonomy(CPTree *tree)
	{
	// calc the taxonomy of the tree
	if (taxonomy == NULL) return FALSE;
	
    CTaxEntry		taxEntries("Tree Taxonomy");	// keeps track of diversity seen in taxonomy of tree sequences
    
	nSeqRemaining = tree->nodeList.size()/2;

	CPNodeList		&nodes = tree->nodeList;
	for (CPNodeListIter iter=nodes.begin() ; iter != nodes.end() ; ++iter)
		{
		CPNode		*node = *iter;

		if (node == NULL) continue;
		if (!node->IsBranchNode()) 
            {
            CSequenceItem *seq = taxonomy->GetSequence(node->title.c_str());
                if (seq != NULL)
                    {
                    string tax = seq->hdr;
                    Trim(tax, "\n\r\t ");
                    // find first \t rest of line is taxonomy
                    int pos = tax.find_first_of('\t');
                    if (pos >= 0)
                        tax = tax.substr(pos+1,tax.length());
                    seq->ReleaseSeqData();
                    taxEntries.Add(tax.c_str());	// update diversity
                    }
            continue;
            }

		CPNodeList		descent;
		node->GetDecendants(descent, TRUE);

		map<string,int>				commonTax;
		map<string,int>::iterator	cIter;
		string						common;
		string						tax;
		int							nDescent = 0;

//        DisplayL("Taxonomy for [%s]\n", node->title.c_str());
		for (CPNodeListIter diter=descent.begin() ; diter != descent.end() ; ++diter)
			{
			CPNode		*des = *diter;
			if (des == NULL) continue;

            tax.clear();
			CSequenceItem *seq = taxonomy->GetSequence(des->title.c_str());
			if (seq != NULL)
				{
				tax = seq->hdr;
				Trim(tax, "\n\r\t ");
				// find first \t rest of line is taxonomy
				int pos = tax.find_first_of('\t');
				if (pos >= 0)
					tax = tax.substr(pos+1,tax.length());
				seq->ReleaseSeqData();
                }

			// add each part of the taxonomy to the common list
			int				p;
			string			entry;
			
			entry = tax;
			p = commonTax[entry];
			commonTax[entry] = p+1;
			
			for (int i=0 ; i < tax.length() ; ++i)
				{
				string			entry = tax;
				if ((tax[i] == ';') || (tax[i] == '/'))
					{
					entry = tax.substr(0, i);
					p = commonTax[entry];
					commonTax[entry] = p+1;
					}
				}
				
			if (!tax.empty())
				++nDescent;
			}

		// find the longest, most trusted value
		int			commonThresh = 50 * nDescent / 100;  // 60% of set descendents are required to select value
		int			commonBest = 0;  // best value seen
		common = "";
		for (cIter=commonTax.begin() ; cIter != commonTax.end() ; ++cIter)
			{
			string		key = (*cIter).first;
			int			v = (*cIter).second;
			if ((v >= commonThresh) && (v > commonBest*90/100) && (common.length() < key.length()))
				{
				common = key;
				commonBest = v;
				}
			}

		// DEBUGGING code for assignment of taxonomy for internal branch
//        DisplayL("\tAssigning: [%s]\n\n", common.c_str());
//        if (common.empty())
//            {
//            DisplayL("Taxonomy for [%s]\n", node->title.c_str());
//            for (CPNodeListIter diter=descent.begin() ; diter != descent.end() ; ++diter)
 //               {
//                CPNode		*des = *diter;
//                if (des == NULL) continue;
//                
//                CSequenceItem *seq = taxonomy->GetSequence(des->title.c_str());
//                if (seq != NULL)
//                    {
//                    tax = seq->hdr;
//                    Trim(tax, "\n\r\t ");
//                    // find first \t rest of line is taxonomy
//                    int pos = tax.find_first_of('\t');
//                    if (pos >= 0)
//                        tax = tax.substr(pos+1,tax.length());
//                    seq->ReleaseSeqData();
//                    }
//                DisplayL("\t[%s]\t[%s]\n", des->title.c_str(), tax.c_str());
//                }
//            }

        node->attrs.Add(ATTR_TAXLABEL, common);
	//	DisplayL("Setting Internal Taxonomy: %s [%s]\n", node->title.c_str(), common.c_str());
		if (--nSeqRemaining%1000 == 0)
			DisplayT("\rSetting Internal Taxonomy %d nodes remaining     ", nSeqRemaining);
		
		LPSTR		data = (LPSTR)malloc(common.length()+1);
		strcpy(data, common.c_str());
		CSequenceItem		*seq = new CSequenceItem(node->title.c_str(), -1, common.length(), data);
		taxonomy->seqlist[seq->name] = seq;
		}

	DisplayT("\n");
    taxEntries.Display(" ", 3);
	return TRUE;
	}


///////////////////////////////////////////////////////////////////////////////////////////////////

void FindInsertLocations_Serial(CSequenceFile *treeSeqfile, CPTree *tree, CSequenceFile *insertSeqfile, CInsertPosArray &inserts)
	{
	// for each node in the tree
	//		for each seq for insertion
	//			calc cost

	time_t				start = time(NULL);

	// build insertion list
	insertSeqfile->ResetSeqIterator();
	CSequenceItem		*seq = insertSeqfile->GetNextSeq();
	while (seq != NULL)
		{
		CSequenceItem	*iseq = insertSeqfile->GetSequence(seq->name.c_str());

		// Check to see if this leaf is in the core tree
		if (coreNames.find(seq->name.c_str()) != coreNames.end())
			{
			DisplayL("**** Sequence is already in core tree: %s\n", seq->name.c_str());
			seq = insertSeqfile->GetNextSeq();		// get next sequence
			continue;
			}

		if (iseq == NULL)
			{
			DisplayL("**** Cannot find sequence for %s\n", seq->name.c_str());
			seq = insertSeqfile->GetNextSeq();		// get next sequence
			continue;
			}

		// check that we have the taxonomy
		string				tax;
		if (taxonomy != NULL)
			tax = GetTaxonomy(seq->name.c_str());
			
		CInsertPos		*item = new CInsertPos();
		item->seq = iseq;
		item->pars.Convert(iseq->GetSeqData());

		for (int i=0 ; i < item->pars.len ; ++i)
			{
			if (item->pars.data[i] != 0)
				{
				if (useMask && (mask[i] != 0)) 
					++item->nSites;
				else if (!useMask)
					++item->nSites;

				item->pars.end = i;
				if (item->pars.start == -1)
					item->pars.start = i;
				}
			}

        item->pars.BuildSegCounts((useMask ? mask : NULL));

		if (item->nSites < 1000)
			TRACE("SHORT");

		DisplayL("Insert [%20s](len=%d) [%5d sites] (%d - %d) [%s]\n", 
						item->seq->name.c_str(), seq->len, 
						item->nSites, 
						item->pars.start, item->pars.end,
						tax.c_str());
		inserts.push_back(item);

		seq = insertSeqfile->GetNextSeq();		// get next sequence
		}

	// For each node in the tree
	//		for each seq in insertion list
	//			score the seq to the node
	//			keep the best scores
	
	int				diffs;
	int				parts;
	int				indel;

	CPNodeList		&nodes = tree->nodeList;
	int				count  = 0;
	DisplayL("Searching for position to insert %d sequences in tree of %d nodes. {Setup: %d sec}\n", inserts.size(), tree->nodeList.size(), time(NULL)-start+1);
	for (CPNodeListIter iter=nodes.begin() ; iter != nodes.end() ; ++iter)
		{
		CPNode		*node = *iter;

		if (++count%1000 == 0)
			DisplayL("[%10d/%d] {%d sec} [%s]\n", count, tree->nodeList.size(), time(NULL)-start+1, node->title.c_str());

		CParsimonySet		*node_pars = GetParsimonySet(node);
		if (node_pars == NULL)
			continue;

		// check match to each of the insert seq
		for (int i=0 ; i < inserts.size() ; ++i)
			{
			CInsertPos		*insert = inserts[i];
		
//			if (insert->nSites < 100) continue;

			int				thresh = insert->best.WorstScore();
            
            int nodecost = CalcCost(insert->pars, *node_pars, diffs, parts, indel, thresh);

			//TRACE("[%10s] %d\n", node->title, nodecost);
			if (nodecost != INT_MAX) 
				{
				// get taxonomy
				string		t = GetNodeName(node);

				string			tax = GetTaxonomy(t.c_str(), node);

				insert->best.Add(nodecost, node, tax.c_str());				
				}
			}
		}

	time_t			finish = time(NULL);
	DisplayL("Insertion of %d sequences in [%d sec]\n", inserts.size(), finish-start+1);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// fast insertion that only processes nodes that have a chance of getting a good score
//
void FindInsertLocations(CSequenceFile *treeSeqfile, CPNode *node, CSequenceItem *seq, CInsertPos *item)
	{
	int 				nodecost = INT_MAX;
	int					diffs    = INT_MAX;
	int					parts;
	int					indel;	
	BOOL				recurse = FALSE;

    --nSeqRemaining;
    CParsimonySet		*node_pars  = GetParsimonySet(node);
	
	if (node_pars == NULL)
		{
		DisplayL("**** Error: could not generate Parsimony Set\n");
		return;
		}
    
    nodecost = CalcCost(item->pars, *node_pars, diffs, parts, indel, INT_MAX);
	//int				per = (item->nSites-(diffs+indel))* 100 / item->nSites;

	//TRACE("[%10s] %d\n", node->title, nodecost);
	if (nodecost < item->best.WorstScore()) // found atleast a partial match
		{
		// add to best match list
		string			name = GetNodeName(node);
		string			tax = GetTaxonomy(name.c_str(), node);

		item->best.Add(nodecost, node, tax.c_str());				
		}
		
	LPCSTR		items = "9";
	if (node->IsBranchNode())
		{
		char		c = node->title[node->title.length()-1];
		if (strchr(items, c) == NULL)
			recurse = TRUE;
		else
			{
		//	DisplayL("Checking '%c'\n", c);
			int				thresh = item->best.WorstScore();
		//	if (thresh < INT_MAX)
		//		thresh = thresh * 120 / 100;
			thresh = min(thresh, item->nSites*fast/100);
	
			// get the max value for matching any other node in subtree
		    CParsimonySet		*union_pars = unionList[node->title];
			CalcCost(item->pars, *union_pars, diffs, parts, indel, thresh);
			
			int				best   = diffs + indel; 
		//	DisplayL("[min(%8d,%8d) > %8d (%8d,%8d)[%s]\n", item->best.WorstScore(), item->nSites*fast/100, diffs,parts,indel, node->title.c_str());
			recurse = (best <= thresh);
		//	if (!recurse)
		//		DisplayL("skipping descendents of [%s]\n", GetNodeName(node));
			}
		}
		
	// Only propogate along this subtree if we can possibly obtain better scores
	if ((fast == 0) || recurse)
		{
		// Now try all the children of this node
		CPNodeList		&children = node->edges;
		for (CPNodeListIter iter=children.begin() ; iter != children.end() ; ++iter)
			{
			CPNode *child = *iter;
			FindInsertLocations(treeSeqfile, child, seq, item);
			}
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// fast insertion that only processes nodes that have a chance of getting a good score
//
void FindInsertLocations(CSequenceFile *treeSeqfile, CPTree *tree, CSequenceFile *insertSeqfile, CInsertPosArray &inserts)
	{
	//	for each seq for insertion
	// 		insert in the tree - starting at root

	CPNode			*root = tree->root;
	int				count  = 0;

	time_t				start = time(NULL);

	// build insertion list
	insertSeqfile->ResetSeqIterator();
	CSequenceItem		*seq = insertSeqfile->GetNextSeq();
	while (seq != NULL)
		{
		CSequenceItem	*iseq = insertSeqfile->GetSequence(seq->name.c_str());

		// Check to see if this leaf is in the core tree
		if (coreNames.find(seq->name.c_str()) != coreNames.end())
			{
			DisplayL("**** Sequence is already in core tree: %s\n", seq->name.c_str());
			seq = insertSeqfile->GetNextSeq();		// get next sequence
			continue;
			}

		if (iseq == NULL)
			{
			DisplayL("**** Cannot find sequence for %s\n", seq->name.c_str());
			seq = insertSeqfile->GetNextSeq();		// get next sequence
			continue;
			}

		// check that we have the taxonomy
		string				tax;
		if (taxonomy != NULL)
			tax = GetTaxonomy(seq->name.c_str());
			
		if (tax.empty() && testing)
			{
			DisplayL("**** Cannot find sequence taxonomy for %s\n", seq->name.c_str());
			seq = insertSeqfile->GetNextSeq();		// get next sequence
			//continue;
			}

		CInsertPos		*item = new CInsertPos();
		item->seq = iseq;
		item->pars.Convert(iseq->GetSeqData());

		for (int i=0 ; i < item->pars.len ; ++i)
			{
			if (item->pars.data[i] != 0)
				{
				if (useMask && (mask[i] != 0)) 
					++item->nSites;
				else if (!useMask)
					++item->nSites;

				item->pars.end = i;
				if (item->pars.start == -1)
					item->pars.start = i;
				}
			}

		if (++count%100 == 0)
			DisplayL("[%10d/%d] {%d sec} [%s]\n", count, insertSeqfile->seqCount, time(NULL)-start+1, seq->name.c_str());
		if (item->nSites < 800)
			TRACE("SHORT");

        nSeqRemaining = tree->nodeList.size();
        time_t      start_item = time(NULL);
		DisplayL("Insert [%20s](len=%d) [%5d sites] (%d - %d) [%s]", 
						item->seq->name.c_str(), item->seq->len, 
						item->nSites, 
						item->pars.start, item->pars.end,
						tax.c_str());
		inserts.push_back(item);

	///	DisplayL("Searching for position to insert %d sequences in tree of %d nodes. {Setup: %d sec}\n", inserts.size(), tree->nodeList.size(), time(NULL)-start+1);
        FindInsertLocations(treeSeqfile, root, seq, item);
        if (nSeqRemaining > 0)
	        DisplayL("(%d sec) [%5d skipped]\n", time(NULL)-start_item, nSeqRemaining);
		else
	        DisplayL("(%d sec)\n", time(NULL)-start_item);

    /// DisplayBestInsertLocations(item, TRUE);
			
		seq = insertSeqfile->GetNextSeq();		// get next sequence
		}

	time_t			finish = time(NULL);
	DisplayL("Insertion of %d sequences in [%d sec]\n", inserts.size(), finish-start+1);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

WORD    *scoreMatrix = NULL;
void BuildScoreMatrix()
    {
    if (scoreMatrix != NULL)
        free(scoreMatrix);
    
    scoreMatrix    = (WORD*)calloc(2, 4096*4096);
    
    for (int i=0 ; i < 4*1024 ; ++i)
        for (int j=0 ; j < 4*1024 ; ++j)
            {
            int         key  = ((i << 12) | j);

            int         diffs = 0;
            int         indel = 0;
            int         parts = 0;
            
            // count number of diffs, indel, partial
            for (int k=0 ; k < 3 ; ++k)
                {
                BYTE v1 = (i >> (k*4)) & 0x0F; // shifts by 0, 4, 8
                BYTE v2 = (j >> (k*4)) & 0x0F;
                if (v1 != v2)
                    {
                    if ((v1 & v2) != 0)
                        parts += 1;
                    else if ((v1 != 0) && (v2 != 0))
                        diffs += 1;
                    else 
                        indel += 1;
                    }
                }
            if ((diffs > 3) || (indel > 3) || (parts > 3) || ((diffs+indel+parts) > 3))
                DisplayL("BAD BAD COMPUTER");
            scoreMatrix[key] = (diffs << 0) | (indel << 3) | (parts << 6);
            }
    }
    
///////////////////////////////////////////////////////////////////////////////////////////////////

int CalcCost_NEW1(CParsimonySet &next_pars, CParsimonySet &node_pars, int &diffs, int &partials, int &indel, int thresh)
    {
	// calc delta to insert next_pars
	int                 len  = min(next_pars.len, node_pars.len);
    register LPCSTR		m    = mask;
    register BYTE		*p1  = next_pars.data;
    register BYTE		*p2  = node_pars.data;
    register DWORD      key;
    
   	++nCalcCost;
    diffs    = 0;
	partials = 0;
	indel    = 0;
	for (int i=0 ; (i < len) ; i+=3)
		{
        key = 0;
        
        for (int j=0 ; j < 3 ; ++j)
            {
            key = key << 4;
            if ((!useMask || *m++) && ((i+j >= next_pars.start) && (i+j <= next_pars.end)))
                {
                key |= ((DWORD)(*p1) << 12) | (*p2&0x0F);
                }
            
            ++p1;
            ++p2;
            }
            
        // look up score for all three values
        register WORD   value = scoreMatrix[key];
        
        if (value  == 0) continue;
        
        diffs += value & 0x07;                
        indel += (value >> 3) & 0x07;
        partials += (value >> 6) & 0x07;

        thresh -= (value & 0x07) + ((value >> 3) & 0x07) + ((value >> 6) & 0x07);
        
        if (thresh < 0)
            return INT_MAX;
        }
    
    ++nCalcCostFull;
	return diffs + indel + partials/4;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CalcCost_NEW4(CParsimonySet &next_pars, CParsimonySet &node_pars, int &diffs, int &partials, int &indel, int thresh)
    {
	// calc delta to insert next_pars
	int			len = min(next_pars.len, node_pars.len);
	LPCSTR		m   = mask;
	
    ++nCalcCost;

    if ((next_pars.segCounts != NULL) && (node_pars.segCounts != NULL))
        {
        // compare the counts, looking for known error
        int count = next_pars.CompareSegments(&node_pars);
        
        // if error is over threshold, skip        
        if (count > thresh)
            return INT_MAX;
        }
                    
    ++nCalcCostFull;
    
	diffs    = 0;
	partials = 0;
	indel    = 0;
    BYTE		*p1  = next_pars.data;
    BYTE		*p2  = node_pars.data;
    
	for (int i=0 ; i < len ; ++i)
		{
        if (!useMask || *(m++)) 
            {
            char		s1 = *p1;
            char		s2 = *p2;
            
            if (s1 != s2)
                {
                // there is a difference, is is partial, diff, or indel
                if ((s1 & s2) != 0) 
                    ++partials;
                else if ((s1 != 0) && (s2 != 0))
                    ++diffs;
                else if ((i >= next_pars.start) && (i <= next_pars.end))  // only count indel if within next_pars area
                    ++indel;

                if ((diffs+indel+partials/4) > thresh)
                    return INT_MAX;
                }
            }
        
        ++p1;
        ++p2;
		}
    
	return diffs + indel + partials/4;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

int CalcCost_ORIG(CParsimonySet &next_pars, CParsimonySet &node_pars, int &diffs, int &partials, int &indel, int thresh)
	{
	// calc delta to insert next_pars
	int			len = min(next_pars.len, node_pars.len);
	register LPCSTR		m   = mask;
	
    ++nCalcCost;
    
	diffs    = 0;
	partials = 0;
	indel    = 0;
    
    register BYTE        *p1  = next_pars.data;
    register BYTE        *p2  = node_pars.data;
    
	for (int i=0 ; i < len ; ++i)
		{
		if ((!useMask || *m++) && ((i >= next_pars.start) && (i <= next_pars.end)))   // only count(indel) if within next_pars area
			{
			register char		s1 = *p1;
			register char		s2 = *p2;
			
			if (s1 != s2)
				{
				// there is a difference, is is partial, diff, or indel
				if ((s1 & s2) != 0) 
					++partials;
				else if ((s1 != 0) && (s2 != 0))
					++diffs;
				else 
					++indel;
				
                if ((diffs+indel+partials/4) > thresh)
                    return INT_MAX;
                }
			}
			
		++p1;
		++p2;
		}

    ++nCalcCostFull;
	return diffs + indel + partials/4;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int CalcCost(CParsimonySet &next_pars, CParsimonySet &node_pars, int &diffs, int &partials, int &indel, int thresh)
    {
    switch (useCalcCost)
        {                        
        case 1:  return CalcCost_NEW1(next_pars, node_pars, diffs, partials, indel, thresh);
        case 4:  return CalcCost_NEW4(next_pars, node_pars, diffs, partials, indel, thresh);
        default: return CalcCost_ORIG(next_pars, node_pars, diffs, partials, indel, thresh);
        }
    }
    
///////////////////////////////////////////////////////////////////////////////////////////////////

CParsimonySet* GetParsimonySet(CPNode *node, BOOL create)
	{
	CParsimonySet		*node_pars = NULL;

	string				name = GetNodeName(node);
	node_pars = (CParsimonySet*)parsimonyList[name];
	
	if ((node_pars == NULL) && create)
		{
		// does not exist, create new entry
    //	DisplayL("Create new CParsimonySet for %s\n", name.c_str());
		
		node_pars = new CParsimonySet();
        
        // get sequence from the sequence file
        CSequenceItem *seq = treeSeqfile->GetSeq(name.c_str());
        if ((seq != NULL) && (seq->len > 0))
            {
            node_pars->Convert(seq->data);
            seq->ReleaseSeqData();
            }

        parsimonyList[name] = node_pars;
                    
	//	if (parsimonyList.size() % 100 == 0)
	//		DisplayL("Created %d parsimony sets\n", parsimonyList.size());
		}
        
	return node_pars;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int	CalcAncestorParsimony(CPTree *tree)
	{
	DisplayT("Calculating Parsimony ...");

	if (treeSeqfile == NULL)
		{
		DisplayT("*** Cannot calc parsimony, no sequence file loaded");
		return 0;
		}

	
	time_t				start = time(NULL);
	int					cost = 0;

	nSeqRemaining = tree->nodeList.size();

	cost = CalcNodeParsimony(tree->root);

    if (fitch)
        {
        DisplayL("Applying parsimonious selections back down the tree.\n");
        int     changes = ForceNodeParsimony(tree->root);
        DisplayL("%d nodes were modified\n", changes);
        }
/***
	// Save the calculated parsimony for quicker startup
	//
	DisplayT("Writing Parsimony ...");

    CSequenceFile   parsimony;
    if (!parsimony.Open("Parsimony.FASTA", "w+"))
        {
            DisplayT("*** Cannot open parsimony file");
            return 0;
        }

    // for each branch in the tree
	//		write sequence to file
    nSeqRemaining = (tree->nodeList.size() + 1) / 2;
    WriteParsimonySet(tree->root, &parsimony);
    
	parsimony.Close();
***/

	time_t				delta = time(NULL) - start + 1;
	double				npersec = tree->nodeList.size() * 1. / delta;
	DisplayT("\nProcessed %d nodes in %d seconds, %.0f nodes per sec\n", tree->nodeList.size(), delta, npersec);
	
	return cost;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

int ForceNodeParsimony(CPNode *node)
    {
	// force specific parsimony from parent onto the children 
    
	if (node == NULL) return 0;
    if (!node->IsBranchNode()) return 0;
    
    CParsimonySet	*node_pars = GetParsimonySet(node, FALSE);
    if (node_pars == NULL) return 0;
    
    //DisplayL("Update node [%s]\n", node->title.c_str());
    
    int     count = 0;
    for (CPNodeListIter iter=node->edges.begin() ; iter != node->edges.end() ; ++iter)
        {
        CPNode *child = (*iter);
        if (child == NULL)
            {
            DisplayL("***** Child not found [%s]\n", node->title.c_str());
            continue;
            }
                    
        CParsimonySet		*childParsSet = GetParsimonySet(child, FALSE);
        if (childParsSet == NULL) continue;
        
        if (childParsSet->Force(node_pars) > 0)
            ++count;
        
        count += ForceNodeParsimony(child);
        }
    
    return count;
    }
    
///////////////////////////////////////////////////////////////////////////////////////////////////

int CalcNodeParsimony(CPNode *node)
	{
	// calc parsimony of the tree
	int					cost = 0;
	string				alignment;

	if (node == NULL) return INT_MAX;

    CParsimonySet	*node_pars = GetParsimonySet(node, TRUE);

	if (node->IsBranchNode())
		{
		CParsimonySet	unionSeq;
		CParsimonySet	interSeq;

		int		count = 0;
		for (CPNodeListIter iter=node->edges.begin() ; iter != node->edges.end() ; ++iter)
			{
			CPNode *child = (*iter);
			if (child == NULL)
				{
				DisplayL("***** Child not found [%s]\n", node->title.c_str());
				continue;
				}

			CParsimonySet		*childParsSet = GetParsimonySet(child, FALSE);

			if ((childParsSet == NULL) || (childParsSet->len == 0))
                {
				cost += CalcNodeParsimony(child);
                childParsSet = GetParsimonySet(child);
                }
                
			unionSeq.Union(childParsSet);
			interSeq.Intersect(childParsSet);

			if (childParsSet->len > 0)
				count++;
			}
        //	unionSeq.Display("UNION");
        //	interSeq.Display("INTER");
        cost += node_pars->Set(&unionSeq, &interSeq);

        // save the union information for scoring subtree
            
        string				name = node->title;
        CParsimonySet		*union_pars = new CParsimonySet();
        union_pars->Union(&unionSeq);
                    
        unionList[name] = union_pars;   
        }
    else 
        {                
        // taxa node
        if (node_pars->len == 0)
            {
            // read the sequence from tree sequence file
            // create the parsimony set from sequence
            string          name = GetNodeName(node);
            CSequenceItem   *seq = treeSeqfile->GetSequence(name.c_str());
            if (seq != NULL) 
                {
                node_pars->Convert(seq->GetSeqData());
                seq->ReleaseSeqData();
                }
            }
        }

	//string		t = node->attrs.GetString(ATTR_NAME);
	//if (t.empty()) t = node->title;
	//TRACE("Setting [%s], Cost = %d\n", t.c_str(), cost);
	node->attrs.Add("PARSIMONY_COST", cost);

	if (--nSeqRemaining%1000 == 0)
		DisplayT("\rCalcParsimony %d nodes remaining     ", nSeqRemaining);

	return cost;
	}
	
///////////////////////////////////////////////////////////////////////////////////////////////////

int WriteParsimonySet(CPNode *node, CSequenceFile *parsimony)
    {
    // Write the node's parsimony value to sequence file
    
    if (!node->IsBranchNode())
        return 0;
    
    // calc the parsimony output string
    CParsimonySet	*node_pars = GetParsimonySet(node, TRUE);
    char			*pars = NULL;
    node_pars->BuildString(pars, -1);
    string				name = GetNodeName(node);

    // add new sequence to the parsSequences
    CSequenceItem *seq =  new CSequenceItem(name.c_str(), -1, strlen(pars), pars);
    parsimony->seqlist[name] = seq;
    parsimony->WriteSequence(seq);
    ///	DisplayL("Parsimony for [%s] from %d decendants: %d,[%100.100s]\n", node->title.c_str(), count, strlen(pars), pars);
    free(pars);
    
    if (--nSeqRemaining%1000 == 0)
        DisplayT("\r %d nodes remaining to be written     ", nSeqRemaining);
        
    return 1;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

LPCSTR GetTaxonomy(LPCSTR name, CPNode *node)
	{
	static string		tax;

	tax.clear();
	if (taxonomy != NULL)
		{
		CSequenceItem *seq = taxonomy->GetSequenceHeader(name);
		if ((seq != NULL) && (seq->hdr != NULL))
			{
			tax = seq->hdr;
			TrimRight(tax, "\n\r\t ");

			int		k = tax.find_first_of("\t ");
			if (k > 0)
				tax = tax.substr(k+1, tax.length());
			Trim(tax, "\n\r\t ");
			}
		}

	return tax.c_str();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

void DisplayInsertLocations(CInsertPosArray &inserts, BOOL all)
{
	// show best for each insert
	int				diffs;
	int				parts;
	int				indel;
    
    //	int				nerror       = 0;
    //	int				exact        = 0;
    //	int				partial_seq  = 0;
    //	int				partial_tree = 0;
    //	int				partials[10] = {10*0};	
    
	for (int i=0 ; i < inserts.size() ; ++i)
		{
            string				name = inserts[i]->seq->name;
            
            string				tax = GetTaxonomy(name.c_str());
            
            DisplayL("======================================================================\n"); 
            if (all)
                DisplayL("\nInsertion of [%s][%5d sites] [%d matches] Taxonomy[%s]:\n", 
                			inserts[i]->seq->name.c_str(), inserts[i]->nSites, 
                			inserts[i]->best.list.size(), tax.c_str());
            
            string					taxFirst;
            string					taxBest;
            CBestList::iterator		iterB = inserts[i]->best.list.begin();
            int						count = 0;
            
            if (inserts[i]->best.list.size() == 0)
                DisplayL("\tNO MATCHES found\n");
            
            for (iterB=inserts[i]->best.list.begin() ; iterB != inserts[i]->best.list.end() ; ++iterB)
                {
                ++count;
                if ((iterB->score > inserts[i]->nSites/5) && (count > 3)) 
                    break;
                
                //	TRACE("Getting taxonomy for [%s]", iterB->node->title);
                ///??			SetMyTaxonomy(iterB->node);
                if (iterB->tax.empty())
                    {
                    iterB->node->GetComment(ATTR_TAXLABEL, taxBest, "");
                    
                    if (taxBest.empty() && (taxonomy != NULL))
                        {
                        string		t = GetNodeName(iterB->node);                                    
                        taxBest = GetTaxonomy(t.c_str(), iterB->node);
                        }
                    
                    iterB->tax = taxBest;
                    }
                //	TRACE("[%s]\n", iterB->tax);
                
                taxBest = iterB->tax;
                if (taxFirst.empty())
                    taxFirst = taxBest;
                
                CParsimonySet		*node_pars = GetParsimonySet(iterB->node);
                char				marker;
                if (node_pars != NULL)
                    {
                    int cost = CalcCost(inserts[i]->pars, *node_pars, diffs, parts, indel);
                    ///		ASSERT(cost==iterB->score);
                    
                    int		per = -1;
                    if (inserts[i]->nSites > 0)
                        per = (inserts[i]->nSites-(diffs+indel))* 100 / inserts[i]->nSites;
                    
                    marker = ' ';
                    if (cost <= inserts[i]->nSites*3/100)
                        marker = '*';
                    else if (cost <= inserts[i]->nSites*6/100)
                        marker = '+';
                    
                    if (all)
                        {
                        DisplayL("%c(%3d%%) %6d (%4d diffs, %4d partials, %4d indels) [%12s][%s]\n", 
                                 marker, per,
                                 iterB->score, diffs, parts, indel,
                                 iterB->node->title.c_str(), taxBest.c_str());
                        
                        if (stats)
                            {
                            string		match;
                            GetTaxonomyMatch(tax, taxBest, match, rankCounts);
                            
                            string		name = "";
                            if (iterB == inserts[i]->best.list.begin())
                                name = inserts[i]->seq->name;
                            
                            fprintf(stats, "%15s\t%s\t%3d%%\t%4d\t%4d\t%4d\t%s\t%s\t%s\n", 
                                    name.c_str(), iterB->node->title.c_str(), per, 
                                    diffs, indel, parts, match.c_str(), taxBest.c_str(), tax.c_str());
                            }
                        }
                    }
                }
            
            
            iterB = inserts[i]->best.list.begin();
            if (iterB != inserts[i]->best.list.end())
                {
                    taxBest = iterB->tax;
                    
                    string		marker = "";
                    string		t = (taxBest.empty() ? taxFirst : taxBest);
                    
                    GetTaxonomyMatch(tax, t, marker, rankCounts);
                    
                    CParsimonySet		*node_pars = GetParsimonySet(iterB->node);
                    
                    int cost = CalcCost(inserts[i]->pars, *node_pars, diffs, parts, indel);
                    DisplayL("\nInsertion [%5d sites] [%20.20s] %6d (%4d diffs, %4d partials, %4d indels)\n", 
                             inserts[i]->nSites, inserts[i]->seq->name.c_str(),
                             cost, diffs, parts, indel,
                             0);
                    
                    string				diff;
                    int					i;
                    for (i=0 ; (i < tax.length()) && (i < taxBest.length()) ; ++i)
                        diff += (tax[i] == taxBest[i] ? ' ' : 'X');
                    for ( ; i < tax.length() ; ++i)
                        diff += 'X';
                    for ( ; i < taxBest.length() ; ++i)
                        diff += 'X';
                    
                    DisplayL("%12.12s\t  Seq[%s]\n", marker.c_str(), tax.c_str());
                    DisplayL("%12.12s\t     [%s]\n", "", diff.c_str());
                    DisplayL("%12.12s\t Best[%s]\n", "", taxBest.c_str());
                    if (taxBest != taxFirst)
                        DisplayL("\tFirst[%s]\n", taxFirst.c_str());
                }
		}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void DisplayBestInsertLocations(CInsertPos *item, BOOL all)
    {
	// show best for each insert
	int				diffs;
	int				parts;
	int				indel;
    
    string				name = item->seq->name;
    string				tax = GetTaxonomy(name.c_str());
    
    DisplayL("======================================================================\n"); 
    if (all)
        DisplayL("\nInsertion of [%s][%5d sites] [%d matches] Taxonomy[%s]:\n", 
        			item->seq->name.c_str(), item->nSites, item->best.list.size(), tax.c_str());
            
    string					taxFirst;
    string					taxBest;
    CBestList::iterator		iterB = item->best.list.begin();
    int						count = 0;
            
    if (item->best.list.size() == 0)
        DisplayL("\tNO MATCHES found\n");
    
    for (iterB=item->best.list.begin() ; iterB != item->best.list.end() ; ++iterB)
        {
        ++count;
        if ((iterB->score > item->nSites/5) && (count > 3)) 
            break;
        
        //	TRACE("Getting taxonomy for [%s]", iterB->node->title);
        ///??			SetMyTaxonomy(iterB->node);
        if (iterB->tax.empty())
            {
            iterB->node->GetComment(ATTR_TAXLABEL, taxBest, "");
            
            if (taxBest.empty() && (taxonomy != NULL))
                {
                string		t = GetNodeName(iterB->node);
                
                taxBest = GetTaxonomy(t.c_str(), iterB->node);
                }
            
            iterB->tax = taxBest;
            }
        //TRACE("[%s]\n", iterB->tax);
        
        taxBest = iterB->tax;
        if (taxFirst.empty())
            taxFirst = taxBest;
        
        CParsimonySet		*node_pars = GetParsimonySet(iterB->node);
        char				marker;
        if (node_pars != NULL)
            {
            int cost = CalcCost(item->pars, *node_pars, diffs, parts, indel);
            ///		ASSERT(cost==iterB->score);
            
            int		per = -1;
            if (item->nSites > 0)
                per = (item->nSites-(diffs+indel))* 100 / item->nSites;
            
            marker = ' ';
            if (cost <= item->nSites*3/100)
                marker = '*';
            else if (cost <= item->nSites*6/100)
                marker = '+';
            
            if (all)
                {
                DisplayL("%c(%3d%%) %6d (%4d diffs, %4d partials, %4d indels) [%12s][%s]\n", 
                         marker, per,
                         iterB->score, diffs, parts, indel,
                         iterB->node->title.c_str(), taxBest.c_str());
                
                if (stats)
                    {
                    string		match;
                    GetTaxonomyMatch(tax, taxBest, match, rankCounts);
                    
                    string		name = "";
                    if (iterB == item->best.list.begin())
                        name = item->seq->name;
                    
                    fprintf(stats, "%15s\t%s\t%3d%%\t%4d\t%4d\t%4d\t%s\t%s\t%s\n", 
                            name.c_str(), iterB->node->title.c_str(), per, 
                            diffs, indel, parts, match.c_str(), taxBest.c_str(), tax.c_str());
                    }
                }
            }
        }
    
    iterB = item->best.list.begin();
    if (iterB != item->best.list.end())
        {
        taxBest = iterB->tax;
        
        string		marker = "";
        string		t = (taxBest.empty() ? taxFirst : taxBest);
        
        GetTaxonomyMatch(tax, t, marker, rankCounts);
        
        CParsimonySet		*node_pars = GetParsimonySet(iterB->node);
        
        int cost = CalcCost(item->pars, *node_pars, diffs, parts, indel);
        DisplayL("\nInsertion [%5d sites] [%20.20s] %6d (%4d diffs, %4d partials, %4d indels)\n", 
                 item->nSites, item->seq->name.c_str(),
                 cost, diffs, parts, indel,
                 0);
        
        string				diff;
        int					i;
        for (i=0 ; (i < tax.length()) && (i < taxBest.length()) ; ++i)
            diff += (tax[i] == taxBest[i] ? ' ' : 'X');
        for ( ; i < tax.length() ; ++i)
            diff += 'X';
        for ( ; i < taxBest.length() ; ++i)
            diff += 'X';
        
        DisplayL("%12.12s\t  Seq[%s]\n", marker.c_str(), tax.c_str());
        DisplayL("%12.12s\t     [%s]\n", "", diff.c_str());
        DisplayL("%12.12s\t Best[%s]\n", "", taxBest.c_str());
        if (taxBest != taxFirst)
            DisplayL("\tFirst[%s]\n", taxFirst.c_str());
        }
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

void DisplayRankCounts(LPCSTR header, int counts[][4])
	{
	if (!testing) return;
	if (counts[0][0] == 0)
        return;
        
    DisplayL("\n%s\n                     _________Precision________    __________Recall__________\n", header);
	for (int r=0 ; r < RANK_N ; ++r)
		{
		if (counts[r][0] == 0) continue;
		double	per  = 100. * counts[r][1] / counts[r][0];
        double  per2 = 100. * counts[r][3] / counts[r][2];
        
        DisplayL("\t%10.10s: %8d %8d (%6.2f%%)     %8d %8d (%6.2f%%)\n", 
        			rankNames[r], counts[r][0], counts[r][1], per, counts[r][3], counts[r][2], per2);
		}
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

// STATS FORMAT: name, insert position title, percentage match, diffs, indel, partials, conf level, best taxonomy [, correct taxonomy]
#define TSTATS_HEADER "#Insertion Name\tSite   \tMatch:%\tdiff\tindel\tpartial\tconf\ttaxonomy | correct taxonomy\n"
#define TSTATS_FORMAT "%15s\t%s\t%3d%%\t%4d\t%4d\t%4d\t%1d\t%s\t|\t%s\n"
#define STATS_HEADER "#Insertion Name\tMatch:%\ttaxonomy\n"
#define STATS_FORMAT "%15s\t%3d%%\t%s\n"
#define STATS_FORMAT2 "%15s\t    \t%s\n"

int			LOWER_SCORE = 0;
int			UPPER_SCORE = 100;
int			LOWER_PER = 90;
int			UPPER_PER = 101;

void OutputInsertLocations(CInsertPosArray &inserts, BOOL all)
    {
	// show best for each insert
	int				diffs;
	int				parts;
	int				indel;
    
	CTaxEntry		taxEntries("root");	// keeps track of diversity seen in taxonomy of inserted sequences
    
    if (stats)
		fprintf(stats, "%s\n", STATS_HEADER);
        
	for (int i=0 ; i < inserts.size() ; ++i)
		{
        string                  name = inserts[i]->seq->name;        
        string					taxCorrect = GetTaxonomy(name.c_str());
        string					taxFirst;
        string					taxBest;
        CBestList::iterator		iterB = inserts[i]->best.list.begin();
        int						count = 0;
        CTaxEntry				taxVotes("Votes:");	// keep track of all high scoring taxonomies for single sequence
        
        DisplayL("======================================================================\n"); 
        DisplayL("\nInsertion: [%20.20s] [%5d sites]\n", name.c_str(), inserts[i]->nSites); 
            
        if (inserts[i]->best.list.size() == 0)
            {
            DisplayL("\tNO MATCHES found\n");
            if (stats)
                {
                char   buf[1024];
                sprintf(buf, "NO MATCHES found (%d sites in sequence)", inserts[i]->nSites); 
                fprintf(stats, STATS_FORMAT, name.c_str(), 0, buf);
                if (!taxCorrect.empty())
                    fprintf(stats, STATS_FORMAT2, "",taxCorrect.c_str());
                fprintf(stats, "\n");
                }

            continue;
            }

        int						bestScore = inserts[i]->best.list.begin()->score;
        int						threshScore = bestScore + (bestScore >= 100 ? 50 : (bestScore > 50 ? 20 : 10));
        
        for (iterB=inserts[i]->best.list.begin() ; iterB != inserts[i]->best.list.end() ; ++iterB)
            {
            ++count;
            if ((iterB->score > inserts[i]->nSites/5) && (count > 3)) // score too low for processing, process at least 3
                break;
            
            string		t = GetNodeName(iterB->node);
            //	find the taxonomy for this insertion point
            if (iterB->tax.empty())
                {
                iterB->node->GetComment(ATTR_TAXLABEL, taxBest, "");
                
                if (taxBest.empty() && (taxonomy != NULL))
                    {
                    SetMyTaxonomy(iterB->node); // set taxonomy at insertion point
                    
                    taxBest = GetTaxonomy(t.c_str(), iterB->node);
                    }
                
                iterB->tax = taxBest;
                }
                
            // find distance from position in full tree
            double         dist = -1;
            if (fullTree)
            	{
                dist = ComparePosition(t.c_str(), iterB->node, fullTree, tree);
	            DisplayL("\t\t%5d [%5f] [%s][%s]\n", iterB->score, dist, iterB->node->title.c_str(), iterB->tax.c_str());
                }
            
            if (iterB->score < threshScore)
                {
                int         votes = ((threshScore-iterB->score)*100/threshScore); // most votes for best scores
                taxVotes.Add(iterB->tax.c_str(), votes);
                }
                
            taxBest = iterB->tax;
            if (taxFirst.empty())
                taxFirst = taxBest;
            
            CParsimonySet		*node_pars = GetParsimonySet(iterB->node);
            int					conf;
            if (node_pars != NULL)
                {
                int cost = CalcCost(inserts[i]->pars, *node_pars, diffs, parts, indel);
                //ASSERT(cost==iterB->score);
                
                int		per = -2;
                if (inserts[i]->nSites > 0)
                    per = (inserts[i]->nSites-(diffs+indel))* 100 / inserts[i]->nSites;

                conf = 0;
                if (cost <= inserts[i]->nSites*3/100)
                    conf = 9;
                else if (cost <= inserts[i]->nSites*6/100)
                    conf = 5;
                
	            DisplayL("\t\t%5d [%3d%%] [%s][%s]\n", iterB->score, per, iterB->node->title.c_str(), iterB->tax.c_str());

                if (all)
                    {
                    DisplayL("%1d (%3d%%) %6d (%4d diffs, %4d partials, %4d indels) [%12s]%c[%s]\n", 
                             conf, per,
                             iterB->score, diffs, parts, indel,
                             iterB->node->title.c_str(),(iterB->node->IsBranchNode()?'^':' '), taxBest.c_str());
                    
                    if (stats)
                        {
                        string		name = "";
                        if (iterB == inserts[i]->best.list.begin())
                            name = inserts[i]->seq->name;
 
                        fprintf(stats, TSTATS_FORMAT, name.c_str(), iterB->node->title.c_str(), per, 
                        			diffs, indel, parts, conf, taxBest.c_str(), taxCorrect.c_str());
                        }
                    }
                }
            }
        
        iterB = inserts[i]->best.list.begin();
        if (iterB != inserts[i]->best.list.end())
            {
            taxBest = iterB->tax;
            
            string		marker = "";
            string		t = (taxBest.empty() ? taxFirst : taxBest);
                     
            if (t.empty())
                {
                // try to find taxonomy from ancestory
                CPNode      *node = iterB->node->parent;
                while ((node != NULL) && t.empty())
                    {
					t = node->attrs.GetString(ATTR_TAXLABEL);
		    	//	DisplayL("\t using ancestor (%s) taxonomy: [%s] \n", GetNodeName(node), t.c_str());
                    node = node->parent;
                    }
                taxBest = t;
                }
                
            CParsimonySet		*node_pars = GetParsimonySet(iterB->node);
            int cost = CalcCost(inserts[i]->pars, *node_pars, diffs, parts, indel);
            DisplayL("\t %6d (%4d diffs, %4d partials, %4d indels) [%s] \n", 
                     cost, diffs, parts, indel, taxBest.c_str());
            
            int			per = -1;
            if (inserts[i]->nSites > 0)
                per = (inserts[i]->nSites-cost)* 100 / inserts[i]->nSites;
            
            if (per >= scoreThresh)
                {
                if (testing > 2)
	                taxVotes.Display("");
                CStringList taxVoteList;
                taxVotes.FindBest(60, taxVoteList);
                string      leader = "Assigned Taxonomy:";
                if (taxVoteList.size() > 0)
                	{
					for (int knox=0; knox < taxVoteList.size() ; ++knox)
						{
						DisplayL("%18s[%s]\n", leader.c_str(), taxVoteList[knox].c_str());
						leader = "";
						}
					if ((taxVoteList.size() > 0) && !taxVoteList[0].empty())
						taxBest = taxVoteList[0];
					}
				else						
					DisplayL("%18s[%s]\n", leader.c_str(), taxBest.c_str());
                	
                taxEntries.Add(taxBest.c_str());	// update diversity

                if (!taxCorrect.empty())
                    {
                    GetTaxonomyMatch(taxCorrect, (taxBest.empty()? taxFirst : taxBest), marker, rankCounts);

                    string      phylum;
                    GetRankFromTax(RANK_PHYLUM, taxCorrect, phylum, FALSE);
                    if (phylumRankCounts.find(phylum.c_str()) == phylumRankCounts.end())
                        {
                        // create ranks
                        phylumRankCounts[phylum.c_str()] = new CRankCounts();
                        }
                    CRankCounts *counts = phylumRankCounts[phylum.c_str()]; 
                    GetTaxonomyMatch(taxCorrect, (taxBest.empty()? taxFirst : taxBest), marker, counts->rankCounts);
                    }
                }    
            else
            	{
                taxBest.clear();
                DisplayL("no assignment (%d%% < %d%%)", per, scoreThresh);
                }
                
            int				conf = 0;
            if (cost <= inserts[i]->nSites*3/100)
                conf = 9;
            else if (cost <= inserts[i]->nSites*6/100)
                conf = 5;
            
            if ( !taxCorrect.empty() && testing )
                DisplayL(" Correct Taxonomy:[%s]\n", taxCorrect.c_str());
            
            if (stats)
                {
                size_t n;
                while ((n=taxBest.find(' ')) != string::npos)
                    taxBest.erase(n,1);
                while ((n=taxCorrect.find(' ')) != string::npos)
                    taxCorrect.erase(n,1);
                    
                if (taxBest.empty())
                    taxBest = "No assignment could be made!";
                fprintf(stats, STATS_FORMAT, name.c_str(),  per, taxBest.c_str());
                if (!taxCorrect.empty())
                    {
                    fprintf(stats, STATS_FORMAT2, "",taxCorrect.c_str());
                    
                    // find first position of difference
                    string diffs;
                    int i = 0;
                    while ((i < taxBest.length()) && (i < taxCorrect.length()) && (taxBest[i] == taxCorrect[i]))
                        {
                        ++i;
                        diffs += " ";
                        }
                    while ((i < taxBest.length()) && (i < taxCorrect.length()))
                        {
                        ++i;
                        if (i < taxBest.length())
                            diffs += "x";
                        else
                            diffs += "o";
                        }
                    
                    //if (diffs.find('x') < diffs.length())
                        fprintf(stats, STATS_FORMAT2, "",diffs.c_str());
                    }
                    
                fprintf(stats, "\n");
                fflush(stats);
                }
            }
		}
    
	DisplayL("Taxonomies List:\n");
	taxEntries.Display(" ");
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

int ParseTaxonomy(string t, CStringList &taxList)
	{
	// find next separator, add to 
	string		tax = t;
	ToLower(tax);
	while (!tax.empty())
		{
		// find separator
		int			sep = 0;
		while ((sep < tax.length()) && (tax[sep] != ';') && (tax[sep] != '/'))
			++sep;

		string		entry = tax.substr(0, sep);
		LPCSTR		unclassified = "unclass";

		string		lower = entry.substr(0, strlen(unclassified));
		for (int i=0 ; i < sizeof(unclassified) ; ++i)
			lower[i] = tolower(lower[i]);
			
		if ((entry.size() > 0) && (strncmp(lower.c_str(), unclassified, strlen(unclassified)) != 0)) // does NOT begin with "Unclassified"
			taxList.push_back(entry);

		if (sep < tax.length())
			tax = tax.substr(sep+1);
		else
			tax.clear();
		}

	return taxList.size();
	}
    
///////////////////////////////////////////////////////////////////////////////////////////////////

void GetRankFromTax(int rank, string &tax, string &rankStr, BOOL only)
    {
    CStringList		taxList;        
    ParseTaxonomy(tax, taxList);
    char            sep = tax.find(';') >= 0 ? ';' : '/';
    
    rankStr.clear();
    
    for (int i=0 ; (i <= rank) && (i < taxList.size()) ; ++i)
        {
        if (!only)
            {
            if (!rankStr.empty())
                rankStr += sep;
            }

        if (!only || (i == rank))
            rankStr += taxList[i];
        }
    }


///////////////////////////////////////////////////////////////////////////////////////////////////

int GetTaxonomyMatch(string &taxA, string &taxB, string &marker, int counts[][4])
	{
    // taxA is the correct solution
    // taxB is our assigned solution
    //
	int		levels = 0;
	BOOL	matching = TRUE;	// TRUE if still matching
	marker = "";

	// parse the entries from each taxonomy
	CStringList		taxListA;
	CStringList		taxListB;

	ParseTaxonomy(taxA, taxListA);
	ParseTaxonomy(taxB, taxListB);

	for (int k=0 ; (k < taxListA.size()) && (k < taxListB.size()) && (k < RANK_N) ; ++k)
		{
		counts[k][0] += 1;  // saw another rank at this level

		if (matching && (k < taxListB.size()))
			{
			LPCSTR		tA = taxListA[k].c_str();
			LPCSTR		tB = taxListB[k].c_str();

			if (strcmp(tA,tB) == 0)
				{
				// matched at this level
				counts[k][1] += 1; // saw another match at this level
				++levels;
				}
			else if (precision) // count all mistakes?
				matching = FALSE;
			}
		}
    
    matching = TRUE;
    for (int k=0 ; (k < taxListA.size()) && (k < RANK_N) ; ++k)
        {
        counts[k][2] += 1;  // saw another rank at this level
        
        if (matching && (k < taxListB.size()))
            {
            LPCSTR		tA = taxListA[k].c_str();
            LPCSTR		tB = taxListB[k].c_str();
            
            if (strcmp(tA,tB) == 0)
                {
                // matched at this level
                counts[k][3] += 1; // saw another match at this level
                }
            else 
                matching = FALSE;
            }
        }
        

	// Create marker for display
	if ((levels < taxListA.size()) && (levels < taxListB.size()))
		{
		// check to see how close we are
		if (levels >= 4)
			{
			marker = ">>>";
			}
		else
			{
			marker = "XXX";
			}
		}
	else if (levels < taxListA.size())
		{
		marker = ">  ";
		}
	else if (levels < taxListB.size())
		{
		marker = "<  ";
		}
	for (int r=levels ; r > 0 ; --r)
		marker += rankNames[r-1][0];

	return levels;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL AddTaxonomyToTree(CPTree *tree, CSequenceFile *taxonomy)
	{
	if (tree == NULL) return FALSE;
	if (taxonomy == NULL) return FALSE;

	CPNodeList		&nodes = tree->nodeList;
	for (CPNodeListIter iter=nodes.begin() ; iter != nodes.end() ; ++iter)
		{
		string				tax;
		CPNode				*node = *iter;

		if (node == NULL) continue;

		CSequenceItem		*seq = taxonomy->GetSequence(node->title.c_str());

		if (seq != NULL) 
            {
            if (seq->data)
                tax = seq->data;
            else if (seq->hdr)
                tax = seq->hdr;
            }

		int		pos = tax.find_first_of('\t');
		if (pos < 20)
			tax = tax.substr(pos+1);
		TrimRight(tax, "\n\r\t ");

		if (!node->IsBranchNode())
			tax = "CORE: " + tax;

		node->attrs.Add(ATTR_COMMENT, tax);
		}

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

// Builds a map of tree node names that have list of new sequences to be inserted as children
CTreeInserts* BuildTreeInsertLists(CInsertPosArray &inserts)
	{
	if (inserts.size() == 0) return NULL;
	CTreeInserts*	treeInserts = new CTreeInserts();

	for (int i=0 ; i < inserts.size() ; ++i)
		{
		CInsertPos		*in = inserts[i];
		if (in->best.list.size() == 0) continue;
		if (in->seq == NULL) continue;

		CPNode					*node   = in->best.list[0].node;
		if (!node->IsBranchNode())
			node = node->parent;

		CTreeInsertList			*iList = (*treeInserts)[node];

		if (iList == NULL)
			{
			iList = new CTreeInsertList();
			(*treeInserts)[node] = iList; 
			}

		iList->push_back(in);
		}

	return treeInserts;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL BuildInsertTree(CTreeInserts* insertLists, CPTree *tree)
	{
	if (tree == NULL) return FALSE;
    if (insertLists == NULL) return FALSE;
	if (insertLists->size() == 0) return FALSE;

	// for each item in the map
	for (CTreeInserts::iterator iter=insertLists->begin() ; iter != insertLists->end() ; ++iter)
		{
		CPNode 			*node    = (*iter).first;
		CTreeInsertList *inserts = (*iter).second;

		CPNode			*child = NULL;
		
		if (inserts->size() == 0) continue;
		
		CInsertPos		*in = (*inserts)[0];
		int				score = in->best.list[0].score * 100 / in->nSites;
		if (score > 100-scoreThresh) 
			{
			DisplayL("[%s] at %d%% does not meet threshold (%d%%)\n", in->seq->name.c_str(), score, 100-scoreThresh);
			continue;
			}
		
		if (inserts->size() == 1)
			{
			// insert as child of the node
			CInsertPos		*in = (*inserts)[0];

			child = AddChild(tree, node, in);
			node->AddEdge(child);
			}
		else
			{
			// Creaate new branch node
			// for each item in the list
			//		add as child of the new branch

			double			shortest = 1000000.;
			CPNode			*branch = new CPNode(tree);
			node->AddEdge(branch);

			for (int i=0 ; i < inserts->size() ; ++i)
				{
				CInsertPos		*in = (*inserts)[i];
				child = AddChild(tree, node, in);
				double		dist = child->attrs.GetDouble(ATTR_DISTANCE);
				if (dist < shortest)
					shortest = dist;
				branch->AddEdge(child);
				}

			if (shortest > 0)
				{
			//	double		dist = branch->attrs.GetDouble(ATTR_DISTANCE);
			//	DisplayL("Shortest distance %f\n", shortest);
				shortest /= 2;
				branch->attrs.Add(ATTR_DISTANCE, shortest);
				// subtract shortest from each of the children
				for (CPNodeListIter iter=branch->edges.begin() ; iter != branch->edges.end() ; ++iter)
					{
					CPNode		*child = (*iter);
					double		dist = child->attrs.GetDouble(ATTR_DISTANCE);
					dist -= shortest;
					child->attrs.Add(ATTR_DISTANCE, dist);
					}
				}
			}

		}

	return TRUE;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

CPNode* AddChild(CPTree *tree, CPNode *node, CInsertPos *in)
	{
	int				diffs;
	int				parts;
	int				indel;

	CPNode			*inserted = new CPNode(tree);

	inserted->attrs.Add(ATTR_NAME, in->seq->name);

	double dist = CalcDist(node, in, diffs, parts, indel);
    
    if (jukes_cantor)
        {
        //        distance = -3/4 * ln(1 - 4/3*p)
        //        where p is the proportion of sites with different nucleotides.       
        double       p = (diffs + indel) * 1. / in->nSites;
 
        dist = -3./4. * log(1. - (4./3.*p));
        }
        
	inserted->attrs.Add(ATTR_DISTANCE, dist);

	// Take parent taxonomy as node's taxonomy
	string 		tax = in->best.list[0].tax;
	
	while (tax.empty() && (node != NULL))
		{
		tax = node->attrs.GetString(ATTR_TAXLABEL);
	//	DisplayL("\t\ttrying ancestor [%s]\n", GetNodeName(node));
		node = node->parent;
		}

	int			per = -1;
	if (in->nSites > 0)
		per = (in->nSites-(diffs+indel))* 100 / in->nSites;
		
	char			comment[2*1024];
	sprintf(comment,"Sites:%d, %d%% diffs[%d,%d,%d], Score:%d, [%s]", in->nSites, per, diffs, parts, indel, in->best.list[0].score, tax.c_str());
	inserted->attrs.Add(ATTR_COMMENT, comment);

	return inserted;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

double CalcDist(CPNode *parent, CInsertPos* in, int &diffs, int &parts, int &indel)
	{
	// calc the distance for tree branch edge
	double				dist = 0.1;
	CParsimonySet		*node_pars = GetParsimonySet(parent);
	if (node_pars != NULL)
		{
		CParsimonySet		pars(&in->pars);
		CalcCost(pars, *node_pars, diffs, parts, indel); // just need the details

		if (in->nSites > 0)
			dist = (diffs+indel)*1. / in->nSites;
		}

	return dist;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL SetMyTaxonomy(CPNode *node)
    {
	// calc the taxonomy of this node
    
	if (taxonomy == NULL) return FALSE;
	if (node == NULL) return FALSE;
	if (!node->IsBranchNode()) return FALSE;
    
	string				common;
	string				tax;
    
	LPCSTR				def = "not set";
	common = node->attrs.GetString(ATTR_TAXLABEL, def);
	if (common != def) return TRUE;
	common.clear();
    
	for (CPNodeListIter iter=node->edges.begin() ; iter != node->edges.end() ; ++iter)
		{
        CPNode		*child = *iter;
        
        if (node->IsBranchNode())
            {
            SetMyTaxonomy(child);
            }
        
        CSequenceItem *seq = taxonomy->GetSequence(child->title.c_str());
        if (seq != NULL)
            {
            tax = seq->hdr;
            int		pos = tax.find_first_of('\t');
            if (pos < 20)
                tax = tax.substr(pos+1,tax.length());
            TrimRight(tax, "\n\r\t ");
            if (tax.empty())
                {
                string		name = GetNodeName(child);
                DisplayL("No Taxnomy for [%s][%s]\n", name.c_str(), child->title.c_str());
                }
            }
        else if ((seq->hdr == NULL) || (seq->hdr))
            {
            string		name = GetNodeName(child);
            DisplayL("No Taxnomy for [%s][%s]\n", name.c_str(), child->title.c_str());
            }
                    
        if (common.empty())
            common = tax;
        else if (!tax.empty())
            {
            // calc how much is the same between two list
            int			i;
            for (i=0 ; (i < common.length()) && (i < tax.length()) ; ++i)
                {
                    if (common[i] != tax[i])
                        break;
                }
            if (i > 0)
                common = common.substr(0,i);
            else
                common.clear();
            }
		}
    
    if (tax.empty())
        {
        string		name = GetNodeName(node);
        DisplayL("No Taxnomy for [%s][%s]\n", name.c_str(), node->title.c_str());
        }
    
    node->attrs.Add(ATTR_TAXLABEL, common);
	CSequenceItem		*seq = new CSequenceItem(node->title.c_str(), -1, 0); //, common.GetLength(), common);
	seq->AllocateSeqHeader(common.length()+1);
	strcpy(seq->hdr, common.c_str());
    
	CSequenceListIter iter = taxonomy->seqlist.find(seq->name);
	if (iter != taxonomy->seqlist.end())
		{
        delete (*iter).second;
		}
	taxonomy->seqlist[seq->name] = seq;
    
	return TRUE;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////
#define CORE_NAME1 "CORE_NAME1"
#define CORE_NAME2 "CORE_NAME2"

// routine to mark the branch nodes with multiple subtrees with core node
int CalcCoreBranches(CPNode *node, CSequenceFile *core_seq)
    {
    int             count = 0;
    
    if (node == NULL) return count;
    
    if (node->IsBranchNode())
        {
        // for each child
        //      check to see if it has a core sequence
        string          core1;
        string          core2;
        CPNodeList		&children = node->edges;
        for (CPNodeListIter iter=children.begin() ; iter != children.end() ; ++iter)
            {
            CPNode *child = *iter;
            if (CalcCoreBranches(child, core_seq) > 0)
                {
                ++count;
                string      name = child->attrs.GetString(CORE_NAME1);
                if (core1.empty())
                    core1 = name;
                else
                    core2 = name;
                }
            }
            
        node->SetAttr(CORE_NAME1, core1.c_str());
        node->SetAttr(CORE_NAME2, core2.c_str());
        }
    else
        {
        // check to see if the node is in the core sequences
        string      name = GetNodeName(node); 
        if (core_seq->GetSeq(name.c_str()) != NULL)
            {
            ++count; // found in core sequence list
            node->SetAttr(CORE_NAME1, name);
            }
        }

    return count;
    }

// routine to find closest ancestor marked with core multiple core branches

// routine to find branch node containing two given core nodes
CPNode *FindBranchContaining(CPTree *tree, LPCSTR one, LPCSTR two)
    {
    if ((one == NULL) || (two == NULL)) return NULL;
    
    // find the node for each name given
    CPNode *node_one = tree->Find(one);
    if (node_one == NULL) return NULL;
    
    CPNode *node_two = tree->Find(two);
    if (node_two == NULL) return NULL;
    
    // Create list of ancestors for each node
    
    CPNodeList ancestors1;
    CPNodeList ancestors2;
    
    node_one->GetAncestors(ancestors1, FALSE);   // ordered from root to direct parent
    node_two->GetAncestors(ancestors2, FALSE);   // ordered from root to direct parent
    
    // Find the last (lowest in tree) common ancestor
    int i = 0;
    while ((i < ancestors1.size()) && (i < ancestors2.size()) && (ancestors1[i] == ancestors2[i]))
        ++i;
    // ith position is first difference

    if (i == 0) return NULL;
    
    return ancestors1[i-1];
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

// routine to calc the distance between two given branch nodes
double DistanceBetween(CPNode *one, CPNode *two)
    {
    double         dist = 0;
    
    if ((one == NULL) || (two == NULL)) return INT_MAX;
    
    // Create list of ancestors for each node

    CPNodeList ancestors1;
    CPNodeList ancestors2;
    
    one->GetAncestors(ancestors1, FALSE);   // ordered from root to direct parent
    two->GetAncestors(ancestors2, FALSE);   // ordered from root to direct parent

    // Find the last (lowest in tree) common ancestor
    int i = 0;
    while ((i < ancestors1.size()) && (i < ancestors2.size()) && (ancestors1[i] == ancestors2[i]))
            ++i;
    // ith position is first difference
    
    // Distance is the number of branches from each node to common ancestor
//    dist += (ancestors1.size() - i); // distance to common ancestor from node one
//    dist += (ancestors2.size() - i); // distance to common ancestor from node two
    while ((i < ancestors1.size()) || (i < ancestors2.size()))
        {
        if (i < ancestors1.size())
            dist += ancestors1[i]->attrs.GetDouble(ATTR_DISTANCE);
        if (i < ancestors2.size())
            dist += ancestors2[i]->attrs.GetDouble(ATTR_DISTANCE);
            
        ++i;
        }
    
    return dist;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

// Compare the position of insertion to the original position in a full tree
//      name - name of inserted node
//      parent - node in coreTree where "name" node has been inserted
//      origTree - original full tree containing correct position for "name" node
//      coreTree - core tree in which parent is a member
//  Returns the distance from current position (parent) to correct position
//
double ComparePosition(LPCSTR name, CPNode *parent, CPTree *origTree, CPTree *coreTree)
    {
    double     dist = -1;
    
    // Find the position in the original tree
    CPNode      *orig = origTree->Find(name);
    
    // Find the first parent with at least 2 branchs with core sequences
    string      core2;
    while (core2.empty() && (orig != NULL))
        {
        orig = orig->parent;
        if (orig != NULL)
            core2 = orig->attrs.GetString(CORE_NAME2);
        }
        
    if (orig == NULL)
        {
        //DisplayL("Could not find original position for [%s]\n", name);
        return -3;    // did not find original position
        }
        
    // Find the two core names associated with core ancestor
    string      core1 = orig->attrs.GetString(CORE_NAME1);
    
    // Find same core ancestor in core tree
    CPNode      *core = FindBranchContaining(coreTree, core1.c_str(), core2.c_str());
    if (core == NULL)
        {
        DisplayL("Could not find core position for [%s][%s]\n", core1.c_str(), core2.c_str());
        return -2;    // did not find core position
        }
        
    // if parent is not a branch node, get parent's parent
    if (!parent->IsBranchNode())
        parent = parent->parent;
        
    // Calc distance between core ancestor and insertion point
    dist = DistanceBetween(parent, core);
    //DisplayL("Distance error [%20.20s] = %d\n", name, dist);
    
    return dist;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////

// Compare the position of insertions to the original position in a full tree
void ComparePositions(CInsertPosArray &inserts, CPTree *origTree, CPTree *coreTree)
    {
        
    for (int i=0 ; i < inserts.size() ; ++i)
        {
        string                  name = inserts[i]->seq->name;        

        DisplayL("======================================================================\n"); 
        DisplayL("\nDistance to Insertion: [%20.20s] [%5d sites]\n", name.c_str(), inserts[i]->best.list.size()); 
        
        if (inserts[i]->best.list.size() == 0)
            {
            DisplayL("\tNO MATCHES found for [%s]\n", name.c_str());
            continue;
            }
                  
        CBestList::iterator		iterB = inserts[i]->best.list.begin();

        for (iterB=inserts[i]->best.list.begin() ; iterB != inserts[i]->best.list.end() ; ++iterB)
            {
            CPNode      *parent = iterB->node;
            string      pname = GetNodeName(parent);
            
            int         dist = ComparePosition(name.c_str(), parent, origTree, coreTree);

            DisplayL("\tDistance error from [%20.20s] = %d\n", pname.c_str(), dist);
            }
        }
                
    }
    
///////////////////////////////////////////////////////////////////////////////////////////////////

BOOL FullTree_Initialization()
    {
    if (tree == NULL) return FALSE;
    if (treeSeqfile == NULL) return FALSE;
    if (fullTreeName == NULL) return FALSE;
    
    // Read in original tree
    DisplayL("Reading Newick Tree: [%s]\n", fullTreeName);
    fullTree = ReadNewickTree(fullTreeName);
    if (fullTree == NULL)
        {
        DisplayL("Error accessing Newick Tree file [%s]\n", fullTreeName);
        return FALSE;
        }
    
    DisplayL("Newick Tree Read Completed: %d taxa\n", (fullTree->nodeList.size()+1)/2);
        
    // Set the core names into ancestors
    CalcCoreBranches(fullTree->root, treeSeqfile);
    
    return TRUE;
    }

void TestInsertPositions()
    {
    if (inserts.size() == 0) return;
    if (fullTree == NULL) 
        {
        if (!FullTree_Initialization())
            return;
        }
        
    ComparePositions(inserts, fullTree , tree);
    }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
