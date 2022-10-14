```
###############################################################################
#
# ParsInsert
#
# ParsInsert efficiently produces both a phylogenetic tree and taxonomic
# classification for sequences for microbial community sequence analysis. 
#
###############################################################################

      Summary: ParsInsert is a C++ implementation of Parsimonious Insertion. 
               The algorithm exploits the knowledge provided by publicly 
               available curated phylogenetic trees to efficiently insert new 
               sequences and infer their taxonomic classification.  The 
               ParsInsert placement of new sequences in the tree is 
               deterministic which allows for distributed processing to quickly 
               handle millions of reads.

 Availability: ParsInsert source code and documentation are available at:
                    http://parsinsert.sourceforge.net
		
      Contact: david.knox@colorado.edu

###############################################################################
##
##  Developer : David Knox (david.knox@colorado.edu) Jan 2011
##  Copyright : Copyright (C) 2007-2011 David Knox
##
##  Web site  : http://parsinsert.sourceforge.net/
##
###############################################################################
##	  This file is part of ParsInsert.
##
##    ParsInsert is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    ParsInsert is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with ParsInsert.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

##########################
Installation Instructions:
##########################

1. Download and unzip the distribution file to a directory of your choice.

2. Compile the application using the following commands from 
   the installation directory:
		make
   this will compile the application using the 'gcc' compiler.
	
3. Test the application.  There is a testing data set in the TestData directory.
		make TEST
	will run the application on the test data and compare the results to preiously 
	run test and results are reported.

##############################
Sample command for ParsInsert:
##############################
	ParsInsert -x rdp.taxonomy -t core_set.tree -s core_set_aligned.fasta unclassified_seq.fasta  


########
Support:
########

	Please send any comments, bug reports or questions to
		david.knox@colorado.edu


###################################
The directory contains these files:
###################################
   AttrList.cpp         - Source to manage sets of key,value pairs
   AttrList.h
   Attrs.h              - Predefined attribute names for Phylogenetic Trees
   GNU_license_agpl.txt - License Agreement (GNU AFFERO GENERAL PUBLIC LICENSE)
   Knox_Stddef.cpp      - Set of routines and old habit data types
   Knox_Stddef.h
   makefile             - rules for making application and performing tests
   ParsimonySet.h       - Handles the parsimony functions
   ParsInsert.cpp       - Source for main algorithm and application
   ParsInsert.h
   PNode.cpp            - Source to manage Phylogenetic Tree Node
   PNode.h
   ReleaseNotes.txt     - Release Notes for each release
   SeqList.cpp          - Source to manage sequence file functions
   SeqList.h
   Taxonomy.cpp         - Source to manage heirarchy of taxonomic classifications
   Taxonomy.h
   
   TestData             - Data used for testing the application installation
      core_rdp.fasta            - sequences for taxa in core tree
      core_rdp.ntree            - core tree
      lanemask.txt              - Lane Mask for important positions in the aligned sequences
      rdp.taxonomy              - taxonomy of taxa in core tree
      set100.fasta              - small test set
      set1000.fasta             - larger test set
      set1000.results           - results from ParsInsert on larger test set
      short1000_NAST.fasta      - sample set of short reads aligned
      short1000_unaligned.fasta - sample set of short reads unaligned
      TaxonomyCounts.txt        - List of the taxonomys in the core tree
   
###################################
###################################
   
```