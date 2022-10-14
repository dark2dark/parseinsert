#
# Makefie for ParsInsert
#
# ParsInsert efficiently produces both a phylogenetic tree and taxonomic classification 
# for sequences for microbial community sequence analysis. 
#
# This is a C++ implementation of the Parsimonious Insertion algorithm.
#
PROGRAM = ParsInsert
OBJECTS = AttrList.o Knox_Stddef.o PNode.o SeqList.o Taxonomy.o ParsInsert.o

CFLAGS   = -g -O3
CC       = g++
INCLUDES =
LIBS     = -lm -lc

.SUFFIXES:	.o .cpp

.cpp.o:
	$(CC) $(CFLAGS) -c -o $@ $<

all:    $(PROGRAM) 

clean: 
	rm -f $(PROGRAM) $(OBJECTS) *.fasta *.log *.results *.tree

TEST_OBJECTS = PNode.o AttrList.o Knox_Stddef.o Test.o
test:  $(TEST_OBJECTS)
	$(CC) $(CFLAG) -o Test $(TEST_OBJECTS) $(LIBS)

Knox_Stddef.o: Knox_Stddef.cpp Knox_Stddef.h

$(PROGRAM):	$(OBJECTS)
	$(CC) $(CFLAG) -o $(PROGRAM) $(OBJECTS) $(LIBS)

TEST_DIR = ./TestData
test_short = short1000_NAST
test_set = set1000

TAXONOMY = $(TEST_DIR)/rdp.taxonomy
TREE     = $(TEST_DIR)/core_rdp.ntree
TREE_SEQ = $(TEST_DIR)/core_rdp.fasta

# create compressed version of the release
RELEASE_VER = 1.04
RELEASE:
	tar -czf ../ParsInsert.$(RELEASE_VER).tgz * 

# use simple scoring function (-d0), threshold cutoff at 80% (-c80), only display best location (-n1)
test_options = -d0 -c80 -n1

TEST: $(PROGRAM)
	rm -f $(test_set).fasta
	cp TestData/$(test_set).fasta .
	./$(PROGRAM)  		$(test_options) \
						-m $(TEST_DIR)/lanemask.txt\
						-x $(TAXONOMY) \
						-l $(test_set).log \
						-o $(test_set).results \
						-t $(TREE) \
						-s $(TREE_SEQ) \
						$(test_set).fasta 
	diff $(test_set).results TestData/$(test_set).results
	echo Differences:
	diff $(test_set).results TestData/$(test_set).results | wc -l

