#
#  Makefile
#  NearTree
#
#  Part if the NearTree API
#  Copyright 2008, 2009 Larry Andrews.  All rights reserved
#
#  Makefile created by Herbert J. Bernstein on 1/08/09
#  with permission from Larry Andrews.
#  Copyright 2008 Larry Andrews and Herbert J. Bernstein. 
#  All rights reserved.
#
#  Revised 30 May 2009, release with full containerization of C++
#                       version and KNear/Far in C++ and C, LCA + HJB
#
#  Revised 31 October 2009, adding triple.h
#
#  Revised 8 October 2012, remove duplicate library load
#
#  Revised 21 January 2022, cleanup for gcc-11


#**********************************************************************
#*                                                                    *
#* YOU MAY REDISTRIBUTE THE NearTree API UNDER THE TERMS OF THE LGPL  *
#*                                                                    *
#**********************************************************************/

#************************* LGPL NOTICES *******************************
#*                                                                    *
#* This library is free software; you can redistribute it and/or      *
#* modify it under the terms of the GNU Lesser General Public         *
#* License as published by the Free Software Foundation; either       *
#* version 2.1 of the License, or (at your option) any later version. *
#*                                                                    *
#* This library is distributed in the hope that it will be useful,    *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
#* Lesser General Public License for more details.                    *
#*                                                                    *
#* You should have received a copy of the GNU Lesser General Public   *
#* License along with this library; if not, write to the Free         *
#* Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
#* MA  02110-1301  USA                                                *
#*                                                                    *
#**********************************************************************/

# Version string
VERSION = 7:0:1
RELEASE = 4.0.1
SHELL = /bin/bash

#
# Compiler and compilation flags
#
CC	?= gcc
CXX     ?= g++
OPTLEVEL = -O3
ifneq ($(MSYSTEM),MINGW32)
CFLAGS  = -g $(OPTLEVEL) -fopenmp  -Wall  -ansi -pedantic -DCNEARTREE_SAFE_TRIANG=1
TIME = time
else
CFLAGS  = -g -fopenmp -O2 -Wall  -DCNEARTREE_SAFE_TRIANG=1
TIME = 
endif
CXXFLAGS = $(CFLAGS) -std=c++11

#
# libtool path if system default is not suitable
#
#LIBTOOL = $(HOME)/bin/libtool
ifndef LIBTOOL
  LIBTOOL = libtool
endif

#
# If local headers must be quoted uncomment the next line
#
#USE_LOCAL_HEADERS = 1

#
CVECTOR_INCLUDE = -I$(NEARTREE_PREFIX)/include
CVECTOR_LIBLOC = -L$(NEARTREE_PREFIX)/lib

#
# Directories
#
ROOT     = .
LIB      = $(ROOT)/lib
BIN      = $(ROOT)/bin
SRC      = $(ROOT)
INC      = $(ROOT)
EXAMPLES = $(ROOT)
TESTDATA = $(ROOT)
NEARTREE_PREFIX ?= /usr/local

#
# Include directories
#
ifdef USE_LOCAL_HEADERS
INCLUDES = -DUSE_LOCAL_HEADERS
else
INCLUDES = -I$(INC) $(CVECTOR_INCLUDE)
endif

#
# C++ Libraries
#
CPPLIBRARIES = -lm


#
# C Libraries
#
CLIBRARIES = $(CVECTOR_LIBLOC) -lCVector -lm

COMPILE_COMMAND        =  $(LIBTOOL) --tag=CC  --mode=compile $(CC) $(CFLAGS) $(INCLUDES) $(WARNINGS) -c
CPPCOMPILE_COMMAND     =  $(LIBTOOL) --tag=CXX --mode=compile $(CXX) $(CFLAGS) $(INCLUDES) $(WARNINGS) -c
LIBRARY_LINK_COMMAND   =  $(LIBTOOL) --tag=CC  --mode=link $(CC) -version-info $(VERSION) -no-undefined -rpath $(NEARTREE_PREFIX)/lib
BUILD_COMMAND_LOCAL    =  $(LIBTOOL) --tag=CC  --mode=link $(CC) $(CFLAGS) $(INCLUDES)
CPPBUILD_COMMAND_LOCAL =  $(LIBTOOL) --tag=CXX --mode=link $(CXX) -no-undefined $(CFLAGS) $(INCLUDES)
BUILD_COMMAND_DYNAMIC  =  $(LIBTOOL) --tag=CC  --mode=link $(CC) -no-undefined $(CFLAGS) -shared -I$(NEARTREE_PREFIX)/include
BUILD_COMMAND_STATIC   =  $(LIBTOOL) --tag=CC  --mode=link $(CC) $(CFLAGS) -static-libtool-libs -I$(NEARTREE_PREFIX)/include
INSTALL_COMMAND        =  $(LIBTOOL) --mode=install cp
INSTALL_FINISH_COMMAND =  $(LIBTOOL) --mode=finish

OBJ_EXT                =  lo
LIB_EXT                =  la

######################################################################
#  You should not need to make modifications below this line         #
######################################################################

#
# Suffixes of files to be used or built
#
.SUFFIXES:	.c .$(OBJ_EXT) .$(LIB_EXT)

#
# Common dependencies
#
COMMONDEP = Makefile

#
# Source files
#
CPPSOURCE =	$(SRC)/main.cpp $(SRC)/CNearTreeTest.cpp $(SRC)/v.cpp
CSOURCE   =	$(SRC)/CNearTree.c $(SRC)/main.c $(SRC)/CNearTreeTest.c
SOURCE   =  $(CPPSOURCE) $(CSOURCE)

#
# Header files
#
HEADERS   =  $(INC)/TNear.h $(INC)/CNearTree.h $(INC)/v.h $(INC)/rhrand.h $(INC)/triple.h

# Default: instructions
#
default:
	@echo ' '
	@echo '***************************************************************'
	@echo ' '
	@echo ' PLEASE READ README_NearTree.txt and lgpl.txt'
	@echo ' '
	@echo ' Before making the NearTree libraries and example programs, check'
	@echo ' that the chosen settings are correct'
	@echo ' '
	@echo ' The current C++ and C compile commands are:'
	@echo ' '
	@echo '   $(CPPCOMPILE_COMMAND)'
	@echo '   $(COMPILE_COMMAND)'
	@echo ' '
	@echo ' The C API, CNearTree.c, depends on the sourceforge project CVector '
ifdef CVECTOR_INCLUDE
	@echo ' The variable CVECTOR_INCLUDE is defined as '
	@echo '   $(CVECTOR_INCLUDE)'
else
	@echo ' You are currently setup to use the system defaults for CVector'
	@echo ' If that is not correct, define the variable CVECTOR_INCLUDE '
endif
	@echo ' '
	@echo ' The current library link command is:'
	@echo ' '
	@echo '   $(LIBRARY_LINK_COMMAND)'
	@echo ' '
	@echo ' The current C++ and C library local, and C dynamic and static build commands are:'
	@echo ' '
	@echo '   $(CPPBUILD_COMMAND_LOCAL)'
	@echo '   $(BUILD_COMMAND_LOCAL)'
	@echo '   $(BUILD_COMMAND_DYNAMIC)'
	@echo '   $(BUILD_COMMAND_STATIC)'
	@echo ' '
	@echo ' Before installing the NearTree library and example programs, check'
	@echo ' that the install directory and install commands are correct:'
	@echo ' '
	@echo ' The current values are :'
	@echo ' '
	@echo '   $(NEARTREE_PREFIX) '	
	@echo '   $(INSTALL_COMMAND) '	
	@echo '   $(INSTALL_FINISH) '	
	@echo ' '
	@echo ' To compile the NearTree library and example programs type:'
	@echo ' '
	@echo '   make clean'
	@echo '   make all'
	@echo ' '
	@echo ' To run a set of tests type:'
	@echo ' '
	@echo '   make tests'
	@echo ' '
	@echo ' To clean up the directories type:'
	@echo ' '
	@echo '   make clean'
	@echo ' '
	@echo ' To install the library and headers type:'
	@echo ' '
	@echo '   make install'
	@echo ' '
	@echo '***************************************************************'
	@echo ' '

#
# Compile the library and examples
#
all:	$(LIB) $(BIN) $(SOURCE) $(HEADERS) \
		$(LIB)/libCNearTree.$(LIB_EXT) \
		$(BIN)/CPPMain $(BIN)/CMain $(BIN)/CPPNearTreeTest $(BIN)/CNearTreeTest $(BIN)/CNearTreeTestInst 
		
install:  all $(NEARTREE_PREFIX) $(NEARTREE_PREFIX)/lib $(NEARTREE_PREFIX)/include \
		  $(INC) $(LIB)/libCNearTree.$(LIB_EXT)  $(INC)/TNear.h $(INC)/CNearTree.h \
		  $(INC)/rhrand.h $(INC)/triple.h
		  $(INSTALL_COMMAND) $(LIB)/libCNearTree.$(LIB_EXT) $(NEARTREE_PREFIX)/lib/libCNearTree.$(LIB_EXT)
		  $(INSTALL_FINISH_COMMAND) $(NEARTREE_PREFIX)/lib/libCNearTree.$(LIB_EXT)
		  -cp $(NEARTREE_PREFIX)/include/TNear.h $(NEARTREE_PREFIX)/include/TNear_old.h
		  -cp $(NEARTREE_PREFIX)/include/CNearTree.h $(NEARTREE_PREFIX)/include/CNearTree_old.h
		  cp $(INC)/TNear.h $(NEARTREE_PREFIX)/include/TNear.h
		  cp $(INC)/rhrand.h $(NEARTREE_PREFIX)/include/rhrand.h
		  cp $(INC)/triple.h $(NEARTREE_PREFIX)/include/triple.h
		  cp $(INC)/CNearTree.h $(NEARTREE_PREFIX)/include/CNearTree.h
		  chmod 644 $(NEARTREE_PREFIX)/include/TNear.h
		  chmod 644 $(NEARTREE_PREFIX)/include/CNearTree.h
		  chmod 644 $(NEARTREE_PREFIX)/include/rhrand.h
		  chmod 644 $(NEARTREE_PREFIX)/include/triple.h
		  echo "Testing final install static"
		  $(BUILD_COMMAND_STATIC) $(INCLUDES) $(EXAMPLES)/CNearTreeTest.c \
		  -L$(NEARTREE_PREFIX)/lib -lCNearTree $(CLIBRARIES) -o $(BIN)/CNearTreeTest_static
		  $(BIN)/CNearTreeTest_static 
		  echo "Testing final install dynamic"
		  $(BUILD_COMMAND_DYNAMIC) $(INCLUDES) $(EXAMPLES)/CNearTreeTest.c \
		  -L$(NEARTREE_PREFIX)/lib -lCNearTree $(CLIBRARIES) -o $(BIN)/CNearTreeTest_dynamic
		  $(BIN)/CNearTreeTest_dynamic 
		
		  
		  
#
# Directories
#
$(NEARTREE_PREFIX):
	mkdir -p $(NEARTREE_PREFIX)

$(NEARTREE_PREFIX)/lib:  $(NEARTREE_PREFIX)
	mkdir -p $(NEARTREE_PREFIX)/lib

$(NEARTREE_PREFIX)/bin:  $(NEARTREE_PREFIX)
	mkdir -p $(NEARTREE_PREFIX)/bin

$(NEARTREE_PREFIX)/include:  $(NEARTREE_PREFIX)
	mkdir -p $(NEARTREE_PREFIX)/include
	

$(LIB):
	mkdir $(LIB)

$(BIN):
	mkdir $(BIN)

#
# CNearTree library
#
$(LIB)/libCNearTree.$(LIB_EXT): $(SRC)/CNearTree.c $(INC)/CNearTree.h $(LIB) $(COMMONDEP)
	$(COMPILE_COMMAND) -c $(SRC)/CNearTree.c
	$(LIBRARY_LINK_COMMAND) -o $(LIB)/libCNearTree.$(LIB_EXT) CNearTree.$(OBJ_EXT)  $(INCLUDES) \
	$(CLIBRARIES)

#
#  example programs
#
$(BIN)/CPPMain: $(INC)/TNear.h $(EXAMPLES)/main.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) $(EXAMPLES)/main.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@

$(BIN)/CPPMain_Flip: $(INC)/TNear.h $(EXAMPLES)/main_flip.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) $(EXAMPLES)/main_flip.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@
		 
$(BIN)/CPPNearTreeTest: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@

$(BIN)/CPPNearTreeTest_2SAMP: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
		$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_DIMSAMPLES=2 $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@

$(BIN)/CPPNearTreeTest_6SAMP: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_DIMSAMPLES=6 $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@


$(BIN)/CPPNearTreeTest_8SAMP: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
		$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_DIMSAMPLES=8 $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@


$(BIN)/CPPNearTreeTestInst: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_INSTRUMENTED $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
	  -o $@

$(BIN)/CPPNearTreeTestInstFlip: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_INSTRUMENTED -DCNEARTREE_FORCEFLIP $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
	  -o $@
		  
$(BIN)/CPPNearTreeTestInstNoDefer: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_INSTRUMENTED -DCNEARTREE_NODEFER $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@

$(BIN)/CPPNearTreeTestInstNoFlip: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_INSTRUMENTED -DCNEARTREE_NOFLIP $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		      -o $@

$(BIN)/CPPNearTreeTestInstNoFlipNoDefer: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) -DCNEARTREE_INSTRUMENTED -DCNEARTREE_NOFLIP -DCNEARTREE_NODEFER $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		  -o $@
		  
$(BIN)/CMain: $(INC)/CNearTree.h $(EXAMPLES)/main.c $(SRC)/CNearTree.c $(BIN) $(LIB)/libCNearTree.$(LIB_EXT)
	$(BUILD_COMMAND_LOCAL) $(EXAMPLES)/main.c  $(LIB)/libCNearTree.$(LIB_EXT) $(CLIBRARIES) \
		  -o $@

$(BIN)/CNearTreeTest: $(INC)/CNearTree.h $(EXAMPLES)/CNearTreeTest.c $(SRC)/CNearTree.c $(BIN) $(LIB)/libCNearTree.$(LIB_EXT)
	$(BUILD_COMMAND_LOCAL) $(EXAMPLES)/CNearTreeTest.c $(LIB)/libCNearTree.$(LIB_EXT) $(CLIBRARIES) \
		  -o $@

$(BIN)/CNearTreeTestInst: $(INC)/CNearTree.h $(EXAMPLES)/CNearTreeTest.c $(SRC)/CNearTree.c $(BIN) $(LIB)/libCNearTree.$(LIB_EXT)
	$(BUILD_COMMAND_LOCAL) -DCNEARTREE_INSTRUMENTED $(EXAMPLES)/CNearTreeTest.c $(LIB)/libCNearTree.$(LIB_EXT) $(CLIBRARIES) \
		  -o $@


#
# Tests
#
tests:		$(LIB) $(BIN) $(BIN)/CPPMain $(BIN)/CPPMain_Flip $(BIN)/CPPNearTreeTest \
		$(BIN)/CPPNearTreeTest_2SAMP \
		$(BIN)/CPPNearTreeTest_6SAMP \
		$(BIN)/CPPNearTreeTest_8SAMP \
		$(BIN)/CPPNearTreeTestInst $(BIN)/CMain $(BIN)/CNearTreeTest \
		$(BIN)/CPPNearTreeTestInstFlip $(BIN)/CPPNearTreeTestInstNoDefer \
		$(BIN)/CPPNearTreeTestInstNoFlip $(BIN)/CPPNearTreeTestInstNoFlipNoDefer \
		$(BIN)/CNearTreeTestInst
		$(TIME) $(BIN)/CPPMain 3458923 > $(TESTDATA)/CPPMain.lst
		$(TIME) $(BIN)/CPPMain_Flip 3458923 > $(TESTDATA)/CPPMain_Flip.lst
		$(TIME) $(BIN)/CMain 3458923 > $(TESTDATA)/CMain.lst
		$(TIME) $(BIN)/CPPNearTreeTest_2SAMP
		$(TIME) $(BIN)/CPPNearTreeTest
		$(TIME) $(BIN)/CPPNearTreeTest_6SAMP
		$(TIME) $(BIN)/CPPNearTreeTest_8SAMP
		$(TIME) $(BIN)/CPPNearTreeTestInst
		$(TIME) $(BIN)/CPPNearTreeTestInstFlip
		$(TIME) $(BIN)/CPPNearTreeTestInstNoDefer
		$(TIME) $(BIN)/CPPNearTreeTestInstNoFlip
		$(TIME) $(BIN)/CNearTreeTest --debug
		$(TIME) $(BIN)/CNearTreeTestInst --debug
ifneq ($(MSYSTEM),MINGW32)
		diff -b -c $(TESTDATA)/CPPMain_orig.lst $(TESTDATA)/CPPMain.lst
		diff -b -c $(TESTDATA)/CPPMain_orig.lst $(TESTDATA)/CPPMain_Flip.lst
		diff -b -c $(TESTDATA)/CMain_orig.lst $(TESTDATA)/CMain.lst
else
		diff -b -c $(TESTDATA)/CPPMain_MINGW_orig.lst $(TESTDATA)/CPPMain.lst
		diff -b -c $(TESTDATA)/CPPMain_MINGW_orig.lst $(TESTDATA)/CPPMain_Flip.lst
		diff -b -c $(TESTDATA)/CMain_MINGW_orig.lst $(TESTDATA)/CMain.lst
endif
		
#
# Remove all non-source files
#
empty:
		  @-rm -rf $(LIB)
		  @-rm -rf $(BIN)
		  
#
# Remove temporary files
#
clean:
	@-rm -f core 
	@-rm -f *.$(OBJ_EXT)
	@-rm -f *.c.*
	@-rm -rf .libs
	@-rm -f CMain.lst
	@-rm -f CPPMain.lst

#
# Restore to distribution state
#
distclean:	clean empty





