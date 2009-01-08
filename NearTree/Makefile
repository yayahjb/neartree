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
#


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
VERSION = 1:0:1
RELEASE = 1.0


#
# Compiler and compilation flags
#
CC	= gcc
CFLAGS  = -g -O2  -Wall -ansi -pedantic

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
# Uncomment the next line if CVector is installed locally in $(HOME)/include
# and $(HOME)/lib
#
#CVECTOR_INCLUDE = -I$(HOME)/include -L$(HOME)/lib

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
#INSTALLDIR = /usr/local
INSTALLDIR  = $(HOME)

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
CLIBRARIES = -lCNearTree -lCVector -lm

COMPILE_COMMAND        =  $(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(INCLUDES) $(WARNINGS) -c
CPPCOMPILE_COMMAND     =  $(LIBTOOL) --mode=compile $(CXX) $(CFLAGS) $(INCLUDES) $(WARNINGS) -c
LIBRARY_LINK_COMMAND   =  $(LIBTOOL) --mode=link  $(CC) -version-info $(VERSION) -release $(RELEASE) -rpath $(INSTALLDIR)/lib
BUILD_COMMAND_LOCAL    =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) $(INCLUDES) -L$(LIB)
CPPBUILD_COMMAND_LOCAL =  $(LIBTOOL) --mode=link $(CXX) $(CFLAGS) $(INCLUDES)
BUILD_COMMAND_DYNAMIC  =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) -dynamic -I $(INSTALLDIR)/include -L$(INSTALLDIR)/lib
BUILD_COMMAND_STATIC   =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) -static -I $(INSTALLDIR)/include -L$(INSTALLDIR)/lib
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
HEADERS   =  $(INC)/TNear.h $(INC)/CNearTree.h $(INC)/v.h

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
	@echo '   $(INSTALLDIR) '	
	@echo '   $(INSTALL_COMMAND) '	
	@echo '   $(INSTALL_FINISH) '	
	@echo ' '
	@echo ' To compile the CVector library and example programs type:'
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
	@echo ' To install the library and binaries type:'
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
		$(BIN)/CPPMain $(BIN)/CMain $(BIN)/CPPNearTreeTest $(BIN)/CNearTreeTest
		
install:  all $(INSTALLDIR) $(INSTALLDIR)/lib $(INSTALLDIR)/include \
          $(INC) $(LIB)/libCNearTree.$(LIB_EXT)  $(INC)/TNear.h $(INC)/CNearTree.h
		  $(INSTALL_COMMAND) $(LIB)/libCNearTree.$(LIB_EXT) $(INSTALLDIR)/lib/libCNearTree.$(LIB_EXT)
		  $(INSTALL_FINISH_COMMAND) $(INSTALLDIR)/lib/libCNearTree.$(LIB_EXT)
		  -cp $(INSTALLDIR)/include/TNear.h $(INSTALLDIR)/include/TNear_old.h
		  -cp $(INSTALLDIR)/include/CNearTree.h $(INSTALLDIR)/include/CNearTree_old.h
		  cp $(INC)/TNear.h $(INSTALLDIR)/include/TNear.h
		  cp $(INC)/CNearTree.h $(INSTALLDIR)/include/CNearTree.h
		  chmod 644 $(INSTALLDIR)/include/TNear.h
		  chmod 644 $(INSTALLDIR)/include/CNearTree.h
		  echo "Testing final install dynamic"
		  $(BUILD_COMMAND_DYNAMIC)  $(EXAMPLES)/CNearTreeTest.c \
		  $(CLIBRARIES) -o $(BIN)/CNearTreeTest_dynamic
		  $(BIN)/CNearTreeTest_dynamic 
		  echo "Testing final install static"
		  $(BUILD_COMMAND_STATIC) $(EXAMPLES)/CNearTreeTest.c \
		  $(CLIBRARIES) -o $(BIN)/CNearTreeTest_static
		  $(BIN)/CNearTreeTest_static 
			
		  
		  
#
# Directories
#
$(INSTALLDIR):
	mkdir -p $(INSTALLDIR)

$(INSTALLDIR)/lib:  $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/lib

$(INSTALLDIR)/bin:  $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/bin

$(INSTALLDIR)/include:  $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/include
	

$(LIB):
	mkdir $(LIB)

$(BIN):
	mkdir $(BIN)

#
# CNearTree library
#
$(LIB)/libCNearTree.$(LIB_EXT): $(SRC)/CNearTree.c $(INC)/CNearTree.h $(LIB) $(COMMONDEP)
	$(COMPILE_COMMAND) -c $(SRC)/CNearTree.c
	$(LIBRARY_LINK_COMMAND) -o $(LIB)/libCNearTree.$(LIB_EXT) CNearTree.$(OBJ_EXT) 

#
#  example programs
#
$(BIN)/CPPMain: $(INC)/TNear.h $(EXAMPLES)/main.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) $(EXAMPLES)/main.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		      -o $@
			  
$(BIN)/CPPNearTreeTest: $(INC)/TNear.h $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(BIN)
	$(CPPBUILD_COMMAND_LOCAL) $(EXAMPLES)/CNearTreeTest.cpp $(EXAMPLES)/v.cpp $(CPPLIBRARIES) \
		      -o $@
			  
$(BIN)/CMain: $(INC)/CNearTree.h $(EXAMPLES)/main.c $(SRC)/CNearTree.c $(BIN) $(LIB)/libCNearTree.$(LIB_EXT)
	$(BUILD_COMMAND_LOCAL) $(EXAMPLES)/main.c $(SRC)/CNearTree.c $(LIB)/libCNearTree.$(LIB_EXT) -lCVector -lm \
		      -o $@

$(BIN)/CNearTreeTest: $(INC)/CNearTree.h $(EXAMPLES)/CNearTreeTest.c $(SRC)/CNearTree.c $(BIN) $(LIB)/libCNearTree.$(LIB_EXT)
	$(BUILD_COMMAND_LOCAL) $(EXAMPLES)/CNearTreeTest.c $(SRC)/CNearTree.c $(LIB)/libCNearTree.$(LIB_EXT) -lCVector -lm \
		      -o $@


#
# Tests
#
tests:		$(LIB) $(BIN) $(BIN)/CPPMain $(BIN)/CPPNearTreeTest \
				$(BIN)/CMain $(BIN)/CNearTreeTest
			$(BIN)/CPPMain 3458923 > $(TESTDATA)/CPPMain.lst
			$(BIN)/CMain 3458923 > $(TESTDATA)/CMain.lst
			$(BIN)/CPPNearTreeTest
			$(BIN)/CNearTreeTest
			diff -b -c $(TESTDATA)/CPPMain_orig.lst $(TESTDATA)/CPPMain.lst
			diff -b -c $(TESTDATA)/CMain_orig.lst $(TESTDATA)/CMain.lst
			
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
	@-rm -f *.o
	@-rm -f *.$(OBJ_EXT)

#
# Restore to distribution state
#
distclean:	clean empty





