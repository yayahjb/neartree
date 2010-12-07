#
#  Makefile
#  DataGen
#
#  Part of the DataGen Package
#  Copyright 2010 Larry Andrews.  All rights reserved
#
#  Makefile created by Herbert J. Bernstein on 11/29/10
#  with permission from Larry Andrews.
#  Copyright 2010 Larry Andrews and Herbert J. Bernstein. 
#  All rights reserved.
#

#######################################################################
#                                                                     #
# YOU MAY REDISTRIBUTE THE DataGen Package UNDER THE TERMS OF THE GPL #
#                                                                     #
# ALTERNATIVELY YOU MAY REDISTRIBUTE DataGen Package as an API UNDER  #
# THE TERMS OF THE LGPL                                               #
#                                                                     #
#######################################################################

########################### GPL NOTICES ##############################
#                                                                    #
# This program is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU General Public License as     #
# published by the Free Software Foundation; either version 2 of     #
# (the License, or (at your option) any later version.               #
#                                                                    #
# This program is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      #
# GNU General Public License for more details.                       #
#                                                                    #
# You should have received a copy of the GNU General Public License  #
# along with this program; if not, write to the Free Software        #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA           #
# 02111-1307  USA                                                    #
#                                                                    #
######################################################################

######################### LGPL NOTICES ###############################
#                                                                    #
# This library is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU Lesser General Public         #
# License as published by the Free Software Foundation; either       #
# version 2.1 of the License, or (at your option) any later version. #
#                                                                    #
# This library is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# Lesser General Public License for more details.                    #
#                                                                    #
# You should have received a copy of the GNU Lesser General Public   #
# License along with this library; if not, write to the Free         #
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    #
# MA  02110-1301  USA                                                #
#                                                                    #
######################################################################


######################################################################
# Programs
#
# DataGen	
# Create data sets
# Creates both 3-D and n-D datasets
#
# $(BIN)/TestMaxDist	
# Analyze data sets by various criteria
# Accepts n-D data files
#
# FixedTestSet	
# Analyzes a dataset using a specified set of probe data
# Accepts n-D data files
#
# NNStat	
# Provides some statistics about a csv file’s data
#
# Scale	
# Change the dimensions of the data
# Accepts and writes 3-D data files
#
# Rotateq	
# Rotate the data set (around the origin). The rotation 
# is specified by a quaternion on the command line. 
# (The quaternion will be normalized to a unit quaternion.)
# Accepts and writes 3-D data files
#
# Rotatea	
# Rotate the data set (around the origin). The rotation is 
# specified by an angle (in degrees) and 3 vector components
#
# BoxIt	
# Rescale and translate a dataset so that it fits within a 
# box of –1 to +1 (for every axis).
# Accepts and writes n-D data files
#
# Transl	
# Translate the points of a dataset (the 3 axial translations 
# are specified on the command line).
# Accepts and writes 3-D data files
#
# Svg	
# Write an SVG file from the input. The first two coordinates 
# (only) are plotted.
# Accepts n-D data files, but only uses the first two
#
######################################################################

#
# Compiler and compilation flags
#
CC	= gcc
ifneq ($(MSYSTEM),MINGW32)
CFLAGS  = -g -O2  -Wall  -ansi -pedantic
else
CFLAGS  = -g -O2 -Wall
endif

#
# libtool path if system default is not suitable
#
#LIBTOOL = $(HOME)/bin/libtool
ifndef LIBTOOL
  LIBTOOL = libtool
endif

#
# Directories
#
ROOT     = .
LIB      = $(ROOT)/lib
BIN      = $(ROOT)/bin
SRC      = $(ROOT)/src
INC      = $(ROOT)/include
EXAMPLES = $(ROOT)/examples
TESTDATA = $(ROOT)/testdata
ifndef INSTALL_PREFIX
INSTALL_PREFIX = /usr/local
#INSTALL_PREFIX  = $(HOME)
endif

OBJ_EXT   =  lo
LIB_EXT   =  la
INCLUDES  =  -I $(SRC)


CPPCOMPILE_COMMAND     =  $(LIBTOOL) --mode=compile $(CXX) $(CFLAGS) $(INCLUDES) $(WARNINGS) -c
CPPBUILD_COMMAND_LOCAL =  $(LIBTOOL) --mode=link $(CXX) -no-undefined $(CFLAGS) $(INCLUDES)
BUILD_COMMAND_LOCAL    =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) $(INCLUDES)
BUILD_COMMAND_LOCAL_STATIC    =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) -static-libtool-libs $(INCLUDES)

examples:	all \
			$(EXAMPLES)/geodb200.svg     \
			$(EXAMPLES)/randb.svg        \
			$(EXAMPLES)/ransp.svg        \
			$(EXAMPLES)/ranusp.svg       \
			$(EXAMPLES)/geosp500.svg     \
			$(EXAMPLES)/geospRot500.svg  \
			$(EXAMPLES)/qrgeor500.svg    \
			$(EXAMPLES)/LAT3-25.svg      \
			$(EXAMPLES)/LAT3Rot-8.svg    \
			$(EXAMPLES)/HAM200.svg       \
			$(EXAMPLES)/HAM200ROT.svg    \
			$(EXAMPLES)/cloud.svg        \
			$(EXAMPLES)/ranbox.svg       \
			$(EXAMPLES)/LAT4.csv         \
			$(EXAMPLES)/HAM.csv          \
			$(EXAMPLES)/CLOUD.csv        \
			$(EXAMPLES)/rot.svg          \
			$(EXAMPLES)/alphaC.csv       \
			$(EXAMPLES)/both1.svg        \
			$(EXAMPLES)/both2.csv        \
			$(EXAMPLES)/both2ROT.svg     \
			$(EXAMPLES)/both2.svg


Data2CSV.$(OBJ_EXT):	$(SRC)/Data2CSV.cpp $(SRC)/Data2CSV.h
	$(CPPCOMPILE_COMMAND) $(SRC)/Data2CSV.cpp

geodesic.$(OBJ_EXT):	$(SRC)/geodesic.cpp $(SRC)/geodesic.h
	$(CPPCOMPILE_COMMAND) $(SRC)/geodesic.cpp
	
hammersley.$(OBJ_EXT):	$(SRC)/hammersley.cpp $(SRC)/hammersley.h
	$(CPPCOMPILE_COMMAND) $(SRC)/hammersley.cpp
	
vector_3d.$(OBJ_EXT):	$(SRC)/vector_3d.cpp $(SRC)/vector_3d.h
	$(CPPCOMPILE_COMMAND) $(SRC)/vector_3d.cpp


all:	$(BIN)/DataGen      \
		$(BIN)/TestMaxDist  \
		$(BIN)/FixedTestSet \
		$(BIN)/NNStat       \
		$(BIN)/Scale        \
		$(BIN)/Rotateq      \
		$(BIN)/Rotatea      \
		$(BIN)/BoxIt        \
		$(BIN)/Transl       \
		$(BIN)/Svg

#  DataGen -- Create data sets

$(BIN)/DataGen:	$(SRC)/DataGen.cpp \
		geodesic.$(OBJ_EXT) \
		vector_3d.$(OBJ_EXT) \
		Data2CSV.$(OBJ_EXT) \
		hammersley.$(OBJ_EXT) \
		vector_3d.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/DataGen.cpp
	$(CPPBUILD_COMMAND_LOCAL) DataGen.$(OBJ_EXT) \
	geodesic.$(OBJ_EXT) \
	hammersley.$(OBJ_EXT) \
	vector_3d.$(OBJ_EXT) \
	Data2CSV.$(OBJ_EXT) -o $@

#  $(BIN)/TestMaxDist -- Analyze data sets by various criteria

$(BIN)/TestMaxDist:	src/TestMaxDist.cpp \
		Data2CSV.$(OBJ_EXT) \
		hammersley.$(OBJ_EXT) \
		vector_3d.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/TestMaxDist.cpp
	$(CPPBUILD_COMMAND_LOCAL) TestMaxDist.$(OBJ_EXT) \
	hammersley.$(OBJ_EXT) \
	vector_3d.$(OBJ_EXT) \
	Data2CSV.$(OBJ_EXT) -o $@
	
#  FixedTestSet -- Analyzes a dataset using a specified set of probe data

$(BIN)/FixedTestSet:	src/FixedTestSet.cpp \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/FixedTestSet.cpp
	$(CPPBUILD_COMMAND_LOCAL) FixedTestSet.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) -o $@

#  NNStat -- Provides some statistics about a csv file’s data

$(BIN)/NNStat:	src/NNStat.cpp \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/NNStat.cpp
	$(CPPBUILD_COMMAND_LOCAL) NNStat.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) -o $@

#  Scale -- Change the dimensions of the data

$(BIN)/Scale:	src/Scale.cpp \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/Scale.cpp
	$(CPPBUILD_COMMAND_LOCAL) Scale.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) -o $@

#  Rotateq -- Rotate the data set (around the origin). The rotation 
#  is specified by a quaternion on the command line. (The quaternion 
#  will be normalized to a unit quaternion.)

$(BIN)/Rotateq:	src/Rotateq.cpp \
		vector_3d.$(OBJ_EXT) \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/Rotateq.cpp
	$(CPPBUILD_COMMAND_LOCAL) Rotateq.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) vector_3d.$(OBJ_EXT) -o $@

#  Rotatea -- Rotate the data set (around the origin). The rotation 
#  is specified by an angle (in degrees) and 3 vector components

$(BIN)/Rotatea:	src/Rotatea.cpp \
		vector_3d.$(OBJ_EXT) \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/Rotatea.cpp
	$(CPPBUILD_COMMAND_LOCAL) Rotatea.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) vector_3d.$(OBJ_EXT) -o $@
	
#  BoxIt -- Rescale and translate a dataset so that it fits within a 
#  box of –1 to +1 (for every axis).

$(BIN)/BoxIt:	src/BoxIt.cpp \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/BoxIt.cpp
	$(CPPBUILD_COMMAND_LOCAL) BoxIt.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) -o $@

#  Transl -- Translate the points of a dataset (the 3 axial 
#  translations are specified on the command line).

$(BIN)/Transl:	src/Transl.cpp \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/Transl.cpp
	$(CPPBUILD_COMMAND_LOCAL) Transl.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) -o $@

#  Svg -- Write an SVG file from the input. The first two coordinates 
#  (only) are plotted.

$(BIN)/Svg:	src/Svg.cpp \
		Data2CSV.$(OBJ_EXT)
	$(CPPCOMPILE_COMMAND) $(SRC)/Svg.cpp
	$(CPPBUILD_COMMAND_LOCAL) Svg.$(OBJ_EXT) Data2CSV.$(OBJ_EXT) -o $@
	
# Examples
	
$(EXAMPLES)/geodb200.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen geodb 200 1 5 | $(BIN)/Svg 500 > $(EXAMPLES)/geodb200.svg
	
$(EXAMPLES)/randb.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen randb 200 1 5 | $(BIN)/Svg 500 > $(EXAMPLES)/randb500.svg

$(EXAMPLES)/ransp.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen ransp 500 | $(BIN)/Svg 500 > $(EXAMPLES)/ransp500.svg

$(EXAMPLES)/ranusp.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen ranusp 500 | $(BIN)/Svg 500 > $(EXAMPLES)/ranusp500.svg

$(EXAMPLES)/geosp500.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen geosp 500 | $(BIN)/Svg 500 > $(EXAMPLES)/geosp500.svg

$(EXAMPLES)/geospRot500.svg: $(BIN)/DataGen $(BIN)/Rotateq $(BIN)/Svg
	$(BIN)/DataGen geosp 500 | $(BIN)/Rotateq .2 1 1 1 | $(BIN)/Svg 500 > $(EXAMPLES)/geospRot500.svg

$(EXAMPLES)/qrgeor500.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen qrgeo 500 | $(BIN)/Svg 500 > $(EXAMPLES)/qrgeo500.svg

$(EXAMPLES)/LAT3-25.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen LAT3 25 | $(BIN)/Svg > $(EXAMPLES)/LAT3-25.svg

$(EXAMPLES)/LAT3Rot-8.svg: $(BIN)/DataGen $(BIN)/Rotateq $(BIN)/Svg
	$(BIN)/DataGen LAT3 8 | $(BIN)/Rotateq .2 1 1 1 | $(BIN)/Svg 500 > $(EXAMPLES)/LAT3Rot-8.svg

$(EXAMPLES)/HAM200.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen HAM 3 200 | $(BIN)/Svg > $(EXAMPLES)/HAM200.svg

$(EXAMPLES)/HAM200ROT.svg: $(BIN)/DataGen $(BIN)/Rotateq $(BIN)/Svg
	$(BIN)/DataGen HAM 3 200 | $(BIN)/Rotateq .2 1 1 1 | $(BIN)/Svg > $(EXAMPLES)/HAM200ROT.svg

$(EXAMPLES)/cloud.svg: $(BIN)/DataGen $(BIN)/Svg
	$(BIN)/DataGen CLOUD 15 2 9 | $(BIN)/Svg > $(EXAMPLES)/cloud.svg

$(EXAMPLES)/ranbox.svg: $(BIN)/DataGen $(BIN)/Rotateq $(BIN)/Svg
	$(BIN)/DataGen ranbox 3 1500 | $(BIN)/Rotateq .95 .1 .1 .05 |$(BIN)/Svg 500 > $(EXAMPLES)/ranbox.svg

$(EXAMPLES)/LAT4.csv: $(BIN)/DataGen $(BIN)/TestMaxDist
	$(BIN)/TestMaxDist ? > $(EXAMPLES)/LAT4.csv
	$(BIN)/DataGen LAT4 15 | $(BIN)/TestMaxDist >> $(EXAMPLES)/LAT4.csv
	$(BIN)/DataGen LAT4 20 | $(BIN)/TestMaxDist >> $(EXAMPLES)/LAT4.csv
	$(BIN)/DataGen LAT4 25 | $(BIN)/TestMaxDist >> $(EXAMPLES)/LAT4.csv
	$(BIN)/DataGen LAT4 30 | $(BIN)/TestMaxDist >> $(EXAMPLES)/LAT4.csv

$(EXAMPLES)/HAM.csv: $(BIN)/DataGen $(BIN)/TestMaxDist
	$(BIN)/TestMaxDist ? > $(EXAMPLES)/HAM.csv
	$(BIN)/DataGen HAM 1 10000 | $(BIN)/TestMaxDist >> $(EXAMPLES)/HAM.csv
	$(BIN)/DataGen HAM 2 10000 | $(BIN)/TestMaxDist >> $(EXAMPLES)/HAM.csv
	$(BIN)/DataGen HAM 3 10000 | $(BIN)/TestMaxDist >> $(EXAMPLES)/HAM.csv
	$(BIN)/DataGen HAM 4 10000 | $(BIN)/TestMaxDist >> $(EXAMPLES)/HAM.csv
	$(BIN)/DataGen HAM 5 10000 | $(BIN)/TestMaxDist >> $(EXAMPLES)/HAM.csv

$(EXAMPLES)/CLOUD.csv: $(BIN)/DataGen $(BIN)/TestMaxDist
	$(BIN)/TestMaxDist ? > $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 3 7 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 3 8 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 3 9 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 4 6 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 4 7 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 4 8 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 5 6 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	$(BIN)/DataGen CLOUD 5 7 2 | $(BIN)/TestMaxDist >> $(EXAMPLES)/CLOUD.csv
	
$(EXAMPLES)/rot.svg: $(BIN)/DataGen $(BIN)/Rotateq $(BIN)/Transl $(BIN)/Svg
	$(BIN)/DataGen lat3 100 | $(BIN)/Rotateq 0.5 1 1 1 > $(EXAMPLES)/out1.txt
	$(BIN)/DataGen lat3 100 | $(BIN)/Transl 0.5 1 1 1 >> $(EXAMPLES)/out1.txt
	$(BIN)/DataGen geodb 5000 1 5 |  $(BIN)/Transl 5 5 0 >> $(EXAMPLES)/out1.txt
	cat $(EXAMPLES)/out1.txt | $(BIN)/Svg > $(EXAMPLES)/rot.svg
	
$(EXAMPLES)/alphaC.csv: $(BIN)/DataGen $(BIN)/Scale $(BIN)/BoxIt \
	$(BIN)/Transl $(BIN)/TestMaxDist $(BIN)/Svg
	$(BIN)/TestMaxDist ? > $(EXAMPLES)/alphaC.csv
	$(BIN)/DataGen cloud 9 5 1 | $(BIN)/Transl 200 0 0 > $(EXAMPLES)/cloud1.txt
	$(BIN)/DataGen cloud 9 6 2 > $(EXAMPLES)/cloud2.txt
	cat $(EXAMPLES)/cloud1.txt $(EXAMPLES)/cloud2.txt > $(EXAMPLES)/cloud.txt
	cat $(EXAMPLES)/cloud.txt | $(BIN)/BoxIt | $(BIN)/TestMaxDist >> $(EXAMPLES)/alphaC.csv
	
$(EXAMPLES)/both1.svg: $(BIN)/DataGen $(BIN)/Transl $(BIN)/Scale $(BIN)/Svg
	$(BIN)/DataGen geodb 400 1 5 | $(BIN)/Transl .5 -.5 0 > $(EXAMPLES)/geodb.csv
	$(BIN)/DataGen ranusp 4000 | $(BIN)/Scale .5 > $(EXAMPLES)/sphere.csv
	cat $(EXAMPLES)/geodb.csv $(EXAMPLES)/sphere.csv > $(EXAMPLES)/both1.csv
	cat $(EXAMPLES)/both1.csv | $(BIN)/Svg 500 > $(EXAMPLES)/both1.svg
	
$(EXAMPLES)/both2.csv: $(BIN)/DataGen $(BIN)/Transl $(BIN)/Rotateq $(BIN)/Scale $(BIN)/Svg
	$(BIN)/DataGen geodb 400 1 5 | $(BIN)/Rotateq .998 0 -.063 0 | \
	$(BIN)/Rotateq .998 0 -.063 0 | $(BIN)/Rotateq .998 -.063 0 0 | \
	$(BIN)/Transl .8 -.5 -1 > $(EXAMPLES)/geodb.csv
	$(BIN)/DataGen ranusp 400 | $(BIN)/Scale .5 > $(EXAMPLES)/sphere.csv
	cat $(EXAMPLES)/geodb.csv $(EXAMPLES)/sphere.csv > $(EXAMPLES)/both2.csv
	
$(EXAMPLES)/both2ROT.svg: $(BIN)/Rotateq $(BIN)/Svg
	cat $(EXAMPLES)/both2.csv | $(BIN)/Rotateq .998 0 -.063 0 |\
		$(BIN)/Rotateq .998 0 -.063 0 | $(BIN)/Svg 500 > $(EXAMPLES)/bothROT.svg

$(EXAMPLES)/both2.svg: $(BIN)/Rotateq $(BIN)/Svg
	cat $(EXAMPLES)/both2.csv | $(BIN)/Rotateq .998 0 -.063 0 |\
		$(BIN)/Svg 500 > $(EXAMPLES)/both2.svg



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
	@-rm -f *.c.*
	@-rm -rf .libs
	@-rm -f $(SRC)/*.o
	@-rm -f $(SRC)/*.$(OBJ_EXT)
	@-rm -r -f $(SRC)/.libs


#
# Restore to distribution state
#
distclean:	clean empty
