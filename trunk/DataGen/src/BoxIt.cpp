// BoxIt.cpp
//
// BoxIt reads from standard input and writes to standard output. It is intended
// to read n-space vectors (as csv) and scale the input into a box bounded
// by -1 to +1 for each of the axes. All directions are scaled equally.
//
// This is intended to be used with other programs for generating test data
// sets for NearTree.
//
#include <iostream>
#include <cmath>

#include "Data2CSV.h"

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int main(int argc, char* argv[])
//---------------------------------------------------------------------
{
    std::vector<vecN> v = ReadGeneralFile( std::cin );
    BoxIt( v );
    WriteVectorFile( v );
	return 0;
}

