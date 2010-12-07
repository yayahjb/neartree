// scale.cpp
//
// scale reads from standard input and writes to standard output. It is intended
// only to read 3-space vectors (as csv) and scale the input by the
// input scale factor (from the command line). All directions are scaled equally.
//
// This is intended to be used with other programs for generating test data
// sets for NearTree.
//
//
#include "vector_3d.h"
#include <iostream>

#include "Data2CSV.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char* argv[])
//---------------------------------------------------------------------
{
	// make sure there is a scale factor on the command line
    if ( argc<2 )
        return( 0 );

    const double scale = atof( argv[1] );

    const std::vector<Vector_3> v = ReadVector_3File( std::cin );

    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        fprintf( stdout, "%f,%f,%f\n", scale*v[i][0], scale*v[i][1], scale*v[i][2] );
    }

    return 0;
}

