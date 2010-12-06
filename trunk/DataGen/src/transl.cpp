// transl.cpp
//
// transl reads from standard input and writes to standard output. It is intended
// only to read 3-space vectors (as csv) and shift the input by the
// 3-vector specified on the command line.
//
// This is intended to be used with other programs for generating test data
// sets for NearTree.
//

#include <iostream>

#include "Data2CSV.h"
#include "Vector_3d.h"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char* argv[])
//---------------------------------------------------------------------
{
	// assure that there are at least 3 translations on the command line
    if ( argc<4 )
        return( 0 );

    const Vector_3 trans = Vector_3( atof( argv[1] ), atof( argv[2] ), atof( argv[3] ) );

    std::vector<Vector_3> v = ReadVector_3File( std::cin );

    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        v[i] += trans;
    }

    WriteVector_3File( v );

    return 0;
}

