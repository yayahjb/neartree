// rotateq.cpp
//
// rotateq reads from standard input and writes to standard output. It is intended
// only to read 3-space vectors (as csv) and rotates the input by the
// quaternion specified on the command line.
//
// This is intended to be used with other programs for generating test data
// sets for NearTree.
//

#include "vector_3d.h"
#include <iostream>

#include "Data2CSV.h"
#include "cqrlib.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char* argv[])
//---------------------------------------------------------------------
{
	// make sure there are at least the 4 quaternion components on the input line
    if ( argc<5 )
        return( 0 );

    const SQR<double> q = SQR<double>( atof( argv[1] ), atof( argv[2] ), atof( argv[3] ), atof( argv[4] ) );
    Matrix_3x3 m = q.Quaternion2Matrix( );

    const std::vector<Vector_3> v = ReadVector_3File( std::cin );

    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        const Vector_3 vr = m.MV( v[i] );
        fprintf( stdout, "%f,%f,%f\n", vr[0], vr[1], vr[2] );
    }

    return 0;
}

