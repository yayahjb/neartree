// Data2CSV.cpp
//
// utility code for the various small programs intended to be
// used with other programs for generating test data
// sets for NearTree. Also included is a limited n-dimensional
// vector class//

#include "Data2CSV.h"

#include <sstream>

#include "rhrand.h"

static RHrand rhr;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
std::vector<std::string> GetNextLine( std::string& line )
//---------------------------------------------------------------------
{
    std::vector<std::string> result;

    std::stringstream lineStream(line);
    std::string       cell;

    while(std::getline(lineStream,cell,','))
    {
        result.push_back(cell);
    }
    return( result );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
std::vector<Vector_3> ReadVector_3File( std::istream& str )
//---------------------------------------------------------------------
{
    std::vector<Vector_3> v;

    while ( ! str.fail() && ! str.eof( ) )
    {
        std::string line;
        getline( str, line, '\n' );
        if ( ! str.fail( ) && ! str.eof( ) )
        {
            std::vector<std::string> vs = GetNextLine( line );
            const Vector_3 v1( atof(vs[0].c_str()), atof(vs[1].c_str()), atof(vs[2].c_str()) );
            v.push_back( v1 );
        }
    }
    return( v );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
std::vector<vecN> ReadGeneralFile( std::istream& str )
//---------------------------------------------------------------------
{
    std::vector<vecN > v;


    while ( ! str.fail() && ! str.eof( ) )
    {
        std::string line;
        getline( str, line, '\n' );
        if ( ! str.fail( ) && ! str.eof( ) )
        {
            if ( ! line.empty( ) )
            {
                const std::vector<std::string> vs = GetNextLine( line );
                std::vector<double> vd(vs.size( ));

                for ( unsigned int i=0; i<vs.size( ); ++i )
                {
                    vd[i] = atof(vs[i].c_str( ));
                }

                const vecN v1( vd );

                v.push_back( v1 );
            }
        }
    }
    return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// WriteVector_3File
//
// prints xyz of Vector_3 objects in CSV format (to stdout)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void WriteVector_3File( const std::vector<Vector_3>& v )
//---------------------------------------------------------------------
{
    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        const Vector_3& p = v[i];
        fprintf( stdout, "%f,%f,%f\n", p[0],p[1],p[2] );
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void WriteVectorFile( const std::vector<vecN>& vin )
//---------------------------------------------------------------------
{
    unsigned int dim = (unsigned int)vin[0].size( );
    for ( unsigned int i=0; i<vin.size( ); ++i )
    {
        const vecN& v = vin[i];
        for ( unsigned int k=0; k<dim; ++k )
        {
            if ( k == dim-1 )
            {
                fprintf( stdout, "%f", v[k] );
            }
            else
            {
                fprintf( stdout, "%f,", v[k] );
            }
        }
        fprintf( stdout, "\n" );
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FindBox
//
// Determines the boundary extents of the vecN object along each of
// the n axes. The maximums and the minimums are returned.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::pair<vecN,vecN> FindBox( const std::vector<vecN>& v )
//---------------------------------------------------------------------
{
    const long dim = (long)v[0].size( );
    vecN mins( dim, DBL_MAX);
    vecN maxs( dim, -DBL_MAX);

    for ( unsigned int i=0; i<v.size(); ++i )
    {
        const vecN& vn = v[i];
        for ( unsigned int k=0; k<vn.size( ); ++k )
        {
            const double d = vn[k];
            if ( d > maxs[k] ) maxs[k] = d;
            if ( d < mins[k] ) mins[k] = d;
        }

    }
    return( std::make_pair( mins, maxs ) );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BoxIt
//
// Converts a set of vecN points so that they fit into the range
// -1 to +1 along the largest axial extent of the set of points
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void BoxIt( std::vector<vecN>& v )
//---------------------------------------------------------------------
{
    const std::pair<vecN,vecN> box = FindBox( v );

    const vecN limits = box.second - box.first;
    const vecN middle = 0.5 * ( box.second + box.first );

    // determine the largest extent along any axis and
    // use that for scaling
    double scale = -DBL_MAX;
    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        for ( unsigned int k=0; k<limits.size(); ++k )
        {
            if ( scale < limits[k] ) scale = limits[k];
        }
    }

    if ( scale == -DBL_MAX  || scale == 0.0)
    {
        scale = 1.0;
    }

    // factor of two because -1 to +1 is length 2
    scale = 2.0 / scale;

    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        vecN& w = v[i];
        w = scale * ( w - middle );
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
vecN operator* ( const double& a, const vecN& v )
//-------------------------------------------------------------------------------------
{
    const long siz = (long)v.pd.size( );
    vecN vout( siz );

    for ( int i=0; i<siz; ++i )
    {
        vout.pd[i] = a * v.pd[i];
    }
    return ( vout );
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
vecN::vecN( const long n ) :
   dim(n), length(0)
//---------------------------------------------------------------------
{
    pd.resize( n );
    for( int i=0; i<dim; ++i )
    {
        pd[i] = rhr.urand( );
    }
    length = Norm( );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
vecN::vecN( const int n ):
   dim(n), length(0)
//---------------------------------------------------------------------
{
    pd.resize( n );
    for( int i=0; i<dim; ++i )
    {
        pd[i] = rhr.urand( );
    }
    length = Norm( );
}
