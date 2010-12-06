// DataGen.cpp : Defines the entry point for the console application.
//
// Generates data sets in csv for nearest neighbor testing. All output
// data sets fit within -1 to +1
// Several data types can be created:
/*
    GEODB",3) );  // geodesic sphere for dumbbell, # pts, sphere rad, separation of centers
    RANDB",3) );  // random sphere for dumbbell, # pts, sphere rad, separation of centers
    RANSP",1) );  // random pts within a unit sphere, n points
    RANUSP",1) ); // random pts on a unit sphere, n points
    GEOSP",1) );  // geodesic sphere, n - min # pts
    QRGEO",1) );  // 4-D quaternion geodesic, min # of points
    LAT2",1) );   // 2-D lattice - n pts in each direction
    LAT3",1) );   // 3-D lattice - n pts in each direction
    LAT4",1) );   // 4-D lattice - n pts in each direction
    HAM",2) );    // dimension of data point, n data points
    CLOUD",3) );  // cloud points per level, n levels, relative cloud spacing
*/
//

// STL libraries
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>

// system includes
#include <cfloat>
#include <cmath>


#include "geodesic.h"
#include "hammersley.h"
#include "Data2CSV.h"
#include "Vector_3D.h"

std::map<std::string,int> theMap;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SVG_Output
//
// produces an svg file of circles, one for each Vector_3 input
// The radius of the PICTURE is circleRadiusInPix
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SVG_Output(
    const long circleRadiusInPix,
    const std::vector<Vector_3>& points )
//---------------------------------------------------------------------
{
    const long sqSide = (long)(2.1 * circleRadiusInPix);
    fprintf( stdout, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
    fprintf( stdout, "<svg width=\"%ld\" height=\"%ld\" version=\"1.1\"\n", sqSide, sqSide);
    fprintf( stdout, "xmlns=\"http://www.w3.org/2000/svg\">\n");
    fprintf( stdout, "<rect width=\"%ld\" height=\"%ld\" style=\"fill:rgb(255,255,255); stroke-width:10;  stroke:rgb(255,255,255)\"/>\n",
        sqSide, sqSide);


    double maxx = -DBL_MAX;
    double minx = DBL_MAX;
    double maxy = -DBL_MAX;
    double miny = DBL_MAX;
    for( unsigned int i=0; i<points.size( ); ++i )
    {
       maxx = MAX(maxx,points[i][0]);
       maxy = MAX(maxy,points[i][1]);
       minx = MIN(minx,points[i][0]);
       miny = MIN(miny,points[i][1]);
    }

    const double scale = MAX( (maxx-minx), (maxy-miny) ) / circleRadiusInPix;

    for( unsigned int i=0; i<points.size( ); ++i )
    {
        const int r = 2;
        const double xd = (points[i][0]-minx)/scale;
        const double yd = (points[i][1]-miny)/scale;
        const int x = int( sign_t(xd)*xd );
        const int y = int( sign_t(yd)*yd );
        fprintf( stdout, "<circle r=\"%d\" cx=\"%d\" cy=\"%d\" stroke=\"black\" stroke-width=\"2\" fill=\"red\"/>\n",
            r, x, y );
    }

    fprintf( stdout, "</svg>\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  FindBox
//
// Determines the enclosing box around a set of Vector_3 points
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::pair<Vector_3,Vector_3> FindBox( const std::vector<Vector_3>& v )
//---------------------------------------------------------------------
{
    double xmax = -DBL_MAX;
    double ymax = -DBL_MAX;
    double zmax = -DBL_MAX;
    double xmin = DBL_MAX;
    double ymin = DBL_MAX;
    double zmin = DBL_MAX;

    for ( unsigned int i=0; i<v.size(); ++i )
    {
        const Vector_3& w = v[i];
        if ( w[0] > xmax ) xmax = w[0];
        if ( w[1] > ymax ) ymax = w[1];
        if ( w[2] > zmax ) zmax = w[2];
        if ( w[0] < xmin ) xmin = w[0];
        if ( w[1] < ymin ) ymin = w[1];
        if ( w[2] < zmin ) zmin = w[2];
    }

    return( std::make_pair( Vector_3( xmin, ymin, zmin ), Vector_3( xmax, ymax, zmax ) ) );
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BoxIt
//
// Converts a set of Vector_3 points so that they fit into the range
// -1 to +1 along the largest axial extent of the set of points
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void BoxIt( std::vector<Vector_3>& v )
//---------------------------------------------------------------------
{
    const std::pair<Vector_3,Vector_3> box = FindBox( v );

    const double dx = box.second[0] - box.first[0];
    const double dy = box.second[1] - box.first[1];
    const double dz = box.second[2] - box.first[2];

    double scale = 1.0;
    if ( scale < dx ) scale = dx;
    if ( scale < dy ) scale = dy;
    if ( scale < dz ) scale = dz;
    scale = 2.0 / scale;

    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        Vector_3& w = v[i];
        w = scale * Vector_3( w[0]-dx/2.0, w[1]-dy/2.0, w[2]-dz/2.0 );
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// HAM_Gen
//
// Uses the Hammersley quasirandom distribution to generate a set of
// n-dimensional (dimension dim) points.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<double> HAM_Gen( const long dim, const long nPoints )
//---------------------------------------------------------------------
{
  std::vector<int> seed( dim, 0 );
  std::vector<int> leap( dim, 1 );
  std::vector<int> base;
  base.push_back( -1000 );
  for ( int i=1; i<dim; ++i )
  {
     base.push_back( prime(i) );
  }
  std::vector<double> r(nPoints*dim);
  hammersley_dim_num_set( dim );
  hammersley_sequence ( nPoints, &r[0] );
  return( r );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GEN_GeoSphere
//
// Generates a geodesic sphere with at least nSpherePoints (but most
// likely more)
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<Vector_3> GEO_GenSphere( const long nSpherePoints )
//---------------------------------------------------------------------
{
    CNearTree<Vector_3> nt;
    const long count = GEO_AtLeastNGeodesic( nt, nSpherePoints );
    count;
    return( nt.GetObjectStore() );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RAN_GenSphere
//
// Generates nSpherePoints of Vector_3 points randomly distributed
// within the unit sphere
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<Vector_3> RAN_GenSphere( const long nSpherePoints )
//---------------------------------------------------------------------
{
    std::vector<Vector_3> v(nSpherePoints);
    int iseed = 19191;

    for ( int i=0; i<nSpherePoints; ++i )
    {
        v[i] = randVector( iseed );
    }

    return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RAN_GenSphere
//
// Generates nSpherePoints of Vector_3 points randomly distributed
// on the surface of the unit sphere
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<Vector_3> RAN_UNIT_GenSphere( const long nSpherePoints )
//---------------------------------------------------------------------
{
    std::vector<Vector_3> v(nSpherePoints);
    v = RAN_GenSphere( nSpherePoints );

    for ( int i=0; i<nSpherePoints; ++i )
    {
        v[i] = v[i].UnitV( );
    }

    return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// SVG_Output
//
// produces an svg file of circles, one for each n-dimensional input
// by plotting only the x and y (first two) values
// The radius of the PICTURE is circleRadiusInPix
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SVG_Output(
    const long circleRadiusInPix,
    const std::vector<std::vector<double> >& points )
//---------------------------------------------------------------------
{
    const long sqSide = (long)(2.1 * circleRadiusInPix);
    fprintf( stdout, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
    fprintf( stdout, "<svg width=\"%ld\" height=\"%ld\" version=\"1.1\"\n", sqSide, sqSide);
    fprintf( stdout, "xmlns=\"http://www.w3.org/2000/svg\">\n");
    fprintf( stdout, "<rect width=\"%ld\" height=\"%ld\" style=\"fill:rgb(255,255,255); stroke-width:10;  stroke:rgb(255,255,255)\"/>\n",
        sqSide, sqSide);

    for( unsigned int i=0; i<points.size( ); ++i )
    {
        const std::vector<double>& v = points[i];
        const int r = 2;
        const double xd = circleRadiusInPix*v[0];
        const double yd = circleRadiusInPix*v[1];
        const int x = int( sign_t(xd)*xd );
        const int y = int( sign_t(yd)*yd );
        fprintf( stdout, "<circle r=\"%d\" cx=\"%d\" cy=\"%d\" stroke=\"black\" stroke-width=\"2\" fill=\"red\"/>\n",
            r, x, y );
    }

    fprintf( stdout, "</svg>\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GEO_Dumbell
//
// Generates a dumbbell like object. A geodesic sphere is constructed.
// Then the right and left (along the x-axis) halves are treated
// in two parts. First, each is translated along x by +/- separation.
// Then the same sets are used to generate a conical set of points
// extending from the two hemispheres to the origin (where the cone
// vertices are).
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<Vector_3> GEO_Dumbell(
                                  const long nSpherePoints,
                                  const double sphereRadius,
                                  const double separation )
//---------------------------------------------------------------------
{
    const double normalsep = separation/sphereRadius;
    std::vector<Vector_3> vOneSphere = GEO_GenSphere( nSpherePoints/2 );

    for ( unsigned int i=0; i<vOneSphere.size( ); ++i )
    {
        vOneSphere[i] = sphereRadius*vOneSphere[i];
    }

    std::vector<Vector_3> dbell;

    for ( unsigned int i=0; i<vOneSphere.size(); ++i )
    {
        const double newx = -0.5*normalsep*(1.0-vOneSphere[i][0]);
        const double newy = vOneSphere[i][1]*(1.0-abs(vOneSphere[i][0]));
        const Vector_3 vTemp = Vector_3( newx, newy, vOneSphere[i][2] );
        if ( vOneSphere[i][0] < 0.0 )
        {
            dbell.push_back( vOneSphere[i] - Vector_3(0.5*normalsep, 0.0, 0.0 ) );
            dbell.push_back( vTemp + Vector_3(1.0*normalsep, 0.0, 0.0 ) );
        }
        else
        {
            dbell.push_back( vOneSphere[i] + Vector_3(0.5*normalsep, 0.0, 0.0 ) );
            dbell.push_back( vTemp + Vector_3(0.0*normalsep, 0.0, 0.0 ) );
        }
    }

    for ( unsigned int i=0; i<dbell.size( ); ++i )
    {
        dbell[i] = sphereRadius * dbell[i];
    }

    BoxIt( dbell );
    return( dbell );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RAN_Dumbell
//
// Generates a dumbbell like object. A random sphere is constructed.
// Then the right and left (along the x-axis) halves are treated
// in two parts. First, each is translated along x by +/- separation.
// Then the same sets are used to generate a conical set of points
// extending from the two hemispheres to the origin (where the cone
// vertices are).
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<Vector_3> RAN_Dumbell(
                                  const long nSpherePoints,
                                  const double sphereRadius,
                                  const double separation )
//---------------------------------------------------------------------
{
    const double normalsep = separation/sphereRadius;
    std::vector<Vector_3> vOneSphere = RAN_GenSphere( nSpherePoints/2 );
    for ( unsigned int i=0; i<vOneSphere.size( ); ++i )
    {
        vOneSphere[i] = sphereRadius*vOneSphere[i];
    }

    std::vector<Vector_3> dbell;

    for ( unsigned int i=0; i<vOneSphere.size(); ++i )
    {
        const double newx = -0.5*normalsep*(1.0-vOneSphere[i][0]);
        const double newy = vOneSphere[i][1]*(1.0-abs(vOneSphere[i][0]));
        const Vector_3 vTemp = Vector_3( newx, newy, vOneSphere[i][2] );
        if ( vOneSphere[i][0] < 0.0 )
        {
            dbell.push_back( vOneSphere[i] - Vector_3(0.5*normalsep, 0.0, 0.0 ) );
            dbell.push_back( vTemp + Vector_3(1.0*normalsep, 0.0, 0.0 ) );
        }
        else
        {
            dbell.push_back( vOneSphere[i] + Vector_3(0.5*normalsep, 0.0, 0.0 ) );
            dbell.push_back( vTemp + Vector_3(0.0*normalsep, 0.0, 0.0 ) );
        }
    }

    for ( unsigned int i=0; i<dbell.size( ); ++i )
    {
        dbell[i] = sphereRadius * dbell[i];
    }

    BoxIt( dbell );
    return( dbell );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LAT4_Gen
//
// Creates a 4-D lattice of points
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<vecN> LAT4_Gen( const int elemCount )
//---------------------------------------------------------------------
{
    std::vector<vecN> v;

    for ( int m=0; m<elemCount; ++m )
    {
        std::vector<double> vTemp(4);
        vTemp[0] = double(m);
        for ( int n=0; n<elemCount; ++n )
        {
            vTemp[1] = double(n);
            for ( int o=0; o<elemCount; ++o )
            {
                vTemp[2] = double(o);
                for ( int p=0; p<elemCount; ++p )
                {
                    vTemp[3] = double(p);
                    v.push_back( vecN( vTemp ) );
                }

            }
        }
    }

    BoxIt( v );
    return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LAT2_Gen
//
// Creates a 2-D lattice of points
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<vecN> LAT2_Gen( const int elemCount )
//---------------------------------------------------------------------
{
    std::vector<vecN> v;

    for ( int i=0; i<elemCount; ++i )
    {
        std::vector<double> vTemp(2);
        vTemp[0] = double(i);
        for ( int k=0; k<elemCount; ++k )
        {
            vTemp[1] = double(k);
            v.push_back( vecN( vTemp ) );
        }
    }

    BoxIt( v );
    return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LAT3_Gen
//
// Creates a 3-D lattice of points
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<Vector_3> LAT3_Gen( const int elemCount )
//---------------------------------------------------------------------
{
    std::vector<Vector_3> v;

    for ( int i=0; i<elemCount; ++i )
    {
        for ( int j=0; j<elemCount; ++j )
        {
            for ( int k=0; k<elemCount; ++k )
            {
                v.push_back( Vector_3( double(i), double(j), double(k) ) );
            }
        }
    }

    BoxIt( v );
    return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// LAT_GetSize
//
// Determines how large a lattice the user has requested.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int LAT_GetSize( const int argc, char* argv[] )
//---------------------------------------------------------------------
{
    int siz = 3;
    if ( argc <= 2 )
    {
        siz = 3;
    }
    else
    {
        int k = atoi( argv[2] );
        if ( k <= 0 )
        {
            siz = 3;
        }
        else
        {
            siz = k;
        }
    }
    return( siz );
}
int iseed = 19191;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CLOUD_Gen(
               std::vector<Vector_3>& v,
               const long levels,
               const long points,
               const Vector_3& newOrigin,
               const double baseSeparation,
               const double separation )
//---------------------------------------------------------------------
{
    if ( levels == 1 )
    {
        //const Vector_3 trans = randVector(iseed)/separation;
        for ( int i=0; i<points; ++i )
        {
            v.push_back( randVector(iseed)/separation + newOrigin );
        }
    }
    else
    {
        //std::vector<Vector_3> vp;
        for ( int i=0; i<points; ++i )
        {
            const Vector_3 vTemp = randVector(iseed)/separation + newOrigin;
            //v.push_back( vTemp );
            CLOUD_Gen( v, levels-1, points, vTemp, baseSeparation, separation*baseSeparation );
        }
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool ParseArgs( const int argc, char* argv[], std::string& sRet, std::vector<double>& vRet  )
//---------------------------------------------------------------------
{
    if ( argc <=1 )
        return( false );


    std::string s(argv[1]);
    if ( theMap[s]>argc-2 )
        return( false );

    std::transform(s.begin(), s.end(), s.begin(), toupper);

    sRet = s;

    for ( int i=2; i<=argc-1; ++i )
    {
        vRet.push_back( atof( argv[i] ) );
    }

    return( true );
}

std::vector<vecN> BOX_Gen( const int dim, const int nPoints )
{
    vecN v(dim);
    std::vector<vecN> vout;
    RHrand rhr;

    for ( int ipt=0; ipt<nPoints; ++ipt )
    {
        for ( int i=0; i<dim; ++i )
        {
            // put it into -1 to +1
            v[i] = 2.0*rhr.urand( ) - 1.0;
        }
        vout.push_back( v );
    }

    return( vout );
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char* argv[])
//---------------------------------------------------------------------
{
    theMap.insert( std::make_pair( "GEODB",3) );  // geodesic sphere for dumbbell, # pts, sphere rad, separation of centers
    theMap.insert( std::make_pair( "RANDB",3) );  // random sphere for dumbbell, # pts, sphere rad, separation of centers
    theMap.insert( std::make_pair( "RANSP",1) );  // random pts within a unit sphere, n points
    theMap.insert( std::make_pair( "RANUSP",1) ); // random pts on a unit sphere, n points
    theMap.insert( std::make_pair( "RANBOX",2) );    // dimension of data point, n data points
    theMap.insert( std::make_pair( "GEOSP",1) );  // geodesic sphere, n - min # pts
    theMap.insert( std::make_pair( "QRGEO",1) );  // 4-D quaternion geodesic, min # of points
    theMap.insert( std::make_pair( "LAT2",1) );   // 2-D lattice - n pts in each direction
    theMap.insert( std::make_pair( "LAT3",1) );   // 3-D lattice - n pts in each direction
    theMap.insert( std::make_pair( "LAT4",1) );   // 4-D lattice - n pts in each direction
    theMap.insert( std::make_pair( "HAM",2) );    // dimension of data point, n data points
    theMap.insert( std::make_pair( "FRAC",1) );    // fractions in 0-1 scaled into -1 to +1, from 1/1 to 1/n
    theMap.insert( std::make_pair( "COMB",1) );    // n integers from 0-1 scaled into -1 to +1
    theMap.insert( std::make_pair( "CLOUD",3) );  // m=cloud points per level, n levels, relative cloud spacing
                                                  // number of points = m^n

    std::string cmd;
    std::vector<double> vArgs;
    ParseArgs( argc, argv, cmd, vArgs );

    if ( theMap.find( cmd ) == theMap.end() )
    {
        std::map<std::string,int>::iterator it;
        for ( it=theMap.begin(); it!=theMap.end( ); ++it )
        {
            fprintf( stdout, "    %s   %d parameters\n", (*it).first.c_str( ), (*it).second );
        }

        fprintf( stdout, "\nAll output should be centered in a box of -1 to +1\n" );
        exit( 1 );
    }

    bool bPrint = true;

    if ( cmd == "GEODB" )
    {
        const std::vector<Vector_3> vGEO = GEO_Dumbell( int(vArgs[0]), vArgs[1], vArgs[2] );
        if ( bPrint )
        {
            WriteVector_3File( vGEO );
        }
        else
        {
            SVG_Output( 500, vGEO );
        }
    }
    else if ( cmd == "RANDB" )
    {
        std::vector<Vector_3> vRAN = RAN_Dumbell( int(vArgs[0]), vArgs[1], vArgs[2] );
        const Vector_3 vtrans( 0.0, 5.0, 0.0 );
        for ( unsigned int i=0; i<vRAN.size( ); ++i )
        {
            vRAN[i] = vRAN[i] + vtrans;
        }

        if ( bPrint )
        {
            WriteVector_3File( vRAN );
        }
        else
        if ( ! bPrint )
        {
            SVG_Output( 500, vRAN );
        }
    }
    else if ( cmd == "RANSP" )
    {
        const std::vector<Vector_3> vRAN = RAN_GenSphere( int(vArgs[0]) );
        if ( bPrint )
        {
            WriteVector_3File( vRAN );
        }
        else
        {
            SVG_Output( 500, vRAN );
        }
    }
    else if ( cmd == "CLOUD" )
    {
        std::string cmd;
        std::vector<double> ret;
        ParseArgs( argc, argv, cmd, ret );
        const long points = long(ret[0]);
        const long levels = long(ret[1]);
        const double separation = ret[2];

        std::vector<Vector_3> v;

        CLOUD_Gen( v, levels, points, Vector_3::GetZeroVector( ), separation, 1.0 );

        BoxIt( v );

        if ( bPrint )
        {
            WriteVector_3File( v );
        }
        else
        {
            SVG_Output( 500, v );
        }
    }
    else if ( cmd == "RANUSP" )
    {
        const std::vector<Vector_3> vRAN = RAN_UNIT_GenSphere( int(vArgs[0]) );
        if ( bPrint )
        {
            WriteVector_3File( vRAN );
        }
        else
        {
            SVG_Output( 500, vRAN );
        }
    }
    else if ( cmd == "GEOSP" )
    {
        const std::vector<Vector_3> vRAN = GEO_GenSphere( int(vArgs[0]) );
        if ( bPrint )
        {
            WriteVector_3File( vRAN );
        }
        else
        {
            SVG_Output( 500, vRAN );
        }
    }
    else if ( cmd == "LAT2" )
    {
        const std::vector<vecN> vLAT = LAT2_Gen( int(vArgs[0]) );
        if ( bPrint )
        {
            WriteVectorFile( vLAT );
        }
        else
        {
            //SVG_Output( 500, vLAT );
        }
    }
    else if ( cmd == "LAT4" )
    {
        const std::vector<vecN> vLAT = LAT4_Gen( int(vArgs[0]) );
        if ( bPrint )
        {
            WriteVectorFile( vLAT );
        }
        else
        {
            //SVG_Output( 500, vLAT );
        }
    }
    else if ( cmd == "LAT3" )
    {
        const std::vector<Vector_3> vLAT = LAT3_Gen( int(vArgs[0]) );
        if ( bPrint )
        {
            WriteVector_3File( vLAT );
        }
        else
        {
            SVG_Output( 500, vLAT );
        }
    }
    else if ( cmd == "HAM" )
    {
        const std::vector<double> vHAM = HAM_Gen( int(vArgs[0]), int(vArgs[1]) );

        for ( unsigned int i=0; i<vHAM.size( ); i+=int(vArgs[0]) )
        {
            for ( int k=0; k<int(vArgs[0]); ++k )
            {
                if ( k == int(vArgs[0])-1 )
                {
                    fprintf( stdout, "%f", vHAM[i+k] );
                }
                else
                {
                    fprintf( stdout, "%f,", vHAM[i+k] );
                }
            }
            fprintf( stdout, "\n" );
        }
    }
    else if ( cmd == "QRGEO" )
    {
        // Geodesic for quaternions

        //long GEO_AtLeastNQGeodesic( CNearTree<SQR<double> >& t, const long minPoints );
        CNearTree<SQR<double> > nt;
        const long nGEOPoints = GEO_AtLeastNQGeodesic( nt, int(vArgs[0]) );
        nGEOPoints;

        std::vector<Vector_3> vq;
        for ( unsigned int i=0; i<nt.size( ); ++i )
        {
            vq.push_back( Vector_3( nt[i][0], nt[i][1], nt[i][2] ) );
        }

        for ( unsigned int i=0; i<nt.size( ); ++i )
        {
            if ( bPrint )
            {
                fprintf( stdout, "%f,%f,%f,%f\n", nt[i][0], nt[i][1], nt[i][2], nt[i][3] );
            }
        }

        if ( bPrint )
        {
            // WriteVector_3File( vq ); // quarternions are printed above
        }
        else
        {
            SVG_Output( 500, vq );
        }
    }
    else if ( cmd == "RANBOX" )
    {  // random set of points in an n-dimensional box
        const std::vector<vecN> vBOX = BOX_Gen( int(vArgs[0]), int(vArgs[1]) );
        WriteVectorFile( vBOX );
    }
    else if ( cmd == "FRAC" )
    {
        for ( unsigned int i=0; i<(unsigned int)int(vArgs[0]); ++i )
        {
            fprintf( stdout, "%f\n", 2.0 / double(i+1) - 1.0 );
        }
    }
    else if ( cmd == "COMB" )
    {
        for ( unsigned int i=0; i<(unsigned int)int(vArgs[0]); ++i )
        {
            fprintf( stdout, "%g\n", double(2*i-1) );
        }
    }
    else
    {
        fprintf( stdout, "The input was not recognized\n" );
    }

    //if ( bPrint )
    //{
    //    const std::vector<Vector_3> vInput = ReadVector_3File( std::cin );
    //    WriteVector_3File( vInput );
    //}

    return 0;
}

