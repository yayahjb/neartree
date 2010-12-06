/*
 *  geodesics.cpp
 *  NearTree
 *
 *  Copyright 2010 Larry Andrews.  All rights reserved
 *
 *
 */

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE CNearTree API UNDER THE TERMS OF THE LGPL *
 *                                                                    *
 **********************************************************************/

/************************* LGPL NOTICES *******************************
 *                                                                    *
 * This library is free software; you can redistribute it and/or      *
 * modify it under the terms of the GNU Lesser General Public         *
 * License as published by the Free Software Foundation; either       *
 * version 2.1 of the License, or (at your option) any later version. *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
 * Lesser General Public License for more details.                    *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License along with this library; if not, write to the Free         *
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
 * MA  02110-1301  USA                                                *
 *                                                                    *
 **********************************************************************/

#include "TNear.h"
#include <cmath>
#include "geodesic.h"
#include <vector>
#include "vector_3d.h"

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AtLeastNGeodesic()
 Description: Computes the points of a particular geodesic sphere (ala Buckminster
              Fuller). Uses GEO_GenerateGeodesic2 to create a relatively uniform set of
              at least minPoints points.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_AtLeastNGeodesic( CNearTree<Vector_3>& t, const long minPoints )
/*-------------------------------------------------------------------------------------*/
{
    double d;
    d = minPoints > 6?((sqrt(((double)minPoints) -2.)-2.)/2.):0.;
    return( GEO_GenerateGeodesic2( t, long(ceil(d) )) );
}

double cbrt( const double d )
{
    return( pow(d,1.0/3.0) );
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AtLeastNQGeodesic()
 Description: Computes the points of a particular geodesic sphere (ala Buckminster
     Fuller) at the equator (w=0) of the unit quaternions and then in steps of the
     spacing of the geodesic computes proportionately small geodesic spheres at
     increasing and decreasing latitudes until it puts a single point at each pole. 
     Uses GenerateGeodesic2 to create a relatively uniform set of points.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_AtLeastNQGeodesic( CNearTree<SQR<double> >& t, const long minPoints )
/*-------------------------------------------------------------------------------------*/
{
    long level;
    level = (minPoints > 44?(((GEOM_SQRlevelforpointcount(minPoints)))):0);
    return( GEO_QGenerateGeodesic2( t,level ) );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AtLeastNHQGeodesic()
 Description: Computes the points of a particular geodesic sphere (ala Buckminster
 Fuller) at the equator (w=0) of the unit quaternions and then in steps of the
 spacing of the geodesic computes proportionately small geodesic spheres at
 increasing latitudes until it puts a single point at the North pole. 
 Uses GenerateGeodesic2 to create a relatively uniform set of points.
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_AtLeastNHQGeodesic( CNearTree<HQR<double> >& t, const long minPoints )
/*-------------------------------------------------------------------------------------*/
{
    long level;
    level = (minPoints > 24?(((GEOM_HQRlevelforpointcount(minPoints)))):0);
    return( GEO_HQGenerateGeodesic2( t, level ) );
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_LevelforPoints()
 Description: Computes a level that will produce the requires number of
 points in a geodesic sphere
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_LevelforPoints( const long minPoints )
/*-------------------------------------------------------------------------------------*/
{
    double d;
    d = minPoints > 6?((sqrt(((double)minPoints) -2.)-2.)/2.):0.;
    return  ((long)ceil(d) );
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AverageDistanceBetweenPoints()
 Description: Computes an approximate distance between point. Since the points are
 not precisely evenly distributed, this is only an average for an 
 environment.  This is an empircal approximation.
 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double GEO_AverageDistanceBetweenPoints( const long level )
/*-------------------------------------------------------------------------------------*/
{
    return ( GEOM_avgdist(level) );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AverageAreaPerPoint()
 Description: Computes an approximate area (steradians) per point. Since the points are
 not precisely evenly distributed, this is only an average.  This is an emperical
 approximation.
 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double GEO_AverageAreaPerPoint( const long level )
/*-------------------------------------------------------------------------------------*/
{
    return (GEOM_avgarea(level));
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_GeodesicPointCount()
 Description: Computes the number of points for a particular geodesic sphere 
              (ala Buckminster Fuller). The returned number of points is the actual 
              number that will be computed. It is the number of points of the smallest
              geodesic that has at least as many points as requested (and very likely
              more).

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_GeodesicPointCount( const long level )
/*-------------------------------------------------------------------------------------*/
{
    const double pointCount = GEOM_pointcount(level) + 0.5;
    return ( long(pointCount) );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_QGeodesicPointCount()
 Description: Computes the number of points for a particular geodesic sphere 
 (ala Buckminster Fuller) at the equator of a unit quaternion sphere, plus
 progressively smaller quaternion spheres at lines of latitude separated
 by the same spacing, plus a point at each pole.  The returned number of 
 points is the actual number that will be computed. 
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_QGeodesicPointCount( const long level )
/*-------------------------------------------------------------------------------------*/
{
    long epointcount;
    long pointcount;
    double pi;
    double omega;
    double target;
    long N;
    long i;
    double spacing;
    pi = 4.*atan2(1.,1.);
    epointcount = (long)(GEOM_pointcount(level) + 0.5);
    N = (long)(.5 + pi/(1.5* GEO_AverageDistanceBetweenPoints(level)));
    omega = pi/((double)(2*N));
    pointcount = 2+epointcount;
    spacing = GEOM_avgdist(level);
    for (i=1; i < N; i++) {
        long glevel;
        target = spacing/sin(((double)i)*omega);
        glevel = (long)(1.5+GEOM_levelformaxdist(target));
        pointcount += 2*GEO_GeodesicPointCount(glevel);
    }
    return ( pointcount );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_HQGeodesicPointCount()
 Description: Computes the number of points for a particular geodesic sphere 
 (ala Buckminster Fuller) at the equator of a unit quaternion sphere, plus
 progressively smaller quaternion spheres at lines of latitude separated
 by the same spacing in the northern hemisphere, plus a point at the
 north pole (w = 1).  The returned number of points is the actual number that 
 will be computed. 
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_HQGeodesicPointCount( const long level )
/*-------------------------------------------------------------------------------------*/
{
    long epointcount;
    long pointcount;
    double pi;
    double omega;
    double target;
    long N;
    long i;
    long spacing;
    pi = 4.*atan2(1.,1.);
    epointcount = (long)(GEOM_pointcount(level) + 0.5);
    N = (long)(.5 + pi/(1.5* GEO_AverageDistanceBetweenPoints(level)));
    omega = pi/((double)(2*N));
    pointcount = 1+epointcount/2;
    spacing = (long)(GEOM_avgdist(level)+0.5);
    for (i=1; i < N; i++) {
        long glevel;
        target = spacing/sin(((double)i)*omega);
        glevel = (long)(1.5+GEOM_levelformaxdist(target));
        pointcount += GEO_GeodesicPointCount(glevel);
    }
    return ( pointcount );
}



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_GenerateGeodesic2()

 EQUATION: N_points = 4*level^2 + 8*level + 6

 Description: Computes the points of a particular geodesic sphere (ala Buckminster
              Fuller). Three points at the corners of an octant of a sphere are chosen
              to start the generation. The arcs from one of the points (001 here) to
              the other two are divided into the input number of segments. The points
              along the arc are connected by great circles, with varying numbers of
              points such that the distribution of the points is approximately even.
              This is the same algorithm used to generate geodesic domes, except that
              they usually start from an icosahedron face. All that is needed to 
              use a different polyhedron for the generator is to define the segment
              that is the seed and define the symmetry operations. NearTree is used
              to prevent duplicate points from being entered.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_GenerateGeodesic2( CNearTree<Vector_3>& t, const long level )
/*-------------------------------------------------------------------------------------*/
{

    t.clear( );

    // start by defining an octant of the sphere
    Vector_3 v0(Vector_3( 0.0, 0.0, 1.0 ) );
    Vector_3 v1(Vector_3( 1.0, 0.0, 0.0 ) );
    Vector_3 v2(Vector_3( 0.0, 1.0, 0.0 ) );

    // computer the step angle from z to the equator 
    const double angle1 = Vector_3::Angle( v0, Vector_3::GetZeroVector(), v1 ) / (level+1);
    const double angle2 = Vector_3::Angle( v0, Vector_3::GetZeroVector(), v2 ) / (level+1);

    // convert the angle steps to rotation matrices
    const Matrix_3x3 m1 = (v0.Cross(v1)).Rotmat( angle1 );
    const Matrix_3x3 m2 = (v0.Cross(v2)).Rotmat( angle2 );

    Vector_3 vStart1 = v0;
    Vector_3 vStart2 = v0;

    // insert the start point (z) and its negative
    t.insert( v0 );
    t.insert( -v0 );

    // march down each great circle (from z to x and from z to y)
    const double minRadius = 0.00001;
    for ( int i=0; i<=level; ++i )
    {
        // one step down the two lines
        vStart1 = m1.MV( vStart1 ); // points on the zx arc
        vStart2 = m2.MV( vStart2 ); // points on the zy arc

            Vector_3 vClose;

            Vector_3 vt;
            const double x = vStart1[0];
            const double y = vStart1[1];
            const double z = vStart1[2];

            // insert all of the points, one per octant (ignoring duplicates)
            vt = Vector_3( x, y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( x, y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( x, -y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( x, -y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, -y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, -y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );

        // find the great circle from vStart1 to vStart2
        const Vector_3 axis = vStart1.Cross( vStart2 );
        // find the angle step from vStart1 to vStart2 along the great circle
        const double angle = Vector_3::Angle( vStart1, Vector_3::GetZeroVector(), vStart2 ) / (i+1);
        const Matrix_3x3 m = axis.Rotmat( angle );

        Vector_3 vTemp = vStart1;

        // march along the great circle from the point on the xz line to the one
        // on the yz line
        for ( int j=0; j<=i; ++j )
        {
            vTemp = m.MV( vTemp );
            const double x = vTemp[0];
            const double y = vTemp[1];
            const double z = vTemp[2];

            // insert all of the points, one per octant (ignoring duplicates)
            vt = Vector_3( x, y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( x, y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( x, -y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( x, -y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, -y, z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
            vt = Vector_3( -x, -y, -z );
            if ( !  t.NearestNeighbor( minRadius, vClose, vt ) )
                t.insert( vt );
        }
        t.CompleteDelayedInsert( );

    }
    return( long(t.size()) );
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_QGenerateSubGeodesic2()
 
 EQUATIONS: 
 Latitudinal geodesic sphere
 N_points = 4*level^2 + 8*level + 6
 at angle omega from north pole
 and, if antipode is set,
 at angle -omega from south pole
 
 Description: Computes the points of a particular geodesic sphere (ala Buckminster
 Fuller) as the equatorial points of a covering of the unit quaternion sphere.
 Three points at the corners of an octant of a 3-sphere are chosen
 to start the generation. The arcs from one of the points (001 here) to
 the other two are divided into the input number of segments. The points
 along the arc are connected by great circles, with varying numbers of
 points such that the distribution of the points is approximately even.
 This is the same algorithm used to generate geodesic domes, except that
 they usually start from an icosahedron face. All that is needed to 
 use a different polyhedron for the generator is to define the segment
 that is the seed and define the symmetry operations. NearTree is used
 to prevent duplicate points from being entered.
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_QGenerateSubGeodesic2( CNearTree<SQR<double> >& t, 
                                const long level,
                                const double omega,
                                const int antipode)
/*-------------------------------------------------------------------------------------*/
{
    /* start by defining an octant of the sphere */
    
    Vector_3 v0(Vector_3( 0.0, 0.0, 1.0 ));
    Vector_3 v1(Vector_3( 1.0, 0.0, 0.0 ));
    Vector_3 v2(Vector_3( 0.0, 1.0, 0.0 ));
    
    /* compute the step angle from z to the equator */
    const double angle1 = Vector_3::Angle( v0, Vector_3::GetZeroVector(), v1 ) / (level+1);
    const double angle2 = Vector_3::Angle( v0, Vector_3::GetZeroVector(), v2 ) / (level+1);

    /* convert the angle steps to rotation matrices */
    const Matrix_3x3 m1 = (v0.Cross(v1)).Rotmat( angle1 );
    const Matrix_3x3 m2 = (v0.Cross(v2)).Rotmat( angle2 );

    int i, j;
    Vector_3 negv, tempv;
    Vector_3 axis, vStart1, vStart2;
    const double minRadius = 0.00001;
    SQR<double>  qClose;
    SQR<double> qn, qs; 
    double cosomega, sinomega;
    long initcount;
    
    cosomega = cos(omega);
    sinomega = sin(omega);
    initcount = (long)t.size();
    
    
    vStart1 = v0;
    vStart2 = v0;
    
    /* insert the start point (z) and its negative */
    qn.Set(cosomega,sinomega*(v0[0]),sinomega*(v0[1]),sinomega*(v0[2]));
    t.insert(qn);
    if (antipode) {
        qs.Set(-cosomega,-sinomega*(v0[0]),-sinomega*(v0[1]),-sinomega*(v0[2]));
        t.insert(qs);
    }

    qn.Set(cosomega,-sinomega*(v0[0]),-sinomega*(v0[1]),-sinomega*(v0[2]));
    t.insert(qn);
    if (antipode) {
        qn.Set(-cosomega,sinomega*(v0[0]),sinomega*(v0[1]),sinomega*(v0[2]));
        t.insert(qs);
    }
    
    /* march down each great circle (from z to x and from z to y) */
    for ( i=0; i<=level; ++i )
    {
        
        /* one step down the two lines */
        vStart1 = m1.MV( vStart1 ); // points on the zx arc
        vStart2 = m2.MV( vStart2 ); // points on the zy arc
        
        Vector_3 vt;
        const double x = vStart1[0];
        const double y = vStart1[1];
        const double z = vStart1[2];
        
        
        /* insert all of the points, one per octant (ignoring duplicates) */
        
        qn.Set(cosomega, sinomega*x,  sinomega*y,  sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );
        
        qn.Set(cosomega, sinomega*x,  sinomega*y, -sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );

        qn.Set(cosomega, sinomega*x, -sinomega*y,  sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );
        
        qn.Set(cosomega, sinomega*x, -sinomega*y, -sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );

	    qn.Set(cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
	    if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	            t.insert( qn );
	    
	    qn.Set(cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
	    if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	            t.insert( qn );

	    qn.Set(cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
	    if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	            t.insert( qn );
	    
	    qn.Set(cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
	    if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	            t.insert( qn );
	    
        if (antipode) {
        
            qs.Set(-cosomega, sinomega*x,  sinomega*y,  sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );
            
            qs.Set(-cosomega, sinomega*x,  sinomega*y, -sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );

            qs.Set(-cosomega, sinomega*x, -sinomega*y,  sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );
            
            qs.Set(-cosomega, sinomega*x, -sinomega*y, -sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );

            qs.Set(-cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );
            
            qs.Set(-cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );

            qs.Set(-cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );
            
            qs.Set(-cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
            if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
            t.insert( qs );
        }
        
        /* find the great circle from vStart1 to vStart2 */
        const Vector_3 axis = vStart1.Cross( vStart2 );
        /* find the angle step from vStart1 to vStart2 along the great circle */
        const double angle = Vector_3::Angle( vStart1, Vector_3::GetZeroVector(), vStart2 ) / (i+1);
        const Matrix_3x3 m = axis.Rotmat( angle );

        Vector_3 vTemp = vStart1;

        
        /* march along the great circle from the point on the xz line to the one
         on the yz line */
        for ( j=0; j<=i; ++j )
        {
            vTemp = m.MV( vTemp );
            const double x = vTemp[0];
            const double y = vTemp[1];
            const double z = vTemp[2];
            
	        qn.Set(cosomega, sinomega*x,  sinomega*y,  sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );
	        
	        qn.Set(cosomega, sinomega*x,  sinomega*y, -sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );

	        qn.Set(cosomega, sinomega*x, -sinomega*y,  sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );
	        
	        qn.Set(cosomega, sinomega*x, -sinomega*y, -sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );

	        qn.Set(cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );
	        
	        qn.Set(cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );

	        qn.Set(cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );
	        
	        qn.Set(cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );
	        
            if (antipode) {
            
                qs.Set(-cosomega, sinomega*x,  sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, sinomega*x,  sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );

                qs.Set(-cosomega, sinomega*x, -sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, sinomega*x, -sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );

                qs.Set(-cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );

                qs.Set(-cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
            }
         }
         t.CompleteDelayedInsert( );
        
    }
    return( long(t.size()) - initcount );
}

long GEO_HQGenerateSubGeodesic2( CNearTree<HQR<double> >& t, 
                                const long level,
                                const double omega,
                                const int antipode)
/*-------------------------------------------------------------------------------------*/
{
    /* start by defining an octant of the sphere */
    
    Vector_3 v0(Vector_3( 0.0, 0.0, 1.0 ));
    Vector_3 v1(Vector_3( 1.0, 0.0, 0.0 ));
    Vector_3 v2(Vector_3( 0.0, 1.0, 0.0 ));
    
    /* compute the step angle from z to the equator */
    const double angle1 = Vector_3::Angle( v0, Vector_3::GetZeroVector(), v1 ) / (level+1);
    const double angle2 = Vector_3::Angle( v0, Vector_3::GetZeroVector(), v2 ) / (level+1);

    /* convert the angle steps to rotation matrices */
    const Matrix_3x3 m1 = (v0.Cross(v1)).Rotmat( angle1 );
    const Matrix_3x3 m2 = (v0.Cross(v2)).Rotmat( angle2 );

    int i, j;
    Vector_3 negv, tempv;
    Vector_3 axis, vStart1, vStart2;
    const double minRadius = 0.00001;
    HQR<double>  qClose;
    HQR<double>  qn, qs; 
    double cosomega, sinomega;
    long initcount;
    int ishemi;
    
    cosomega = cos(omega);
    sinomega = sin(omega);
    initcount = (long)t.size();
    ishemi = true;
    
    
    if (fabs(cosomega) > minRadius) ishemi = false;

    
    vStart1 = v0;
    vStart2 = v0;
    
    /* insert the start point (z) and its negative */
    qn.Set(cosomega,sinomega*(v0[0]),sinomega*(v0[1]),sinomega*(v0[2]));
    t.insert(qn);
    if (antipode) {
        qs.Set(-cosomega,-sinomega*(v0[0]),-sinomega*(v0[1]),-sinomega*(v0[2]));
        t.insert(qs);
    }
    if (fabs(cosomega) > minRadius) {

        qn.Set(cosomega,-sinomega*(v0[0]),-sinomega*(v0[1]),-sinomega*(v0[2]));
        t.insert(qn);
        if (antipode) {
            qn.Set(-cosomega,sinomega*(v0[0]),sinomega*(v0[1]),sinomega*(v0[2]));
            t.insert(qs);
        }
    }
    
    /* march down each great circle (from z to x and from z to y) */
    for ( i=0; i<=level; ++i )
    {
        
        /* one step down the two lines */
        vStart1 = m1.MV( vStart1 ); // points on the zx arc
        vStart2 = m2.MV( vStart2 ); // points on the zy arc
        
        Vector_3 vt;
        const double x = vStart1[0];
        const double y = vStart1[1];
        const double z = vStart1[2];
        
        
        /* insert all of the points, one per octant (ignoring duplicates) */
        
        qn.Set(cosomega, sinomega*x,  sinomega*y,  sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );
        
        qn.Set(cosomega, sinomega*x,  sinomega*y, -sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );

        qn.Set(cosomega, sinomega*x, -sinomega*y,  sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );
        
        qn.Set(cosomega, sinomega*x, -sinomega*y, -sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );

        if (fabs(cosomega) > minRadius) {

        qn.Set(cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );
        
        qn.Set(cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );

        qn.Set(cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );
        
        qn.Set(cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
                t.insert( qn );
        
            if (antipode) {
            
                qs.Set(-cosomega, sinomega*x,  sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, sinomega*x,  sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );

                qs.Set(-cosomega, sinomega*x, -sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, sinomega*x, -sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );

                qs.Set(-cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );

                qs.Set(-cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
                
                qs.Set(-cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
                t.insert( qs );
            }
        }
        
        /* find the great circle from vStart1 to vStart2 */
        const Vector_3 axis = vStart1.Cross( vStart2 );
        /* find the angle step from vStart1 to vStart2 along the great circle */
        const double angle = Vector_3::Angle( vStart1, Vector_3::GetZeroVector(), vStart2 ) / (i+1);
        const Matrix_3x3 m = axis.Rotmat( angle );

        Vector_3 vTemp = vStart1;

        
        /* march along the great circle from the point on the xz line to the one
         on the yz line */
        for ( j=0; j<=i; ++j )
        {
            vTemp = m.MV( vTemp );
            const double x = vTemp[0];
            const double y = vTemp[1];
            const double z = vTemp[2];
            
	        qn.Set(cosomega, sinomega*x,  sinomega*y,  sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );
	        
	        qn.Set(cosomega, sinomega*x,  sinomega*y, -sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );

	        qn.Set(cosomega, sinomega*x, -sinomega*y,  sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );
	        
	        qn.Set(cosomega, sinomega*x, -sinomega*y, -sinomega*z);
	        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
	                t.insert( qn );

	        if (fabs(cosomega) > minRadius) {

		        qn.Set(cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
		        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
		                t.insert( qn );
		        
		        qn.Set(cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
		        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
		                t.insert( qn );

		        qn.Set(cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
		        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
		                t.insert( qn );
		        
		        qn.Set(cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
		        if ( !  t.NearestNeighbor( minRadius, qClose, qn ) )
		                t.insert( qn );
		        
	            if (antipode) {
	            
	                qs.Set(-cosomega, sinomega*x,  sinomega*y,  sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );
	                
	                qs.Set(-cosomega, sinomega*x,  sinomega*y, -sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );

	                qs.Set(-cosomega, sinomega*x, -sinomega*y,  sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );
	                
	                qs.Set(-cosomega, sinomega*x, -sinomega*y, -sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );

	                qs.Set(-cosomega, -sinomega*x,  sinomega*y,  sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );
	                
	                qs.Set(-cosomega, -sinomega*x,  sinomega*y, -sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );

	                qs.Set(-cosomega, -sinomega*x, -sinomega*y,  sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );
	                
	                qs.Set(-cosomega, -sinomega*x, -sinomega*y, -sinomega*z);
	                if ( !  t.NearestNeighbor( minRadius, qClose, qs ) )
	                t.insert( qs );
	            }
	        }
         }
         t.CompleteDelayedInsert( );
        
    }
    return( long(t.size()) - initcount );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_QGenerateGeodesic2()
 
 Description: Covers a quaternion sphere with latitudunal bands of
 geodesic spheres setting the equatorial geodesic sphere to level
 level
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_QGenerateGeodesic2( CNearTree<SQR<double> >& t, 
                             const long level)
/*-------------------------------------------------------------------------------------*/
{
    long epointcount;
    long pointcount;
    long xpointcount;
    double pi;
    double omega;
    double target;
    long N;
    long i;
    SQR<double> q = SQR<double>(1.,0.,0.,0.);
    double spacing;
    pi = 4.*atan2(1.,1.);
    epointcount = (long)(GEOM_pointcount(level) + 0.5);
    N = (long)(.5 + pi/(1.5* GEO_AverageDistanceBetweenPoints(level)));
    omega = pi/((double)(2*N));
    pointcount = 2+epointcount;

    t.clear();
    t.insert(q);
    q = SQR<double>(-1.,0.,0.,0.);
    t.insert(q);
    xpointcount = 2 + GEO_QGenerateSubGeodesic2( t,  level, pi/2.,0);
    spacing = GEOM_avgdist(level);
    for (i=1; i < N; i++) {
        long glevel;
        target = spacing/sin(((double)i)*omega);
        glevel = (long)(1.5+GEOM_levelformaxdist(target));
        xpointcount += GEO_QGenerateSubGeodesic2( t,  glevel, ((double)i)*omega, 1);
        pointcount += 2*GEO_GeodesicPointCount(glevel);
    }
    return ( xpointcount );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_HQGenerateGeodesic2()
 
 Description: Covers the northern hemisphere of a  quaternion sphere with latitudunal 
 bands of geodesic spheres setting the equatorial geodesic sphere to level
 level
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_HQGenerateGeodesic2( CNearTree<HQR<double> >& t, 
                              const long level)
/*-------------------------------------------------------------------------------------*/
{
    long epointcount;
    long pointcount;
    long xpointcount;
    double pi;
    double omega;
    double target;
    double spacing;
    long N;
    long i;
    HQR<double> q = HQR<double>(1.,0.,0.,0.);
    pi = 4.*atan2(1.,1.);
    epointcount = (long)(GEOM_pointcount(level) + 0.5);
    N = (long)(.5 + pi/(1.5* GEO_AverageDistanceBetweenPoints(level)));
    omega = pi/((double)(2*N));
    pointcount = 1+epointcount/2;
    
    t.clear();
    t.insert(q);
    xpointcount = 1 + GEO_HQGenerateSubGeodesic2( t,  level, pi/2.,0);
    spacing = GEOM_avgdist(level);
    for (i=1; i < N; i++) {
        long glevel;
        target = spacing/sin(((double)i)*omega);
        glevel = (long)(1.5+GEOM_levelformaxdist(target));
        xpointcount += GEO_HQGenerateSubGeodesic2( t,  glevel, ((double)i)*omega, 0);
        pointcount += 2*GEO_GeodesicPointCount(glevel);
    }
    return ( xpointcount );
}


