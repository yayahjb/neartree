/*
 *  geodesic.h
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
#ifndef GEODESIC_H
#define GEODESIC_H

#include "TNear.h"
#include "vector_3d.h"
#include <cmath>

#ifndef LOCAL_CBRT
#define CBRT cbrt
#else
#define CBRT local_cbrt
double local_cbrt(const double);
#endif

#define GEOM_avgdist(level)  1./(.8033489316118848*(double)(level) + .8491319884243533)
#define GEOM_avgarea(level)  1./(.6453695059219567*((double)(level))*((double)(level))  + 1.364298551396359*((double)(level)) + .7210251337654962)
#define GEOM_levelforavgdist(avgdist)  1.24478910800758/(avgdist)-1.056990250451453
#define GEOM_levelformaxdist(maxdist)  1.247376673501387/(maxdist)-1.005871648960158
#define GEOM_pointcount(level) 4.0*((double)(level))*((double)(level)) + 8.0*((double)(level)) + 6.0
#define GEOM_levelforpointcount(points) (points>6)?((long)floor(.5+(sqrt(((double)(points))-2.)/2. -1.))):0
#define GEOM_SQRpointcount(level) sqrt(GEOM_pointcount(level)*GEOM_pointcount(level)*GEOM_pointcount(level))
#define GEOM_SQRlevelforpointcount(points) (points>44)?((long)floor(sqrt(CBRT((((double)(points))*((double)(points))))-2.)/2. -1.+CBRT((double)(points))/50.)):0
#define GEOM_HQRpointcount(level) 0.5*(GEOM_SQRpointcount(level))
#define GEOM_HQRlevelforpointcount(points) ((points>22)?(GEOM_SQRlevelforpointcount(2.*(double)(points))):0)

long GEO_GenerateGeodesic2( CNearTree<Vector_3>& t, const long level );

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AtLeastNGeodesic()
 Description: Computes the points of a particular geodesic sphere (ala Buckminster
              Fuller). Uses GEO_GenerateGeodesic2 to create a relatively uniform set of
              at least minPoints points.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_AtLeastNGeodesic( CNearTree<Vector_3>& t, const long minPoints );
/*-------------------------------------------------------------------------------------*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AtLeastNQGeodesic()
 Description: Computes the points of a particular geodesic sphere (ala Buckminster
     Fuller) at the equator (w=0) of the unit quaternions and then in steps of the
     spacing of the geodesic computes proportionately small geodesic spheres at
     increasing and decreasing latitudes until it puts a single point at each pole. 
     Uses GenerateGeodesic2 to create a relatively uniform set of points.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_AtLeastNQGeodesic( CNearTree<SQR<double> >& t, const long minPoints );
/*-------------------------------------------------------------------------------------*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_AtLeastNHQGeodesic()
 Description: Computes the points of a particular geodesic sphere (ala Buckminster
 Fuller) at the equator (w=0) of the unit quaternions and then in steps of the
 spacing of the geodesic computes proportionately small geodesic spheres at
 increasing latitudes until it puts a single point at the North pole. 
 Uses GenerateGeodesic2 to create a relatively uniform set of points.
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_AtLeastNHQGeodesic( CNearTree<HQR<double> >& t, const long minPoints );
/*-------------------------------------------------------------------------------------*/


double GEO_AverageDistanceBetweenPoints( const long minPoints );

double GEO_AverageAreaPerPoint( const long minPoints );

long GEO_GeodesicPointCount( const long minPoints );

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_QGeodesicPointCount()
 Description: Computes the number of points for a particular geodesic sphere 
 (ala Buckminster Fuller) at the equator of a unit quaternion sphere, plus
 progressively smaller quaternion spheres at lines of latitude separated
 by the same spacing, plus a point at each pole.  The returned number of 
 points is the actual number that will be computed. 
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_QGeodesicPointCount( const long level );
/*-------------------------------------------------------------------------------------*/


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_HQGeodesicPointCount()
 Description: Computes the number of points for a particular geodesic sphere 
 (ala Buckminster Fuller) at the equator of a unit quaternion sphere, plus
 progressively smaller quaternion spheres at lines of latitude separated
 by the same spacing in the northern hemisphere, plus a point at the
 north pole (w = 1).  The returned number of points is the actual number that 
 will be computed. 
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_HQGeodesicPointCount( const long level );
/*-------------------------------------------------------------------------------------*/


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
long GEO_GenerateGeodesic2( CNearTree<Vector_3>& t, const long level );
/*-------------------------------------------------------------------------------------*/


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
                                const int antipode);
/*-------------------------------------------------------------------------------------*/

long GEO_HQGenerateSubGeodesic2( CNearTree<HQR<double> >& t, 
                                const long level,
                                const double omega,
                                const int antipode);
/*-------------------------------------------------------------------------------------*/


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_QGenerateGeodesic2()
 
 Description: Covers a quaternion sphere with latitudunal bands of
 geodesic spheres setting the equatorial geodesic sphere to level
 level
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_QGenerateGeodesic2( CNearTree<SQR<double> >& t, 
                             const long level);
/*-------------------------------------------------------------------------------------*/


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Name: GEO_HQGenerateGeodesic2()
 
 Description: Covers the northern hemisphere of a  quaternion sphere with latitudunal 
 bands of geodesic spheres setting the equatorial geodesic sphere to level
 level
 
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
long GEO_HQGenerateGeodesic2( CNearTree<HQR<double> >& t, 
                              const long level);
/*-------------------------------------------------------------------------------------*/


#endif // GEODESIC_H
