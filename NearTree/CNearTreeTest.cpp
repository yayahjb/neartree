//*
//*  CNearTreeTest.cpp
//*  NearTree
//*
//*  test harness for the templated neartree implementation, TNear.h
//*  Copyright 2008 Larry Andrews.  All rights reserved
//*

//**********************************************************************
//*                                                                    *
//* YOU MAY REDISTRIBUTE NearTree UNDER THE TERMS OF THE LGPL          *
//*                                                                    *
//**********************************************************************/

//************************* LGPL NOTICES *******************************
//*                                                                    *
//* This library is free software; you can redistribute it and/or      *
//* modify it under the terms of the GNU Lesser General Public         *
//* License as published by the Free Software Foundation; either       *
//* version 2.1 of the License, or (at your option) any later version. *
//*                                                                    *
//* This library is distributed in the hope that it will be useful,    *
//* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
//* Lesser General Public License for more details.                    *
//*                                                                    *
//* You should have received a copy of the GNU Lesser General Public   *
//* License along with this library; if not, write to the Free         *
//* Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
//* MA  02110-1301  USA                                                *
//*                                                                    *
//**********************************************************************/

/*
 This is a test harness for the templated neartree implementation, CNearTree.
 
 */


#include <float.h>
#include <limits.h>
#include <math.h>
#ifndef USE_LOCAL_HEADERS
#include <TNear.h>
#else
#include "TNear.h"
#endif

void testEmptyTree( void );
void testLinearTree( const int n );
void testFindFirstObject( void );
void testFindLastObject( void );
void testFindInSphereFromBottom( void );
void testFindInSphereFromTop( void );

void testRandomTree( const int n );
void testBigVector( );

long g_errorCount;

#define MYRAND_MAX  32767

/*=======================================================================*/
int main(int argc, char* argv[])
{
    
    /* test the interface with an empty tree */
    testEmptyTree( );
    
    /* test the interface with trees with varying content, one entry and several */
    for( int i=1; i<10; ++i )
    {
        testLinearTree( i );
    }
    
    testFindFirstObject( );
    testFindLastObject( );
    testFindInSphereFromBottom( );
    testFindInSphereFromTop( );
    testRandomTree( 10000 );
    testBigVector( );
    
    if( g_errorCount == 0 )
    {
        printf( "No errors were detected while testing CNearTree\n" );
    }
    
    return g_errorCount;
}

/*
 For an empty tree of int's, test the public interface for CNearTree.
 */
/*=======================================================================*/
void testEmptyTree( void )
{
    CNearTree<int> tree;
    int close;
    int nFar;
    std::vector<int> v;
    bool bTreeEmpty;
    bool bTreeHasNearest;
    bool bTreeHasFarthest;
    long lFoundPointsInSphere;
    
    bTreeEmpty = tree.empty( );
    if( ! bTreeEmpty )
    {
        ++g_errorCount;
        printf( ".empty incorrect for empty tree\n" );
    }
    
    bTreeHasNearest  = tree.NearestNeighbor( 0.0, close, 1 );
    if( bTreeHasNearest )
    {
        ++g_errorCount;
        printf( "NearestNeighbor incorrect for empty tree\n" );
    }
    
    bTreeHasFarthest = tree.FarthestNeighbor( nFar, 0 );
    if( bTreeHasFarthest )
    {
        ++g_errorCount;
        printf( "FarthestNeighbor incorrect for empty tree\n" );
    }
    
    lFoundPointsInSphere = tree.FindInSphere( 1000.0, v, 1 );
    if( lFoundPointsInSphere != 0 )
    {
        ++g_errorCount;
        printf( "FindInSphere incorrect for empty tree\n" );
    }
}

/*
 Test CNearTree using a tree built only using int's. This is not a robust test
 in some ways because the tree is only composed of right branches and left leaves.
 Perform some general tests.
 */
/*=======================================================================*/
void testLinearTree( const int n )
{
    CNearTree<int> tree;
    
    /* generate an unbalanced tree*/
    for( int i=1; i<=n; ++i )
    {
        tree.Insert( i );
    }
    
    /*
     Search for the nearest value using a probe point that is larger than 
     the largest value that was input. The returned values should be the
     last value entered into the tree. 
     */
    int closest=-1;   
    const bool bClose = tree.NearestNeighbor( 22, closest, 2*n );
    if( ! bClose || closest != n )
    {
        ++g_errorCount;
        printf( "NearestNeighbor failed in testLinearTree, got %d, should be %d\n", closest, n );
    }
    
    /*
     Search for the farthest value using a probe point that is larger than 
     the largest value that was input. The returned values should be the
     first value entered into the tree. 
     */
    int farthest;
    const bool bFar = tree.FarthestNeighbor( farthest, 2*n );
    if( ! bFar || farthest != 1 )
    {
        ++g_errorCount;
        printf( "FarthestNeighbor failed in testLinearTree, got %d, should be %d\n", farthest, n );
    }
    
    /*
     Find all of the points in the tree using a negative radius, using the first
     input point as the probe point. Nothing should be found.
     */
    std::vector<int> v;
    if( tree.FindInSphere( -100.0, v, 1 ) != 0 )
    {
        ++g_errorCount;
        printf( "FindInSphere found points for negative radius\n" );
    }
    v.clear( );
    
    /*
     Find all of the points in the tree using a small radius. Do separate
     searches for every point entered. In every case, only a single point
     should be found.
     */
   long localErrorCount = 0;
   long localErrorMax = 0;
    for( int i=1; i<=n; ++i )
    {
        const int found = tree.FindInSphere( 0.1, v, i );
        if( found != 1 )
        {
         ++localErrorCount;
         if( found > localErrorMax ) 
         {
            localErrorMax = found;
         }
      }
   }

   if( localErrorCount != 0 )
   {
            ++g_errorCount;
      printf( "FindInSphere found too many points (as many as %ld) %ld times\n", localErrorMax, localErrorCount );
    }
    v.clear( );
    
    /*
     Find all of the points in the tree that are within a large radius, using 
     the first input point as the probe point. All of the input points should 
     be found within the radius.
     */
    if( tree.FindInSphere( (double)(10*n), v, 0 ) != n )
    {
        ++g_errorCount;
        printf( "FindInSphere did not find all the points, found %ld\n", tree.FindInSphere( 1000, v, 0 ) );
    }
    v.clear( );
    
}

/*
 Perform general tests using floating point numbers. Two test sets are
 included, one for float and one for double.
 The values are computed starting from some initial value, and
 each succeeding value is one half of the previous until zero is
 computed (the zero is NOT inserted into the tree). The tree will consist
 of only right branches and left leaves.
 */
/*=======================================================================*/
void testFindFirstObject( void )
{
    {
        double fFinal = FLT_MAX; /* just initialization */
        CNearTree<float> tree;
        long count = 0;
        
        /* build the float tree starting with 1.0 */
        float f = 1.0;
        while( f > 0.0 )
        {
            tree.Insert( f );
            fFinal = f;
            f /= 2.0;
            ++count;
        }
        
        if( tree.empty( ) )
        {
            ++g_errorCount;
            printf( "testFindFirstObject incorrectly found empty tree for float\n" );
        }
        
        /*
         Search for the value closest to zero. It should be a very small, probably
         denormalized number.
         */
        float closest = 0.0;
        const bool bReturnNear = tree.NearestNeighbor( 1.0e-10, closest, 0.0 );
        if( ! bReturnNear || closest != fFinal )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for float, got %f\n", closest );
        }
        
        /*
         Search for the value farthest from a large number. It should be a
         very small, probably denormalized number.
         */
        float farthest = 0.0;
        const bool bReturnFar = tree.FarthestNeighbor( farthest, 100.0 );
        if( ! bReturnFar || farthest != fFinal )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for float, got %f\n", farthest );
        }
        
        /*
         Determine if FindInSphere can find all of the input data.
         */
        std::vector<float> v;
        const long lFound = tree.FindInSphere( 100.0, v, 1.0 );
        if( lFound != count )
        {
            ++g_errorCount;
            printf( "testFindFirstObject: found wrong count for FindInSphere for float, should be%ld, got %ld\n", count, lFound );
        }
    }
    
    {
        double dFinal = DBL_MAX; /* just initialization */
        CNearTree<double> tree;
        long count = 0;
        
        /* build the double tree starting with 1.0 */
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 )
        {
            tree.Insert( f );
            dFinal = f;
            f /= 2.0;
            ++count;
        }
        
        if( tree.empty( ) )
        {
            ++g_errorCount;
            printf( "testFindFirstObject incorrectly found empty tree for double\n" );
        }
        
        /*
         Search for the value closest to zero. It should be a very small, probably
         denormalized number.
         */
        double closest = 0.0;
        const bool bReturnNear = tree.NearestNeighbor( 1.0e-10, closest, 0.0 );
        if( ! bReturnNear || closest != dFinal )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for double, got %f\n", closest );
        }
        
        /*
         Search for the value farthest from a large number. It should be a
         very small, probably denormalized number.
         */
        double farthest = 0.0;
        const bool bReturn = tree.FarthestNeighbor( farthest, 100.0 );
        if( ! bReturn || farthest != dFinal )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for double, got %f\n", farthest );
        }
        
        /*
         Determine if FindInSphere can find all of the input data.
         */
        std::vector<double> v;
        const long lFound = tree.FindInSphere( 100.0, v, 1.0 );
        if( lFound != count )
        {
            ++g_errorCount;
            printf( "testFindFirstObject: found wrong count for FindInSphere for double, should be%ld, got %ld\n", count, lFound );
        }
    }
}

/*
 Perform general tests using floating point numbers.
 Build a tree of floats. The values are computed starting from some initial 
 value, and each succeeding value is one half of the previous until zero
 is computed (the zero is NOT inserted into the tree). The tree will consist
 of only right branches and left leaves.
 */
/*=======================================================================*/
void testFindLastObject( void )
{
    {
        double fFinal = FLT_MAX;
        CNearTree<float> tree;
        long count = 0;
        
        float f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 )
        {
            tree.Insert( f );
            fFinal = f;
            f /= 2.0;
            ++count;
        }
        
        float closest = 0.0;
        const bool bReturn = tree.NearestNeighbor( 100.0, closest, 11.0 );
        if( ! bReturn || closest != 1.0 )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for float, got %f\n", closest );
        }
    }
    
    {
        double dFinal = DBL_MAX;
        CNearTree<double> tree;
        long count = 0;
        
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 )
        {
            tree.Insert( f );
            dFinal = f;
            f /= 2.0;
            ++count;
        }
        
        double closest = 0.0;
        const bool bReturn = tree.NearestNeighbor( 100.0, closest, 11.0 );
        if( ! bReturn || closest != 1.0 )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for float, got %f\n", closest );
        }
    }
    {
        double fFinal = FLT_MAX;
        CNearTree<float> tree;
        long count = 0;
        
        float f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 )
        {
            tree.Insert( f );
            fFinal = f;
            f /= 2.0;
            ++count;
        }
        
        float closest = 0.0;
        const bool bReturn = tree.FarthestNeighbor( closest, -11.0 );
        if( ! bReturn || closest != 1.0 )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for float, got %f\n", closest );
        }
    }
    
    {
        double dFinal = DBL_MAX;
        CNearTree<double> tree;
        long count = 0;
        
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 )
        {
            tree.Insert( f );
            dFinal = f;
            f /= 2.0;
            ++count;
        }
        
        double closest = 0.0;
        const bool bReturn = tree.FarthestNeighbor( closest, -11.0 );
        if( ! bReturn || closest != 1.0 )
        {
            ++g_errorCount;
            printf( "testFindFirstObject failed for float, got %f\n", closest );
        }
    }
}

/*=======================================================================*/
void testFindInSphereFromBottom( void )
{
    const int nmax = 100;
    CNearTree<double> tree;
    
    for( int i=1; i<=nmax; ++i )
    {
        tree.Insert( (double)i );
    }
    
    std::vector<double> v;
    /* generate an unbalanced tree*/
    for( int i=1; i<=nmax; ++i )
    {
        v.clear( );
        const double radius = 0.05 + (double)i;
        const long lReturned = tree.FindInSphere( radius, v, 0.9 );
        if( lReturned != (long)i )
        {
            ++g_errorCount;
            printf( "FindInSphere failed in testFindInSphereFromBottom for i=%d\n", i );
        }
    }
}

/*=======================================================================*/
void testFindInSphereFromTop( void )
{
    const int nmax = 100;
    CNearTree<double> tree;
    
    for( int i=1; i<=nmax; ++i )
    {
        tree.Insert( (double)i );
    }
    
    std::vector<double> v;
    /* generate an unbalanced tree*/
    for( int i=1; i<=nmax; ++i )
    {
        v.clear( );
        const double radius = 0.05 + (double)i;
        const long lReturned = tree.FindInSphere( radius, v, 0.1 + (double)nmax );
        if( lReturned != (long)i )
        {
            ++g_errorCount;
            printf( "FindInSphere failed in testFindInSphereFromTop for i=%d\n", i );
        }
    }
}

/*
 Test CNearTree with a bunch of random numbers (integers).
 */
/*=======================================================================*/
void testRandomTree( const int nRequestedRandoms )
{
   const int n = nRequestedRandoms;
    CNearTree<int> tree;
    int nmax = INT_MIN;
    int nmin = INT_MAX;
    
    /* Build the tree with n random numbers. Remember the largest and smallest values. */
    for( int i=0; i<n; ++i )
    {
        const int next = rand( )%MYRAND_MAX;
        tree.Insert( next );
        if( next > nmax ) nmax = next;
        if( next < nmin ) nmin = next;
    }
    
    {
        /*verify that the correct extremal point is detected (from below)*/
        int closest=INT_MAX;
        const bool bNear = tree.NearestNeighbor( (double)LONG_MAX, closest, INT_MIN/2 );
        if( ! bNear || closest != nmin )
        {
            ++g_errorCount;
            printf( "NearestNeighbor failed in testRandomTree for min\n" );
        }
        
        int farthest;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MIN/2 );
        if( !bFar || farthest != nmax )
        {
            ++g_errorCount;
            printf( "FarthestNeighbor failed in testRandomTree for min\n" );
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from above)*/
        int closest=INT_MIN;
        const bool bNear = tree.NearestNeighbor( (double)LONG_MAX, closest, INT_MAX/2 );
        if( ! bNear || closest != nmax )
        {
            ++g_errorCount;
            printf( "NearestNeighbor failed in testRandomTree for max\n" );
        }
        
        int farthest;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MAX/2 );
        if( !bFar || farthest != nmin )
        {
            ++g_errorCount;
            printf( "FarthestNeighbor failed in testRandomTree for max\n" );
        }
    }
    
    {
        /*verify that for very large radius, every point is detected (from below)*/
        std::vector<int> v;
        const double radius = DBL_MAX;
        const long lReturn = tree.FindInSphere( radius, v, INT_MIN/2 );
        if( lReturn != n )
        {
            ++g_errorCount;
            printf( "FindInSphere failed in testRandomTree, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*verify that we find NO points if we are below the lowest and with too small radius*/
        std::vector<int> v;
        const double radius = .5;
        const long lReturn = tree.FindInSphere( radius, v, nmin-1 );
        if( lReturn != 0 )
        {
            ++g_errorCount;
            printf( "FindInSphere failed in testRandomTree found points incorrectly, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*test that the number of found points in a sphere is non-deceasing with increasing radius*/
        std::vector<int> v;
        long lReturn = 0;
        int lastFoundCount = 0;
        int cycleCount = 0;
        
        double radius = 0.00001; /* start with a very small radius (remember these are int's) */
        while( radius < 5*(nmax-nmin) )
        {
            ++cycleCount;
            lReturn = tree.FindInSphere( radius, v, nmin-1 );
            if( lReturn < lastFoundCount )
            {
                ++g_errorCount;
                printf( "FindInSphere in testRandomTree found DECREASING count on increasing radius for radius=%f\n", radius );
                break;
            }
            else
            {
                lastFoundCount = lReturn;
                radius *= 1.414;
            }
        }
        
        /* verify that all points were included in the final check (beyond all reasonable distances) */
        if( lReturn != n )
        {
            ++g_errorCount;
            printf( "FindInSphere in testRandomTree did not find all the points\n" );
        }
    }
}

/*=======================================================================*/
/* make a 17-dimension vector class for testing */
class vec17
    {
    public:
        double pd[17];
        int dim;
        vec17( )
        {
            dim = 17;
            for( int i=0; i<dim; ++i )
            {
                pd[i] = (double)(rand( )%MYRAND_MAX);
            }
        }
        vec17( const double d )
        {
            dim = 17;
            for( int i=0; i<dim; ++i )
            {
                pd[i] = d;
            }
        }
        ~vec17( void )
        {
        }
        operator double( void ) const
        {
            double dtemp = 0.0;
            for( int i=0; i<dim; ++i )
            {
                dtemp = pd[i]*pd[i];
            }
            return( ::sqrt( dtemp ) );
        }
        vec17 operator-( const vec17& v ) const /* USERS: be sure to make both const */
        {
            vec17 vtemp;
            for( int i=0; i<dim; ++i )
            {
                vtemp.pd[i] = this->pd[i] - v.pd[i];
            }
            return( vtemp );
        }
    };


/*=======================================================================*/
void testBigVector(  )
{
    
    CNearTree<vec17> tree;
    vec17 vAll[1000]; /* keep a list of all of the input so we can find particular entries */
    double rmax = -DBL_MAX;
    double rmin =  DBL_MAX;
    vec17 v17min; /* to be the point nearest to the origin */
    vec17 v17max; /* to be the point farthest from the origin */
    
    /* All of the coordinate values will be in the range 0-MYRAND_MAX. In other words,
     all of the data points will be within a 17-D cube that has a corner at
     the origin of the space.
     */
    for( int i=0; i<1000; ++i )
    {
        vec17 v;
        if( double( v ) < rmin )
        {
            rmin = double( v );
            v17min = v;
        }
        
        if( double( v ) > rmax )
        {
            rmax = double( v );
            v17max = v;
        }
        
        tree.Insert( v );
        vAll[i] = v;
    }
    
    {
        /* Find the point farthest from the point that was nearest the origin. */
        vec17 vFarthest;
        tree.FarthestNeighbor( vFarthest, v17min );
        
        /* Brute force search for the farthest */
        vec17 vSearch;
        double dmax = -DBL_MAX;
        for( int i=0; i<1000; ++i )
        {
            if( double( vAll[i] - v17min ) > dmax )
            {
                dmax = double( vAll[i] - v17min );
                vSearch = vAll[i];
            }
        }
        
        if( double(vSearch-v17max) > DBL_MIN )
        {
            ++g_errorCount;
            printf( "in testBigVector, apparently FarthestNeighbor has failed\n" );
        }
    }
    
    {
        /* somewhere in the middle, find a point and its nearest neighbor */
        /* make sure that each includes the other in sphere search */
        
        std::vector<vec17> v;
        vec17 vCenter( (double)(MYRAND_MAX/2) );
        vec17 vNearCenter;
        vec17 vCloseToNearCenter;
        tree.NearestNeighbor( 1000.0, vNearCenter, vCenter );
        long iFoundNearCenter = tree.FindInSphere( 100.0, v, vNearCenter );
        
        /* Brute force search for the point closest to the point closest to the center */
        double dmin = DBL_MAX;
        for( int i=0; i<iFoundNearCenter; ++i )
        {
            if( vNearCenter != v[i] && double( vNearCenter-v[i] ) < dmin )
            {
                dmin = double( vNearCenter-v[i] );
                vCloseToNearCenter = v[i];
            }
        }
        
        if( dmin == DBL_MAX ) 
        {
            ++g_errorCount;
            printf( "testBigVector: apparently FindInSphere failed\n" );
        }
        
        {
            /* Using zero radius, check that only one point is found when a point is searched
             with FindInSphere */
            v.clear( );
            long iFound = tree.FindInSphere( 0.0, v, vNearCenter );
            if( iFound < 1 )
            {
                ++g_errorCount;
                printf( "testBigVector: FindInSphere found no points using zero radius\n" );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                printf( "testBigVector: FindInSphere found more than %ld points using zero radius\n", iFound );
            }
        }
        
        {
            /* Using minimal radius, check that at least 2 points are found when a point is searched
             with FindInSphere */
            v.clear( );
            long iFound = tree.FindInSphere( (double)(vCloseToNearCenter-vNearCenter), v, vNearCenter );
            if( iFound < 2 )
            {
                ++g_errorCount;
                printf( "testBigVector: FindInSphere found only 1 point\n" );
            }
        }
        
        {
            /* Using small radius, check that only one point is found when a point is searched
             with FindInSphere */
            v.clear( );
            long iFound = tree.FindInSphere( (double)(vCloseToNearCenter-vNearCenter)*0.9, v, vNearCenter );
            if( iFound < 1 )
            {
                ++g_errorCount;
                printf( "testBigVector: FindInSphere found no points using %f radius\n", (double)(vCloseToNearCenter-vNearCenter)*0.9 );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                printf( "testBigVector: FindInSphere found %ld points using %f radius\n", iFound, (double)(vCloseToNearCenter-vNearCenter)*0.9 );
            }
        }
        
        /* Just make sure that FarthestNeighbor works for objects */
        vec17 vExtreme;
        tree.FarthestNeighbor( vExtreme, vNearCenter );
    }
}



