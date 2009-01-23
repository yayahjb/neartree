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
#include <stdio.h>
#include <stdlib.h>
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
void testBackwardForward( void );
void testFindInSphereFromTop( void );
void testDelayedInsertion( void );
void testIterators( void );

void testRandomTree( const int n );
void testBigVector( void );

long g_errorCount;

int debug;

/*=======================================================================*/
int main(int argc, char* argv[])
{
    
    g_errorCount = 0;
    debug = 0;
    if (argc > 1 && !strcmp(argv[1],"--debug")) debug = 1;
    
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
    testBackwardForward( );
    testDelayedInsertion( );
    testIterators( );
    
    if( g_errorCount == 0 )
    {
        fprintf(stdout,  "No errors were detected while testing CNearTree\n" );
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
        fprintf(stdout, "testEmptyTree: ..empty incorrect for empty tree\n" );
    }
    if( tree.size( ) != 0 )
    {
        ++g_errorCount;
        fprintf(stdout, "testEmptyTree: .size incorrect for empty tree\n" );
    }
    
    bTreeHasNearest  = tree.NearestNeighbor( 0.0, close, 1 );
    if( bTreeHasNearest )
    {
        ++g_errorCount;
        fprintf(stdout, "testEmptyTree: .NearestNeighbor incorrect for empty tree\n" );
    }
    
    bTreeHasFarthest = tree.FarthestNeighbor( nFar, 0 );
    if( bTreeHasFarthest )
    {
        ++g_errorCount;
        fprintf(stdout, "testEmptyTree: .FarthestNeighbor incorrect for empty tree\n" );
    }
    
    lFoundPointsInSphere = tree.FindInSphere( 1000.0, v, 1 );
    if( lFoundPointsInSphere != 0 )
    {
        ++g_errorCount;
        fprintf(stdout, "testEmptyTree: .FindInSphere incorrect for empty tree\n" );
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
    
    if( tree.GetDepth( ) > (size_t)(n+1)/2 )
    {
        ++g_errorCount;
        fprintf(stdout, "testDelayedInsertion: tree depth is too large, %lu is greater than %d\n", (unsigned long)tree.GetDepth( ), (n+1)/2 );
    }
    /*
     Search for the nearest value using a probe point that is larger than 
     the largest value that was input. The returned values should be the
     last value entered into the tree. 
     */
    int closest=-1;   
    const bool bClose = tree.NearestNeighbor( 22.0, closest, 2*n );
    if( ! bClose )
    {
        ++g_errorCount;
        fprintf(stdout, "testLinearTree: NearestNeighbor failed to find anything\n" );
    }
    else if( closest != n )
    {
        ++g_errorCount;
        fprintf(stdout, "NearestNeighbor failed in testLinearTree, got %d, should be %d\n", closest, n );
    }
    
    /*
     Search for the farthest value using a probe point that is larger than 
     the largest value that was input. The returned values should be the
     first value entered into the tree. 
     */
    int farthest;
    const bool bFar = tree.FarthestNeighbor( farthest, 2*n );
    if( ! bFar )
    {
        ++g_errorCount;
        fprintf(stdout, "testLinearTree: FarthestNeighbor failed to find anything\n" );
    }
    else if( farthest != 1 )
    {
        ++g_errorCount;
        fprintf(stdout, "FarthestNeighbor failed in testLinearTree, got %d, should be %d\n", farthest, n );
    }
    
    /*
     Find all of the points in the tree using a negative radius, using the first
     input point as the probe point. Nothing should be found.
     */
    std::vector<int> v;
    if( tree.FindInSphere( -100.0, v, 1 ) != 0 )
    {
        ++g_errorCount;
        fprintf(stdout, "FindInSphere found points for negative radius\n" );
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
        fprintf(stdout, "FindInSphere found too many points (as many as %ld) %ld times\n", localErrorMax, localErrorCount );
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
        fprintf(stdout, "FindInSphere did not find all the points, found %ld\n", tree.FindInSphere( (double)(10*n), v, 0 ) );
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
            fprintf(stdout, "testFindFirstObject incorrectly found empty tree for float\n" );
        }
        
        //   /*
        //   Search for the value closest to zero. It should be a very small, probably
        //   denormalized number.
        //   */
        float closest = 0.0;
        const bool bReturnNear = tree.NearestNeighbor( 1.0e-10, closest, 0.0 );
        if( ! bReturnNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject: Near failed to find anything\n" );
        }
        else if( closest != fFinal )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject: Near failed for float, got %f\n", closest );
        }
        
        //   /*
        //   Search for the value farthest from a large number. It should be a
        //   very small, probably denormalized number.
        //   */
        float farthest = -1000.0;
        const bool bReturnFar = tree.FarthestNeighbor( farthest, 100.0 );
        if( ! bReturnFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject Far failed to any anything for floatn" );
        }
        else if( farthest != fFinal )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject Far failed for float, got %g\n", farthest );
        }
        
        //   /*
        //   Determine if FindInSphere can find all of the input data.
        //   */
        std::vector<float> v;
        const long lFound = tree.FindInSphere( 100.0, v, 1.0 );
        if( lFound != count )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject: found wrong count for FindInSphere for float, should be%ld, got %ld\n", count, lFound );
        }
        }
        
    {
        double dFinal = DBL_MAX; /* just initialization */
        CNearTree<double> tree;
        long count = 0;
        
        /* build the double tree starting with 1.0 */
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 && f >= DBL_MIN)
        {
            tree.Insert( f );
            dFinal = f;
            f /= 2.0;
            ++count;
        }
        
        if( tree.empty( ) )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject incorrectly found empty tree for double\n" );
        }
        
        /*
         Search for the value closest to zero. It should be a very small, probably
         denormalized number.
         */
        double closest = 0.0;
        const bool bReturnNear = tree.NearestNeighbor( 1.0e-10, closest, 0.0 );
        if( ! bReturnNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject: Near failed to find anything\n" );
        }
        else if( closest != dFinal )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject: Near failed for double, got %f\n", closest );
        }
        
        /*
         Search for the value farthest from a large number. It should be a
         very small, probably denormalized number.
         */
        double farthest = 10000.0;
        const bool bReturn = tree.FarthestNeighbor( farthest, 100.0 );
        if( ! bReturn )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject: Far failed to find anything for double\n"  );
        }
        else if( farthest != dFinal )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject Far failed for double, got %g\n", farthest );
        }
        
        /*
         Determine if FindInSphere can find all of the input data.
         */
        std::vector<double> v;
        const long lFound = tree.FindInSphere( 100.0, v, 1.0 );
        if( lFound != count )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindFirstObject: found wrong count for FindInSphere for double, should be%ld, got %ld\n", count, lFound );
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
        CNearTree<float> tree;
        
        float f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0  && f >= FLT_MIN )
        {
            tree.Insert( f );
            f /= 2.0;
        }
        
        float closest = 0.0;
        const bool bReturn = tree.NearestNeighbor( 100.0, closest, 11.0 );
        if( ! bReturn )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Near failed to find anything\n" );
        }
        else if( closest != 1.0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Near failed for float, got %g\n", closest );
        }
    }
    
    {
        CNearTree<double> tree;
        
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 && f >= DBL_MIN)
        {
            tree.Insert( f );
            f /= 2.0;
        }
        
        double closest = 0.0;
        const bool bReturn = tree.NearestNeighbor( 100.0, closest, 11.0 );
        if( ! bReturn )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Near failed to find anything\n" );
        }
        else if( closest != 1.0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Near failed for float, got %g\n", closest );
        }
    }
    
    {
        CNearTree<float> tree;
        
        float f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 && f >= FLT_MIN)
        {
            tree.Insert( f );
            f /= 2.0;
        }
        
        float closest = 0.0;
        const bool bReturn = tree.FarthestNeighbor( closest, -11.0 );
        if( ! bReturn )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Far failed to find anything\n" );
        }
        else if( closest != 1.0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Far failed for float, got %f\n", closest );
        }
        }
        
    {
        CNearTree<double> tree;
        
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 && f >= DBL_MIN)
        {
            tree.Insert( f );
            f /= 2.0;
        }
        
        double closest = 0.0;
        const bool bReturn = tree.FarthestNeighbor( closest, -11.0 );
        if( ! bReturn )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Far failed to find anything\n" );
        }
        else if( closest != 1.0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject: Far failed for float, got %f\n", closest );
        }
        }
        
    {
        double dFinal = DBL_MAX;
        CNearTree<double> tree;
        
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 && f >= DBL_MIN)
        {
            dFinal = f;
            f /= 2.0;
        }
        
        f = DBL_MIN;
        int count = 0;
        while( f < 1.01 )
        {
            tree.Insert( f );
            dFinal = f;
            f *= 2.0;
            ++count;
        }
        
        double farthest = 0.0;
        const bool bReturnNear = tree.NearestNeighbor( 1000.0, farthest, 1000.0 );
        if( ! bReturnNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject Near failed to find anything for double\n" );
        }
        else if( farthest != 1.0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject Near failed for double, reverse tree, got %f\n", farthest );
        }
        
        double closest = 0.0;
        const bool bReturnFar = tree.FarthestNeighbor( closest, -11.0 );
        if( ! bReturnFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject Far failed to find anything for double\n" );
        }
        else if( closest != 1.0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindLastObject Far failed for double, reverse tree, got %f\n", closest );
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
            fprintf(stdout, "FindInSphere failed in testFindInSphereFromBottom for i=%d\n", i );
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
            fprintf(stdout, "FindInSphere failed in testFindInSphereFromTop for i=%d\n", i );
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
        if( ! bNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: NearestNeighbor failed to find anything\n" );
        }
        else if( closest != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: NearestNeighbor failed for min\n" );
        }
        
        int farthest=INT_MIN;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MIN/2 );
        if( !bFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: FarthestNeighbor failed to find anything\n" );
        }
        else if( farthest != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: FarthestNeighbor failed in testRandomTree for min\n" );
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from above)*/
        int closest=INT_MIN;
        const bool bNear = tree.NearestNeighbor( (double)LONG_MAX, closest, INT_MAX/2 );
        if( ! bNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: NearestNeighbor failed to find anything for max\n" );
        }
        else if( closest != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: NearestNeighbor failed for max\n" );
        }
        
        int farthest;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MAX/2 );
        if( !bFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: FarthestNeighbor failed to find anything for min\n" );
        }
        else if( farthest != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree: FarthestNeighbor failed for min\n" );
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
            fprintf(stdout, "FindInSphere failed in testRandomTree, n=%d, lReturn=%ld\n", n, lReturn );
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
            fprintf(stdout, "FindInSphere failed in testRandomTree found points incorrectly, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*test that the number of found points in a sphere is non-deceasing with increasing radius*/
        std::vector<int> v;
        long lReturn = 0;
        int lastFoundCount = 0;
        
        double radius = 0.00001; /* start with a very small radius (remember these are int's) */
        while( radius < (double)(5*(nmax-nmin)) )
        {
            lReturn = tree.FindInSphere( radius, v, nmin-1 );
            if( lReturn < lastFoundCount )
            {
                ++g_errorCount;
                fprintf(stdout, "FindInSphere in testRandomTree found DECREASING count on increasing radius for radius=%f\n", radius );
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
            fprintf(stdout, "FindInSphere in testRandomTree did not find all the points\n" );
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
        explicit vec17( const double d )
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
            fprintf(stdout, "in testBigVector, apparently FarthestNeighbor has failed\n" );
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
        unsigned long iFoundNearCenter = (unsigned long)tree.FindInSphere( 100.0, v, vNearCenter );
        
        /* Brute force search for the point closest to the point closest to the center */
        double dmin = DBL_MAX;
        for( unsigned long i=0; i<iFoundNearCenter; ++i )
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
            fprintf(stdout, "testBigVector: apparently FindInSphere failed\n" );
        }
        
        {
            /* Using zero radius, check that only one point is found when a point is searched
             with FindInSphere */
            v.clear( );
            long iFound = tree.FindInSphere( 0.0, v, vNearCenter );
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigVector: FindInSphere found no points using zero radius\n" );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigVector: FindInSphere found more than %ld points using zero radius\n", iFound );
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
                fprintf(stdout, "testBigVector: FindInSphere found only 1 point\n" );
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
                fprintf(stdout, "testBigVector: FindInSphere found no points using %f radius\n", (double)(vCloseToNearCenter-vNearCenter)*0.9 );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigVector: FindInSphere found %ld points using %f radius\n", iFound, (double)(vCloseToNearCenter-vNearCenter)*0.9 );
            }
        }
        
        /* Just make sure that FarthestNeighbor works for objects */
        vec17 vExtreme;
        tree.FarthestNeighbor( vExtreme, vNearCenter );
    }
    }
    
/*=======================================================================*/
void testBackwardForward( void )
{
    CNearTree<double> tree;
    const int nMax = 1000;
    
    for( int i=0; i<nMax; ++i )
    {
        tree.Insert( (double)i );
        tree.Insert( (double)(nMax-i) );
        tree.Insert( (double)i + 0.25 );
        tree.Insert( (double)(nMax-i) + 0.75 );
    }
    
    for( int i=100; i<300; ++i )
    {
        double closest;
        const bool bReturn = tree.NearestNeighbor( 1000.0, closest, (double)i+0.25 );
        if( ! bReturn )
        {
            ++g_errorCount;
            fprintf(stdout, "testBackwardForward: NearestNeighbor failed to find anything\n" );
        }
        else if( ::fabs( closest-((double)i+0.25) ) > DBL_MIN )
        {
            fprintf(stdout, "testBackwardForward::NearestNeighbor failed, closest=%f\n", closest );
        }
    }
}

/*=======================================================================*/
void testDelayedInsertion( void )
{
    
    {
        // make sure that CompleteDelayedInsert flushes the delayed data
        const long nmax = 10001;
        CNearTree<double> tree;
        
        for( int i=1; i<=nmax; ++i )
        {
            const double insertValue = (double)(i);
            tree.DelayedInsert( insertValue );
        }
        
        tree.CompleteDelayedInsert( );
        const size_t depth = tree.GetDepth( );
        const size_t insertedSize = tree.GetTotalSize( );
        const size_t delayed = tree.GetDeferredSize( );
        
        if( delayed != 0 || (long)insertedSize != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testDelayedInsertion: CompleteDelayedInsert completion is incorrect\n" );
        }
        else if( depth >=nmax/2 )
        {
            ++g_errorCount;
            fprintf(stdout, "testDelayedInsertion: tree depth is too large, %lu is greater than %ld\n", (unsigned long)depth, nmax/2 );
        }
    }
    
    {
        // make sure that FindInSphere flushes the delayed data
        const long nmax = 100;
        CNearTree<double> tree;
        
        for( int i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 )
            {
                tree.Insert( (double)i );
            }
            else
            {
                tree.DelayedInsert( (double)i );
            }
        }
        
        std::vector<double> v;
        const double radius = 1000.0;
        const long lReturned = tree.FindInSphere( radius, v, 0.9 );
        
        if( lReturned != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testDelayedInsertion: FindInSphere failed for nmax=%ld, found %ld points\n", nmax, lReturned );
        }
    }
    
    {
        // make sure that NearestNeighbor flushes the delayed data
        const long nmax = 100;
        CNearTree<double> tree;
        double fFinal = DBL_MAX;
        
        for( int i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 && i!=nmax) // ensure that the last one is delayed
            {
                tree.Insert( (double)i );
            }
            else
            {
                tree.DelayedInsert( (double)i );
                fFinal = (double)i;
            }
        }
        
        double closest;
        const double radius = 0.1;
        const bool bReturned = tree.NearestNeighbor( radius, closest, fFinal );
        
        if( ! bReturned )
        {
            ++g_errorCount;
            fprintf(stdout, "testDelayedInsertion: NearestNeighbor failed\n" );
        }
        else if( closest != fFinal )
        {
            ++g_errorCount;
            fprintf(stdout, "testDelayedInsertion: NearestNeighbor failed to find the data\n" );
        }
    }
    
    {
        // make sure that FarthestNeighbor flushes the delayed data
        const long nmax = 100;
        CNearTree<double> tree;
        double fFinal = DBL_MAX;
        
        for( int i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 && i!=nmax) // ensure that the last one is delayed
            {
                tree.Insert( (double)i );
            }
            else
            {
                tree.DelayedInsert( (double)i );
                fFinal = (double)i;
            }
        }
        
        double farthest;
        const bool bReturned = tree.FarthestNeighbor( farthest, -100.0 );
        
        if( ! bReturned )
        {
            ++g_errorCount;
            fprintf(stdout, "testDelayedInsertion: FarthestNeighbor failed\n" );
        }
        else if( farthest != fFinal )
        {
            ++g_errorCount;
            fprintf(stdout, "testDelayedInsertion: FarthestNeighbor failed to find the data\n" );
        }
    }
    
}

/*=======================================================================*/
void testIterators( void )
{
    CNearTree<int> tree;
    CNearTree<int>::iterator itEmpty = tree.back( );
    if( itEmpty != tree.end( ) )
    {
        ++g_errorCount;
        fprintf(stdout, "testIterators: back failed for empty tree\n" );
    }
    
    if( itEmpty != tree.end( ) )
    {
        ++g_errorCount;
        fprintf(stdout, "testIterators: back failed for empty tree\n" );
    }
    
    const int nMax = 1000;
    
    for( int i=0; i<nMax; ++i )
    {
        tree.Insert( (int)i );
        if( i == 1 )
        {
            const CNearTree<int>::iterator itSingle = tree.back( );
            if( (*itSingle) != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testIterators: size failed for size=2 tree\n" );
            }
        }
    }
    
    {
        CNearTree<int>::iterator it;
        it = --tree.end( );
        if( *it != nMax-1 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test01 failed\n" );
        }
        
        it = tree.begin( );
        if( *it != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test02 failed\n" );
        }
        
        it = it+1;
        if( *it != 1 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test03 failed\n" );
        }
        
        it = it-1;
        if( *it != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test04 failed\n" );
        }
        
        const int i0 = *it;
        if( i0 != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test05 failed\n" );
        }
        
        ++it;
        if( *it != 1 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test06 failed\n" );
        }
        
        --it;
        if( *it != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test07 failed\n" );
        }
        
        const CNearTree<int>::iterator itPlus = it++;
        if( *it != 1 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test08 failed\n" );
        }
        
        it--;
        if( i0 != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test09 failed\n" );
        }
        
        ++it;
        if( *it != 1 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: test10 failed\n" );
        }
    }
    
    {
        CNearTree<int>::iterator it;
        
        const int searchValue = 14;
        const CNearTree<int>::iterator itEnd = tree.end( );
        it = tree.NearestNeighbor( 0.1, searchValue );
        if( it == itEnd )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: NearestNeighbor failed with iterator\n" );
        }
        else if( *it != searchValue )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: NearestNeighbor found wrong value with iterator, %d\n", *it );
        }
        
        it = tree.FarthestNeighbor( searchValue );
        if( it == tree.end( ) )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators:FarthestNeighbor failed with iterator\n" );
        }
        else if( *it != nMax-1 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: FarthestNeighbor found wrong value with iterator, %d\n", *it );
        }
    }
    
    {
        vec17 v;
        CNearTree<vec17> nt;
        nt.Insert( v );
        CNearTree<vec17>::iterator itv = nt.begin( );
        const int n = itv->dim;
        if( n != 17 )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: operator-> got wrong value for dim, %d\n", n );
        }
        const double d = itv->pd[0];
        if( d != v.pd[0] )
        {
            ++g_errorCount;
            fprintf(stdout, "testIterators: operator-> got wrong value for pd[0], %g\n", d );
        }
    }
}