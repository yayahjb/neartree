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

#include <algorithm>
#include <utility>

#include <float.h>
#include <limits.h>
#include <cmath>
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
void testOutSphere( void );
void testDelayedInsertion( void );
void testIterators( void );
void testFindInAnnulus( void );
void testMisc( void );

void testRandomTree1( const int n );
void testRandomTree2( const int n );
void testBigVector( void );
void testIntegerReturn( void );
void testBigIntVec( void );
void testSTLContainerInput( void );

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
    fprintf( stdout, "testEmptyTree\n" );
    
    /* test the interface with trees with varying content, one entry and several */
    for( int i=1; i<10; ++i )
    {
        testLinearTree( i );
    fprintf( stdout, "testLinearTree %d\n", i );
    }
    
    testFindFirstObject( );
    fprintf( stdout, "testFindFirstObject\n" );
    testFindLastObject( );
    fprintf( stdout, "testFindLastObject\n" );
    testFindInSphereFromBottom( );
    fprintf( stdout, "testFindInSphereFromBottom\n" );
    testFindInSphereFromTop( );
    fprintf( stdout, "testFindInSphereFromTop\n" );
    testOutSphere( );
    fprintf( stdout, "testOutSphere\n" );
    testFindInAnnulus( );
    fprintf( stdout, "testFindInAnnulus\n" );
    testRandomTree1( 10000 );
    fprintf( stdout, "testRandomTree1\n" );
    testRandomTree2( 1000 );
    fprintf( stdout, "testRandomTree2\n" );
    testBigVector( );
    fprintf( stdout, "testBigVector\n" );
    testBigIntVec( );
    fprintf( stdout, "testBigIntVec\n" );
    testBackwardForward( );
    fprintf( stdout, "testBackwardForward\n" );
    testDelayedInsertion( );
    fprintf( stdout, "testDelayedInsertion\n" );
    testIterators( );
    fprintf( stdout, "testIterators\n" );
    testIntegerReturn( );
    fprintf( stdout, "testIntegerReturn\n" );
    testMisc( );
    fprintf( stdout, "testMisc\n" );
    testSTLContainerInput( );
    fprintf( stdout, "testSTLContainerInput\n" );
    
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
    

    CNearTree<int> ntReturn;
    lFoundPointsInSphere = tree.FindInSphere( 1000.0, ntReturn, 1 );
    if( lFoundPointsInSphere != 0 )
    {
        ++g_errorCount;
        fprintf(stdout, "testEmptyTree: .FindInSphere incorrect for empty tree\n" );
    }

    tree.clear( );
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
        tree.insert( i );
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
    CNearTree<int> ntReturn;
    const long lFoundPointsInSphere = tree.FindInSphere( -100.0, ntReturn, 1 );
    if( lFoundPointsInSphere != 0 )
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
        ntReturn.clear( );
        const long found = tree.FindInSphere( 0.1, ntReturn, 1 );
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
    if( tree.FindInSphere( (double)(10*n), ntReturn, 0 ) != n )
    {
        ++g_errorCount;
        fprintf(stdout, "FindInSphere did not find all the points, found %ld\n", tree.FindInSphere( (double)(10*n), ntReturn, 0 ) );
    }

    /* test operator[ ] */
    for( size_t i=1; i<=(size_t)n; ++i )
    {
       if( i != (size_t)tree[i-1] )
    {
        ++g_errorCount;
        fprintf(stdout, "testLinearTree: found tree[i]!=i for tree[i-1]=%d, i=%ld\n", tree[i-1], (long)i );
       }
    }
    
    v.clear( );

    tree.clear( );
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
            tree.insert( f );
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
        CNearTree<float> sphereReturn;
        const long lFound = tree.FindInSphere( 100.0, sphereReturn, 1.0 );
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
            tree.insert( f );
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
        CNearTree<double> sphereReturn;
        const long lFound = tree.FindInSphere( 100.0, sphereReturn, 1.0 );
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
            tree.insert( f );
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
            tree.insert( f );
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
            tree.insert( f );
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
            tree.insert( f );
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
        CNearTree<double> tree;
        
        double f = 1.0;
        /* generate an unbalanced tree*/
        while( f > 0.0 && f >= DBL_MIN)
        {
            f /= 2.0;
        }
        
        f = DBL_MIN;
        while( f < 1.01 )
        {
            tree.insert( f );
            f *= 2.0;
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
        tree.insert( (double)i );
    }
    
    std::vector<double> v;
    /* generate an unbalanced tree*/
    for( int i=1; i<=nmax; ++i )
    {
        v.clear( );
        const double radius = 0.05 + (double)i;
        CNearTree<double> sphereReturn;
        const long lReturned = tree.FindInSphere( radius, sphereReturn, 0.9 );
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
        tree.insert( (double)i );
    }
    
    /* generate an unbalanced tree*/
    for( int i=1; i<=nmax; ++i )
    {
        const double radius = 0.05 + (double)i;
        CNearTree<double> sphereReturn;
        const long lReturned = tree.FindInSphere( radius, sphereReturn, 0.1 + (double)nmax );
        if( lReturned != (long)i )
        {
            ++g_errorCount;
            fprintf(stdout, "FindInSphere failed in testFindInSphereFromTop for i=%d\n", i );
        }
    }
    
    std::vector<double> v;
    /* generate an unbalanced tree*/
    for( int i=1; i<=nmax; ++i )
    {
        v.clear( );
        const double radius = 0.05 + (double)i;
        CNearTree<double> sphereReturn;
        const long lReturned = tree.FindInSphere( radius, v, 0.1 + (double)nmax );
        if( lReturned != (long)i )
        {
            ++g_errorCount;
            fprintf(stdout, "FindInSphere failed for vector in testFindInSphereFromTop for i=%d\n", i );
        }
    }
}

/*=======================================================================*/
void testOutSphere( void )
{
    const int nmax = 100;
    CNearTree<double> tree;
    
    /* generate an unbalanced tree*/
    for( int i=1; i<=nmax; ++i )
    {
        tree.insert( (double)i );
    }
    
    // now test the NearTree return
    for( int i=1; i<=nmax; ++i )
    {
        CNearTree<double> sphereReturn;
        const double radius = 0.05 + (double)i;
        const long lReturned = tree.FindOutSphere( radius, sphereReturn, 0.1 + (double)nmax );
        if( lReturned != 100-(long)i )
        {
            ++g_errorCount;
            fprintf(stdout, "testOutSphere: FindOutSphere failed for i=%d\n", i );
        }
    }
    
    // now test the vector return
    for( int i=1; i<=nmax; ++i )
    {
        std::vector<double> v;
        v.clear( );
        const double radius = 0.05 + (double)i;
        const long lReturned = tree.FindOutSphere( radius, v, 0.1 + (double)nmax );
        if( lReturned != 100-(long)i )
        {
            ++g_errorCount;
            fprintf(stdout, "testOutSphere: FindOutSphere failed for for vector for i=%d\n", i );
        }
    }
}

/*
 Test CNearTree with a bunch of random numbers (integers).
 */
/*=======================================================================*/
void testRandomTree1( const int nRequestedRandoms )
{
    const int n = nRequestedRandoms;
    CNearTree<int> tree;
    int nmax = INT_MIN;
    int nmin = INT_MAX;
    
    /* Build the tree with n random numbers. Remember the largest and smallest values. */
    for( int i=0; i<n; ++i )
    {
        const int next = random( )%MYRAND_MAX;
        tree.insert( next );
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
            fprintf(stdout, "testRandomTree1: NearestNeighbor failed to find anything\n" );
        }
        else if( closest != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree1: NearestNeighbor failed for min\n" );
        }
        
        int farthest=INT_MIN;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MIN/2 );
        if( !bFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree1: FarthestNeighbor failed to find anything\n" );
        }
        else if( farthest != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree1: FarthestNeighbor failed in testRandomTree1 for min\n" );
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from above)*/
        int closest=INT_MIN;
        const bool bNear = tree.NearestNeighbor( (double)LONG_MAX, closest, INT_MAX/2 );
        if( ! bNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree1: NearestNeighbor failed to find anything for max\n" );
        }
        else if( closest != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree1: NearestNeighbor failed for max\n" );
        }
        
        int farthest;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MAX/2 );
        if( !bFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree1: FarthestNeighbor failed to find anything for min\n" );
        }
        else if( farthest != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree1: FarthestNeighbor failed for min\n" );
        }
    }
    
    {
        /*verify that for very large radius, every point is detected (from below)*/
        std::vector<int> v;
        const double radius = DBL_MAX;
        CNearTree<int> sphereReturn;
        const long lReturn = tree.FindInSphere( radius, sphereReturn, INT_MIN/2 );
        if( lReturn != n )
        {
            ++g_errorCount;
            fprintf(stdout, "FindInSphere failed in testRandomTree1, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*verify that we find NO points if we are below the lowest and with too small radius*/
        const double radius = .5;
        CNearTree<int> sphereReturn;
        const long lReturn = tree.FindInSphere( radius, sphereReturn, nmin-1 );
        if( lReturn != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "FindInSphere failed in testRandomTree1 found points incorrectly, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*test that the number of found points in a sphere is non-deceasing with increasing radius*/
        long lReturn = 0;
        int lastFoundCount = 0;
        
        double radius = 0.00001; /* start with a very small radius (remember these are int's) */
        while( radius < (double)(5*(nmax-nmin)) )
        {
            CNearTree<int> sphereReturn;
            lReturn = tree.FindInSphere( radius, sphereReturn, nmin-1 );
            if( lReturn < lastFoundCount )
            {
                ++g_errorCount;
                fprintf(stdout, "FindInSphere in testRandomTree1 found DECREASING count on increasing radius for radius=%f\n", radius );
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
            fprintf(stdout, "FindInSphere in testRandomTree1 did not find all the points\n" );
        }
    }
} // end testRandomTree1
    
/*
 Test CNearTree with a bunch of random numbers (integers).
 */
/*=======================================================================*/
void testRandomTree2( const int nRequestedRandoms )
{
    const int n = nRequestedRandoms;
    CNearTree<int, double, -11> tree;
    int nmax = INT_MIN;
    int nmin = INT_MAX;
    
    /* Build the tree with n random numbers. Remember the largest and smallest values. */
    for( int i=0; i<n; ++i )
    {
        const int next = random( )%MYRAND_MAX;
        tree.insert( next );
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
            fprintf(stdout, "testRandomTree2: NearestNeighbor failed to find anything\n" );
        }
        else if( closest != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree2: NearestNeighbor failed for min\n" );
        }
        
        int farthest=INT_MIN;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MIN/2 );
        if( !bFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree2: FarthestNeighbor failed to find anything\n" );
        }
        else if( farthest != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree2: FarthestNeighbor failed in testRandomTree2 for min\n" );
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from above)*/
        int closest=INT_MIN;
        const bool bNear = tree.NearestNeighbor( (double)LONG_MAX, closest, INT_MAX/2 );
        if( ! bNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree2: NearestNeighbor failed to find anything for max\n" );
        }
        else if( closest != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree2: NearestNeighbor failed for max\n" );
        }
        
        int farthest;
        const bool bFar = tree.FarthestNeighbor( farthest, INT_MAX/2 );
        if( !bFar )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree2: FarthestNeighbor failed to find anything for min\n" );
        }
        else if( farthest != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "testRandomTree2: FarthestNeighbor failed for min\n" );
        }
    }
    
    {
        /*verify that for very large radius, every point is detected (from below)*/
        std::vector<int> v;
        const double radius = DBL_MAX;
        CNearTree<int,double,-11> sphereReturn;
        const long lReturn = tree.FindInSphere( radius, sphereReturn, INT_MIN/2 );
        if( lReturn != n )
        {
            ++g_errorCount;
            fprintf(stdout, "FindInSphere failed in testRandomTree2, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*verify that we find NO points if we are below the lowest and with too small radius*/
        const double radius = .5;
        CNearTree<int,double,-11> sphereReturn;
        const long lReturn = tree.FindInSphere( radius, sphereReturn, nmin-1 );
        if( lReturn != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "FindInSphere failed in testRandomTree2 found points incorrectly, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*test that the number of found points in a sphere is non-deceasing with increasing radius*/
        long lReturn = 0;
        int lastFoundCount = 0;
        
        double radius = 0.00001; /* start with a very small radius (remember these are int's) */
        while( radius < (double)(5*(nmax-nmin)) )
        {
            CNearTree<int,double,-11> sphereReturn;
            lReturn = tree.FindInSphere( radius, sphereReturn, nmin-1 );
            if( lReturn < lastFoundCount )
            {
                ++g_errorCount;
                fprintf(stdout, "FindInSphere in testRandomTree2 found DECREASING count on increasing radius for radius=%f\n", radius );
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
            fprintf(stdout, "FindInSphere in testRandomTree2 did not find all the points\n" );
        }
    }
    }
    
/*=======================================================================*/
/* make a 17-dimension vector class for testing */
        class vec17;
class vec17
    {
    public:
        double pd[17];
        int dim;
        double length;  // just for a signature for debugging
        vec17( ) :
            dim(17), length(0)
        {
            for( int i=0; i<dim; ++i )
            {
                pd[i] = (double)(random( )%MYRAND_MAX);
            }
            length = Norm( );
        }
        explicit vec17( const double d ):
            dim(17), length(0)
        {
            for( int i=0; i<dim; ++i )
            {
                pd[i] = d;
            }
            length = Norm( );
        }
        explicit vec17( const int n ):
            dim(17), length(0)
        {
            for( int i=0; i<dim; ++i )
            {
                pd[i] = double(n);
            }
            length = Norm( );
        }
        ~vec17( void )
        {
        } 
        
        vec17 operator-( const vec17& v ) const /* USERS: be sure to make both const */
        {
            vec17 vtemp;
            for( int i=0; i<dim; ++i )
            {
                vtemp.pd[i] = this->pd[i] - v.pd[i];
            }
            vtemp.length = Norm(  );
            return( vtemp );
        }
        double Norm( void ) const
        {
            long double dtemp = 0.0;
            for( int i=0; i<dim; ++i )
            {
                dtemp += pd[i]*pd[i];  //  L2 measure here
            }
            return( double(::sqrtl( dtemp )) );
        }
    };


/*=======================================================================*/
void testBigVector(  )
{
    const int vectorsize = 1000;   
    CNearTree<vec17> tree;
    vec17 vAll[vectorsize]; /* keep a list of all of the input so we can find particular entries */
    double rmax = -DBL_MAX;
    double rmin =  DBL_MAX;
    vec17 v17min; /* to be the point nearest to the origin */
    vec17 v17max; /* to be the point farthest from the origin */
    
    /* All of the coordinate values will be in the range 0-MYRAND_MAX. In other words,
     all of the data points will be within a 17-D cube that has a corner at
     the origin of the space.
     */
    for( int i=0; i<vectorsize; ++i )
    {
        vec17 v;
        if( v.Norm( ) < rmin )
        {
            rmin = v.Norm( );
            v17min = v;
        }
        
        if( v.Norm( ) > rmax )
        {
            rmax = v.Norm( );
            v17max = v;
        }
        
        tree.insert( v );
        vAll[i] = v;
    }
    
    {
        /* Find the point farthest from the point that was nearest the origin. */
        vec17 vFarthest;
        tree.FarthestNeighbor( vFarthest, v17min );
        
        /* Brute force search for the farthest */
        vec17 vSearch;
        double dmax = -DBL_MAX;
        for( int i=0; i<vectorsize; ++i )
        {
            if( ( vAll[i] - v17min ).Norm( ) > dmax )
            {
                dmax = ( vAll[i] - v17min ).Norm( );
                vSearch = vAll[i];
            }
        }
        
        double distdiff = ( (vSearch-v17min) - (vFarthest-v17min) ).Norm( );
        distdiff = (distdiff<0)?-distdiff:distdiff;
        if( distdiff > DBL_MIN )
        {
            ++g_errorCount;
            fprintf(stdout, "in testBigVector, apparently FarthestNeighbor has failed\n" );
        }
    }
    
    {
        /* somewhere in the middle, find a point and its nearest neighbor */
        /* make sure that each includes the other in sphere search */
        
        const vec17 vBox17Center( (double)(MYRAND_MAX/2) );
        vec17 vNearCenter;
        vec17 vCloseToNearCenter;
        tree.NearestNeighbor( double(MYRAND_MAX/2)*sqrt(17.), vNearCenter, vBox17Center );
        CNearTree<vec17> sphereReturn;
        unsigned long iFoundNearCenter = (unsigned long)tree.FindInSphere( double(MYRAND_MAX/2)*sqrt(17.)/2., sphereReturn, vNearCenter );
        
        /* Brute force search for the point closest to the point closest to the center */
        double dmin = DBL_MAX;
        for( unsigned long i=0; i<iFoundNearCenter; ++i )
        {
            if( (vNearCenter - sphereReturn[i]).Norm( ) > DBL_MIN && ( vNearCenter-sphereReturn[i] ).Norm( ) < dmin )
            {
                dmin = ( vNearCenter-sphereReturn[i] ).Norm( );
                vCloseToNearCenter = sphereReturn[i];
            }
        }

        {
            //const double radius = ( vCloseToNearCenter - vNearCenter ).Norm( );
            const double radius = MYRAND_MAX*::sqrt(17.0);
            unsigned long iSphereFoundNearCenter = (unsigned long)tree.FindInSphere( radius, sphereReturn, vNearCenter );

            double searchRadius = radius/2;
            double delta        = searchRadius;
            int count = 0;
            while( iSphereFoundNearCenter != 2 && count < 100 )
            {
                iSphereFoundNearCenter = (unsigned long)tree.FindInSphere( searchRadius, sphereReturn, vNearCenter );
                if( iSphereFoundNearCenter > 2 )
                {
                   searchRadius = searchRadius - delta/2;
                }
                else if( iSphereFoundNearCenter < 2 )
                {
                   searchRadius = searchRadius + delta/2;
                }
                delta /= 2;
                ++count;
            }  // end while

            if( iSphereFoundNearCenter != 2 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigVector: FindInSphere failed to find only 2 points\n" );
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
            sphereReturn.clear( );
            const long iFound = tree.FindInSphere( 0.0, sphereReturn, vNearCenter );
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
            sphereReturn.clear( );
            const long iFound = tree.FindInSphere( (vCloseToNearCenter-vNearCenter).Norm( ), sphereReturn, vNearCenter );
            if( iFound < 2 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigVector: FindInSphere found only 1 point\n" );
            }
        }
        
        {
            /* Using small radius, check that only one point is found when a point is searched
             with FindInSphere */
            sphereReturn.clear( );
            const long iFound = tree.FindInSphere( (vCloseToNearCenter-vNearCenter).Norm( )*0.9, sphereReturn, vNearCenter );
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigVector: FindInSphere found no points using %f radius\n", (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigVector: FindInSphere found %ld points using %f radius\n", iFound, (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
        }
    }
}  // testBigVector
    
/*=======================================================================*/
void testBackwardForward( void )
{
    CNearTree<double> tree;
    const int nMax = 1000;
    
    for( int i=0; i<nMax; ++i )
    {
        tree.insert( (double)i );
        tree.insert( (double)(nMax-i) );
        tree.insert( (double)i + 0.25 );
        tree.insert( (double)(nMax-i) + 0.75 );
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
                tree.insert( (double)i );
            }
            else
            {
                tree.DelayedInsert( (double)i );
            }
        }
        
        CNearTree<double> sphereReturn;
        const double radius = 1000.0;
        const long lReturned = tree.FindInSphere( radius, sphereReturn, 0.9 );
        
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
                tree.insert( (double)i );
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
                tree.insert( (double)i );
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
        tree.insert( (int)i );
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
       // test InSphere
       CNearTree<int> sphereTree;
       const double radius = 5.0;
       const long sizeReturned = tree.FindInSphere( radius, sphereTree, 50 );
       CNearTree<int>::iterator itSphere = sphereTree.begin( );
       const size_t sphereSizeIn = sphereTree.size( );

       if( (long)sphereSizeIn != sizeReturned )
       {
            ++g_errorCount;
            fprintf(stdout, "testIterators: FindInSphere and size found different values, %ld and %ld\n", sizeReturned, (long)sphereSizeIn );
       }
       const long sizeReturned2In = tree.FindInSphere( radius, sphereTree, 50 );

       if( sizeReturned != sizeReturned2In )
       {
            ++g_errorCount;
            fprintf(stdout, "testIterators: 2nd FindInSphere found different value, %ld and %ld\n", sizeReturned, sizeReturned2In );
       }

       // test OutSphere
       const long sizeReturnedOut = tree.FindOutSphere( radius, sphereTree, 50 );
       itSphere = sphereTree.begin( );
       const size_t sphereSizeOut = sphereTree.size( );

       if( (long)sphereSizeOut != sizeReturnedOut )
       {
            ++g_errorCount;
            fprintf(stdout, "testIterators: FindOutSphere and size found different values, %ld and %ld\n", sizeReturned, (long)sphereSizeOut );
       }
       const long sizeReturned2Out = tree.FindOutSphere( radius, sphereTree, 50 );

       if( sizeReturned2Out != sizeReturnedOut )
       {
            ++g_errorCount;
            fprintf(stdout, "testIterators: 2nd FindOutSphere found different value, %ld and %ld\n", sizeReturned, sizeReturned2Out );
       }

       if( sizeReturned2In+sizeReturned2Out != nMax+2 )
       {
           // because there are 2 ints on the boundary, the count should be nMax+2
            ++g_errorCount;
            fprintf(stdout, "testIterators: total count from InSphere and OutSphere is wrong, %ld, should be%d\n",
                sizeReturned2Out+sizeReturned2Out, nMax+2 );
       }

       sphereTree.insert( 999999 );
       CNearTree<int>::iterator it = sphereTree.end( );
       --it;

       if( *(it) != 999999 )
       {
            ++g_errorCount;
            fprintf(stdout, "testIterators: Error after incrementing iterator at end(1)\n" );
       }

       it = --sphereTree.end( );
       if( *(it) != 999999 )
       {
            ++g_errorCount;
            fprintf(stdout, "testIterators: Error after incrementing iterator at end(2)\n" );
       }

    }
    
    {
        vec17 v;
        CNearTree<vec17> sphereReturn;
        sphereReturn.insert( v );
        CNearTree<vec17>::iterator itv = sphereReturn.begin( );
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

/*=======================================================================*/
    class v2
    {
        int m,n;
    public:
        v2( const int mm, const int nn ) :
            m(mm), n(nn)
        {
        }

        v2 operator- ( const v2& v ) const
        {
            return( v2( m-v.m, n-v.n ) );
        }
        double Norm( void ) const
        {
            return( ::sqrt( (double)( n*n + m*m ) ) );  //  L2 measure
        }

    }; // end of v2



/*=======================================================================*/
void testFindInAnnulus( void )
{
    {
        CNearTree<int> tree;
        CNearTree<int>::iterator itEmpty = tree.back( ); // make sure back works

        const int nMax = 1000;

        for( int i=0; i<nMax; ++i )
        {
            tree.insert( (int)i );
        }
        const double r1 = 100.1;
        const double r2 = 299.9;

        {

            CNearTree<int> annulusTree;
            const long lInAnnulus = tree.FindInAnnulus( r1, r2, annulusTree, 0 );

            if( lInAnnulus != (299-101+1) )
            {
                ++g_errorCount;
                fprintf(stdout, "testFindInAnnulus: wrong number of objects found\n" );
            }

            CNearTree<int>::iterator itNear = annulusTree.NearestNeighbor( 1000.0, 0 );

            if( *itNear != 101 )
            {
                ++g_errorCount;
                fprintf(stdout, "testFindInAnnulus: lowest value (%d) is incorrect\n", *itNear );
            }

            CNearTree<int>::iterator itFar = annulusTree.FarthestNeighbor( 0 );
            if( *itFar != 299 )
            {
                ++g_errorCount;
                fprintf(stdout, "testFindInAnnulus: highest value (%d) is incorrect\n", *itFar );
            }
        }

        {
            std::vector<int> annulusvector;
            const long lInAnnulus = tree.FindInAnnulus( r1, r2, annulusvector, 0 );

            if( lInAnnulus != (299-101+1) )
            {
                ++g_errorCount;
                fprintf(stdout, "testFindInAnnulus: for vector, wrong number of objects found\n" );
            }
        }
    }
    {
        // test annulus in the case of more than 1-D using CNearTree return
        CNearTree< v2 > tree;

        for( int i=0; i<30; ++i )
        {
            for( int j=0; j<30; ++j )
            {
                tree.insert( v2(i,j) );
            }
        }

        const double r1 = 4.0;
        const double r2 = 8.0;

        const v2 probe( -1, -1 );

        // use FindInAnnulus to find the points in the annulus
        CNearTree< v2 > annulusResult;
        const long searchFoundInAnnulus = tree.FindInAnnulus( r2, r1, annulusResult, probe );

        // now count all the points in the annulus by brute force

        int count = 0;
        const double radius = 0.001; // specify a very small radius so integer grid will only detect single points
        for( int i=0; i<30; ++i )
        {
            for( int j=0; j<30; ++j )
            {
                const v2 probe2( i, j );
                v2 vReturn( 0, 0 ); // dummy coords for constructor
                if( annulusResult.NearestNeighbor( radius, vReturn, probe2 ) )
                {
                    ++count;
                }
            }
        }

        if( searchFoundInAnnulus != count )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindInAnnulus: wrong number of points found in annulus by NearTree, %ld, direct count=%d\n", searchFoundInAnnulus, count );
        }
    }

}

/*=======================================================================*/
void testMisc( void )
{
    CNearTree<int, double, -3> nt;
    for( int i=0; i<20; ++i )
    {
        nt.insert( i );
    }

    const std::vector<int> v = nt;
    if( v.size( ) != nt.size( ) )
    {
            ++g_errorCount;
            fprintf(stdout, "testMisc: wrong number of points %ld in vector\n", (long)v.size( ) );
    }

    const size_t indexToRetrieve = 3;
    if( nt.at( indexToRetrieve ) != (int)indexToRetrieve )
    {
            ++g_errorCount;
            fprintf(stdout, "testMisc: 'at' returned wrong value, %d\n",nt.at( 3 ) );
    }
}

/*=======================================================================*/
void testIntegerReturn( void )
{
    {
        CNearTree<int, int> nt;

        const int nmax = 20;
        for( int i=0; i<nmax; ++i )
        {
            nt.insert( i );
        }

        int closest=-1;
        if( ! nt.NearestNeighbor( 1000, closest, 100 ) )
        {
                ++g_errorCount;
                fprintf(stdout, "testIntegerReturn: NearestNeighbor failed\n" );
        }
        else if( closest != nmax-1 )
        {
                ++g_errorCount;
                fprintf(stdout, "testIntegerReturn: NearestNeighbor return wrong value, %d\n", closest );
        }
    }

    {
        CNearTree<int, long> nt;

        const int nmax = 20;
        for( int i=-nmax; i<nmax; ++i )
        {
            nt.insert( i );
        }

        int closest=-nmax-1;
        if( ! nt.NearestNeighbor( 1000, closest, 100 ) )
        {
                ++g_errorCount;
                fprintf(stdout, "testIntegerReturn: NearestNeighbor failed\n" );
        }
        else if( closest != nmax-1 )
        {
                ++g_errorCount;
                fprintf(stdout, "testIntegerReturn: NearestNeighbor return wrong value, %d\n", closest );
        }
    }
}

    
/*=======================================================================*/
/* make a 17-dimension all-integer vector class for testing */
class intVec17
{
public:
    int pd[17];
    int dim;
    double length;  // just for a signature for debugging
    intVec17( ) :
        dim(17)
    {
        for( int i=0; i<dim; ++i )
        {
            pd[i] = random( )%MYRAND_MAX;
        }
        length = this->Norm( );
    }
    explicit intVec17( const int d )
    {
        dim = 17;
        for( int i=0; i<dim; ++i )
        {
            pd[i] = d;
        }
        length = this->Norm( );
    }
    ~intVec17( void )
    {
    }
    int Norm( void ) const
    {
        long dtemp = 0;
        for( int i=0; i<dim; ++i )
        {
            dtemp += abs(pd[i]);  // L1 measure here
        }
        return( dtemp );
    } 
    
    intVec17 operator-( const intVec17& v ) const /* USERS: be sure to make both const */
    {
        intVec17 vtemp;
        for( int i=0; i<dim; ++i )
        {
            vtemp.pd[i] = pd[i] -v.pd[i];
        }
        vtemp.length = vtemp.Norm( );
        return( vtemp );
    }
}; // intVec17

/*=======================================================================*/

void testBigIntVec( void )
{
    const int vectorsize = 1000;   
    CNearTree<intVec17, int, -2> tree;
    intVec17 vAll[vectorsize]; /* keep a list of all of the input so we can find particular entries */
    double rmax = -DBL_MAX;
    double rmin =  DBL_MAX;
    intVec17 v17min; /* to be the point nearest to the origin */
    intVec17 v17max; /* to be the point farthest from the origin */
    
    /* All of the coordinate values will be in the range 0-MYRAND_MAX. In other words,
     all of the data points will be within a 17-D cube that has a corner at
     the origin of the space.
     */
    for( int i=0; i<vectorsize; ++i )
    {
        intVec17 v;
        if( v.Norm( ) < rmin )
        {
            rmin = v.Norm( );
            v17min = v;
        }
        
        if( v.Norm( ) > rmax )
        {
            rmax = v.Norm( );
            v17max = v;
        }
        
        tree.insert( v );
        vAll[i] = v;
    }
    
    {
        /* Find the point farthest from the point that was nearest the origin. */
        intVec17 vFarthest;
        tree.FarthestNeighbor( vFarthest, v17min );
        
        /* Brute force search for the farthest */
        intVec17 vSearch;
        double dmax = -DBL_MAX;
        for( int i=0; i<vectorsize; ++i )
        {
            if( ( vAll[i] - v17min ).Norm( ) > dmax )
            {
                dmax = ( vAll[i] - v17min ).Norm( );
                vSearch = vAll[i];
            }
        }
        
        int distdiff = (int)((vSearch-v17min).Norm( )-(vFarthest-v17min).Norm( ));
        distdiff = (distdiff<0)?-distdiff:distdiff;
        if( distdiff > DBL_MIN )
        {
            ++g_errorCount;
            fprintf(stdout, "in testBigIntVector, apparently FarthestNeighbor has failed\n" );
        }
    }
    
    {
        /* somewhere in the middle, find a point and its nearest neighbor */
        /* make sure that each includes the other in sphere search */
        
        const intVec17 vBox17Center( MYRAND_MAX/2 );
        intVec17 vNearCenter;
        intVec17 vCloseToNearCenter;
        const bool bFoundNear = tree.NearestNeighbor( (int)MYRAND_MAX*17, vNearCenter, vBox17Center );

        if( ! bFoundNear )
        {
            ++g_errorCount;
            fprintf(stdout, "testBigIntVector: NearestNeighbor found no points\n" );
        }

        CNearTree<intVec17, int, -2> sphereReturn;
        // really, the only way to get the next section perfectly right is to do a 
        // binary search until exactly 2 points are found. (or do a linear search thru the points)
        const int radius = 2*(int)(MYRAND_MAX/2*sqrt(17.)); // hand-tuned value
        unsigned long iFoundNearCenter = (unsigned long)tree.FindInSphere( radius, sphereReturn, vNearCenter );

        int searchRadius = radius/2;
        int delta        = searchRadius;
        int count = 0;
        while( iFoundNearCenter != 2 && count < 100 )
        {
            iFoundNearCenter = (unsigned long)tree.FindInSphere( searchRadius, sphereReturn, vNearCenter );
            if( iFoundNearCenter > 2 )
            {
               searchRadius = searchRadius - delta/2;
            }
            else if( iFoundNearCenter < 2 )
            {
               searchRadius = searchRadius + delta/2;
            }
            delta /= 2;
            ++count;
        }  // end while

        if( iFoundNearCenter < 2 )
        {
            ++g_errorCount;
            fprintf(stdout, "testBigIntVector: FindInSphere found %d points, instead of 2\n", iFoundNearCenter );
        }
        
        /* Brute force search for the point closest to the point closest to the center */
        double dmin = DBL_MAX;
        for( unsigned long i=0; i<iFoundNearCenter; ++i )
        {
            if( (vNearCenter - sphereReturn[i]).Norm( )!=0 && ( vNearCenter-sphereReturn[i] ).Norm( ) < dmin )
            {
                dmin = ( vNearCenter-sphereReturn[i] ).Norm( );
                vCloseToNearCenter = sphereReturn[i];
            }
        }

        {
            sphereReturn.clear( );
            const int radius2 = (int)( vCloseToNearCenter - vNearCenter ).Norm( );
            const unsigned long iSphereFoundNearCenter = (unsigned long)tree.FindInSphere( radius2, sphereReturn, vNearCenter );
            if( iSphereFoundNearCenter != 2 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FindInSphere failed to find only 2 points\n" );
            }
        }
        
        if( dmin == DBL_MAX ) 
        {
            ++g_errorCount;
            fprintf(stdout, "testBigIntVector: apparently FindInSphere failed\n" );
        }
        
        {
            /* Using zero radius, check that only one point is found when a point is searched
             with FindInSphere */
            sphereReturn.clear( );
            const long iFound = tree.FindInSphere( 0, sphereReturn, vNearCenter );
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FindInSphere found no points using zero radius\n" );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FindInSphere found more than %ld points using zero radius\n", iFound );
        }
    }
        
        {
            /* Using minimal radius, check that at least 2 points are found when a point is searched
             with FindInSphere */
            sphereReturn.clear( );
            const long iFound = tree.FindInSphere( (vCloseToNearCenter-vNearCenter).Norm( ), sphereReturn, vNearCenter );
            if( iFound < 2 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FindInSphere found only 1 point\n" );
            }
        }
        
        {
            /* Using small radius, check that only one point is found when a point is searched
             with FindInSphere */
            sphereReturn.clear( );
            const long iFound = tree.FindInSphere( (int)((vCloseToNearCenter-vNearCenter).Norm( )*0.9), sphereReturn, vNearCenter );
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FindInSphere found no points using %f radius\n", (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FindInSphere found %ld points using %f radius\n", iFound, (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
        }
         
        {
            const intVec17 vOrigin( 0 );
            const intVec17 vTopCorner( MYRAND_MAX );
            tree.insert( vOrigin );
            tree.insert( vTopCorner );
            sphereReturn.clear( );
            intVec17 vReturn;
            const long iFound = tree.FarthestNeighbor( vReturn, vOrigin );
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FarthestNeighbor from origin found no points using %f radius\n", (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FarthestNeighbor from origin found %ld points using %f radius\n", iFound, (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
            else if( (vReturn - vTopCorner ).Norm( ) != 0 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FarthestNeighbor from origin found wrong top point\n" );
            }
        }
         
        {
            const intVec17 vOrigin( 0 );
            const intVec17 vTopCorner( MYRAND_MAX );
            tree.insert( vOrigin );
            tree.insert( vTopCorner );
            sphereReturn.clear( );
            intVec17 vReturn;
            const long iFound = tree.FarthestNeighbor( vReturn, vTopCorner );
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FarthestNeighbor from top found no points using %f radius\n", (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FarthestNeighbor from top found %ld points using %f radius\n", iFound, (vCloseToNearCenter-vNearCenter).Norm( )*0.9 );
            }
            else if( (vReturn - vOrigin ).Norm( ) != 0 )
            {
                ++g_errorCount;
                fprintf(stdout, "testBigIntVector: FarthestNeighbor from top found wrong top point\n" );
            }
        }
    }
} // testBigIntVec

/*=======================================================================*/

void testSTLContainerInput( void )
{
   std::vector<int> v;
   std::list<int>   l;
   std::set<int>    s;

   for( int i=0; i<10; ++i )
   {
      v.push_back( i );
      l.push_back( i );
      s.insert( i );
   }

   {  // test the constructors
      CNearTree<int> tv(v);
      CNearTree<int> tl(l);
      CNearTree<int> ts(s);

      if( tv.size() != 10 )
      {
          ++g_errorCount;
          fprintf(stdout, "testSTLContainerInput: vector constructor has wrong count\n" );
      }

      if( tl.size() != 10 )
      {
          ++g_errorCount;
          fprintf(stdout, "testSTLContainerInput: list constructor has wrong count\n" );
      }

      if( ts.size() != 10 )
      {
          ++g_errorCount;
          fprintf(stdout, "testSTLContainerInput: set constructor has wrong count\n" );
      }
   }

   {  // test the inserters
      CNearTree<int> tv;
      CNearTree<int> tl;
      CNearTree<int> ts;

      tv.insert( v );
      ts.insert( s );
      tl.insert( l );

      if( tv.size() != 10 )
      {
          ++g_errorCount;
          fprintf(stdout, "testSTLContainerInput: vector insertion has wrong count\n" );
      }

      if( tl.size() != 10 )
      {
          ++g_errorCount;
          fprintf(stdout, "testSTLContainerInput: list insertion has wrong count\n" );
      }

      if( ts.size() != 10 )
      {
          ++g_errorCount;
          fprintf(stdout, "testSTLContainerInput: set insertion has wrong count\n" );
      }
   }

} // testSTLContainerInput

