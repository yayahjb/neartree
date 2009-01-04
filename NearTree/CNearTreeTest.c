/*
 *  CNearTreeTest.c
 *  NearTree
 *
 *  Based on CNearTreeTest.cpp
 *  Copyright 2008 Larry Andrews.  All rights reserved
 *
 *  C Version created by Herbert J. Bernstein on 1/1/09
 *  with permission from Larry Andrews.
 *  Copyright 2008 Larry Andrews and Herbert J. Bernstein. 
 *  All rights reserved.
 *
 */

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE NearTree UNDER THE TERMS OF THE LGPL          *
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

/*
 This is a test harness for the C version neartree API CNearTree.c.
 
 */

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <CNearTree.h>

void testEmptyTree( void );
void testLinearTree( const int n );
void testFindFirstObject( void );
void testFindLastObject( void );
void testFindInSphereFromBottom( void );
void testFindInSphereFromTop( void );

void testRandomTree( const int n );
void testBigVector( );

long g_errorCount;

#define bool int    /* for older C compilers */

#define MYRAND_MAX  32767

/*=======================================================================*/
int main(int argc, char** argv)
{
    int i;
    
    /* test the interface with an empty tree */
    testEmptyTree( );
    
    /* test the C API with trees with varying content, one entry and several */
    for( i=1; i<10; ++i )
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
        fprintf(stdout, "No errors were detected while testing CNearTree\n" );
    }
    
    return g_errorCount;
}

/*
 For an empty tree of int's, test the public interface for CNearTree.
 */
/*=======================================================================*/
void testEmptyTree( void )
{
    CNearTreeHandle tree;
    CVectorHandle v;
    int FAR * close;
    int FAR * nFar;
    int probe[1];
    bool bTreeEmpty;
    bool bTreeHasNearest;
    bool bTreeHasFarthest;
    bool bInSphere;
    size_t lFoundPointsInSphere;
    
    if (CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER))
        fprintf(stdout,"CNearTreeTest: testEmptyTree: CNearTreeCreate failed.\n");
    
    bTreeEmpty = !CNearTreeZeroIfEmpty(tree);
    if( ! bTreeEmpty )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: CNearTreeZeroIfEmpty not zero for empty tree\n" );
    }
    
    probe[0] = 1;
    bTreeHasNearest =  !CNearTreeNearestNeighbor(tree, 0.0, (void FAR * FAR *) &close, NULL, probe); 
    if( bTreeHasNearest )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: NearestNeighbor incorrect for empty tree\n" );
    }
    
    probe[0] = 0;
    bTreeHasFarthest = !CNearTreeFarthestNeighbor(tree, (void FAR * FAR *) &nFar, NULL, probe); 
    if( bTreeHasFarthest )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: FarthestNeighbor incorrect for empty tree\n" );
    }
    
    probe[0] = 1;
    if (CVectorCreate(&v,sizeof(int FAR *),10))
        fprintf(stdout,"CNearTreeTest: testEmptyTree: CVectorCreate failed\n" );
    
    bInSphere = !CNearTreeFindInSphere( tree, 1000.0, v, NULL, probe, 1);
    lFoundPointsInSphere = CVectorSize(v);
    if( bInSphere || lFoundPointsInSphere != 0 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: FindInSphere incorrect for empty tree\n" );
    }
    
    if (CVectorFree(&v))
        fprintf(stdout, "CNearTreeTest: testEmptyTree: CNearVectorFree failed\n");
    
    if (CNearTreeFree(&tree))
        fprintf(stdout, "CNearTreeTest: testEmptyTree: CNearTreeFree failed\n");
}

/*
 Test CNearTree using a tree built only using int's. This is not a robust test
 in some ways because the tree is only composed of right branches and left leaves.
 Perform some general tests.
 */
/*=======================================================================*/
void testLinearTree( const int n )
{
    CNearTreeHandle tree;
    CVectorHandle v;
    bool bInsert, bClose, bFar, bInSphere, bResult;
    int FAR * closest;   
    int FAR * farthest;
    int probe[1];
    size_t lFoundPointsInSphere;
    int i;
    size_t found;
    long localErrorCount = 0;
    long localErrorMax = 0;
    
    if (CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER))
        fprintf(stdout,"CNearTreeTest: testLinearTree: CNearTreeCreate failed.\n");
    
    /* generate an unbalanced tree*/
    for( i=1; i<=n; ++i )
    {
        probe[0] = i;
        bInsert = !CNearTreeInsert(tree, probe, NULL);
        if (!bInsert)
            fprintf(stdout, "CNearTreeTest: testLinearTree: CNearTreeInsert failed.\n");
    }
    
    /*
     Search for the nearest value using a probe point that is larger than 
     the largest value that was input. The returned values should be the
     last value entered into the tree. 
     */
    probe[0] = 2*n;
    bClose = !CNearTreeNearestNeighbor(tree,22.,(void FAR * FAR *)&closest,NULL,probe);
    if( ! bClose || closest[0] != n )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testLinearTree: NearestNeighbor failed, got %d, should be %d\n", closest[0], n );
    }
    
    /*
     Search for the farthest value using a probe point that is larger than 
     the largest value that was input. The returned values should be the
     first value entered into the tree. 
     */
    
    probe[0] = 2*n;
    bFar = !CNearTreeFarthestNeighbor(tree,(void FAR * FAR *)&farthest,NULL,probe);
    if( ! bFar || farthest[0] != 1 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testLinearTree: FarthestNeighbor failed, got %d, should be %d\n", farthest[0], 1 );
    }
    
    /*
     Find all of the points in the tree using a negative radius, using the first
     input point as the probe point. Nothing should be found.
     */
    if (CVectorCreate(&v,sizeof(int FAR *),10))
        fprintf(stdout,"CNearTreeTest: testLinearTree: CVectorCreate failed\n" );
    
    probe[0] = 1;
    bInSphere = !CNearTreeFindInSphere( tree, -100., v, NULL, probe, 1);
    lFoundPointsInSphere = CVectorSize(v);
    
    if( bInSphere || lFoundPointsInSphere != 0 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testLinearTree: FindInSphere found points for negative radius\n" );
    }
    CVectorClear(v);
    
    /*
     Find all of the points in the tree using a small radius. Do separate
     searches for every point entered. In every case, only a single point
     should be found.
     */
    for( i=1; i<=n; ++i )
    {
        size_t found;
        
        probe[0] = i;
        bInSphere = !CNearTreeFindInSphere( tree, 0.1, v, NULL, probe, 1);
        found = CVectorSize(v);
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
        fprintf(stdout, "CNearTreeTest: testLinearTree: FindInSphere found too many points (as many as %ld) %ld times\n", localErrorMax, localErrorCount );
    }
    CVectorClear(v);
    
    /*
     Find all of the points in the tree that are within a large radius, using 
     the first input point as the probe point. All of the input points should 
     be found within the radius.
     */
    
    probe[0] = 0;
    bInSphere = !CNearTreeFindInSphere( tree,(double)(10*n), v, NULL, probe, 1);
    found = CVectorSize(v);
    if( found != n )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testLinearTree: FindInSphere did not find all the points, found %ld\n", found );
    }

    bResult = !CVectorFree(&v);
    if (!bResult)
        fprintf(stdout, "CNearTreeTest: testLinearTree: CVectorFree has failed\n" );
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
        fprintf(stdout, "CNearTreeTest: testLinearTree: CNearTreeFree has failed\n" );
    
}

/*
 Perform general tests using floating point numbers. For the C vesion
 only one set of tests is performed, for double.
 The values are computed starting from some initial value, and
 each succeeding value is one half of the previous until zero is
 computed (the zero is NOT inserted into the tree). The tree will consist
 of only right branches and left leaves.
 */
/*=======================================================================*/
void testFindFirstObject( void )
{
    double fFinal = DBL_MAX; /* just initialization */
    CNearTreeHandle tree;
    CVectorHandle v;
    bool bReturn, bResult;
    bool bReturnNear;
    bool bReturnFar;
    long count = 0;
    double f[1];
    double FAR * closest;
    double FAR * farthest;
    size_t lFound;
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CNearTreeCreate failed\n" );
    
    /* build the double tree starting with 1.0 */
    f[0] = 1.0;
    while( f[0] > 0.0 )
    {
        bReturn = !CNearTreeInsert(tree,f, NULL);
        fFinal = f[0];
        f[0] /= 2.0;
        ++count;
    }
    
    if( (bReturn = !CNearTreeZeroIfEmpty(tree) ) )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: incorrectly found empty tree for double\n" );
    }
    
    /*
     Search for the value closest to zero. It should be a very small, probably
     denormalized number.
     */
    f[0] = 0.0;
    bReturnNear = !CNearTreeNearestNeighbor(tree, 1.0e-10, (void FAR * FAR *)&closest, NULL, f);
    if( ! bReturnNear || closest[0] != fFinal )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: failed for double, got %f\n", closest[0] );
    }
    
    /*
     Search for the value farthest from a large number. It should be a
     very small, probably denormalized number.
     */
    f[0] = 100.0;
    bReturnFar = !CNearTreeFarthestNeighbor(tree, (void FAR * FAR *)&farthest, NULL,f);
    if( ! bReturnFar || farthest[0] != fFinal )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: failed for double, got %f\n", farthest[0] );
    }
    
    /*
     Determine if FindInSphere can find all of the input data.
     */
    f[0] = 1.0;
    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CVectorCreate failed\n" );
    bReturn = !CNearTreeFindInSphere( tree, 100.0, v, NULL, f, 1);
    lFound = CVectorSize(v);
    if( lFound != count )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: found wrong count for FindInSphere for float, should be%ld, got %ld\n", count, lFound );
    }
 
    bResult = !CVectorFree(&v);
    if (!bResult)
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CVectorFree has failed\n" );
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CNearTreeFree has failed\n" );
    
    
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
    double fFinal = DBL_MAX;
    CNearTreeHandle tree;
    bool bReturn, bResult;
    long count = 0;
    double f[1];
    double FAR * closest;
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeCreate failed\n" );
    
    f[0] = 1.0;
    /* generate an unbalanced tree*/
    while( f[0] > 2.0*sqrt(DBL_MIN) )
    {
        bReturn = !CNearTreeInsert(tree,f, NULL);
        fFinal = f[0];
        f[0] /= 2.0;
        ++count;
    }
    
    
    f[0] = 11.0;
    bReturn = !CNearTreeNearestNeighbor(tree, 100.0, (void FAR * FAR *)&closest, NULL, f);
    if( ! bReturn || closest[0] != 1.0 )
    {
        ++g_errorCount;
      fprintf(stdout, "CNearTreeTest: testFindLastObject: Near(1) failed for double, got %g\n", closest[0] );
    }
    
    f[0] = -11.0;
    bReturn = !CNearTreeNearestNeighbor(tree, 100.0, (void FAR * FAR *)&closest, NULL, f);
    if( ! bReturn || closest[0] != 1.0 )
    {
        ++g_errorCount;
      fprintf(stdout, "CNearTreeTest: testFindLastObject: Near(2) failed for double, got %g\n", closest[0] );
    }
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeFree has failed\n" );
    
}

/*=======================================================================*/
void testFindInSphereFromBottom( void )
{
    const int nmax = 100;
    CNearTreeHandle tree;
    CVectorHandle v;
    bool bReturn;
    double f[1];
    double radius;
    size_t lReturned;
    int i;
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromBottom: CNearTreeCreate failed\n" );

    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromBottom: CVectorCreate failed\n" );
    
    for( i=1; i<=nmax; ++i )
    {
        f[0] = (double)i;
        bReturn = !CNearTreeInsert(tree,f, NULL);
        if (!bReturn)
            fprintf(stdout, "CNearTreeTest: testFindInSphereFromBottom: CNearTreeInsert failed\n" );
    }
    
    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    /* generate an unbalanced tree*/
    for( i=1; i<=nmax; ++i )
    {
        bReturn = !CVectorClear(v);
        radius = 0.05 + (double)i;
        f[0] = 0.9;
        bReturn = !CNearTreeFindInSphere( tree, radius, v, NULL, f, 1);
        lReturned = CVectorSize(v);
        if( lReturned != (size_t)i )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: FindInSphere failed in testFindInSphereFromBottom for i=%d\n", i );
        }
    }
}

/*=======================================================================*/
void testFindInSphereFromTop( void )
{
    const int nmax = 100;
    CNearTreeHandle tree;
    CVectorHandle v;
    double f[1];
    int i;
    bool bReturn;
    double radius;
    size_t lReturned;
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CNearTreeCreate failed\n" );
    
    for( i=1; i<=nmax; ++i )
    {
        f[0] = (double)i;
        bReturn = !CNearTreeInsert(tree,f, NULL);
    }
   
    
    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CVectorCreate failed\n" );
    
    for( i=1; i<=nmax; ++i )
    {
        CVectorClear(v);
        f[0] = 0.1 + (double)nmax;
        radius = 0.05 + (double)i;
        bReturn = !CNearTreeFindInSphere( tree, radius, v, NULL, f, 1);
        lReturned = CVectorSize(v);
        if( lReturned != (size_t)i )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CNearTreeFindInSphere failed in testFindInSphereFromTop for i=%d\n", i );
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
    CNearTreeHandle tree;
    CVectorHandle v;
    int nmax = INT_MIN;
    int nmin = INT_MAX;
    int i;
    int next[1];
    int FAR * closest;
    int FAR * farthest;
    bool bReturn, bNear, bFar;
    double radius;
    size_t lReturn;
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER);
    if (!bReturn)
        fprintf(stdout, "CNearTreeTest: testRandomTree: CNearTreeCreate failed\n" );
    
    
    /* Build the tree with n random numbers. Remember the largest and smallest values. */
    for( i=0; i<n; ++i )
    {
        next[0] = rand( )%MYRAND_MAX;
        bReturn = !CNearTreeInsert(tree,next,NULL);
        if( next[0] > nmax ) nmax = next[0];
        if( next[0] < nmin ) nmin = next[0];
    }
    
    {
        /*verify that the correct extremal point is detected (from below)*/
        next[0] = INT_MIN/2;
        bNear = !CNearTreeNearestNeighbor(tree, (double)LONG_MAX, (void FAR * FAR *)&closest, NULL, next);
        if( ! bNear || closest[0] != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree: NearestNeighbor failed in testRandomTree for min\n" );
        }
        
        next[0] = INT_MIN/2;
        bFar = !CNearTreeFarthestNeighbor(tree, (void FAR * FAR *)&farthest, NULL, next);
        if( !bFar || farthest[0] != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree: FarthestNeighbor failed in testRandomTree for min\n" );
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from above)*/
        next[0] = INT_MAX/2;
        bNear = !CNearTreeNearestNeighbor(tree, (double)LONG_MAX, (void FAR * FAR *)&closest, NULL, next);
        if( ! bNear || closest[0] != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree: NearestNeighbor failed for max\n" );
        }
        
        next[0] = INT_MAX/2;
        bFar = !CNearTreeFarthestNeighbor(tree, (void FAR * FAR *)&farthest, NULL, next);
        if( !bFar || farthest[0] != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree: FarthestNeighbor failed for max\n" );
        }
    }
    
    {
        /*verify that for very large radius, every point is detected (from below)*/
        bReturn = !CVectorCreate(&v,sizeof(int FAR *), 10);
        if (!bReturn)
            fprintf(stdout, "CNearTreeTest: testRandomTree: CVectorCreate failed\n" );
        radius = DBL_MAX;
        next[0] = 1;
        bReturn = !CNearTreeFindInSphere( tree, radius, v, NULL, next, 1);
        lReturn = CVectorSize(v);
        if( lReturn != n )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree: FindInSphere failed, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*verify that we find NO points if we are below the lowest and with too small radius*/
        bReturn = !CVectorClear(v);
        radius = .5;
        next[0] = nmin-1;
        bReturn = !CNearTreeFindInSphere( tree, radius, v, NULL, next, 1);
        lReturn = CVectorSize(v);
        if( lReturn != 0 )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree: FindInSphere failed found points incorrectly, n=%d, lReturn=%ld\n", n, lReturn );
        }
    }
    
    {
        /*test that the number of found points in a sphere is non-deceasing with increasing radius*/
        int lastFoundCount = 0;
        int cycleCount = 0;
        bReturn = !CVectorClear(v);
       
        radius = 0.00001; /* start with a very small radius (remember these are int's) */
        while( radius < 5*(nmax-nmin) )
        {
            ++cycleCount;
            next[0] = nmin-1;
            bReturn = !CNearTreeFindInSphere( tree, radius, v, NULL, next, 1);
            lReturn = CVectorSize(v);
            if( lReturn < lastFoundCount )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testRandomTree: FindInSphere found DECREASING count on increasing radius for radius=%f\n", radius );
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

/*  Utility routines for 17D double vectors */

int LoadVec17(double vec[17]) {
    size_t i;
    for (i = 0; i < 17; i++) {
        vec[i] = (double)(rand( )%MYRAND_MAX);
    }
    return 0;
}

int ConstVec17(double vec[17], double d) {
    size_t i;
    for (i = 0; i < 17; i++) {
        vec[i] = d;
    }
    return 0;
}
double NormVec17(double vec[17]) {
    size_t i;
    double sum;
    sum = vec[0]*vec[0];
    for (i = 1; i < 17; i++) {
        sum += vec[i]*vec[i];
    }
    return sqrt(sum);
}

int CopyVec17(double dst[17], double src[17]) {
    size_t i;
    
    for (i = 0; i < 17; i++) {
      dst[i] = src[i];
    }
    return 0;
}



/*=======================================================================*/
void testBigVector(  )
{
    
    CNearTreeHandle tree;
    CVectorHandle vhand;
    double vAll[1000][17]; /* keep a list of all of the input so we can find particular entries */
    double v17min[17];     /* to be the point nearest to the origin */
    double v17max[17];     /* to be the point farthest from the origin */
    double v[17], vSearch[17], vCenter[17], vCloseToNearCenter[17];
    double FAR * vv;
    double FAR * vFarthest;
    double FAR * vNearCenter;
    double FAR * vExtreme;
    
    double rmax = -DBL_MAX;
    double rmin =  DBL_MAX;
    double dmin, dmax;
    
    size_t iFoundNearCenter, iFound, i;
    
    bool bResult;
    
    bResult = !CNearTreeCreate(&tree,17,CNEARTREE_TYPE_DOUBLE);
    if (!bResult)
        fprintf(stdout,"CNearTreeTest: testBigVector: CNearTreeCreate failed\n");
    
    /* All of the coordinate values will be in the range 0-MYRAND_MAX. In other words,
     all of the data points will be within a 17-D cube that has a corner at
     the origin of the space.
     */
    
    for( i=0; i<1000; ++i )
    {
        LoadVec17(v);
        if( NormVec17( v ) < rmin )
        {
            rmin = NormVec17( v );
            CopyVec17(v17min, v);
        }
        
        if( NormVec17( v )  > rmax )
        {
            rmax = NormVec17( v );
            CopyVec17(v17max, v);
        }
        
        bResult = !CNearTreeInsert(tree,v,NULL);
        if (!bResult)
            fprintf(stdout,"CNearTreeTest: testBigVector: CNearTreeInsert failed\n");
        CopyVec17(vAll[i],v);
    }
    
    {        
        /* Find the point farthest from the point that was nearest the origin. */
        bResult = !CNearTreeFarthestNeighbor(tree,(void FAR * FAR *)&vFarthest,NULL,v17min);
        if (!bResult)
            fprintf(stdout,"CNearTreeTest: testBigVector: CNearTreeFarthestNeighbor failed\n");
        
        /* Brute force search for the farthest */
        dmax = -DBL_MAX;
        for( i=0; i<1000; ++i )
        {
            if (sqrt(CNearTreeDistsq(vAll[i],v17min,17,CNEARTREE_TYPE_DOUBLE)) > dmax)
            {
                dmax = sqrt(CNearTreeDistsq(vAll[i],v17min,17,CNEARTREE_TYPE_DOUBLE));
                CopyVec17(vSearch, vAll[i]);
            }
        }
        
      if( sqrt(CNearTreeDistsq(vSearch,vFarthest,17,CNEARTREE_TYPE_DOUBLE)) > DBL_MIN )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: FarthestNeighbor has failed\n" );
        }
    }
    
    {
        /* somewhere in the middle, find a point and its nearest neighbor */
        /* make sure that each includes the other in sphere search */
        
        ConstVec17(vCenter, (double)(MYRAND_MAX/2) );
        CVectorCreate(&vhand,sizeof(double *),10);

      bResult = !CNearTreeNearestNeighbor( tree, 100000.0, (void FAR * FAR *)&vNearCenter, NULL, vCenter );
        if (!bResult)
          fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeNearestNeighbor has failed\n" );
        bResult = !CNearTreeFindInSphere(tree, 100.0, vhand, NULL, vNearCenter,1);
        if (!bResult)
            fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFindInSphere has failed\n" );
        iFoundNearCenter = CVectorSize(vhand);
        
        /* Brute force search for the point closest to the point closest to the center */
        dmin = DBL_MAX;
      for( i=0; i<1000; ++i )
      {
         const int iGetElementReturn = CVectorGetElement(vhand,&vv,i);
         if( iGetElementReturn == CVECTOR_NOT_FOUND )
         {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: CVectorGetElement failed\n" );
         }
         else
         {
            if( CNearTreeDistsq(vNearCenter,vv,17,CNEARTREE_TYPE_DOUBLE)!=0. && 
                sqrt(CNearTreeDistsq(vNearCenter,vv,17,CNEARTREE_TYPE_DOUBLE))< dmin)
            {
                dmin = sqrt(CNearTreeDistsq(vNearCenter,vv,17,CNEARTREE_TYPE_DOUBLE));
                CopyVec17(vCloseToNearCenter,vv);
            }
            }
        }
        
        if( dmin == DBL_MAX ) 
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere failed\n" );
        }
      else
      {
        {
            /* Using zero radius, check that only one point is found when a point is searched
             with FindInSphere */
            CVectorClear(vhand);
            bResult = !CNearTreeFindInSphere(tree, 0.0, vhand, NULL, vNearCenter,1);
            if (!bResult)
                fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFindInSphere has failed\n" );
            iFound = CVectorSize(vhand);
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere found no points using zero radius\n" );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest:testBigVector: FindInSphere found more than %ld points using zero radius\n", iFound );
            }
        }
        
        {
            /* Using minimal radius, check that at least 2 points are found when a point is searched
             with FindInSphere */
            CVectorClear(vhand);
            bResult = !CNearTreeFindInSphere(tree, 
                        sqrt(CNearTreeDistsq(vCloseToNearCenter,vNearCenter,17,CNEARTREE_TYPE_DOUBLE)),
                        vhand, NULL, vNearCenter,1);
            if (!bResult)
                fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFindInSphere has failed\n" );
            iFound = CVectorSize(vhand);
            if( iFound < 2 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere found only 1 point\n" );
            }
        }
        
        {
            /* Using small radius, check that only one point is found when a point is searched
             with FindInSphere */
            CVectorClear(vhand);
            bResult = !CNearTreeFindInSphere(tree, 
                         0.9* sqrt(CNearTreeDistsq(vCloseToNearCenter,vNearCenter,17,CNEARTREE_TYPE_DOUBLE)),
                                            vhand, NULL, vNearCenter,1);
            if (!bResult)
                fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFindInSphere has failed\n" );
            iFound = CVectorSize(vhand);
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere found no points using %f radius\n", 
                        0.9* sqrt(CNearTreeDistsq(vCloseToNearCenter,vNearCenter,17,CNEARTREE_TYPE_DOUBLE)) );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere found %ld points using %f radius\n", 
                       iFound, 
                       0.9* sqrt(CNearTreeDistsq(vCloseToNearCenter,vNearCenter,17,CNEARTREE_TYPE_DOUBLE)) );
            }
            }
        }
        
        /* Just make sure that FarthestNeighbor works for objects */
         bResult = !CNearTreeFarthestNeighbor(tree, (void FAR * FAR *)&vExtreme, NULL, vNearCenter );
         if (!bResult)
             fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFarthestNeighbor has failed\n" );
           
         bResult = !CVectorFree(&vhand);
         if (!bResult)
             fprintf(stdout, "CNearTreeTest: testBigVector: CVectorFree has failed\n" );
           
         bResult = !CNearTreeFree(&tree);
         if (!bResult)
            fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFree has failed\n" );
    }
}



