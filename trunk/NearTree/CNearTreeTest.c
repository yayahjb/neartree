/*
 *  CNearTreeTest.c
 *  NearTree
 *
 *  Based on CNearTreeTest.cpp
 *  Copyright 2008 Larry Andrews.  All rights reserved
 *
 *  C Version created by Herbert J. Bernstein on 1/1/09
 *  with permission from Larry Andrews.
 *
 *  Rev 18 Jan 2009, Tests of DelayedInsert
 *
 *  Copyright 2008, 2009 Larry Andrews and Herbert J. Bernstein. 
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
#include <stdlib.h>
#ifndef USE_LOCAL_HEADERS
#include <CNearTree.h>
#else
#include "CNearTree.h"
#endif

#ifdef USE_MINGW_RAND
#define random(x) rand(x)
#define srandom(x) srand(x)
#endif

void testEmptyTree( void );
void testLinearTree( const int n );
void testFindFirstObject( void );
void testFindLastObject( void );
void testFindInSphereFromBottom( void );
void testBackwardForward( void );
void testFindInSphereFromTop( void );
void testFindOutSphere( void );
void testDelayedInsertion( void );
void testIterators( void );
void testFindInAnnulus( void );
void testMisc( void );
void testRandomTree1( const int n );
void testBigVector( void );
void testDelayedInsertion( void );
void testKNearFar( void );

long g_errorCount;
int dbgflg;

#define bool int    /* for older C compilers */

#define MYRAND_MAX  32767

/*=======================================================================*/
int main(int argc, char** argv)
{
    int i;
    dbgflg = 0;
    
    if (argc > 1 && !strncmp(argv[1],"--debug",7)) {
        dbgflg = 1;
    }
    
    if (dbgflg) fprintf(stderr,"Debug enabled\n");
    
    /* test the interface with an empty tree */
    testEmptyTree( );
    
    /* test the C API with trees with varying content, one entry and several */
    for( i=1; i<10; ++i )
    {
        testLinearTree( i );
    }
    
    testFindFirstObject( );
    fprintf( stdout, "testFindFirstObject\n" );
    testFindLastObject( );
    fprintf( stdout, "testFindLastObject\n" );
    testFindInSphereFromBottom( );
    fprintf( stdout, "testFindInSphereFromBottom\n" );
    testFindInSphereFromTop( );
    fprintf( stdout, "testFindInSphereFromTop\n" );
    testFindOutSphere();
    fprintf( stdout, "testFindOutSphere\n" );
    testFindInAnnulus( );
    fprintf( stdout, "testFindInAnnulus\n" );
    testRandomTree1( 10000 );
    fprintf( stdout, "testRandomTree1\n" );
    testBackwardForward( );
    fprintf( stdout, "testBackwardForward\n" );
    testBigVector( );
    fprintf( stdout, "testBigVector\n" );
    testDelayedInsertion( );
    fprintf( stdout, "testDelayedInsertion\n" );
    testKNearFar();
    fprintf( stdout, "testKNearFar\n" );
    
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
    void FAR * vclose;
    int FAR * nFar;
    void FAR * vnFar;
    int probe[1];
    bool bTreeEmpty;
    bool bTreeHasNearest;
    bool bTreeHasFarthest;
    bool bInSphere;
    size_t lFoundPointsInSphere;
    
    if (CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER))
    {
        ++g_errorCount;
        fprintf(stdout,"CNearTreeTest: testEmptyTree: CNearTreeCreate failed.\n");
    }
    
    bTreeEmpty = !CNearTreeZeroIfEmpty(tree);
    if( ! bTreeEmpty )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: CNearTreeZeroIfEmpty not zero for empty tree\n" );
    }
    
    probe[0] = 1;
    bTreeHasNearest =  !CNearTreeNearestNeighbor(tree, 0.0, &vclose, NULL, probe); 
    close = (int FAR *)vclose;
    if( bTreeHasNearest )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: NearestNeighbor incorrect for empty tree\n" );
    }
    
    probe[0] = 0;
    bTreeHasFarthest = !CNearTreeFarthestNeighbor(tree, &vnFar, NULL, probe);
    nFar = (int FAR *)vnFar;
    if( bTreeHasFarthest )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: FarthestNeighbor incorrect for empty tree\n" );
    }
    
    probe[0] = 1;
    if (CVectorCreate(&v,sizeof(int FAR *),10))
    {
        ++g_errorCount;
        fprintf(stdout,"CNearTreeTest: testEmptyTree: CVectorCreate failed\n" );
    }
    
    bInSphere = !CNearTreeFindInSphere( tree, 1000.0, v, NULL, probe, 1);
    lFoundPointsInSphere = CVectorSize(v);
    if( bInSphere || lFoundPointsInSphere != 0 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: FindInSphere incorrect for empty tree\n" );
    }
    
    if (CVectorFree(&v))
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: CNearVectorFree failed\n");
    }
    
    if (CNearTreeFree(&tree))
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testEmptyTree: CNearTreeFree failed\n");
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
    CNearTreeHandle tree;
    CVectorHandle v;
    bool bInsert, bClose, bFar, bInSphere, bResult;
    int FAR * closest;   
    void FAR * vclosest;
    int FAR * farthest;
    void FAR * vfarthest;
    int probe[1];
    size_t lFoundPointsInSphere;
    int i;
    size_t found;
    long localErrorCount = 0;
    long localErrorMax = 0;
    
    if (CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER))
    {
        ++g_errorCount;
        fprintf(stdout,"CNearTreeTest: testLinearTree: CNearTreeCreate failed.\n");
    }
    
    /* generate an unbalanced tree*/
    for( i=1; i<=n; ++i )
    {
        probe[0] = i;
        bInsert = !CNearTreeInsert(tree, probe, NULL);
        if (!bInsert)
        {
           ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testLinearTree: CNearTreeInsert failed.\n");
        }
    }
    
    /*
     Search for the nearest value using a probe point that is larger than 
     the largest value that was input. The returned values should be the
     last value entered into the tree. 
     */
    probe[0] = 2*n;
    bClose = !CNearTreeNearestNeighbor(tree,22.,&vclosest,NULL,probe);
    closest = (int FAR *)vclosest;
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
    bFar = !CNearTreeFarthestNeighbor(tree,&vfarthest,NULL,probe);
    farthest = (int FAR *)vfarthest;
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
    {
        ++g_errorCount;
        fprintf(stdout,"CNearTreeTest: testLinearTree: CVectorCreate failed\n" );
    }
    
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
            if( found > (size_t)localErrorMax ) 
            {
                localErrorMax = (long)found;
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
    if( found != (size_t)n )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testLinearTree: FindInSphere did not find all the points, found %lu\n", (long unsigned int)found );
    }

    bResult = !CVectorFree(&v);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testLinearTree: CVectorFree has failed\n" );
    }
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testLinearTree: CNearTreeFree has failed\n" );
    }
    
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
    void FAR * vclosest;
    double FAR * farthest;
    void FAR * vfarthest;
    size_t lFound;
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CNearTreeCreate failed\n" );
    }
    
    /* build the double tree starting with 1.0 */
    f[0] = 1.0;
    while( f[0]*f[0] > 0.0 )
    {
        bReturn = !CNearTreeInsert(tree,f, NULL);
        fFinal = f[0];
        f[0] /= 2.0;
        ++count;
    }
    
    if( (!CNearTreeZeroIfEmpty(tree) ) )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: incorrectly found empty tree for double\n" );
    }
    
    /*
     Search for the value closest to zero. It should be a very small, probably
     denormalized number.
     */
    f[0] = 0.0;
    bReturnNear = !CNearTreeNearestNeighbor(tree, 1.0e-10, &vclosest, NULL, f);
    closest = (double FAR *)vclosest;
    if( ! bReturnNear )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: failed for double\n" );
    }
    
    if( bReturnNear && closest[0] != fFinal )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: failed for double, got %g\n", closest[0] );
    }
    
    /*
     Search for the value farthest from a large number. It should be a
     very small, probably denormalized number.
     */
    f[0] = 100.0;
    bReturnFar = !CNearTreeFarthestNeighbor(tree, &vfarthest, NULL,f);
    farthest = (double FAR *)vfarthest;
    if( ! bReturnFar || farthest[0] != fFinal )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: failed for double, got %g\n", farthest[0] );
    }
    
    /*
     Determine if FindInSphere can find all of the input data.
     */
    f[0] = 1.0;
    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CVectorCreate failed\n" );
    }
    bReturn = !CNearTreeFindInSphere( tree, 100.0, v, NULL, f, 1);
    lFound = CVectorSize(v);
    if( lFound != (size_t)count )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: found wrong count for FindInSphere for float, should be%ld, got %lu\n", count, (long unsigned int)lFound );
    }
 
    bResult = !CVectorFree(&v);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CVectorFree has failed\n" );
    }
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindFirstObject: CNearTreeFree has failed\n" );
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
    double fFinal;
    CNearTreeHandle tree;
    bool bReturn, bResult;
    double f[1];
    double FAR * closest;
    void FAR * vclosest;
    int count;
    
    count = 0;
    fFinal = DBL_MAX;
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeCreate failed\n" );
    }
    
    f[0] = 1.0;
    /* generate an unbalanced tree*/
    while( f[0] > 2.0*sqrt(DBL_MIN) )
    {
        bReturn = !CNearTreeInsert(tree,f, NULL);
        fFinal = f[0];
        f[0] /= 2.0;
        ++count;
    }

    if (dbgflg) {
        size_t depth, size;
        bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
        if (bResult) {
            fprintf(stderr,"CNearTreeTest:  testFindLastObject: depth=%lu, size=%lu\n",
                    (unsigned long)depth,(unsigned long)size);
        } else {
            fprintf(stderr,"CNearTreeTest:  testFindLastObject: get size and depth failed\n");
        }
    }
    
    
    f[0] = 11.0;
    bReturn = !CNearTreeNearestNeighbor(tree, 100.0, &vclosest, NULL, f);
    closest = (double FAR *)vclosest;
    if( ! bReturn || closest[0] != 1.0 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeNearestNeighbor failed for double, got %g\n", closest[0] );
    }
    
    f[0] = -11.0;
    bReturn = !CNearTreeFarthestNeighbor(tree, &vclosest, NULL, f);
    closest = (double FAR *)vclosest;
    if( ! bReturn || closest[0] != 1.0 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeFarthestNeighbor failed for double, got %g\n", closest[0] );
    }
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeFree has failed\n" );
    }
  
    count = 0;
    fFinal = DBL_MAX;
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE|CNEARTREE_FLIP);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeCreate with flip failed\n" );
    }
    
    f[0] = 1.0;
    /* generate an unbalanced tree*/
    while( f[0] > 2.0*sqrt(DBL_MIN) )
    {
        bReturn = !CNearTreeInsert(tree,f, NULL);
        fFinal = f[0];
        f[0] /= 2.0;
        ++count;
    }

    if (dbgflg) {
        size_t depth, size;
        bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
        if (bResult) {
            fprintf(stderr,"CNearTreeTest:  testFindLastObject: with flip depth=%lu, size=%lu\n",
                    (unsigned long)depth,(unsigned long)size);
        } else {
            fprintf(stderr,"CNearTreeTest:  testFindLastObject: with flip get size and depth failed\n");
        }
    }
    
    
    f[0] = 11.0;
    bReturn = !CNearTreeNearestNeighbor(tree, 100.0, &vclosest, NULL, f);
    closest = (double FAR *)vclosest;
    if( ! bReturn || closest[0] != 1.0 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeNearestNeighbor with flip failed for double, got %g\n", closest[0] );
    }
    
    f[0] = -11.0;
    bReturn = !CNearTreeFarthestNeighbor(tree, &vclosest, NULL, f);
    closest = (double FAR *)vclosest;
    if( ! bReturn || closest[0] != 1.0 )
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeFarthestNeighbor with flip failed for double, got %g\n", closest[0] );
    }
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindLastObject: CNearTreeFree with flip has failed\n" );
    }
    
    
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
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromBottom: CNearTreeCreate failed\n" );
    }

    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromBottom: CVectorCreate failed\n" );
    }
    
    for( i=1; i<=nmax; ++i )
    {
        f[0] = (double)i;
        bReturn = !CNearTreeInsert(tree,f, NULL);
        if (!bReturn)
        {
           ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindInSphereFromBottom: CNearTreeInsert failed\n" );
        }
    }
    
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
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CNearTreeCreate failed\n" );
    }
    
    for( i=1; i<=nmax; ++i )
    {
        f[0] = (double)i;
        bReturn = !CNearTreeInsert(tree,f, NULL);
    }
   
    
    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CVectorCreate failed\n" );
    }
    
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
            fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CNearTreeFindInSphere failed for i=%d\n", i );
        }
    }
    
    bReturn = !CVectorFree(&v);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CVectorFree failed\n");
    }
    bReturn = !CNearTreeFree(&tree);
    if (!bReturn) {
            ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInSphereFromTop: CNearTreeFree failed\n");
    }
    
}


/*=======================================================================*/
void testFindOutSphere( void )
{
    const int nmax = 100;
    CNearTreeHandle tree;
    CNearTreeHandle sphereReturn;
    CVectorHandle v;
    bool bReturn;
    double f[1];
    double radius;
    size_t lReturned;
    int i;

    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CNearTreeCreate failed\n" );
    }
    
    bReturn = !CNearTreeCreate(&sphereReturn,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CNearTreeCreate failed\n" );
    }

    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CVectorCreate failed\n" );
    }
    
    for( i=1; i<=nmax; ++i )
    {
        f[0] = (double)i;
        bReturn = !CNearTreeInsert(tree,f, NULL);
        if (!bReturn)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindOutSphere: CNearTreeInsert failed\n" );
        }
    }
    
    /* now test the NearTree return */
    for( i=1; i<=nmax; ++i )
    {
        radius = 0.05 + (double)i;
        f[0] = 0.1 + (double)nmax;
        bReturn = !CNearTreeFindTreeOutSphere( tree, radius, sphereReturn, f, 1);
        lReturned =  CNearTreeSize(sphereReturn);
        if( lReturned != 100-(size_t)i )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindOutSphere: FindOutSphere failed for i=%d\n", i );
        }
    }
    
    /* now test the vector return */
    for( i=1; i<=nmax; ++i )
    {
        bReturn = !CVectorClear(v);
        radius = 0.05 + (double)i;
        f[0] = 0.1 + (double)nmax;
        bReturn = !CNearTreeFindOutSphere( tree, radius, v, NULL, f, 1);
        lReturned = CVectorSize(v);
        if( lReturned != 100-(size_t)i )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindOutSphere: FindOutSphere failed for for vector for i=%d\n", i );
        }
    }
    
    bReturn = !CVectorFree(&v);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CVectorFree failed\n");
    }
    
    bReturn = !CNearTreeFree(&sphereReturn);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CNearTreeFree failed\n");
    }

    bReturn = !CNearTreeFree(&tree);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CNearTreeFree failed\n");
    }
    
}

/*=======================================================================*/
void testFindInAnnulus( void )
{
    const int nmax = 1000;
    CNearTreeHandle tree;
    CNearTreeHandle annulusReturn;
    CVectorHandle v;
    bool bReturn;
    double f[1];
    double FAR * nearest;
    void FAR * vnearest;
    double radiusin, radiusout;
    size_t lReturned;
    int i;
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInAnnulus: CNearTreeCreate failed\n" );
    }
    
    bReturn = !CNearTreeCreate(&annulusReturn,1,CNEARTREE_TYPE_DOUBLE);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInAnnulus: CNearTreeCreate failed\n" );
    }
    
    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindInAnnulus: CVectorCreate failed\n" );
    }
    
    for( i=1; i<=nmax; ++i )
    {
        f[0] = (double)i;
        bReturn = !CNearTreeInsert(tree,f, NULL);
        if (!bReturn)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindInAnnulus: CNearTreeInsert failed\n" );
        }
    }
    
    radiusin = 100.1;
    radiusout = 299.9;

    /* now test the NearTree return */
    {
        f[0] = 0.;
        bReturn = !CNearTreeFindTreeInAnnulus( tree, radiusin, radiusout, annulusReturn, f, 1);
        lReturned =  CNearTreeSize(annulusReturn);
        if( lReturned != (299-101+1) )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindInAnnulus: wrong number of objects found\n" );
        }
        
        f[0] = 0.;
        bReturn = !CNearTreeNearestNeighbor(annulusReturn,1000.,&vnearest,NULL,f);
	nearest = (double FAR *)vnearest;
        if (!bReturn) {
            ++g_errorCount;
            fprintf(stdout, "testFindInAnnulus: no lowest value found\n" );
        } else if( nearest[0] != 101. )
        {
            ++g_errorCount;
            fprintf(stdout, "testFindInAnnulus: lowest value (%g) is incorrect\n", nearest[0] );
        }
        
    }
    
    /* now test the vector return */
    {
        f[0] = 0.;
        bReturn = !CNearTreeFindInAnnulus( tree, radiusin, radiusout, v, NULL, f, 1);
        lReturned = CVectorSize(v);
        if( lReturned != (299-101+1) )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testFindInAnnulus: for vector, wrong number of objects found\n" );
        }
    }
    
    bReturn = !CVectorFree(&v);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CVectorFree failed\n");
    }
    
    bReturn = !CNearTreeFree(&annulusReturn);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CNearTreeFree failed\n");
    }
    
    bReturn = !CNearTreeFree(&tree);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testFindOutSphere: CNearTreeFree failed\n");
    }
    
}


/*
 Test CNearTree with a bunch of random numbers (integers).
 */
/*=======================================================================*/
void testRandomTree1( const int nRequestedRandoms )
{
    const int n = nRequestedRandoms;
    CNearTreeHandle tree;
    CVectorHandle v;
    int nmax;
    int nmin;
    int i;
    int next[1];
    int FAR * closest;
    void FAR * vclosest;
    int FAR * farthest;
    void FAR * vfarthest;
    bool bReturn, bResult, bNear, bFar;
    double radius;
    size_t lReturn;
    
    nmax = INT_MIN;
    nmin = INT_MAX;

    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testRandomTree1: CNearTreeCreate failed\n" );
    }
    
    /* Build the tree with n random numbers. Remember the largest and smallest values. */
    for( i=0; i<n; ++i )
    {
        next[0] = random( )%MYRAND_MAX;
        bReturn = !CNearTreeInsert(tree,next,NULL);
        if( next[0] > nmax ) nmax = next[0];
        if( next[0] < nmin ) nmin = next[0];
    }

    if (dbgflg) {
        size_t depth, size;
        bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
        if (bResult) {
            fprintf(stderr,"CNearTreeTest:  testRandomTree1: depth=%lu, size=%lu\n",
                    (unsigned long)depth,(unsigned long)size);
        } else {
            fprintf(stderr,"CNearTreeTest:  testRandomTree1: get size and depth failed\n");
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from below)*/
        next[0] = INT_MIN/2;
        bNear = !CNearTreeNearestNeighbor(tree, (double)LONG_MAX, &vclosest, NULL, next);
        closest = (int FAR *)vclosest;
        if( ! bNear || closest[0] != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: NearestNeighbor failed in testRandomTree for min\n" );
        }
        
        next[0] = INT_MIN/2;
        bFar = !CNearTreeFarthestNeighbor(tree, &vfarthest, NULL, next);
        farthest = (int FAR *)vfarthest;
        if( !bFar || farthest[0] != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FarthestNeighbor failed in testRandomTree for min\n" );
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from above)*/
        next[0] = INT_MAX/2;
        bNear = !CNearTreeNearestNeighbor(tree, (double)LONG_MAX, &vclosest, NULL, next);
        closest = (int FAR *)vclosest;
        if( ! bNear || closest[0] != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: NearestNeighbor failed for max\n" );
        }
        
        next[0] = INT_MAX/2;
        bFar = !CNearTreeFarthestNeighbor(tree, &vfarthest, NULL, next);
        farthest = (int FAR *)vfarthest;
        if( !bFar || farthest[0] != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FarthestNeighbor failed for max\n" );
        }
    }
    
    {
        /*verify that for very large radius, every point is detected (from below)*/
        bReturn = !CVectorCreate(&v,sizeof(int FAR *), 10);
        if (!bReturn)
        {
           ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: CVectorCreate failed\n" );
        }
        radius = DBL_MAX;
        next[0] = 1;
        bReturn = !CNearTreeFindInSphere( tree, radius, v, NULL, next, 1);
        lReturn = CVectorSize(v);
        if( lReturn != (size_t)n )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere failed, n=%d, lReturn=%lu\n", n, (long unsigned int)lReturn );
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
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere failed found points incorrectly, n=%d, lReturn=%lu\n", n, (long unsigned int)lReturn );
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
            if( lReturn < (size_t)lastFoundCount )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere found DECREASING count on increasing radius for radius=%g\n", radius );
                break;
            }
            else
            {
                lastFoundCount = (int)lReturn;
                radius *= 1.414;
            }
        }
        
        /* verify that all points were included in the final check (beyond all reasonable distances) */
        if( lReturn != (size_t)n )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere did not find all the points\n" );
        }
    }

    nmax = INT_MIN;
    nmin = INT_MAX;

    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER|CNEARTREE_FLIP);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testRandomTree1: CNearTreeCreate with flip failed\n" );
    }
    
    
    /* Build the tree with n random numbers. Remember the largest and smallest values. */
    for( i=0; i<n; ++i )
    {
        next[0] = random( )%MYRAND_MAX;
        bReturn = !CNearTreeInsert(tree,next,NULL);
        if( next[0] > nmax ) nmax = next[0];
        if( next[0] < nmin ) nmin = next[0];
    }
 
    if (dbgflg) {
        size_t depth, size;
        bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
        if (bResult) {
            fprintf(stderr,"CNearTreeTest:  testRandomTree1: with flip depth=%lu, size=%lu\n",
                    (unsigned long)depth,(unsigned long)size);
        } else {
            fprintf(stderr,"CNearTreeTest:  testRandomTree1: with flip get size and depth failed\n");
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from below)*/
        next[0] = INT_MIN/2;
        bNear = !CNearTreeNearestNeighbor(tree, (double)LONG_MAX, &vclosest, NULL, next);
        closest = (int FAR *)vclosest;
        if( ! bNear || closest[0] != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: NearestNeighbor with flip failed in testRandomTree for min\n" );
        }
        
        next[0] = INT_MIN/2;
        bFar = !CNearTreeFarthestNeighbor(tree, &vfarthest, NULL, next);
        farthest = (int FAR *)vfarthest;
        if( !bFar || farthest[0] != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FarthestNeighbor with flip failed in testRandomTree for min\n" );
        }
    }
    
    {
        /*verify that the correct extremal point is detected (from above)*/
        next[0] = INT_MAX/2;
        bNear = !CNearTreeNearestNeighbor(tree, (double)LONG_MAX, &vclosest, NULL, next);
        closest = (int FAR *)vclosest;
        if( ! bNear || closest[0] != nmax )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: NearestNeighbor with flip failed for max\n" );
        }
        
        next[0] = INT_MAX/2;
        bFar = !CNearTreeFarthestNeighbor(tree, &vfarthest, NULL, next);
        farthest = (int FAR *)vfarthest;
        if( !bFar || farthest[0] != nmin )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FarthestNeighbor with flip failed for max\n" );
        }
    }
    
    {
        /*verify that for very large radius, every point is detected (from below)*/
        bReturn = !CVectorCreate(&v,sizeof(int FAR *), 10);
        if (!bReturn)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: CVectorCreate failed\n" );
        }
        radius = DBL_MAX;
        next[0] = 1;
        bReturn = !CNearTreeFindInSphere( tree, radius, v, NULL, next, 1);
        lReturn = CVectorSize(v);
        if( lReturn != (size_t)n )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere with flip failed, n=%d, lReturn=%lu\n", n, (long unsigned int)lReturn );
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
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere with flip failed found points incorrectly, n=%d, lReturn=%lu\n", n, (long unsigned int)lReturn );
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
            if( lReturn < (size_t)lastFoundCount )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere with flip found DECREASING count on increasing radius for radius=%g\n", radius );
                break;
            }
            else
            {
                lastFoundCount = (int)lReturn;
                radius *= 1.414;
            }
        }
        
        /* verify that all points were included in the final check (beyond all reasonable distances) */
        if( lReturn != (size_t)n )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testRandomTree1: FindInSphere with flip  did not find all the points\n" );
        }
    }
} /* end testRandomTree1 */

/*  Utility routines for 17D double vectors */

int LoadVec17(double vec[17]) {
    size_t i;
    for (i = 0; i < 17; i++) {
        vec[i] = (double)(random( )%MYRAND_MAX);
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
    void FAR * vvFarthest;
    double FAR * vNearCenter;
    void FAR * vvNearCenter;
    double FAR * vExtreme;
    void FAR * vvExtreme;
    
    double rmax = -DBL_MAX;
    double rmin =  DBL_MAX;
    double dmin, dmax;
    
    size_t iFoundNearCenter, iFound, i;
    
    bool bResult;
    
    bResult = !CNearTreeCreate(&tree,17,CNEARTREE_TYPE_DOUBLE);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout,"CNearTreeTest: testBigVector: CNearTreeCreate failed\n");
    }
    
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
        {
           ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testBigVector: CNearTreeInsert failed\n");
        }
        CopyVec17(vAll[i],v);
    }
    
    {        
        /* Find the point farthest from the point that was nearest the origin. */
        bResult = !CNearTreeFarthestNeighbor(tree,&vvFarthest,NULL,v17min);
        vFarthest = (double FAR *)vvFarthest;
        if (!bResult)
        {
           ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testBigVector: CNearTreeFarthestNeighbor failed\n");
        }
        
        /* Brute force search for the farthest */
        dmax = -DBL_MAX;
        for( i=0; i<1000; ++i )
        {
            if (CNearTreeDist(tree,vAll[i],v17min) > dmax)
            {
                dmax = CNearTreeDist(tree,vAll[i],v17min);
                CopyVec17(vSearch, vAll[i]);
            }
        }
        
      if( CNearTreeDist(tree,vSearch,vFarthest) > DBL_MIN )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: FarthestNeighbor has failed\n" );
        }
    }
    
    {
        /* somewhere in the middle, find a point and its nearest neighbor */
        /* make sure that each includes the other in sphere search */
        
        ConstVec17(vCenter, (double)(MYRAND_MAX/2) );
        CVectorCreate(&vhand,sizeof(double FAR *),10);

        bResult = !CNearTreeNearestNeighbor( tree, 100000.0, &vvNearCenter, NULL, vCenter );
        vNearCenter = (double FAR *)vvNearCenter;
        if (!bResult)
        {
           ++g_errorCount;
          fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeNearestNeighbor has failed\n" );
        }
        bResult = !CNearTreeFindInSphere(tree, 1000000.0, vhand, NULL, vNearCenter,1);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFindInSphere has failed\n" );
        }
        iFoundNearCenter = CVectorSize(vhand);
        
        /* Brute force search for the point closest to the point closest to the center */
        dmin = DBL_MAX;
        for( i=0; i<iFoundNearCenter; ++i )
        {
         const int iGetElementReturn = CVectorGetElement(vhand,&vv,i);
         if( iGetElementReturn == CVECTOR_NOT_FOUND )
         {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: CVectorGetElement failed\n" );
         }
         else
         {
            if( CNearTreeDist(tree,vNearCenter,vv)!=0. && 
                CNearTreeDist(tree,vNearCenter,vv)< dmin)
            {
                dmin = CNearTreeDist(tree,vNearCenter,vv);
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
            iFound = CVectorSize(vhand);
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere found no points using zero radius\n" );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest:testBigVector: FindInSphere found more than 1 point using zero radius, found=%lu\n", (long unsigned int)iFound );
            }
        }
        
        {
            /* Using minimal radius, check that at least 2 points are found when a point is searched
             with FindInSphere */
            CVectorClear(vhand);
            bResult = !CNearTreeFindInSphere(tree, 
                        1.001*CNearTreeDist(tree,vCloseToNearCenter,vNearCenter),
                        vhand, NULL, vNearCenter,1);
            if (!bResult)
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFindInSphere has failed\n" );
            }
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
                         0.9* CNearTreeDist(tree,vCloseToNearCenter,vNearCenter),
                                            vhand, NULL, vNearCenter,1);
            if (!bResult)
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFindInSphere has failed\n" );
            }
            iFound = CVectorSize(vhand);
            if( iFound < 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere found no points using %g radius\n", 
                        0.9* CNearTreeDist(tree,vCloseToNearCenter,vNearCenter) );
            }
            else if( iFound != 1 )
            {
                ++g_errorCount;
                fprintf(stdout, "CNearTreeTest: testBigVector: FindInSphere found %lu points using %g radius\n", 
                       (long unsigned int)iFound, 
                       0.9* CNearTreeDist(tree,vCloseToNearCenter,vNearCenter) );
            }
            }
        }
        
        /* Just make sure that FarthestNeighbor works for objects */
         bResult = !CNearTreeFarthestNeighbor(tree, &vvExtreme, NULL, vNearCenter );
         vExtreme = (double FAR *)vvExtreme;
         if (!bResult)
            {
                ++g_errorCount;
             fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFarthestNeighbor has failed\n" );
         }
           
         bResult = !CVectorFree(&vhand);
         if (!bResult)
         {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: CVectorFree has failed\n" );
         }
           
         bResult = !CNearTreeFree(&tree);
         if (!bResult)
         {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBigVector: CNearTreeFree has failed\n" );
         }
    }
}

/*=======================================================================*/
void testBackwardForward( void )
{
    CNearTreeHandle tree;
    const int nMax = 1000;
    double data[1];
    double probe[1];
    double FAR * closest;
    bool bResult;
    bool bResult1, bResult2, bResult3, bResult4;
    int i;
    
    bResult = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);    
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout,"CNearTreeTest: testBackwardForward: CNearTreeCreate failed\n");
    }
    
    for( i=0; i<nMax; ++i )
    {
        data[0]=(double)i;        bResult1 = ! CNearTreeInsert(tree,data,NULL); 
        data[0]=(double)(nMax-i); bResult2 = ! CNearTreeInsert(tree,data,NULL); 
        data[0]=(double)i + 0.25; bResult3 = ! CNearTreeInsert(tree,data,NULL); 
        data[0]=(double)(nMax-i) + 0.75; 
                                  bResult4 = ! CNearTreeInsert(tree,data,NULL);
        if (!bResult1 || !bResult2 || !bResult3 || !bResult4) {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testBackwardForward: CNearTreeInserts failed for i = %d\n", i);
        }
    }
 
    if (dbgflg) {
        size_t depth, size;
        bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
        if (bResult) {
            fprintf(stderr,"CNearTreeTest:  testBackwardForward: depth=%lu, size=%lu\n",
                    (unsigned long)depth,(unsigned long)size);
        } else {
            fprintf(stderr,"CNearTreeTest:  testBackwardForward: get size and depth failed\n");
        }
        
    }
    
    for( i=100; i<300; ++i )
    {   probe[0] = (double)i+0.25;
        bResult =!CNearTreeNearestNeighbor(tree,1000.0,(void FAR *)&closest,NULL,probe);
        if( !bResult || fabs( closest[0]-((double)i+0.25) ) > DBL_MIN )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBackForward: NearestNeighbor failed to find %g, found %g instead\n", 
                   (double)i+0.25, closest[0] );      
        }
    }

    bResult = !CNearTreeFree(&tree);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testBackwardForward: CNearTreeFree has failed\n" );
    }
    
    bResult = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE|CNEARTREE_FLIP);    
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout,"CNearTreeTest: testBackwardForward: CNearTreeCreate failed\n");
    }
    
    for( i=0; i<nMax; ++i )
    {
        data[0]=(double)i;        bResult1 = ! CNearTreeInsert(tree,data,NULL); 
        data[0]=(double)(nMax-i); bResult2 = ! CNearTreeInsert(tree,data,NULL); 
        data[0]=(double)i + 0.25; bResult3 = ! CNearTreeInsert(tree,data,NULL); 
        data[0]=(double)(nMax-i) + 0.75; 
        bResult4 = ! CNearTreeInsert(tree,data,NULL);
        if (!bResult1 || !bResult2 || !bResult3 || !bResult4) {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testBackwardForward: CNearTreeInserts failed for i = %d\n", i);
        }
    }

    if (dbgflg) {
        size_t depth, size;
        bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
        if (bResult) {
            fprintf(stderr,"CNearTreeTest:  testBackwardForward: with flip depth=%lu, size=%lu\n",
                    (unsigned long)depth,(unsigned long)size);
        } else {
            fprintf(stderr,"CNearTreeTest:  testBackwardForward: with flip get size and depth failed\n");
        }
        
    }
    
    for( i=100; i<300; ++i )
    {   probe[0] = (double)i+0.25;
        bResult =!CNearTreeNearestNeighbor(tree,1000.0,(void FAR *)&closest,NULL,probe);
        if( !bResult || fabs( closest[0]-((double)i+0.25) ) > DBL_MIN )
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testBackForward: NearestNeighbor failed to find %g, found %g instead\n", 
                    (double)i+0.25, closest[0] );      
        }
    }
    
    bResult = !CNearTreeFree(&tree);
    if (!bResult)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testBackwardForward: CNearTreeFree has failed\n" );
    }
    
    
}

/*=======================================================================*/
void testDelayedInsertion( void )
{
    CNearTreeHandle tree;
    CVectorHandle v;
    const long nmax = 100;
    double data[1];
    double probe[1];
    double radius;
    double fFinal;
    double FAR * closest;
    double FAR * farthest;
    size_t lReturned;
    bool bResult, bResult1, bResult2;
    int i;
    
    {
        /* make sure that CompleteDelayInsert works - note there is no way in the
           final test to differentiate that from FindInSphere completing the
           flush. */
        
        bResult = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);    
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeCreate (first block) failed\n");
        }
        
        for( i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 )
            {
                data[0]=(double)i;        bResult1 = ! CNearTreeInsert(tree,data,NULL); 
                if (!bResult1) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeInsert failed for i = %d\n", i);
                }
            }
            else
            {
                data[0]=(double)i;        bResult2 = ! CNearTreeDelayedInsert(tree,data,NULL); 
                if (!bResult2) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeDelayedInsert failed for i = %d\n", i);
                }
            }
        }
        
        bResult = !CVectorCreate(&v,sizeof(double FAR *),10);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CVectorCreate failed\n");
        }
        radius = 1000.0;
        if (dbgflg) {
            size_t depth, size;
            bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
            if (bResult) {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: before CNearTreeCompleteDelayedInsert depth=%lu, size=%lu\n",
                        (unsigned long)depth,(unsigned long)size);
            } else {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: before CNearTreeCompleteDelayedInsert get size and depth failed\n");
            }
            
        }
        bResult = !CNearTreeCompleteDelayedInsert(tree);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeCompleteDelayedInsert failed\n");
        }
        if (dbgflg) {
            size_t depth, size;
            bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
            if (bResult) {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: after CNearTreeCompleteDelayedInsert depth=%lu, size=%lu\n",
                        (unsigned long)depth,(unsigned long)size);
            } else {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: after CNearTreeCompleteDelayedInsert get size and depth failed\n");
            }
            
        }
        probe[0] = 0.9;
        bResult = !CNearTreeFindInSphere(tree, radius,v,NULL,probe,1);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeFindInSphere failed\n");
        } else {
            lReturned = CVectorSize(v);
            if( lReturned != nmax )
            {
                ++g_errorCount;
                printf( "CNearTreeTest:  testDelayedInsertion: CNearTreeFindInSphere failed for nmax=%ld, found %lu points\n", nmax, (unsigned long)lReturned );
            }
        }
        bResult = !CVectorFree(&v);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest:  testDelayedInsertion: CVectorFree (first block) has failed\n" );
        }
        bResult = !CNearTreeFree(&tree);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testDelayedInsertion: CNearTreeFree (first block) has failed\n" );
        }
        
    }

    {
        /* make sure that CompleteDelayInsert works - note there is no way in the
         final test to differentiate that from FindInSphere completing the
         flush. */
        
        bResult = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE|CNEARTREE_FLIP);    
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeCreate (first block) with flip failed\n");
        }
        
        for( i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 )
            {
                data[0]=(double)i;        bResult1 = ! CNearTreeInsert(tree,data,NULL); 
                if (!bResult1) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeInsert with flip failed for i = %d\n", i);
                }
            }
            else
            {
                data[0]=(double)i;        bResult2 = ! CNearTreeDelayedInsert(tree,data,NULL); 
                if (!bResult2) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeDelayedInsert with flip failed for i = %d\n", i);
                }
            }
        }
        
        bResult = !CVectorCreate(&v,sizeof(double FAR *),10);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CVectorCreate failed\n");
        }
        radius = 1000.0;
        if (dbgflg) {
            size_t depth, size;
            bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
            if (bResult) {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: before CNearTreeCompleteDelayedInsert with flip depth=%lu, size=%lu\n",
                        (unsigned long)depth,(unsigned long)size);
            } else {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: before CNearTreeCompleteDelayedInsert with flip get size and depth failed\n");
            }
            
        }
        bResult = !CNearTreeCompleteDelayedInsert(tree);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeCompleteDelayedInsert failed\n");
        }
        if (dbgflg) {
            size_t depth, size;
            bResult = !(CNearTreeGetDepth(tree,&depth)) && !(CNearTreeGetSize(tree,&size));
            if (bResult) {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: after CNearTreeCompleteDelayedInsert depth=%lu, size=%lu\n",
                        (unsigned long)depth,(unsigned long)size);
            } else {
                fprintf(stderr,"CNearTreeTest:  testDelayedInsertion: after CNearTreeCompleteDelayedInsert get size and depth failed\n");
            }
        }
        probe[0] = 0.9;
        bResult = !CNearTreeFindInSphere(tree, radius,v,NULL,probe,1);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeFindInSphere failed\n");
        } else {
            lReturned = CVectorSize(v);
            if( lReturned != nmax )
            {
                ++g_errorCount;
                printf( "CNearTreeTest:  testDelayedInsertion: CNearTreeFindInSphere failed for nmax=%ld, found %lu points\n", nmax, (unsigned long)lReturned );
            }
        }
        bResult = !CVectorFree(&v);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest:  testDelayedInsertion: CVectorFree (first block) has failed\n" );
        }
        bResult = !CNearTreeFree(&tree);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testDelayedInsertion: CNearTreeFree (first block) has failed\n" );
        }
        
    }
    
    
    {
        /* make sure that FindInSphere flushes the delayed data */
        bResult = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);    
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeCreate (second block) failed\n");
        }
        
        for( i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 )
            {
                data[0]=(double)i;        bResult1 = ! CNearTreeInsert(tree,data,NULL); 
                if (!bResult1) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeInsert (second block) failed for i = %d\n", i);
                }
            }
            else
            {
                data[0]=(double)i;        bResult2 = ! CNearTreeDelayedInsert(tree,data,NULL); 
                if (!bResult2) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeDelayedInsert (second block) failed for i = %d\n", i);
                }
            }
        }
        
        bResult = !CVectorCreate(&v,sizeof(double FAR *),10);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CVectorCreate (second block) failed\n");
        }
        
        radius = 1000.0;
        probe[0] = 0.9;
        bResult = !CNearTreeFindInSphere(tree, radius,v,NULL,probe,1);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeFindInSphere (second block) failed\n");
        } else {
            lReturned = CVectorSize(v);
            if( lReturned != nmax )
            {
                ++g_errorCount;
                printf( "CNearTreeTest:  testDelayedInsertion: CNearTreeFindInSphere (second block) failed for nmax=%ld, found %lu points\n", nmax, (unsigned long)lReturned );
            }
        }
        bResult = !CVectorFree(&v);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest:  testDelayedInsertion: CVectorFree (second block) has failed\n" );
        }
        bResult = !CNearTreeFree(&tree);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testDelayedInsertion: CNearTreeFree (second block) has failed\n" );
        }
        
    }
    
    {
        /* make sure that NearestNeighbor flushes the delayed data */
        fFinal = DBL_MAX;

        bResult = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);    
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeCreate (third block) failed\n");
        }
        
        for( i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 && i!=nmax) /* ensure that the last one is delayed */
            {
                data[0]=(double)i;        bResult1 = ! CNearTreeInsert(tree,data,NULL); 
                if (!bResult1) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeInsert failed for i = %d\n", i);
                }
            }
            else
            {
                data[0]=(double)i;        bResult2 = ! CNearTreeDelayedInsert(tree,data,NULL); 
                if (!bResult2) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeDelayedInsert (third block) failed for i = %d\n", i);
                }
                fFinal = (double)i;
            }
        }
        
        radius = 0.1;
        probe[0] = fFinal;
        bResult = !CNearTreeNearestNeighbor( tree, radius, (void FAR *)&closest, NULL, probe );
        if( ! bResult )
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeNearestNeighbor failed \n" );
        }
        else if( closest[0] != fFinal )
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testDelayedInsertion: NearestNeighbor failed to find the data\n" );
        }
        bResult = !CNearTreeFree(&tree);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testDelayedInsertion: CNearTreeFree (third block) has failed\n" );
        }
        
    }
    
    {
        /* make sure that FarthestNeighbor flushes the delayed data */
        bResult = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_DOUBLE);    
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest:  testDelayedInsertion: CNearTreeCreate (fourth block) failed\n");
        }
        
        fFinal = DBL_MAX;
        
        for( i=1; i<=nmax; ++i )
        {
            if( (i%2) == 0 && i!=nmax) /* ensure that the last one is delayed */
            {
                data[0]=(double)i;        bResult1 = ! CNearTreeInsert(tree,data,NULL); 
                if (!bResult1) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeInsert (fourth block) failed for i = %d\n", i);
                }
            }
            else
            {
                data[0]=(double)i;        bResult2 = ! CNearTreeDelayedInsert(tree,data,NULL); 
                if (!bResult2) {
                    ++g_errorCount;
                    fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeDelayedInsert (fourth block) failed for i = %d\n", i);
                }
                fFinal = (double)i;
            }
        }
        
        probe[0] = -100;
        bResult = !CNearTreeFarthestNeighbor( tree,  (void FAR *)&farthest, NULL, probe );
        if( ! bResult )
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testDelayedInsertion: CNearTreeFarthestNeighbor failed \n" );
        } else if( farthest[0] != fFinal )
        {
            ++g_errorCount;
            fprintf(stdout,"CNearTreeTest: testDelayedInsertion: FarthestNeighbor failed to find the data\n" );
        }
        bResult = !CNearTreeFree(&tree);
        if (!bResult)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testDelayedInsertion: CNearTreeFree (third block) has failed\n" );
        }
        
    }
}

/*=======================================================================*/
void testKNearFar( void )
{
    CNearTreeHandle tree;
    CNearTreeHandle outTree;
    CVectorHandle v;
    bool bReturn;
    int searchPoint[1];
    double radius;
    size_t lReturned;
    int i;
    
    
    bReturn = !CNearTreeCreate(&tree,1,CNEARTREE_TYPE_INTEGER);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: CNearTreeCreate failed\n" );
    }
    
    bReturn = !CNearTreeCreate(&outTree,1,CNEARTREE_TYPE_INTEGER);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: CNearTreeCreate failed\n" );
    }
    
    
    bReturn = !CVectorCreate(&v,sizeof(double FAR *),10);
    if (!bReturn)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: CVectorCreate failed\n" );
    }
    
    searchPoint[0] = 50;
    radius = 1000.;
    bReturn = !CNearTreeFindKTreeNearest(tree,0,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 0)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong count #0\n" );   	
    }
    
    for( i=1; i<=100; ++i )
    {
        searchPoint[0] = i;
        bReturn = !CNearTreeInsert(tree,searchPoint, NULL);
        if (!bReturn)
        {
            ++g_errorCount;
            fprintf(stdout, "CNearTreeTest: testKNearFar: CNearTreeInsert failed\n" );
        }
    }
    
    searchPoint[0] = 50;
    radius = 3.5;
    bReturn = !CNearTreeFindKTreeNearest(tree,13,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong tree count #1\n" );   	
    }
    bReturn = !CNearTreeFindKNearest(tree,13,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong vector count #1\n" );   	
    }
    
    searchPoint[0] = 98;
    radius = 3.5;
    bReturn = !CNearTreeFindKTreeNearest(tree,13,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 6)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong tree count #2\n" );   	
    }
    bReturn = !CNearTreeFindKNearest(tree,13,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 6)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong vector count #2\n" );   	
    }
    
    
    searchPoint[0] = 50;
    radius = 12.;
    bReturn = !CNearTreeFindKTreeNearest(tree,7,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong tree count #3\n" );   	
    }
    bReturn = !CNearTreeFindKNearest(tree,7,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong vector count #3\n" );   	
    }
    
    
    searchPoint[0] = 2;
    radius = 12.;
    bReturn = !CNearTreeFindKTreeNearest(tree,7,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong tree count #4\n" );   	
    }
    bReturn = !CNearTreeFindKNearest(tree,7,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong vector count #4\n" );   	
    }
    
    searchPoint[0] = 2;
    radius = 3.5;
    bReturn = !CNearTreeFindKTreeNearest(tree,7,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 5)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong tree count #5\n" );   	
    }
    bReturn = !CNearTreeFindKNearest(tree,7,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 5)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Near: found wrong vector count #5\n" );   	
    }
    
    
    searchPoint[0] = -1;
    radius = 0.;
    bReturn = !CNearTreeFindKTreeFarthest(tree,13,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 13)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong tree count #1\n" );   	
    }
    bReturn = !CNearTreeFindKFarthest(tree,13,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 13)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong vector count #1\n" );   	
    }
    
    searchPoint[0] = 50;
    radius = 0.;
    bReturn = !CNearTreeFindKTreeFarthest(tree,13,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 13)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong tree count #2\n" );   	
    }
    bReturn = !CNearTreeFindKFarthest(tree,13,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 13)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong vector count #2\n" );   	
    }
    
    searchPoint[0] = 150;
    radius = 0.;
    bReturn = !CNearTreeFindKTreeFarthest(tree,7,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong tree count #3\n" );   	
    }
    bReturn = !CNearTreeFindKFarthest(tree,7,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong vector count #3\n" );   	
    }
    
    searchPoint[0] = 46;
    radius = 0.;
    bReturn = !CNearTreeFindKTreeFarthest(tree,12,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 12)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong tree count #4\n" );   	
    }
    bReturn = !CNearTreeFindKFarthest(tree,12,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 12)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong vector count #4\n" );   	
    }
    
    
    searchPoint[0] = 2;
    radius = 0.;
    bReturn = !CNearTreeFindKTreeFarthest(tree,7,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong tree count #5\n" );   	
    }
    bReturn = !CNearTreeFindKFarthest(tree,7,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong vector count #5\n" );   	
    }
    
    searchPoint[0] = 200;
    radius = 0.;
    bReturn = !CNearTreeFindKTreeFarthest(tree,7,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong tree count #6\n" );   	
    }
    bReturn = !CNearTreeFindKFarthest(tree,7,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 7)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong vector count #6\n" );   	
    }
    
    searchPoint[0] = 2;
    radius = 95.;
    bReturn = !CNearTreeFindKTreeFarthest(tree,7,radius,outTree,searchPoint,1);
    lReturned = CNearTreeSize(outTree);
    if (lReturned != 4)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong tree count #7\n" );   	
    }
    bReturn = !CNearTreeFindKFarthest(tree,7,radius,v,NULL,searchPoint,1);
    lReturned = CVectorSize(v);
    if (lReturned != 4)
    {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: Far: found wrong vector count #7\n" );   	
    }
    
    bReturn = !CVectorFree(&v);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: CVectorFree failed\n");
    }
    
    bReturn = !CNearTreeFree(&outTree);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: CNearTreeFree failed\n");
    }
    
    bReturn = !CNearTreeFree(&tree);
    if (!bReturn) {
        ++g_errorCount;
        fprintf(stdout, "CNearTreeTest: testKNearFar: CNearTreeFree failed\n");
    }
    
    
    
    
}




