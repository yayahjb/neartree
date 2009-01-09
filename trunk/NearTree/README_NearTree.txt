
                                   NearTree

                                  Release 1.0
8 January 2009 (c) Copyright 2001, 2008, 2009 Larry Andrews. All rights reserved
                                    based on
         Larry Andrews, "A template for the nearest neighbor problem",
   C/C++ Users Journal, Volume 19 , Issue 11 (November 2001), 40 - 49 (2001),
                                ISSN:1075-2838,
                        www.ddj.com/architect/184401449

   Revised 12 Dec 2008, 8 January 2009 for sourceforge release, Larry Andrews
                            and Herbert J. Bernstein

    YOU MAY REDISTRIBUTE NearTree UNDER THE TERMS OF THE LGPL

+------------------------------------------------------------------------------+
|                                 LGPL NOTICES                                 |
|                                                                              |
| This library is free software; you can redistribute it and/or modify it      |
| under the terms of the GNU Lesser General Public License as published by the |
| Free Software Foundation; either version 2.1 of the License, or (at your     |
| option) any later version.                                                   |
|                                                                              |
| This library is distributed in the hope that it will be useful, but WITHOUT  |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        |
| FITNESS FOR A PARTICULAR PURPOSE. See the GNU* Lesser General Public License |
| for more details.                                                            |
|                                                                              |
| You should have received a copy of the GNU Lesser General Public License     |
| along with this library; if not, write to the Free Software Foundation,      |
| Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA                 |
+------------------------------------------------------------------------------+

   This is a release of an API for finding nearest neighbors among points in
   spaces of arbtrary dimensions. This release provides a C++ template,
   TNear.h, and a C library, CNearTree.c, with example/test programs.

     * Installation
     * The C++ template, TNear.h
     * The C NearTree API CNearTree.c

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

    Installation

   The NearTree package is available at
   www.sourceforge.net/projects/neartree. A source tarball is available at
   downloads.sourceforge.net/neartree/NearTree-1.0.0.tar.gz. Later tarballs
   may be available.

   When the source tarball is dowloaded and unpacked, you should have a
   directory NearTree-1.0.0. To see the current settings for a build execute

   make

   which should give the following information:

  PLEASE READ README_NearTree.txt and lgpl.txt
 
  Before making the NearTree libraries and example programs, check
  that the chosen settings are correct
 
  The current C++ and C compile commands are:
 
    libtool --mode=compile g++ -g -O2  -Wall -ansi -pedantic -I.   -c
    libtool --mode=compile gcc -g -O2  -Wall -ansi -pedantic -I.   -c
 
  The C API, CNearTree.c, depends on the sourceforge project CVector
  You are currently setup to use the system defaults for CVector
  If that is not correct, define the variable CVECTOR_INCLUDE
 
    libtool --mode=compile g++ -g -O2  -Wall -ansi -pedantic -I.   -c
    libtool --mode=compile gcc -g -O2  -Wall -ansi -pedantic -I.   -c
 
  The current library link command is:
 
    libtool --mode=link  gcc -version-info 1:0:1 -release 1.0 -rpath \
        /usr/local/lib
 
  The current C++ and C library local, and C dynamic and static build 
  commands are:
 
    libtool --mode=link g++ -g -O2  -Wall -ansi -pedantic -I.
    libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic -I.  -L./lib
    libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic -dynamic \
       -I /usr/local/include -L/usr/local/lib
    libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic -static \
       -I /usr/local/include -L/usr/local/lib
 
  Before installing the NearTree library and example programs, check
  that the install directory and install commands are correct:
 
  The current values are :
 
    /usr/local
    libtool --mode=install cp
    
 
  To compile the NearTree library and example programs type:
 
    make clean
    make all
 
  To run a set of tests type:
 
    make tests
 
  To clean up the directories type:
 
    make clean
 
  To install the library and headers type:
 
    make install


   If these settings need to be changed, edit Makefile. On some systems, e.g.
   Mac OS X, the default libtool is not appropriate. In that case you should
   install a recent version of libtool. The CVector kit has been tested with
   libtool versions 1.3.5 and 1.5.4. If the system libtool is not to be used,
   define the variable LIBTOOL to be the path to the libtaol executable, e.g.
   in bash

   export LIBTOOL=$HOME/bin/libtool

   or in the Makefie

   LIBTOOL = $(HOME)/bin/libtool

   If you need to include local header files using #include "..." instead of
   #include <...>, define the variable USE_LOCAL_HEADERS

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

    The C++ template, TNear.h

   This is a revised release of

     template <typename T> class CNearTree;

   Nearest Neighbor algorithm after Kalantari and McDonald, (IEEE
   Transactions on Software Engineering, v. SE-9, pp. 631-634,1983) modified
   to use recursion instead of a double-linked tree and simplified so that it
   does a bit less checking for things like is the distance to the right less
   than the distance to the left; it was found that these checks little to no
   difference.

   The user of this class needs to provide at least the following
   functionality for the template to work. For the built-in numerics of C++,
   they are provided by the system.

operator double( ); // conversion constructor from the templated class to double 
                    (usually will return a "length")                             
operator- ( );      // geometrical (vector) difference of two objects            
                    // a copy constructor                                        
                    // a constructor would be nice                               
                    // a destructor would be nice                                

   The provided interface is:

     #include <TNear.h>

     CNearTree( void )   // constructor
        instantiated by something like:      CNearTree <v> vTree;
        for some type v

     ~CNearTree( void )  // destructor
        invoked by  vTree.CNeartree::~CNearTree
        for an object vTree of some type v

     bool empty( ) const
        test for an empty NearTree

     void Insert( T& t )
        where t is an object of the type v

     bool NearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
        dRadius is the largest radius within which to search; make it
           very large if you want to include every point that was loaded; dRadius
           is returned as the closest distance to the probe (or the search radius
           if nothing is found)
        tClosest is returned as the object that was found closest to the probe
           point (if any were within radius dRadius of the probe)
        t is the probe point, used to search in the group of points Insert'ed
        return value is true if some object was found within the search radius, false otherwise

     bool FarthestNeighbor ( T& tFarthest,   const T& t ) const
        tFarthest is returned as the object that was found farthest to the probe
           point
        t is the probe point, used to search in the group of points Insert'ed
        return value is true if some object was found, false otherwise

     long FindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
        dRadius is the radius within which to search; make it very large if you want to
            include every point that was loaded;
        tClosest is returned as the vector of objects that were found within a radius dRadius
           of the probe point
        t is the probe point, used to search in the group of points Insert'ed
        return value is the number of objects found within the search radius
    
     T InSphere ( const double& dRadius,  const T& t )
       dRadius is the radius within which to search; make it very large if you want to
            include every point that was loaded;
       t is the probe point, used to search in the group of points Insert'ed
       return value is an object of the type v which is closest within the sphere to
            the probe point t
      

   So a complete program is:

  #include "TNear.h"
  #include <cstdio>
  void main()
  {
    CNearTree< double > dT;
    double dNear;
    dT.Insert( 1.5 );
    if ( dT.NearestNeighbor( 10000.0,   dNear,  2.0 )) printf( "%f\n",dRad );
  }

   and it should print 0.5 (that's how for 2.0 is from 1.5). For more
   examples of the use of TNear.h, see main.cpp and CNearTreeTest.cpp.

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

    The C NearTree API CNearTree.c

    Synopsis

   #include <CNearTree.h>

   double CNearTreeDist (CNearTreeHandle treehandle, void * coord1, void *
   coord2);

   int CNearTreeCreate (CNearTreeHandle * treehandle, size_t treedim, int
   treetype);

   int CNearTreeFree(CNearTreeHandle * treehandle);

   int CNearTreeInsert( CNearTreeHandle treehandle, const void * coord, const
   void * obj );

   int CNearTreeZeroIfEmpty (CNearTreeHandle treehandle);

   int CNearTreeNearestNeighbor (CNearTreeHandle treehandle, const double
   dRadius, void * * coordClosest, void * * objClosest, const void * coord );

   int CNearTreeFarthestNeighbor (CNearTreeHandle treehandle, void * *
   coordFarthest, void * * objFarthest, const void * coord );

   int CNearTreeFindInSphere ( CNearTreeHandle treehandle, const double
   dRadius, CVectorHandle coords, CVectorHandle objs, const void * coord, int
   resetcount);

   int CNearTreeNearest ( CNearTreeHandle treehandle, double * dRadius, void
   * * coordClosest, void * * objClosest, const void * coord );

   int CNearTreeFindFarthest ( CNearTreeHandle treehandle, double * dRadius,
   void * * coordFarthest, void * * objFarthest, const void * coord );

   The NearTree API works with coordinate vectors in an arbitrary number of
   dimensions. Each neartree is accessed by a pointer of type CNearTreeHandle
   which points to a struct of type CNearTree:

     typedef struct _CNearTree {
         void           * m_coordLeft;    /* coords of left object stored in this node   */
         void           * m_coordRight;   /* coords of right object stored in this node  */
         double           m_dMaxLeft;     /* longest distance from the left object
                                             to anything below it in the tree            */
         double           m_dMaxRight;    /* longest distance from the right object
                                             to anything below it in the tree            */
         struct _CNearTree FAR * m_pLeftBranch; 
                                          /* tree descending from the left object        */
         struct _CNearTree FAR * m_pRightBranch;
                                          /* tree descending from the right object       */
         void FAR *       m_ptobjLeft;    /* pointer to the left object                  */
         void FAR *       m_ptobjRight;   /* pointer to the right object                 */
         size_t           m_dimension;    /* dimension of the coordinates                */
         int              m_flags;        /* flags                                       */
         /* NOTE:  in actual blocks the coordinates will be stored here after the flags         */
     } CNearTree;


   The internal operation of the API depends on the function CNearTreeDist
   that returns the distance (L1, L2 or L-infinity) between two coordinate
   vectors as a double according to the parameters of the given tree. Note
   that the tree may store the coordinates as integers or as doubles, but the
   distance is always computed as a double. If this function is replaced by a
   user function, it is important that the replacement obey the triangle
   inequality.

   A neartree is created by CNearTreeCreate and freed by CNearTreeFree.
   treedim is the dimension of the coordinate vectors and treetype is one of
   the two predefined constants CNEARTREE_TYPE_DOUBLE for double or
   CNEARTREE_TYPE_INTEGER for integer, optionally ORed with
   CNEARTREE_NORM_L1, CNEARTREE_NORM_L2 or CNEARTREE_NORM_LINF for L1, L2 or
   L-infinity norms. When first created, a neartree has no right or left node
   and with the dMax-below set to negative values so that any match found
   will be stored since it will greater than the negative value. The tree is
   then populated by calls to CNearTreeInsert, with each call providing a
   coordinate vector coord and an optional object pointer obj. The API copies
   the coordinate vector, but does not copy the object. The tree is populated
   in the order left, right and then the nearer child.

   The neartree is searched for the nearest or farthest coordinate vector in
   the neartree to a given probe coordinate vector coord by
   CNearTreeNearestNeighbor and CNearTreeFarthestNeighbor, respectively. The
   given radius contains the search to a sphere around the probe. If more
   than a single extremal coordinate point is needed, CNearTreeFindInSphere
   can be used to obtain a CVector result vector of all the coordinate
   vectors that satisfy the constraint. The vectors themselves are not copied
   into the result vector. If the parameter resetcount is true (non zero) the
   result vector is cleared before the search. A parallel CVector result
   vector of the matching object pointers is returned is objs is not NULL.
   The functions CNearTreeNearest and CNearTreeFindFarthest implement
   CNearTreeNearestNeighbor and CNearTreeFarthestNeighbor, respectively,
   adjusting the radius of the search while the searh is in progress, and are
   not normall used by users

    Returns

   If CNearTreeDist fails, it returns -1. Except for CNearTreeDist, all the
   functions in the API return 0 (CNEARTREE_SUCCESS) for success. If dynamic
   memory allocation fails, CNEARTREE_MALLOC_FAILED is returned. If a call is
   made with an improper argument, CNEARTREE_BAD_ARGUMENT is returned. If a
   search fails to find anything, CNEARTREE_NOT_FOUND is returned.

    Examples

   To create a neartree for 3-dimensional vectors of doubles:

 
 #include <CNearTree.h>
 CNearTreeHandle treehandle;
 int bReturn;
 
  ...
 
  bReturn = !CNearTreeCreate(&treehandle,3,CNEARTREE_TYPE_DOUBLE);

   To insert a copy of a 3-dimensional vector of doubles into this tree, with
   no associated object:

 
     double v[3];
 
  ...
     
     v[0] = 1.; v[1] = 2.; v[2] = 3.;
     bReturn = !CNearTreeInsert(treehandle,&v[0],NULL);

   To search for the nearest neighbor to a probe vector vSearch in a radius
   of 3., returning a pointer to the resulting vector in vBest:

 
     double * vBest;
     void * vvBest;
     double vSearch[3];
     double   dRad = =3.;
    
   ...
  
     if ( !CNearTreeNearestNeighbor(treehandle,dRad,&vvBest,NULL,vSearch))
         {   vBest = (double *)vvBest; }

   Note the use of a seaparate void * vvBest instead of a cast of &vBest to
   avoid compiler type punning warnings.

   For more examples of the use of CNearTree.c, see main.c and
   CNearTreeTest.c in the release kit.
