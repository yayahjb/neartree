                                    NearTree

                                 Release 2.1.1
                                  4 June 2009
       (c) Copyright 2001, 2008, 2009 Larry Andrews. All rights reserved
                                    based on
         Larry Andrews, "A template for the nearest neighbor problem",
   C/C++ Users Journal, Volume 19 , Issue 11 (November 2001), 40 - 49 (2001),
                                ISSN:1075-2838,
                        www.ddj.com/architect/184401449

   Revised 12 Dec 2008, for sourceforge release, Larry Andrews and Herbert J.
                                   Bernstein
                       8 Jan 2009 Release 1.0 LCA and HJB
                     11 Jan 2009 Release 1.0.1 LCA and HJB
                     21 March 2009 Release 2.0 LCA and HJB
                      30 May 2009 Release 2.1 LCA and HJB
                     4 June 2009 Release 2.1.1 LCA and HJB

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

   This is a cleanup update to the 2.1 release of 30 May 2009 to increase
   portability, dealing with the following issues:
     * Convert to use of a self-contained portable random-number generator
       from Rob Harrison
     * Ensure wider use of const where appropriate

   The 2.1 release was a minor update to the 2.0 release of 21 March 2009 to
   deal with the following issues:
     * Make delayed insertion the default
     * Complete the containerization of TNear.h
     * Add code for K-nearest/farthest in TNear.h and in CNearTree.c
     * Correct the InAnnulus search filter

   Rlease 2.0 was a major update to the 1.0 release of 8 January 2009 to deal
   with the following issues:
     * Replace use recursion with a stack, except in insertion logic
     * Replace use of double with templated DistanceType (usually double)
     * Provide constuctors to build NearTree from vectors, lists or sets
     * Change "Insert" to "insert" for consistentcy with other containers
     * Add access function "at" or array type references [], and provide
       contents of a neartree as a vector
     * Add iterator support
     * Provide delayed insertion logic
     * Functions added for searches outside of a sphere or in an annular
       region

   Our thanks to Nicolas Brodu for suggesting the more general handling of
   the distance type.

   Note: As Nicolas Brodu has noted, CNearTree is paticularly well-suited to
   multi-threaded applications. However, if the same CNearTree is to be
   searched in multiple threads, it is important to complete all insertions
   and/or delayed insertions before parallel execution of parallel searches.
   Delayed insertions are completed at the beginning if all searches, so for
   normal use, this will not be a problem.

     ----------------------------------------------------------------------

    Contents

     * Installation
     * The C++ template: TNear.h
     * The C NearTree API: CNearTree.c
     * A Portable pseudo-random number generator: rhrand.h

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

    Installation

   The NearTree package is available at
   www.sourceforge.net/projects/neartree. A source tarball is available at
   downloads.sourceforge.net/neartree/NearTree-2.1.1.tar.gz. Later tarballs
   may be available.

   When the source tarball is dowloaded and unpacked, you should have a
   directory NearTree-2.1.1. To see the current settings for a build execute

   make

   which should give the following information:

  PLEASE READ README_NearTree.txt and lgpl.txt
 
  Before making the NearTree libraries and example programs, check
  that the chosen settings are correct
 
  The current C++ and C compile commands are:
 
    libtool --mode=compile g++ -g -O2  -Wall -ansi -pedantic  \
       -DCNEARTREE_SAFE_TRIANG=1 -I.   -c
    libtool --mode=compile gcc -g -O2  -Wall -ansi -pedantic  \
       -DCNEARTREE_SAFE_TRIANG=1 -I.   -c
 
  The C API, CNearTree.c, depends on the sourceforge project CVector
  You are currently setup to use the system defaults for CVector
  If that is not correct, define the variable CVECTOR_INCLUDE

  The current library link command is:
 
    libtool --mode=link  gcc -version-info 2:0:2 -release 2.0.0 \
      -no-undefined -rpath /usr/local/lib
 
  The current C++ and C library local, and C dynamic and static build commands are:
 
    libtool --mode=link g++ -no-undefined -g -O2  -Wall -ansi -pedantic \
       -DCNEARTREE_SAFE_TRIANG=1 -I.
    libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic \
       -DCNEARTREE_SAFE_TRIANG=1 -I.
    libtool --mode=link gcc -no-undefined -g -O2  -Wall -ansi -pedantic \
       -DCNEARTREE_SAFE_TRIANG=1 -shared -I/usr/local/include
    libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic  \
       -DCNEARTREE_SAFE_TRIANG=1 -static-libtool-libs -I/usr/local/include
 
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
   libtool versions 1.3.5 and 1.5.4. For MINGW, libtool version 2.2.6 and gcc
   version 4 are needed to work with shared libraries (DLLs). If the system
   libtool is not to be used, define the variable LIBTOOL to be the path to
   the libtool executable, e.g. in bash

   export LIBTOOL=$HOME/bin/libtool

   or in the Makefie

   LIBTOOL = $(HOME)/bin/libtool

   If you need to include local header files using #include "..." instead of
   #include <...>, define the variable USE_LOCAL_HEADERS

   There is an optimization error in some compilers that miscompiles the
   triangle inequality. The default in this API is to do the triangle
   inequality three different ways under the control of CNEARTREE_SAFE_TRIANG

 #ifdef CNEARTREE_SAFE_TRIANG
 #define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0) \
                       || ((b)-((a)-(c)) >= 0) \
                       || ((c)-((a)-(b)) >= 0))   
 #else
 #define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0))
 #endif

   Problems with the unsafe definition of TRIANG have been seen in Linux
   under gcc version 4 and in MS Window under VS 2003. There is a slight
   performance hit from the triple test. If maximal speed is critical and
   misidentification of nearest points by relative distance errors of about 1
   part in 10**15 is not a serious problem, the definition of
   -DCNEARTREE_SAFE_TRIANG=1 can be removed from the definition of CFLAGS in
   the Makefile.

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

    The C++ template: TNear.h

   This is a revised release of

     template <typename T, typename DistanceType=double, int distMinValue=-1
     > class CNearTree;

   implementing the Nearest Neighbor algorithm after Kalantari and McDonald,
   (IEEE Transactions on Software Engineering, v. SE-9, pp. 631-634,1983)
   modified to use recursion for insertions and recursion (original version)
   or a stack (current version) for searches instead of a double-linked tree
   and simplified so that it does a bit less checking for things like is the
   distance to the right less than the distance to the left; it was found
   that these checks made little to no difference.

   This template is used to contain a collection of objects. After the
   collection has been loaded into this structure, it can be quickly queried
   for which object is "closest" to some probe object of the same type. The
   major restriction on applicability of the near-tree is that the algorithm
   only works if the objects obey the triangle inequality. The triangle rule
   states that the length of any side of a triangle cannot exceed the sum of
   the lengths of the other two sides.

   CNearTree is the root class for the neartree. The actual data of the tree
   is stored in NearTreeNode objects descending from a CNearTree.

   The types of objects that can be stored in the tree is quite broad. The
   biggest limitation is that the objects must reside in some sort of metric
   space and must obey the triangle rule. They must also be all of the same
   size because they are stored in an std::vector. If your application
   requires object of varying storage, then your best way to use this code is
   to store pointers or handles and to write your own distance functions.

   The type of the objects to be stored is the only required template
   argument. The type of the distance measure (DistanceType) defaults to
   double. If your application is for an integer type then the type for
   DistanceType can be your integer type. This has the potential for speeding
   the calculations by avoiding FP computation. Other general types can be
   used if desired, but you may need to also input a value of distMinValue.

   The template argument distMinValue must be something that your class will
   understand as a negative number. The default input is negative one.
   Internally, that is cast to DistanceType. Since most uses will be for
   DistanceType to be double, that is a simple conversion. Obviously, for
   integer types, there is no problem either. The need for this value is to
   have something internally that is recognizable as smaller than the
   smallest "distance" that can exist between any two objects in your type.
   For most users, there is no need to input anything other than the default,
   -1. -1 must be castable to DistanceType. It seems unlikely that anyone
   would actually need this optional parameter, but it is here for
   completeness.

   It is a design decision that this class cannot work for unsigned types. It
   is hard to see how to verify the triangle rule for unsigned types, and
   distance computations become more complex. Sorry, unsigned types are left
   as an exercise for the reader.

   The user of this class needs to provide at least the following
   functionality for the template to work. For the built-in numerics of C++,
   they are provided here (or else you should create them).

                      // conversion constructor from the templated class to      
DistanceType Norm( ); DistanceType                                               
                      // (usually will return a "length" of type double)         
operator- ( );        // geometrical (vector) difference of two objects          
                      // a copy constructor would be nice                        
                      // a constructor would be nice                             
                      // a destructor would be nice                              

   The provided interface is:

     #include <TNear.h>

     CNearTree( void )   // constructor
        instantiated by something like:      CNearTree <T> vTree;
        for some type T
       
     CNearTree( const ContainerType<T> )   // constructor from containers, e.g. ...

     CNearTree( const std::vector<T> )     // constructor
     CNearTree( const std::list<T> )       // constructor
     CNearTree( const std::set<T> )        // constructor
    
     ~CNearTree( void )  // destructor

     void Insert( const T& t )
        where t is an object of the type T

             all inserts are delayed until a search is performed or
             until an explicit call to CompleteDelayedInsertions
             is called or a search is called. The purpose is to distribute
             the objects a bit more  randomly. Excessively ordered objects
             leads to less than optimal trees.
      
             Starting with the 2.1 release, places objects in a queue for
             insertion later when  CompleteDelayInsert is called.  In
             earlier releasesm the default was immediate insertion.
            
             The following additional convenience insert template
             allow insertion of containers of objects


     void Insert( ContainerType ) 
             for containers, std::vector, ..., or CNearTree
             all inserts are delayed until a search is performed or
             until an explicit call to CompleteDelayedInsertions.
             examples follow ...


     void insert( const std::vector<T> )
     void insert( const std::list<T> )
     void insert( const std::set<T> )

     bool NearestNeighbor ( const DistanceType& dRadius,  T& tClosest,   const T& t ) const
        dRadius is the largest radius within which to search; make it
           very large if you want to include every point that was loaded.
           dRadius is returned as the closest distance to the probe (or the search radius
           if nothing is found)
        tClosest is returned as the object that was found closest to the probe
           point (if any were within radius dRadius of the probe)
        t is the probe point, used to search in the group of points insert'ed

        return value is true if some object was found within the search radius, false otherwise.
            If false is returned, tClosest is invalid (at best).

     iterator NearestNeighbor ( const DistanceType & dRadius, const T& t ) const
        returns an iterator to the nearest point to the probe point t or end() if there is none

     bool FarthestNeighbor ( T& tFarthest,   const T& t ) const
        tFarthest is returned as the object that was found farthest to the probe
           point
        t is the probe point, used to search in the group of points Insert'ed
        return value is true if some object was found, false otherwise
           If false is returned, tFarthest is invalid (at best).

     iterator FarthestNeighbor ( T& const T& t ) const
        returns an iterator to the nearest point to the probe pointt or end() if there is none


     The following functions (FindInSphere, FindOutSphere, and FindInAnnulus) all return a container
     (ContainerType) that can be any standard library container (such as std::vector< T >) or CNearTree.

     long FindInSphere ( const DistanceType& dRadius,  ContainerType& tInside,   const T& t ) const
        dRadius is the radius within which to search; make it very large if you want to
            include every point that was loaded;
        tInside is returned as the container of objects that were found within a radius dRadius
           of the probe point
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the search radius

     long FindInSphere ( const DistanceType& dRadius,  CNearTree<  T >& tInside,   const T& t ) const
        dRadius is the radius within which to search; make it very large if you want to
            include every point that was loaded;
        tInside is returned as the CNearTree of objects that were found within a radius dRadius
           of the probe point
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the search radius

     long FindOutSphere ( const DistanceType& dRadius,  ContainerType& tOutside,   const T& t ) const
        dRadius is the radius outside of which to search
        tOutside is returned as the container of objects that were found at or outside of radius dRadius
           of the probe point
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the search radius

     long FindOutSphere ( const DistanceType& dRadius,  CNearTree<  T >& tOutside,   const T& t ) const
        dRadius is the radius outside of which to search.
        tClosest is returned as the CNearTree of objects that were found at or outside of a radius dRadius
           of the probe point
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the search radius
    
     long FindInAnnulus ( const DistanceType& dRadius1, const DistanceType& dRadius2, ContainerType& tInRing,   const T& t ) const
        dRadius1 and  dRadius2 are the two radii between which to find data points
        tInRing is returned as the container of objects that were found at or outside of a radius dRadius1
           and at or inside of radius dRadius2 of the probe point
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the annulus

     long FindInAnnulus ( const DistanceType& dRadius1, const DistanceType& dRadius2,  CNearTree<  T >& tInRing,   const T& t ) const
        dRadius1 and  dRadius2 are the two radii between which to find data points
        tInRing is returned as the CNearTree of objects that were found at or outside of a radius dRadius1
           and at or inside of radius dRadius2 of the probe point
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the annulus
    
     long FindK_NearestNeighbors ( const size_t k, const DistanceType& dRadius, ContainerType& tClosest, const T& t )
        k is the maximum number of nearest neighbors to return. Finds this many if possible
        dRadius within a sphere defined by dRadius, search for the k-nearest-neighbors
        tClosest is returned ContainerType of the objects found within the sphere
        t is the probe point, used to search in the group of points insert'ed

        return value is the count of the number of points found within the sphere

     long FindK_NearestNeighbors ( const size_t k, const DistanceType& dRadius, CNearTree<  T >& tClosest, const T& t )
        k is the maximum number of nearest neighbors to return. Finds this many if possible
        dRadius within a sphere defined by dRadius, search for the k-nearest-neighbors
       tClosest is returned as the CNearTree of objects found within the sphere
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the sphere

     long FindK_FarthestNeighbors ( const size_t k, ContainerType& tFarthest, const T& t )
        k is the maximum number of farthest neighbors to return. Finds this many if possible
        tFarthest is returned ContainerType of the objects found
        t is the probe point, used to search in the group of points insert'ed

        return value is the count of the number of points found within the sphere

     long FindK_NearestNeighbors ( const size_t k, CNearTree<  T >& tFarthest, const T& t )
        k is the maximum number of farthest neighbors to return. Finds this many if possible
        tFarthest is returned as the CNearTree of objects found
        t is the probe point, used to search in the group of points Insert'ed

        return value is the count of the number of points found within the sphere
   

     ----------------------------------------------------------------------


 Access Functions:

      T at ( const size_t n ) const
         returns the n'th item of the internal data store

      T operator[] ( const size_t n )
         returns the n'th item of the internal data store

      operator ContainerType ( void ) const
         returns all of the inserted objects in the tree in a container of type ContainerType.
         ContainerType can be std::vector<T>, etc, or other containers.
         The returned vector contents are not guaranteed to be returned in the order loaded.

      iterator begin ( void ) const
         returns an iterator to the beginning of the internal data store

      iterator end ( void ) const
         returns an iterator to the end of the data store (one beyond the last item)

      iterator back ( void ) const
         returns an iterator to the last data item of the internal data store


     ----------------------------------------------------------------------


 Information and special operation functions:

      void ImmediateInsert( void ) Places objects immediately into the tree. The usual insert function
          delays insertions, allowing them to be inserted into the tree in a more random order. The delay
          can improve the structure of the tree and speed searches.

      void CompleteDelayedInsert ( void )
         calls insert for all delayed objects. sqrt(n) of them are inserted by random choice.
         The rest are inserted in linear order as originally queued. CompleteDelayedInsert
         is invoked at the beginning of all searches, so the average user will never need
         to call it.

      size_t GetDeferredSize ( void )
         returns the number of delayed objects that have not yet been insert'ed. This is
         mainly for information about details of the tree.

      size_t GetTotalSize ( void )
         returns the number of objects that have been insert'ed plus those DelayInsert'ed

      size_t size ( void )
         identical to GetTotalSize

      size_t GetDepth ( void )
         returns the maximum tree layers from the root.  This is mainly for information about
         details of the tree.

      bool empty ( void )
         returns true if the tree is empty, otherwise false
        

     ----------------------------------------------------------------------


 Iterators:

      Random access iterators are provided for accessing the data in a CNearTree. The most important
      expected use is to retrieve the objects returned from one of the sphere search functions that
      return a CNearTree. However, they can be used with any CNearTree.
     
      They should function in a fashion essentially the same as STL iterators. There is no assurance
      that data will be returned in the order it was loaded, just that it is accessible.  The same set is
      provided for const_iterator.

      iterator ( void ) { }; // constructor

      iterator& operator=   ( const iterator& s )     
      iterator  operator++  ( const int n )           
      iterator  operator--  ( const int n )           
      iterator& operator++  ( void )                  
      iterator& operator--  ( void )                  
      iterator  operator+   ( const long n ) const    
      iterator  operator-   ( const long n ) const    
      iterator& operator+=  ( const long n )          
      iterator& operator-=  ( const long n )          
      T         operator*   ( void )        const    

      bool      operator==  ( const iterator& t ) const
      bool      operator!=  ( const iterator& t ) const

      const T * const operator->  ( void )   const


     ----------------------------------------------------------------------

  

   So a complete program is:

  #include "TNear.h"
  #include <cstdio>
  void main()
  {
    CNearTree< double > dT;
    double dNear;
    dT.Insert( 1.5 );
    if ( dT.NearestNeighbor( 10000.0,   dNear,  2.0 )) printf( "%f\n",double(dNear-2.0) );
  }

   and it should print 0.5 (that's how for 2.0 is from 1.5). For more
   examples of the use of TNear.h, see main.cpp and CNearTreeTest.cpp.

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

    The C NearTree API: CNearTree.c

    Synopsis

     #include <CNearTree.h>

     double CNearTreeDist ( CNearTreeHandle treehandle, void * coord1, void *
     coord2 );

     int CNearTreeSetNorm ( const CNearTreeHandle treehandle, int treenorm );

     int CNearTreeNodeCreate ( const CNearTreeHandle treehandle,
     CNearTreeNodeHandle * treenodehandle )

     int CNearTreeCreate ( CNearTreeHandle * treehandle, size_t treedim, int
     treetype );

     int CNearTreeFree ( const CNearTreeHandle treehandle );

     int CNearTreeClear ( CNearTreeHandle * treehandle );

     int CNearTreeNodeFree ( CNearTreeNodeHandle * treenodehandle );

     int CNearTreeInsert( const CNearTreeHandle treehandle, const void *
     coord, const void * obj );

     int CNearTreeNodeInsert ( const CNearTreeHandle treehandle,
     CNearTreeNodeHandle treenodehandle, size_t index; size_t * depth );

     int CNearTreeImmediateInsert ( const CNearTreeHandle treehandle, const
     void * coord, const void * obj );

     int CNearTreeDelayedInsert ( const CNearTreeHandle treehandle, const
     void * coord, const void * obj ); /* ***DEPRECATED*** */

     int CNearTreeCompleteDelayedInsert ( const CNearTreeHandle treehandle )
     int CNearTreeZeroIfEmpty ( const CNearTreeHandle treehandle );

     int CNearTreeGetSize ( const CNearTreeHandle treehandle, size_t * size
     );

     int CNearTreeGetDelayedSize ( const CNearTreeHandle treehandle, size_t *
     size );

     int CNearTreeGetDepth ( const CNearTreeHandle treehandle, size_t * depth
     )

     int CNearTreeNearestNeighbor ( const CNearTreeHandle treehandle, const
     double dRadius, void * * coordClosest, void * * objClosest, const void *
     coord );

     int CNearTreeFarthestNeighbor ( const CNearTreeHandle treehandle, void *
     * coordFarthest, void * * objFarthest, const void * coord );

     int CNearTreeFindInSphere ( const CNearTreeHandle treehandle, const
     double dRadius, CVectorHandle coordInside, CVectorHandle objInside,
     const void * coord, int resetcount );

     int CNearTreeFindTreeInSphere ( const CNearTreeHandle treehandle, const
     double dRadius, CNearTreeHandle foundInside, const void * coord, int
     resetcount )

     int CNearTreeFindOutSphere ( const CNearTreeHandle treehandle, const
     double dRadius, CVectorHandle coordOutside, CVectorHandle objOutside,
     const void * coord, int resetcount );

     int CNearTreeFindTreeOutSphere ( const CNearTreeHandle treehandle, const
     double dRadius, CNearTreeHandle foundOutside, const void * coord, int
     resetcount )

     int CNearTreeFindInAnnulus ( const CNearTreeHandle treehandle, const
     double dRadiusInner, const double dRadiusOuter, CVectorHandle
     coordInRing, CVectorHandle objInRing, const void * coord, int resetcount
     );

     int CNearTreeFindTreeInAnnulus ( const CNearTreeHandle treehandle, const
     double dRadiusInner, const double dRadiusOuter, CNearTreeHandle
     foundInRing, const void * coord, int resetcount )

     int CNearTreeFindKNearest ( const CNearTreeHandle treehandle, const
     size_t k, const double dRadius, CVectorHandle coordClosest,
     CVectorHandle objClosest, const void * coord, int resetcount );

     int CNearTreeFindKTreeNearest ( const CNearTreeHandle treehandle, const
     size_t k, const double dRadius, CNearTreeHandle foundClosest, const void
     * coord, int resetcount )

     int CNearTreeFindKFarthest ( const CNearTreeHandle treehandle, const
     size_t k, const double dRadius, CVectorHandle coordFarthest,
     CVectorHandle objFarthest, const void * coord, int resetcount );

     int CNearTreeFindKTreeFarthest ( const CNearTreeHandle treehandle, const
     size_t k, const double dRadius, CNearTreeHandle foundFarthest, const
     void * coord, int resetcount )

     int CNearTreeNearest ( const CNearTreeHandle treehandle, double *
     dRadius, void * * coordClosest, void * * objClosest, const void * coord
     );

     int CNearTreeFindFarthest ( const CNearTreeHandle treehandle, double *
     dRadius, void * * coordFarthest, void * * objFarthest, const void *
     coord );

     int CNearTreeObjects ( const CNearTreeHandle treehandle, CVectorHandle *
     vectorhandle );

     void * CNearTreeObjectAt ( const CNearTreeHandle treehandle, size_t
     index );

     int CNearTreeCoords ( const CNearTreeHandle treehandle, CVectorHandle *
     vectorhandle );

     void * CNearTreeCoordAt ( const CNearTreeHandle treehandle, size_t index
     );

   The NearTree API works with coordinate vectors in an arbitrary number of
   dimensions. Each neartree is accessed by a pointer of type CNearTreeHandle
   which points to a struct of type CNearTree, which points to a tree of
   nodes of type CNearTreeNode:

     typedef struct _CNearTreeNode {
         size_t           m_indexLeft;    /* index of left coords in m_CoordStore 
                                             and of left object in m_ObjectStore     */
         size_t           m_indexRight;   /* index of right coords in m_CoordStore
                                             and of right object in m_ObjectStore     */
         double           m_dMaxLeft;     /* longest distance from the left object
                                             to anything below it in the tree            */
         double           m_dMaxRight;    /* longest distance from the right object
                                             to anything below it in the tree            */
         struct _CNearTreeNode * m_pLeftBranch; 
                                          /* tree descending from the left object        */
         struct _CNearTreeNode * m_pRightBranch;
                                          /* tree descending from the right object       */
         int              m_iflags;       /* flags                                       */
     } CNearTreeNode;
    
    
     typedef CNearTreeNode * CNearTreeNodeHandle;  
    
     typedef struct {
         CNearTreeNodeHandle m_ptTree;     /* pointer to the actual tree                  */
         size_t           m_szdimension;   /* dimension of the coordinates                */
         size_t           m_szsize;        /* size of this tree                           */
         size_t           m_szdepth;       /* depth of this tree                          */
         int              m_iflags;        /* flags                                       */
         CVectorHandle    m_ObjectStore;   /* all inserted objects                        */
         CVectorHandle    m_CoordStore;    /* all inserted coordinates                    */
         CVectorHandle    m_DelayedIndices;/* objects queued for insertion                */
     } CNearTree;
    
     typedef CNearTree     FAR * CNearTreeHandle;


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
   L-infinity norms. Starting with release 2.1, all insertions are delayed by
   default, unless the insertions is done by a call to
   CNearTreeImmediateInsert. The insertions that have been queued are
   completed by a call to CNearTreeCompleteDelayedInsert or by any search.
   The insertions are actually done in a randomized order, either for an
   initial block of sqrt(#queue) by default, or for the entire queue if the
   flag CNEARTREE_DEFER_ALL is ored with treetype. In addition, if the flag
   CNEARTREE_DEFER_ALL is set, all inserts will be deferred, even in calls to
   CNearTreeImmediateInsert. Alteratively, if minimal reorganization of the
   order of the tree is desired on each insertion, CNEARTREE_FLIP should be
   ored with treetype. When first created, a neartree has no right or left
   node and with the dMax-below set to negative values so that any match
   found will be stored since it will greater than the negative value. The
   tree is then populated by calls to CNearTreeInsert, with each call
   providing a coordinate vector coord and an optional object pointer obj.
   The API copies the coordinate vector, but does not copy the object. Later,
   when a search is requested or an explicit call to
   CNearTreeCompleteDelayedInsert is made, the tree is populated in the order
   left, right and then the nearer child, working from a randomized selection
   from the items queued for insertion.

   Optionally, the actual insertions may done immediately by calling
   CNearTreeImmediateInsert instead of CNearTreeInsert. For upwards
   compatability of the library for existing code, the deprecated
   CNearTreeDelayedInsert is provided as an deprecated alternate call to
   CNearTreeInsert.

   The neartree is searched for the nearest or farthest coordinate vector in
   the neartree to a given probe coordinate vector coord by
   CNearTreeNearestNeighbor and CNearTreeFarthestNeighbor, respectively. The
   given radius confines the search to a sphere around the probe. If more
   than a single extremal coordinate point is needed, CNearTreeFindInSphere
   can be used to obtain a CVector result vector of all the coordinate
   vectors that satisfy the constraint of being within a specified radius, or
   CNearTreeFindOutSphere can be used to obtain a CVector result vector of
   all the coordinates that satisfy the constraint of being outside a
   specified radius. CNearTreeFindIn Annulus can be used to obtain a CVector
   result vector of all the coordinates that satisfy the constraint of being
   between two specified radii from the probe. CNearTreeFindKNearest can be
   used to obtain a CVector result vector of the k coordinates closest to the
   probe point such that all results are within the specified radius of the
   probe point, or CNearTreeFindKFarthest to obtain a CVector result vector
   of the k coordinates farthest to the probe point such that all results are
   at or outside the specified radius of the probe point. The vectors
   themselves are not copied into the result vector. If the parameter
   resetcount is true (non zero) the result vector is cleared before the
   search. A parallel CVector result vector of the matching object pointers
   is returned if objs is not NULL. Aternatives the forms
   CNearTreeFindTreeInSphere, CNearTreeFindTreeOutSphere,
   CNearTreeFindTreeInAnnulus, CNearTreeFindKTreeNearest,
   CNearTreeFindKTreeFarthest can be used to obtain CNearTrees rather than
   CVectors of results. The functions CNearTreeNearest and
   CNearTreeFindFarthest implement CNearTreeNearestNeighbor and
   CNearTreeFarthestNeighbor, respectively, adjusting the radius of the
   search while the search is in progress, and are not normally used by users

    Returns

   If CNearTreeDist fails, it returns -1. Except for CNearTreeDist, all the
   functions in the API return 0 ( CNEARTREE_SUCCESS ) for success. If
   dynamic memory allocation fails, CNEARTREE_MALLOC_FAILED is returned. If a
   call is made with an improper argument, CNEARTREE_BAD_ARGUMENT is
   returned. If a search fails to find anything, CNEARTREE_NOT_FOUND is
   returned. If there is a failure in an attempt to free a CNearTree,
   CNEARTREE_FREE_FAILED is returned. If any of the internal call to CVector
   fail, CNEARTREE_CVECTOR_FAILED is returned. For convenience in debugging,
   the formerly negative values of this returns are now positive.

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

   Note the use of a separate void * vvBest instead of a cast of &vBest to
   avoid compiler type punning warnings.

   For more examples of the use of CNearTree.c, see main.c and
   CNearTreeTest.c in the release kit.

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

    A Portable pseudo-random number generator: rhrand.h

   rhrand.h is a portable pseudo-random number generator based one by Rob
   Harrison, derived from "one in J.M.Hammersley and D.C. Handscomb, 'Monte
   Carlo Methods,' Methuen & Co., London and Wiley & Sons, New York, 1964,
   p47". See also, D. E. Knuth "The Art of Computer Programming", Volume 2,
   "Seminumerical Alogorithms, Third Edition, Addison-Wesley, Reading MA,
   1997.

   rhrand.h is a header file in which a C++ class, RHrand, is defined, and a
   C struct typedef CRHrand is defined.

   The C++ interface is

     static const int RHRAND_MAX = 32767;  /* the integer range accessible as RHrand::RHRAND_MAX */
    
     RHrand(void)                          /* the default constructor */
    
     RHrand( const int iseed )             /* a constructor to start with the given seed */
    
     ~RHrand( void)                        /* a destructor */
    
     void srandom( const int iseed)        /* reset the generator based on the given seed */
    
     double urand( void )                  /* return a random double uniformly distributed in [0,1) */
    
     int random ( void )                   /* return a random integer uniformly distributed in [0, RHRAND_MAX-1] */


   In C++ code, typical use is

 #include <rhhand.h>
     RHrand rhr;

 ...

     x = rhr.urand();

   The C interface is suppressed in RHRAND_NOCCODE is defined. Otherwise the
   C interface is based on defining a struct of type CRHRrand and calling
   macros that refer to a handle of type RCRHrandHandle.

     typedef struct CRHrand_ {           /* the struct used in random number generattion */
         double buffer[55];
         int indx;
         int jndx;
         int kndx;
         double dTemp;
     } CRHrand;

     typedef CRHrand * CRHrandHandle;     /* the type to be used in maro calls */

     #define CRHRAND_MAX 32767            /* the integer range */
    
     #define CRHrandSrandom(randhandle,iseed) ... 
                                          /* a macro to call to initialize CHRrandHandle randhandle
                                             using see int iseed */
                                            
     #define CRHrandUrand(randhandle) ... /* a macro to return a random double uniformly distributed in [0,1) */
    
     #define CRHrandRandom(randhandle) ((int)(CRHrandUrand(randhandle)*(double)CRHRAND_MAX))
                                          /* a macro to return a random integer uniformly distributed in
                                             [0, CRHRAND_MAX-1] */

   Typical use is

 #include <rhhand.h>
     CRHrand rhr;

 ...

     CRHrandSrandom(&rhr, 0 );

 ...

     x = CRHrandUrand(&rhr);

     ----------------------------------------------------------------------

     ----------------------------------------------------------------------

   Updated 4 June 2009
   yaya@bernstein-plus-sons.com
