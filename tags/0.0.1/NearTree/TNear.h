//*
//*  TNear.h
//*  NearTree
//*
//*  Copyright 2001, 2008 Larry Andrews.  All rights reserved
//*  Revised 12 Dec 2008 for sourceforge release -- H. J. Bernstein


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

//  This is a revised release of 
//  template <typename T> class CNearTree;
//
// Nearest Neighbor algorithm after Kalantari and McDonald,
// (IEEE Transactions on Software Engineering, v. SE-9, pp.
//    631-634,1983)
//  modified to use recursion instead of a double-linked tree
//  and simplified so that it does a bit less checking for
//  things like is the distance to the right less than the
//  distance to the left; it was found that these checks little
//  to no difference.



// The user of this class needs to provide at least the following
// functionality for the template to work. For the built-in
// numerics of C++, they are provided by the system.

//    operator double( );   // conversion constructor from the templated class to double
//                            (usually will return a "length")
//    operator- ( );        // geometrical (vector) difference of two objects
//    a copy constructor
//    a constructor would be nice
//    a destructor would be nice

// The provided interface is:
//
//    #include "TNear.h"
//
//    CNearTree( void )   // constructor
//       instantiated by something like:      CNearTree <v> vTree;
//       for some type v
//
//    void m_fnInsert( T& t )
//       where t is an object of the type v
//
//    bool m_bfnNearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//       dRadius is the largest radius within which to search; make it
//          very large if you want to include every point that was loaded; dRadius
//          is returned as the closest distance to the probe (or the search radius
//          if nothing is found)
//       tClosest is returned as the object that was found closest to the probe
//          point (if any were within radius dRadius of the probe)
//       t is the probe point, used to search in the group of points m_fnInsert'ed
//       return value is true if some object was found within the search radius, false otherwise
//
//    bool m_bfnFarthestNeighbor ( T& tFarthest,   const T& t ) const
//       tFarthest is returned as the object that was found farthest to the probe
//          point
//       t is the probe point, used to search in the group of points m_fnInsert'ed
//       return value is true if some object was found, false otherwise
//
//    long m_lfnFindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//       dRadius is the radius within which to search; make it very large if you want to
//           include every point that was loaded;
//       tClosest is returned as the vector of objects that were found within a radius dRadius
//          of the probe point
//       t is the probe point, used to search in the group of points m_fnInsert'ed
//       return value is the number of objects found within the search radius
//
//    ~CNearTree( void )  // destructor
//       invoked by  vTree.CNeartree<v>::~CNearTree
//       for an object vTree of some type v

// So a complete program is:
//
// #include "TNear.h"
// #include <cstdio>
// void main()
// {
//   CNearTree< double > dT;
//   double dNear;
//   dT.m_fnInsert( 1.5 );
//   if ( dT.m_bfnNearestNeighbor( 10000.0,   dNear,  2.0 )) printf( "%f\n",dRad );
// }
//
// and it should print 0.5 (that's how for 2.0 is from 1.5)
//
//
//-------------------------------------------------------------------------


#if !defined(TNEAR_H_INCLUDED)
#define TNEAR_H_INCLUDED


#include <limits.h>
#include <float.h>
#include <math.h>
#include <vector>

template <typename T> class CNearTree
{
   // m_fnInsert copies the input objects into a binary NEAR tree. When a node has
   // two entries, a descending node is used or created. The current datum is
   // put into the branch descending from the nearer of the two
   // objects in the current node.

   // m_bfnNearestNeighbor retrieves the object nearest to some probe by descending
   // the tree to search out the appropriate object. Speed is gained
   // by pruning the tree if there can be no data below that are
   // nearer than the best so far found.

   // The tree is built in time O(n log n), and retrievals take place in
   // time O(log n).


    T *           m_ptLeft;         // left object (of type T) stored in this node
    T *           m_ptRight;        // right object (of type T) stored in this node
    double        m_dMaxLeft;       // longest distance from the left object to
                                    // anything below it in the tree
    double        m_dMaxRight;      // longest distance from the right object to
                                    // anything below it in the tree
    CNearTree *   m_pLeftBranch;    // tree descending from the left object
    CNearTree *   m_pRightBranch;   // tree descending from the right object

public:

//=======================================================================
//  CNearTree ( )
//
//  Default constructor for class CNearTree
//  creates an empty tree with no right or left node and with the dMax-below
//  set to negative values so that any match found will be stored since it will
//  greater than the negative value
//
//=======================================================================
   CNearTree(void)  // constructor
   {
      m_ptLeft       = 0;
      m_ptRight      = 0;
      m_pLeftBranch  = 0;
      m_pRightBranch = 0;
      m_dMaxLeft     = DBL_MIN;
      m_dMaxRight    = DBL_MIN;
   }  //  CNearTree constructor

//=======================================================================
//  ~CNearTree ( )
//
//  Destructor for class CNearTree
//
//=======================================================================
   ~CNearTree(void)  // destructor
   {
      delete m_pLeftBranch  ;  m_pLeftBranch  =0;
      delete m_pRightBranch ;  m_pRightBranch =0;
      delete m_ptLeft       ;  m_ptLeft       =0;
      delete m_ptRight      ;  m_ptRight      =0;

      m_dMaxLeft     = DBL_MIN;
      m_dMaxRight    = DBL_MIN;
   }  //  ~CNearTree

//=======================================================================
//  void m_fnInsert ( const T& t )
//
//  Function to insert some "point" as an object into a CNearTree for
//  later searching
//
//     t is an object of the templated type which is to be inserted into a
//     Neartree
//
//  Three possibilities exist: put the datum into the left
//  postion (first test),into the right position, or else
//  into a node descending from the nearer of those positions
//  when they are both already used.
//
//=======================================================================
   void m_fnInsert( const T& t )
   {
      // do a bit of precomputing if it is possible so that we can
      // reduce the number of calls to operator 'double' as much as possible;
      // 'double' might use square roots in some cases
      double dTempRight =  0;
      double dTempLeft  =  0;

      if ( m_ptRight  != 0 )
      {
         dTempRight  = fabs( double( t - *m_ptRight ));
         dTempLeft   = fabs( double( t - *m_ptLeft  ));
      }

      if ( m_ptLeft == 0 )
      {
         m_ptLeft = new T( t );
      }
      else if ( m_ptRight == 0 )
      {
         m_ptRight   = new T( t );
      }
      else if ( dTempLeft > dTempRight )
      {
         if ( m_pRightBranch == 0 ) m_pRightBranch = new CNearTree;
         // note that the next line assumes that m_dMaxRight is negative for a new node
         if ( m_dMaxRight < dTempRight ) m_dMaxRight = dTempRight;
         m_pRightBranch->m_fnInsert( t );
      }
      else  // ((double)(t - *m_tLeft) <= (double)(t - *m_tRight) )
      {
         if ( m_pLeftBranch  == 0 ) m_pLeftBranch  = new CNearTree;
         // note that the next line assumes that m_dMaxLeft is negative for a new node
         if ( m_dMaxLeft < dTempLeft ) m_dMaxLeft  = dTempLeft;
         m_pLeftBranch->m_fnInsert( t );
      }

   }  //  m_fnInsert

//=======================================================================
//  bool m_bfnNearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the object closest to some probe point, t. This function
//  is only here so that the function m_bfnNearest can be called without having dRadius const
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    tClosest is an object of the templated type and is the returned nearest point
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found
//
//=======================================================================
   bool m_bfnNearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
   {
      double dSearchRadius = dRadius;
      return ( m_bfnNearest ( dSearchRadius, tClosest, t ) );
   }  //  m_bfnNearestNeighbor

//=======================================================================
//  bool m_bfnFarthestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the object closest to some probe point, t. This function
//  is only here so that the function m_bfnFarthestNeighbor can be called without the user
//  having to input a search radius and so the search radius can be guaranteed to be
//  negative at the start.
//
//    tFarthest is an object of the templated type and is the returned farthest point
//             from the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found (should only be false for
//             an empty tree)
//
//=======================================================================
   bool m_bfnFarthestNeighbor ( T& tFarthest,   const T& t ) const
   {
      double dSearchRadius = DBL_MIN;
      return ( m_bfnFindFarthest ( dSearchRadius, tFarthest, t ) );
   }  //  m_bfnFarthestNeighbor

//=======================================================================
//  long m_lfnFindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the set of objects closer to some probe point, t,
//  than dRadius. This is only here so that tClosest can be cleared before starting the work.
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    tClosest is a vector of objects of the templated type and is the returned set of nearest points
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//    return value is the number of points found within dRadius of the probe point
//
//=======================================================================
   long m_lfnFindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
   {
      // clear the contents of the return vector so that things don't accidentally accumulate
      tClosest.clear( );
      return ( m_lfnInSphere( dRadius, tClosest, t ) );
   }  //  m_lfnFindInSphere

private:

//=======================================================================
//  long m_lfnInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//
//  Private function to search a Neartree for the object closest to some probe point, t.
//  This function is only called by m_lfnFindInSphere.
//
//    dRadius is the search radius
//    tClosest is a vector of objects of the templated type found within dRadius of the
//         probe point
//    t  is the probe point
//    the return value is the number of points found within dRadius of the probe
//
//=======================================================================
   long m_lfnInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
   {
      long lReturn = 0;
      // first test each of the left and right positions to see if
      // one holds a point nearer than the search radius.
      if (( m_ptLeft !=0 ) && (( fabs( double( t - *m_ptLeft  ))) <= dRadius ))
      {
         tClosest.push_back( *m_ptLeft );
         lReturn++ ;
      }
      if (( m_ptRight!=0 ) && (( fabs( double( t - *m_ptRight ))) <= dRadius ))
      {
         tClosest.push_back( *m_ptRight );
         lReturn++ ;
      }
      //
      // Now we test to see if the branches below might hold an object
      // nearer than the search radius. The triangle rule is used
      // to test whether it's even necessary to descend.
      //
      if (( m_pLeftBranch  != 0 )  && (( dRadius + m_dMaxLeft  ) >= fabs( double( t - *m_ptLeft  ))))
      {
         lReturn += m_pLeftBranch->m_lfnInSphere( dRadius, tClosest, t );
      }

      if (( m_pRightBranch != 0 )  && (( dRadius + m_dMaxRight ) >= fabs( double( t - *m_ptRight ))))
      {
         lReturn += m_pRightBranch->m_lfnInSphere( dRadius, tClosest, t );
      }

      return ( lReturn );
   }  //  m_lfnInSphere

//=======================================================================
//  bool m_bfnNearest ( double& dRadius,  T& tClosest,   const T& t ) const
//
//  Private function to search a Neartree for the object closest to some probe point, t.
//  This function is only called by m_bfnNearestNeighbor.
//
//    dRadius is the smallest currently known distance of an object from the probe point.
//    tClosest is an object of the templated type and is the returned closest point
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found within dRadius
//
//=======================================================================
   bool m_bfnNearest ( double& dRadius,  T& tClosest,   const T& t ) const
   {
      double   dTempRadius;
      bool  bRet = false;

      // first test each of the left and right positions to see if
      // one holds a point nearer than the nearest so far discovered.
      if (( m_ptLeft!=0 ) && (( dTempRadius = fabs( double( t - *m_ptLeft ))) <= dRadius ))
      {
         dRadius  = dTempRadius;
         tClosest = *m_ptLeft;
         bRet     = true;
      }
      if (( m_ptRight!=0 ) && (( dTempRadius = fabs( double( t - *m_ptRight))) <= dRadius ))
      {
         dRadius  = dTempRadius;
         tClosest = *m_ptRight;
         bRet     = true;
      }

      //
      // Now we test to see if the branches below might hold an object
      // nearer than the best so far found. The triangle rule is used
      // to test whether it's even necessary to descend.
      //
      if (( m_pLeftBranch  != 0 )  && (( dRadius + m_dMaxLeft  ) >= fabs( double( t - *m_ptLeft  ))))
      {
         bRet |= m_pLeftBranch->m_bfnNearest( dRadius, tClosest, t );
      }

      if (( m_pRightBranch != 0 )  && (( dRadius + m_dMaxRight ) >= fabs( double( t - *m_ptRight ))))
      {
         bRet |= m_pRightBranch->m_bfnNearest( dRadius, tClosest, t );
      }

      return ( bRet );
   };   // m_bfnNearest

//=======================================================================
//  bool m_bfnFindFarthest ( double& dRadius,  T& tFarthest,   const T& t ) const
//
//  Private function to search a Neartree for the object farthest from some probe point, t.
//  This function is only called by m_bfnFarthestNeighbor.
//
//    dRadius is the largest currently known distance of an object from the probe point.
//    tFarthest is an object of the templated type and is the returned farthest point
//             from the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found (should only be false for
//             an empty tree)
//
//=======================================================================
   bool m_bfnFindFarthest ( double& dRadius,  T& tFarthest,   const T& t ) const
   {
      double   dTempRadius;
      bool  bRet     = false;

      // first test each of the left and right positions to see if
      // one holds a point farther than the farthest so far discovered.
      // the calling function is presumed initially to have set dRadius to a
      // negative value before the recursive calls to m_bfnFindFarthestNeighbor

      if (( m_ptLeft!=0  ) && (( dTempRadius = fabs( double( t - *m_ptLeft ))) >= dRadius ))
      {
         dRadius   = dTempRadius;
         tFarthest = *m_ptLeft;
         bRet      = true;
      }
      if (( m_ptRight!=0 ) && (( dTempRadius = fabs( double( t - *m_ptRight))) >= dRadius ))
      {
         dRadius   = dTempRadius;
         tFarthest = *m_ptRight;
         bRet      = true;
      }
      //
      // Now we test to see if the branches below might hold an object
      // farther than the best so far found. The triangle rule is used
      // to test whether it's even necessary to descend.
      //
      if (( m_pLeftBranch  != 0 )  && (( dRadius - m_dMaxLeft  ) <= fabs( double( t - *m_ptLeft  ))))
      {
         bRet |=  m_pLeftBranch->m_bfnFindFarthest( dRadius, tFarthest, t );
      }

      if (( m_pRightBranch != 0 )  && (( dRadius - m_dMaxRight ) <= fabs( double( t - *m_ptRight ))))
      {
         bRet |=  m_pRightBranch->m_bfnFindFarthest( dRadius, tFarthest, t );
      }

      return ( bRet );
   };   // m_bfnFindFarthest


}; // template class TNear

#endif // !defined(TNEAR_H_INCLUDED)
