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
//  distance to the left; it was found that these checks made little
//  to no difference in timing.


// This template is used to contain a collection of objects. After the
// collection has been loaded into this structure, it can be quickly
// queried for which object is "closest" to some probe object of the
// same type. The major restriction on applicability of the near-tree
// is that the algorithm only works if the objects obey the triangle
// inequality. The triangle rule states that the length of any side of
// a triangle cannot exceed the sum of the lengths of the other two sides.


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
//    void Insert( T& t )
//       where t is an object of the type v
//
//    bool NearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//       dRadius is the largest radius within which to search; make it
//          very large if you want to include every point that was loaded; dRadius
//          is returned as the closest distance to the probe (or the search radius
//          if nothing is found)
//       tClosest is returned as the object that was found closest to the probe
//          point (if any were within radius dRadius of the probe)
//       t is the probe point, used to search in the group of points Insert'ed
//       return value is true if some object was found within the search radius, false otherwise
//
//    bool FarthestNeighbor ( T& tFarthest,   const T& t ) const
//       tFarthest is returned as the object that was found farthest to the probe
//          point
//       t is the probe point, used to search in the group of points Insert'ed
//       return value is true if some object was found, false otherwise
//
//    long FindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//       dRadius is the radius within which to search; make it very large if you want to
//           include every point that was loaded;
//       tClosest is returned as the vector of objects that were found within a radius dRadius
//          of the probe point
//       t is the probe point, used to search in the group of points Insert'ed
//       return value is the number of objects found within the search radius
//
//    ~CNearTree( void )  // destructor
//       invoked by  vTree.CNeartree<v>::~CNearTree
//       for an object vTree of some type v
//
//    void DelayedInsert( void ) Holds objects in a queue for insertion later when CompleteDelayInsert
//       is called or a search is called.
//
//    void CompleteDelayedInsert( void ) Calls Insert for all delayed objects. sqrt(n) are inserted
//       by random choice. The rest are inserted in linear order as originally queued.
//
//    size_t GetDeferredSize( void ) Returns the number of delayed objects that have not
//       yet been Insert'ed
//
//    size_t GetTotalSize( void ) Returns the number of objects that have been Insert'ed plus
//       those DelayInsert'ed
//
//    size_t size( void ) identical to GetTotalSize
//
//    size_t GetDepth( void ) Returns the maximum tree layers from the root
//
//    bool empty( void )  returns true if the tree is empty
//
//    iterator NearestNeighbor( const double radius, const T& probe ); returns an iterator 
//       to the nearest point to the probe point or end() if there is none
//
//    iterator FarthestNeighbor( const T& probe ); returns an iterator 
//       to the farthest point to the probe point or end() if there is none
//
// So a complete program is:
//
// #include "TNear.h"
// #include <cstdio>
// void main()
// {
//   CNearTree< double > dT;
//   double dNear;
//   dT.Insert( 1.5 );
//   if ( dT.NearestNeighbor( 10000.0,   dNear,  2.0 )) printf( "%f\n",double(dNear-2.0) );
// }
//
// and it should print 0.5 (that's how for 2.0 is from 1.5)
//
//
//-------------------------------------------------------------------------


#if !defined(TNEAR_H_INCLUDED)
#define TNEAR_H_INCLUDED

#include <stdlib.h>
#ifdef USE_MINGW_RAND
#define random(x) rand(x)
#define srandom(x) srand(x)
#endif


#define MYRAND_MAX  32767

#include <limits.h>
#include <float.h>
#include <math.h>
#include <vector>

#ifdef CNEARTREE_SAFE_TRIANG
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0) \
|| ((b)-((a)-(c)) >= 0) \
|| ((c)-((a)-(b)) >= 0))    
#else
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0))
#endif



template <typename T> class CNearTree
{

    // Insert copies the input objects into a binary NEAR tree. When a node has
    // two entries, a descending node is used or created. The current datum is
    // put into the branch descending from the nearer of the two
    // objects in the current node.
    
    // NearestNeighbor retrieves the object nearest to some probe by descending
    // the tree to search out the appropriate object. Speed is gained
    // by pruning the tree if there can be no data below that are
    // nearer than the best so far found.
    
    // The tree is built in time O(n log n), and retrievals take place in
    // time O(log n).
    
    
    size_t            m_ptLeft;            // index of left object (of type T) stored in this node
    size_t            m_ptRight;           // index of right object (of type T) stored in this node
    double            m_dMaxLeft;          // longest distance from the left object to
                                           // anything below it in the tree
    double            m_dMaxRight;         // longest distance from the right object to
                                           // anything below it in the tree
    CNearTree *       m_pLeftBranch;       // tree descending from the left object
    CNearTree *       m_pRightBranch;      // tree descending from the right object
    std::vector<long> m_DelayedIndices;    // objects queued for insertion, possibly in random order
    std::vector<T>    m_DelayedInsertions; // objects queued for insertion, possibly in random order
    std::vector<T>    m_ObjectStore;       // all inserted objects go here
    size_t            m_DeepestDepth;      // maximum diameter of the tree
    size_t            m_InsertedSize;      // objects actually inserted (queued ones are not counted)
    
public:

   // Forward declaration for the nested class. Friend is necessary
   // for the access to the appropriate data elements
   class iterator;
   friend class iterator;

//=======================================================================
//  CNearTree ( )
//
//  Default constructor for class CNearTree
//  creates an empty tree with no right or left node and with the dMax-below
//  set to negative values so that any match found will be stored since it will
//  greater than the negative value
//
//=======================================================================

   CNearTree(void) :  // constructor
      m_ptLeft            ( ULONG_MAX ),
      m_ptRight           ( ULONG_MAX ),
      m_dMaxLeft          ( DBL_MIN ),
      m_dMaxRight         ( DBL_MIN ),
      m_pLeftBranch       ( 0 ),
      m_pRightBranch      ( 0 ),
      m_DelayedIndices    ( 0 ),
      m_DelayedInsertions (   ),
      m_ObjectStore       (   ),
      m_DeepestDepth      ( 0 ),
      m_InsertedSize      ( 0 )

   {
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
        
        m_dMaxLeft     = DBL_MIN;
        m_dMaxRight    = DBL_MIN;
        m_DeepestDepth      = 0;
        m_InsertedSize      = 0;
    }  //  ~CNearTree
    
//=======================================================================
//  empty ( )
//
//  Test for an empty CNearTree
//
//=======================================================================
    
    bool empty( ) const
    {
        return( m_ptLeft == ULONG_MAX && m_DelayedInsertions.empty( ) );
    }
    
//=======================================================================
//  void Insert ( const T& t )
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
    void Insert( const T& t )
    {
        size_t localDepth = 0;
       Inserter( t, localDepth, m_ObjectStore );
        m_DeepestDepth = std::max( localDepth, m_DeepestDepth );
        ++m_InsertedSize;
    }
    
//=======================================================================
//  iterator NearestNeighbor( const double radius, const T& t ) const
//
//  Function to search a Neartree for the object closest to some probe point, t. This function
//  is only here so that the function Nearest can be called without having the radius const.
//  This was necessary because Nearest is recursive, but needs to keep the current radius.
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    t  is the probe point
//
//    the return is an iterator to the templated type and is the returned nearest point
//             to the probe point (t) that can be found in the Neartree
//             or iterator::end if no point was found
//
//=======================================================================
    iterator NearestNeighbor( const double radius, const T& t ) const
    {
       T closest;
       size_t index = ULONG_MAX;
       double tempRadius = radius;
       if( Nearest( tempRadius, closest, t, index ) )
       {
          return( iterator( (long)index, this ) );
       }
       else
       {
          return( end( ) );
       }
    }

//=======================================================================
//  bool NearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the object closest to some probe point, t. This function
//  is only here so that the function Nearest can be called without having the radius const.
//  This was necessary because Nearest is recursive, but needs to keep the current radius.
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    tClosest is an object of the templated type and is the returned nearest point
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//
//    the return value is true only if a point was found
//
//=======================================================================
    bool NearestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
    {
        if( dRadius < 0.0 ) 
        {
            return( false );
        }
        else if( this->empty( ) )
        {
            return( false );
        }
        else
        {
            const_cast<CNearTree*>(this)->CompleteDelayedInsert( );
            double dSearchRadius = dRadius;
         size_t index = ULONG_MAX;
         return ( const_cast<CNearTree*>(this)->Nearest ( dSearchRadius, tClosest, t, index ) );
        }
    }  //  NearestNeighbor
    

//=======================================================================
//  iterator FarthestNeighbor( const double radius, const T& t ) const
//
//  Function to search a Neartree for the object farthest from some probe point, t. This function
//  is only here so that the function Farthest can be called without having the radius const.
//  This was necessary because Farthest is recursive, but needs to keep the current radius.
//
//    t  is the probe point
//
//    the return is an iterator to the templated type and is the returned farthest point
//             to the probe point (t) that can be found in the Neartree
//             or iterator::end if no point was found
//
//=======================================================================
    iterator FarthestNeighbor( const T& t ) const
    {
       T farthest;
       size_t index = ULONG_MAX;
       double radius = DBL_MIN;
       if( Farthest( radius, farthest, t, index ) )
       {
          return( iterator( (long)index, this ) );
       }
       else
       {
          return( end( ) );
       }
    }
//=======================================================================
//  bool FarthestNeighbor ( const double& dRadius,  T& tClosest,   const T& t ) const
//
//  Function to search a Neartree for the object closest to some probe point, t. This function
//  is only here so that the function FarthestNeighbor can be called without the user
//  having to input a search radius and so the search radius can be guaranteed to be
//  negative at the start.
//
//    tFarthest is an object of the templated type and is the returned farthest point
//             from the probe point that can be found in the Neartree
//    t  is the probe point
//
//    the return value is true only if a point was found (should only be false for
//             an empty tree)
//
//=======================================================================
    bool FarthestNeighbor ( T& tFarthest,   const T& t ) const
    {
        if( this->empty( ) )
        {
            return( false );
        }
        else
        {
            const_cast<CNearTree*>(this)->CompleteDelayedInsert( );
            double dSearchRadius = DBL_MIN;
         size_t index = ULONG_MAX;
         return ( const_cast<CNearTree*>(this)->Farthest ( dSearchRadius, tFarthest, t, index ) );
        }
    }  //  FarthestNeighbor
    
//=======================================================================
//  long FindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
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
    long FindInSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
    {
        // clear the contents of the return vector so that things don't accidentally accumulate
        tClosest.clear( );
        const_cast<CNearTree*>(this)->CompleteDelayedInsert( );
        return ( const_cast<CNearTree*>(this)->InSphere( dRadius, tClosest, t ) );
    }  //  FindInSphere
    
    
//=======================================================================
//  void DelayedInsert ( const T& t )
//
//  Function to insert some "point" as an object into a CNearTree for
//  later searching
//
//     t is an object of the templated type which is to be inserted into a
//     Neartree
//
//  The function Insert immediately inserts the object into the tree. 
//  DelayedInsert keeps the object in an internal store, but does not 
//  immediately insert it. The object in the internal store are only inserted
//  when CompleteDelayedInsert is called or when one of the search functions
//  is invoked (they call CompleteDelayedInsert). When that is called, all
//  of the stored objects are then inserted into the list in a way designed
//  to give a relatively balanced tree even if the data are strongly sorted.
//
//=======================================================================
    void DelayedInsert( const T& t )
    {
        m_DelayedInsertions.push_back( t );
        m_DelayedIndices   .push_back( (long)m_DelayedInsertions.size( ) - 1 );
    };
    
//=======================================================================
//  void CompleteDelayedInsert( void )
//
//  When CompleteDelayedInsert is invoked, if there are any objects in the 
//  internal store they are then inserted into the neartree. CompleteDelayedInsert
//  randomly selects enough objects to half fill a well-balanced tree. That is,
//  if there are n objects in the internal store, it randomly selects and inserts
//  sqrt(n) objects. After that, the remaining objects are inserted in a linear
//  sequence as they were entered.
//
//=======================================================================
    void CompleteDelayedInsert( void )
    {
        if( m_DelayedInsertions.empty( ) )
        {
            return;
        }
        
        // Insert a random selection of the objects
        const size_t vectorSize = m_DelayedInsertions.size( );
      const size_t toRandomlyInsert = (size_t)::sqrt( (double)vectorSize );
        for( size_t i=0; i<toRandomlyInsert; ++i )
        {
            size_t n = (size_t)((double)(vectorSize-1u) * (double)(random( )%MYRAND_MAX) / (double) MYRAND_MAX);
            random( ); random( );
            random( ); random( );
            // Find the next pointer that hasn't already had its object "Insert"ed
            // We can do this blindly since sqrt(n)<=n for all cases. n=1 would be the only 
            // bad case here, and that will not trigger the later loop.
            while( m_DelayedIndices[n] == -1 )
            {
                ++n;
                n = n% vectorSize;
            }
            Insert( m_DelayedInsertions[n] );
            m_DelayedIndices[n] = -1;         
        }
        
        // finish by inserting all the remaining objects
        for( size_t i=0; i<vectorSize; ++i )
        {
            if( m_DelayedIndices[i] != -1 )
            {
                Insert( m_DelayedInsertions[i] );
            }
        }
        
        // now get rid of the temporary storage that was used for delayed 
        // insertions (fast way, faster than clear() )
        std::vector<long> DelayedPointersTemp;   
        std::vector<T>  DelayedInsertionsTemp;
        DelayedPointersTemp  .swap( m_DelayedIndices );
        DelayedInsertionsTemp.swap( m_DelayedInsertions );
    };
    
//=======================================================================
//  size_t GetDeferredSize(  void )
//
//  The number of objects currently queued for insertion.
//
//=======================================================================
   size_t GetDeferredSize( void ) const
    {
        return( m_DelayedInsertions.size( ) );
    };
    
//=======================================================================
//  size_t GetTotalSize(  void )
//
//  The total number of objects that have been inserted plus those
//  queued for insertion.
//
//=======================================================================
   size_t GetTotalSize( void ) const
    {
        return( m_InsertedSize + m_DelayedInsertions.size( ) );
    };

//=======================================================================
//  size_t size(  void )
//
//  The total number of objects that have been inserted plus those
//  queued for insertion.
//
//=======================================================================
   size_t size( void ) const
   {
      return( GetTotalSize( ) );
   };

//=======================================================================
//  size_t GetDepth(  void ) const
//
//  The greatest depth of the tree (1-based) from the root.
//
//=======================================================================
   size_t GetDepth( void ) const
    {
        return( m_DeepestDepth );
    };
    
private:
//=======================================================================
//  void Inserter ( const T& t( const T& t, size_t& localDepth )
//
//  Function to insert some "point" as an object into a CNearTree for
//  later searching
//
//     t is an object of the templated type which is to be inserted into a
//     Neartree
//
//     localDepth is the returned deepest tree level reached for the current insert
//
//  Three possibilities exist: put the datum into the left
//  postion (first test),into the right position, or else
//  into a node descending from the nearer of those positions
//  when they are both already used.
//
//=======================================================================
   void Inserter( const T& t, size_t& localDepth, std::vector<T>& ObjectStore )
    {
        // do a bit of precomputing if it is possible so that we can
        // reduce the number of calls to operator 'double' as much as possible;
        // 'double' might use square roots in some cases
        double dTempRight =  0;
        double dTempLeft  =  0;
        ++localDepth;
        
      if ( m_ptRight  != ULONG_MAX )
      {
         dTempRight  = ::fabs( double( t - ObjectStore[m_ptRight] ));
         dTempLeft   = ::fabs( double( t - ObjectStore[m_ptLeft]  ));
        }
        
      if ( m_ptLeft == ULONG_MAX )
      {
         m_ptLeft = ObjectStore.size( );
         ObjectStore.push_back( t );
      }
      else if ( m_ptRight == ULONG_MAX )
      {
         m_ptRight = ObjectStore.size( );
         ObjectStore.push_back( t );
        }
        else if ( dTempLeft > dTempRight )
        {
            if ( m_pRightBranch == 0 ) m_pRightBranch = new CNearTree;
            // note that the next line assumes that m_dMaxRight is negative for a new node
            if ( m_dMaxRight < dTempRight ) m_dMaxRight = dTempRight;
         m_pRightBranch->Inserter( t, localDepth, ObjectStore );
        }
        else  // ((double)(t - *m_tLeft) <= (double)(t - *m_tRight) )
        {
            if ( m_pLeftBranch  == 0 ) m_pLeftBranch  = new CNearTree;
            // note that the next line assumes that m_dMaxLeft is negative for a new node
            if ( m_dMaxLeft < dTempLeft ) m_dMaxLeft  = dTempLeft;
         m_pLeftBranch->Inserter( t, localDepth, ObjectStore );
        }
    }  //  Inserter

   //=======================================================================
//  long InSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
//
//  Private function to search a Neartree for the object closest to some probe point, t.
//  This function is only called by FindInSphere.
//
//    dRadius is the search radius
//    tClosest is a vector of objects of the templated type found within dRadius of the
//         probe point
//    t  is the probe point
//    the return value is the number of points found within dRadius of the probe
//
//=======================================================================
    long InSphere ( const double& dRadius,  std::vector<  T >& tClosest,   const T& t ) const
    {
        std::vector <CNearTree<T>* > sStack;
        long lReturn = 0;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
        CNearTree* pt = const_cast<CNearTree*>(this);
      if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
            const double dDR = ::fabs(double( t - m_ObjectStore[pt->m_ptRight] ));
                if ( dDR <= dRadius )
                {
                    ++lReturn;
               tClosest.push_back( m_ObjectStore[pt->m_ptRight]);
            }
            if ( pt->m_pRightBranch != 0 && pt->m_dMaxRight+dRadius >= dDR )
                { // we did the left and now we finished the right, go down
                    pt = pt->m_pRightBranch;
                    eDir = left;
                }
                else
                {
                    eDir = end;
                }
            }
            if ( eDir == left )
            {
            const double dDL = ::fabs(double( t - m_ObjectStore[pt->m_ptLeft]) );
                if ( dDL <= dRadius )
                {
                    ++lReturn;
               tClosest.push_back( m_ObjectStore[pt->m_ptLeft]);
            }
            if ( pt->m_ptRight != ULONG_MAX ) // only stack if there's a right object
                {
                    sStack.push_back( pt );
                }
            if ( pt->m_pLeftBranch != 0 &&  pt->m_dMaxLeft+dRadius >= dDL )
                { // we did the left, go down
                    pt = pt->m_pLeftBranch;
                }
                else
                {
                    eDir = end;
                }
            }
            
            if ( eDir == end && !sStack.empty( ) )
            {
                pt = sStack.back( );
                sStack.pop_back( );
                eDir = right;
            }
        }

        return ( lReturn );
    }  //  InSphere

//=======================================================================
//  bool Nearest ( double& dRadius,  T& tClosest,   const T& t ) const
//
//  Private function to search a Neartree for the object closest to some probe point, t.
//  This function is only called by NearestNeighbor.
//
//    dRadius is the smallest currently known distance of an object from the probe point.
//    tClosest is an object of the templated type and is the returned closest point
//             to the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found within dRadius
//
//=======================================================================
   bool Nearest ( double& dRadius, T& tClosest, const T& t, size_t& pClosest ) const
    {
        std::vector <CNearTree<T>* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
        CNearTree* pt = const_cast<CNearTree*>(this);
      pClosest = ULONG_MAX;
      if ( pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
            const double dDR = ::fabs(double( t - m_ObjectStore[pt->m_ptRight] ));
                if ( dDR < dRadius )
                {
                    dRadius = dDR;
                    pClosest = pt->m_ptRight;
                }
                if ( pt->m_pRightBranch != 0 && (TRIANG(dDR,pt->m_dMaxRight,dRadius)))
                { // we did the left and now we finished the right, go down
                    pt = pt->m_pRightBranch;
                    eDir = left;
                }
                else
                {
                    eDir = end;
                }
            }
            if ( eDir == left )
            {
            const double dDL = ::fabs(double( t - m_ObjectStore[pt->m_ptLeft] ));
                if ( dDL < dRadius )
                {
                    dRadius = dDL;
                    pClosest = pt->m_ptLeft;
                }
            if ( pt->m_ptRight != ULONG_MAX ) // only stack if there's a right object
                {
                    sStack.push_back( pt );
                }
                if ( pt->m_pLeftBranch != 0 && (TRIANG(dDL,pt->m_dMaxLeft,dRadius)))
                { // we did the left, go down
                    pt = pt->m_pLeftBranch;
                }
                else
                {
                    eDir = end;
                }
            }
            
            if ( eDir == end && !sStack.empty( ) )
            {
                pt = sStack.back( );
                sStack.pop_back( );
                eDir = right;
            }
        }
        while ( !sStack.empty( ) ) // for safety !!!
            sStack.pop_back( );
      if ( pClosest != ULONG_MAX )
         tClosest = m_ObjectStore[pClosest];
      return ( pClosest != ULONG_MAX );
    };   // Nearest

//=======================================================================
//  bool Farthest ( double& dRadius,  T& tFarthest,   const T& t ) const
//
//  Private function to search a Neartree for the object farthest from some probe point, t.
//  This function is only called by FarthestNeighbor.
//
//    dRadius is the largest currently known distance of an object from the probe point.
//    tFarthest is an object of the templated type and is the returned farthest point
//             from the probe point that can be found in the Neartree
//    t  is the probe point
//    the return value is true only if a point was found (should only be false for
//             an empty tree)
//
//=======================================================================
   bool Farthest ( double& dRadius,  T& tFarthest,   const T& t, size_t& pFarthest ) const
    {
        std::vector <CNearTree<T>* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
      CNearTree* pt = const_cast<CNearTree*>(this);
      pFarthest = ULONG_MAX;
      if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
            const double dDR = ::fabs(double( t - m_ObjectStore[pt->m_ptRight]) );
                if ( dDR >= dRadius )
                {
                    dRadius = dDR;
                    pFarthest = pt->m_ptRight;
                }
                if ( pt->m_pRightBranch != 0 && TRIANG(dRadius,dDR,pt->m_dMaxRight))
                { // we did the left and now we finished the right, go down
                    pt = pt->m_pRightBranch;
                    eDir = left;
                }
                else
                {
                    eDir = end;
                }
            }
            if ( eDir == left )
            {
            const double dDL = ::fabs(double( t - m_ObjectStore[pt->m_ptLeft] ) );
                if ( dDL >= dRadius )
                {
                    dRadius = dDL;
                    pFarthest = pt->m_ptLeft;
                }
            if ( pt->m_ptRight != ULONG_MAX ) // only stack if there's a right object
                {
                    sStack.push_back( pt );
                }
                if ( pt->m_pLeftBranch != 0 && TRIANG(dRadius,dDL,pt->m_dMaxLeft) )
                { // we did the left, go down
                    pt = pt->m_pLeftBranch;
                }
                else
                {
                    eDir = end;
                }
            }
            
            if ( eDir == end && !sStack.empty( ) )
            {
                pt = sStack.back( );
                sStack.pop_back( );
                eDir = right;
            }
        }
        while ( !sStack.empty( ) ) // for safety !!!
            sStack.pop_back( );
      if ( pFarthest != ULONG_MAX )
         tFarthest = m_ObjectStore[pFarthest];
      return ( pFarthest != ULONG_MAX );
    };   // Farthest

   public:
         iterator begin( void ) const { return( iterator( 0, this ) ); };
         iterator end  ( void ) const { return( iterator(m_ObjectStore.empty()?1:(long)m_ObjectStore.size( ), this ) ); };
         iterator back ( void ) const { return( iterator( (m_ObjectStore.empty())? 1 :(long)m_ObjectStore.size( )-1, this ) ); };
         T operator[] ( const size_t position ) const { return( m_ObjectStore[position] ); };



      //====================================================================================
      class iterator
      {
         friend class CNearTree<T>;
         long position;
         const CNearTree<T>* parent;
      public:
         iterator( void ) { }; // constructor
         explicit iterator( const long s ) : position(s), parent((CNearTree<T>*)this) { }; // constructor
         explicit iterator( const int  s ) : position(s), parent((CNearTree<T>*)this) { }; // constructor

         iterator& operator=  ( const iterator& s )      { position = s.position; parent = s.parent; return( *this ); };
         iterator  operator++ ( const int n )            { iterator it(*this); position+=1+n; return( it ); };
         iterator  operator-- ( const int n )            { iterator it(*this); position-=1+n; return( it ); };
         iterator& operator++ ( void )                   { ++position; return( *this ); };
         iterator& operator-- ( void )                   { --position; return( *this ); };
         iterator  operator+  ( const long n ) const     { iterator it( position+n, parent); return( it ); };
         iterator  operator-  ( const long n ) const     { iterator it( position-n, parent); return( it ); };
         iterator& operator+= ( const long n )           { position += n; return( *this ); };
         iterator& operator-= ( const long n )           { position -= n; return( *this ); };
         T         operator*  ( void )         const     { return( parent->m_ObjectStore[position] ); };
         T*        operator-> ( void )         const     { return( &(const_cast<CNearTree<T>*>(parent)->m_ObjectStore[position]) ); };

         bool      operator== ( const iterator t ) const { return( t.position==(parent->m_ObjectStore.empty()?1:position) && t.parent==parent ); };
         bool      operator!= ( const iterator t ) const { return( ! (*this==t )); };

      private:
         iterator( const long s, const CNearTree<T>* nt ) { position = s; parent = nt; }; // constructor
      }; // class iterator



}; // template class TNear

#endif // !defined(TNEAR_H_INCLUDED)
