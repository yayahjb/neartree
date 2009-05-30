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
//  Later revisions have replaced the use of recursion with a stack,
//  except for the case of inserting data into the tree.


// This template is used to contain a collection of objects. After the
// collection has been loaded into this structure, it can be quickly
// queried for which object is "closest" to some probe object of the
// same type. The major restriction on applicability of the near-tree
// is that the algorithm only works if the objects obey the triangle
// inequality. The triangle rule states that the length of any side of
// a triangle cannot exceed the sum of the lengths of the other two sides.


// The user of this class needs to provide at least the following
// functionality for the template to work. For the built-in
// numerics of C++, they are provided here (or else you should create them).

//    DistanceType Norm( );   // conversion constructor from the templated class to DistanceType
//                                (usually will return a "length" of type double)
//    operator- ( );          // geometrical (vector) difference of two objects
//    a copy constructor would be nice
//    a constructor would be nice
//    a destructor would be nice

// The provided interface is:
//
//    #include "TNear.h"
//
//    CNearTree( void )   // constructor
//       instantiated by something like:      CNearTree <T> vTree;
//       for some type T
//       the following additional convenience constructors are available
//
//    CNearTree( const ContainerType<T> )   // constructor from containers, std::vector, ..., or CNearTree
//
//    void insert( const T& t )
//       where t is an object of the type T
//       the following additional convenience insert template available
//       all inserts are delayed until a search is performed or until an explicit call to CompleteDelayedInsertions
//       is called or a search is called. The purpose is to distribute the objects a bit more 
//       randomly. Excessively ordered objects leads to less than optimal trees.
//       Places objects in a queue for insertion later when CompleteDelayInsert
//
//    void insert( ContainerType ) // for containers, std::vector, ..., or CNearTree
//       all inserts are delayed until a search is performed or until an explicit call to CompleteDelayedInsertions
//
//    bool NearestNeighbor ( const DistanceType dRadius,  T& tClosest,   const T& t ) const
//       dRadius is the largest radius within which to search; make it
//          very large if you want to include every point that was loaded; dRadius
//          is returned as the closest distance to the probe (or the search radius
//          if nothing is found)
//       tClosest is returned as the object that was found closest to the probe
//          point (if any were within radius dRadius of the probe)
//       t is the probe point, used to search in the group of points insert'ed

//       return value is true if some object was found within the search radius, false otherwise
//
//    iterator NearestNeighbor( const DistanceType radius, const T& probe ); returns an iterator 
//       to the nearest point to the probe point or end() if there is none
//
//    bool FarthestNeighbor ( T& tFarthest,   const T& t ) const
//       tFarthest is returned as the object that was found farthest to the probe
//          point
//       t is the probe point, used to search in the group of points insert'ed
//       return value is true if some object was found, false otherwise
//
//    iterator FarthestNeighbor( const T& probe ); returns an iterator 
//       to the farthest point to the probe point or end() if there is none
//
//
//    the following functions (FindInSphere, FindOutSphere, and FindInAnnulus) all return a container 
//    (ContainerType) that can be any standard library container (such as std::vector< T >) or CNearTree.
//
//    long FindInSphere ( const DistanceType dRadius,  ContainerType& tClosest, const T& t ) const
//       dRadius is the radius within which to search; make it very large if you want to
//           include every point that was loaded;
//       tClosest is returned as the ContainerType of objects that were found within a radius dRadius
//          of the probe point
//       t is the probe point, used to search in the group of points insert'ed
//       return value is the number of objects found within the search radius
//
//    long FindOutSphere ( const DistanceType dRadius,  ContainerType& tClosest, const T& t ) const
//       dRadius is the radius outside which to search; make it very small if you want to
//           include every point that was loaded;
//       tClosest is returned as the ContainerType of objects that were found within a radius dRadius
//          of the probe point
//       t is the probe point, used to search in the group of points insert'ed
//       return value is the number of objects found within the search radius
//
//    long FindInAnnulus (const DistanceType dRadius1, const DistanceType dRadius2, ContainerType& tClosest,   const T& t ) const
//       dRadius1 and dRadius2 are the two radii between which to find  data points
//       tClosest is returned ContainerType of the objects found in the annulus
//       t is the probe point, used to search in the group of points insert'ed
//       return value is the number of objects found within the search radius
//
//    long FindK_NearestNeighbors ( const size_t k, const DistanceType& radius,  OutputContainerType& tClosest,   const T& t )
//       k is the maximum number of nearest neighbors to return. Finds this many if possible
//       radius Within a sphere defined by radius, search for the k-nearest-neighbors
//       tClosest is returned ContainerType of the objects found within the sphere
//       t is the probe point, used to search in the group of points insert'ed
//
//    long FindK_FarthestNeighbors ( const size_t k, OutputContainerType& tClosest,   const T& t )
//       k is the maximum number of farthest neighbors to return. Finds this many if possible
//       tClosest is returned ContainerType of the objects found
//       t is the probe point, used to search in the group of points insert'ed
//
//    ~CNearTree( void )  // destructor
//
// =====================================================================================================
// access functions
//
// T at( const size_t n ) const
//     returns the n'th item of the internal data store
//
// T operator[] ( const size_t n)
//     returns the n'th item of the internal data store
//
// operator ContainerType( void ) const 
//     returns all of the inserted objects in the tree in a container of type ContainerType.
//     ContainerType can be std::vector<T>, etc, or other containers.
//     The returned vector contents are not guaranteed to be returned in the order loaded.
//
// iterator begin ( void ) const
//     returns an iterator to the beginning of the internal data store
//
// iterator end ( void ) const
//     returns an iterator to the end of the data store (one beyond the last item)
//
// iterator back ( void ) const
//     returns an iterator to the last data item of the internal data store
//
// =====================================================================================================
// information and special operation functions
// =====================================================================================================
//
//    void ImmediateInsert( void ) Places objects immediately into the tree. The usual insert function
//       delays insertions, allowing them to be inserted into the tree in a more random order. The delay
//       can improve the structure of the tree and speed searches.
//
//    void CompleteDelayedInsert( void ) Calls insert for all delayed objects. sqrt(n) of them are inserted
//       by random choice. The rest are inserted in linear order as originally queued. CompleteDelayedInsert
//       is invoked at the beginning of all searches, so the average user will never need
//       to call it.
//
//    size_t GetDeferredSize( void ) Returns the number of delayed objects that have not
//       yet been insert'ed. This is mainly for information about details of the tree.
//
//    size_t GetTotalSize( void ) Returns the number of objects that have been insert'ed plus
//       those DelayInsert'ed
//
//    size_t size( void ) identical to GetTotalSize
//
//    size_t GetDepth( void ) Returns the maximum tree layers from the root.  This is 
//       mainly for information about details of the tree.
//
//    bool empty( void )  returns true if the tree is empty, otherwise false
//
// =====================================================================================================
// iterators
//     Random access iterators are provided for accessing the data in a CNearTree. The most important
//     expected use is to retrieve the objects returned from one of the sphere search functions that
//     return a CNearTree. However, they can be used with any CNearTree.
//     They should function in a fashion essentially the same as STL iterators. There is no assurance
//     that data will be returned in the order it was loaded, just that it is accessible. The same set is 
//     provided for const_iterator.
// =====================================================================================================
//      iterator( void ) { }; // constructor
//
//      iterator& operator=  ( const iterator& s )      
//      iterator  operator++ ( const int n )            
//      iterator  operator-- ( const int n )            
//      iterator& operator++ ( void )                   
//      iterator& operator-- ( void )                   
//      iterator  operator+  ( const long n ) const     
//      iterator  operator-  ( const long n ) const     
//      iterator& operator+= ( const long n )           
//      iterator& operator-= ( const long n )           
//      T         operator*  ( void )         const     
//
//      bool      operator== ( const iterator& t ) const 
//      bool      operator!= ( const iterator& t ) const 
//
//      const T * const operator-> ( void )   const
//
// =====================================================================================================
//
// So a complete program is:
//
// #include "TNear.h"
// #include <cstdio>
// void main()
// {
//   CNearTree< double > dT;
//   double dNear;
//   dT.insert( 1.5 );
//   if ( dT.FindNearestNeighbor( 10000.0,   dNear,  2.0 )) printf( "%f\n",DistanceType(dNear-2.0) );
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

#include <cmath>

#define MYRAND_MAX  32767

#include <list>
#include <set>
#include <vector>

#ifdef CNEARTREE_SAFE_TRIANG
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0) \
|| ((b)-((a)-(c)) >= 0) \
|| ((c)-((a)-(b)) >= 0))    
#else
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0))
#endif


//=======================================================================
// CNearTree is the root class for the neartree. The actual data of the
// tree is stored in NearTreeNode objects descending from a CNearTree.
//=======================================================================

template <typename T, typename DistanceType=double, int distMinValue=-1 > class CNearTree
{
    //=======================================================================
    //   NOTES:
    //
    // The types of objects that can be stored in the tree is quite broad. The 
    // biggest limitation is that the objects must reside in some sort of metric
    // space and must obey the triangle rule. They must also be all of the same
    // size because they are stored in an std::vector. If your application
    // requires object of varying storage, then your best way to use this
    // code is to store pointers or handles and to write your own distance functions.
    //
    // The type of the objects to be stored is the only _required_ template argument.
    // The type of the distance measure (DistanceType) defaults to double. If your 
    // applications is for an integer type then the type for DistanceType can be your 
    // integer type. This has the potential for speeding the calculations by
    // avoiding FP computation. Other general types can be used if desired, but you
    // may need to also input a value of distMinValue.
    //
    // The template argument distMinValue must be something that your class will
    // understand as a negative number. The default input is negative one. Internally,
    // that is cast to DistanceType. Since most uses will be for DistanceType
    // to be double, that is a simple conversion. Obviously, for integer types,
    // there is no problem either. The need for this value is to have something
    // internally that is recognizable as smaller than the smallest "distance"
    // that can exist between any two objects in your type. For most users,
    // there is no need to input anything other than the default, -1. -1 must 
    // be castable to DistanceType. It seems unlikely that anyone would actually
    // need this optional parameter, but it is here for completeness.
    //
    //  It is a design decision that this class cannot work for unsigned types.
    //  It is hard to see how to verify the triangle rule for unsigned types,
    //  and distance computations become more complex. Sorry, unsigned types
    //  are left as an exercise for the reader.
    //
    //=======================================================================

    // insert copies the input objects into a binary NEAR tree. When a node has
    // two entries, a descending node is used or created. The current datum is
    // put into the branch descending from the nearer of the two
    // objects in the current node.
    
    // NearestNeighbor retrieves the object nearest to some probe by descending
    // the tree to search out the appropriate object. Speed is gained
    // by pruning the tree if there can be no data below that are
    // nearer than the best so far found.
    
    // The tree is built in time O(n log n), and retrievals take place in
    // average time O(log n). However, worst case is O(n).

public:
    // DistanceBetween
    // template function for calculating the "distance" between two objects.
    // The specific functions for the built-in types must be here also. For
    // the common types (int, float, ...) they are provided.
    template <typename TT>
    static inline DistanceType DistanceBetween( const TT& t1, const TT& t2 )
    {
        DistanceType d = ( t1-t2 ).Norm( );
        return( d>0?d:-d );  // apparent compiler error makes this necessary
    }

    static inline DistanceType DistanceBetween( const double t1, const double t2 )
    {
        return( (DistanceType)fabs( t1-t2 ) ); // encourage the compiler to get the correct abs
    }

    //static inline DistanceType DistanceBetween( const long double t1, const long double t2 )
    //{
    //    return(  (DistanceType)fabsl(t1-t2) ); // encourage the compiler to get the correct abs
    //}

    static inline DistanceType DistanceBetween( const float t1, const float t2 )
    {
        return( (DistanceType)fabsf( t1-t2 )); // encourage the compiler to get the correct abs
    }

    static inline DistanceType DistanceBetween( const int t1, const int t2 )
    {
        return( (DistanceType)abs(t1-t2) ); // encourage the compiler to get the correct abs
    }

    static inline DistanceType DistanceBetween( const long t1, const long t2 )
    {
        return( (DistanceType)labs(t1-t2) ); // encourage the compiler to get the correct abs
    }

    //static inline DistanceType DistanceBetween( const long long t1, const long long t2 )
    //{
    //    return( (DistanceType)llabs(t1-t2) ); // encourage the compiler to get the correct abs
    //}

    static inline DistanceType DistanceBetween( const short t1, const short t2 )
    {
        return( (DistanceType)abs(t1-t2) ); // encourage the compiler to get the correct abs
    }
    
    
private:
    // forward declaration of nested class NearTreeNode
    template <typename TNode, typename DistanceTypeNode, int distMinValueNode > 
    class NearTreeNode;
public:
    // Forward declaration for the nested classes, iterator and const_itertor. Friend is necessary
    // for the access to the appropriate data elements
    class iterator;
    friend class iterator;
    class const_iterator;
    friend class const_iterator;
    
private: // start of real definition of CNearTree
    std::vector<long> m_DelayedIndices;    // objects queued for insertion, possibly in random order
    std::vector<T>    m_ObjectStore;       // all inserted objects go here
    size_t            m_DeepestDepth;      // maximum diameter of the tree

    NearTreeNode<T, DistanceType, distMinValue>      m_BaseNode; // the tree's data is stored down from here
    
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
   CNearTree ( void ) :  // constructor
       m_DelayedIndices (   ),
       m_ObjectStore    (   ),
       m_DeepestDepth   ( 0 ),
       m_BaseNode       (   )

   {
   }  //  CNearTree constructor

//=======================================================================
      // CNearTree ( const InputContainer& o )
      //
      // templated constructor for class CNearTree for input of containers.
      // The containers can be standard library containers or a CNearTree.
      //
//=======================================================================
      template<typename InputContainer> 
      CNearTree ( const InputContainer& o ) :  // constructor
           m_DelayedIndices (   ),
           m_ObjectStore    (   ),
           m_DeepestDepth   ( 0 ),
           m_BaseNode       (   )

   {
          typename InputContainer::const_iterator it;
       for( it=o.begin(); it!=o.end(); ++it )
       {
              insert( *it);
       }
    }  //  CNearTree constructor
    
//=======================================================================
//  ~CNearTree ( )
//
//  Destructor for class CNearTree
//
//=======================================================================
    ~CNearTree ( void )  // destructor
    {
        clear ( );
    }  //  ~CNearTree
    
      //=======================================================================
      // clear( void )
      //
      // removes all content from a tree
      //
//=======================================================================
    void clear ( void )
    {
        m_DeepestDepth = 0;

        if ( ! m_DelayedIndices.empty( ) )
        {
           std::vector<long> vtempLong;
           m_DelayedIndices.swap( vtempLong );  // release any delayed indices list
        }

        std::vector<T> vtempT;
        m_ObjectStore.swap( vtempT );  // release the object store

        this->m_BaseNode .clear( ); // clear the nodes of the tree
    }
    
//=======================================================================
//  empty ( )
//
//  Test for an empty CNearTree
//
//=======================================================================
    bool empty ( void ) const
    {
        return ( m_ObjectStore.empty( ) );
    }
    
//=======================================================================
//  void insert ( const T& t )
//
//  Function to insert some "point" as an object into a CNearTree for
//  later searching
//
//     t is an object of the templated type which is to be inserted into a
//     NearTree
//
      //  The function ImmediateInsert immediately inserts the object into the tree. 
      //  insert keeps the object in an internal store, but does not 
      //  immediately insert it. The object in the internal store are only inserted
      //  when CompleteDelayedInsert is called or when one of the search functions
      //  is invoked (they call CompleteDelayedInsert). When that is called, all
      //  of the stored objects are then inserted into the list in a way designed
      //  to give a relatively more balanced tree even if the data are strongly sorted.
//
//=======================================================================
    void insert ( const T& t )
    {
          m_ObjectStore    .push_back( t );
          m_DelayedIndices .push_back( (long)m_ObjectStore.size( ) - 1 );
      };

      //=======================================================================
      //  insert( const iterator& i, const T& t )
      //
      //  dummy here just for compatibility with vector and list, etc.
      //=======================================================================
      void insert( const iterator& /*i*/, const T& t )
      {
          insert( t );
    }

    //=======================================================================
      // insert ( const InputContainer& o )
      //
      // Function to insert a containerful for data into a CNearTree. Standard 
      // Library containers and CNearTree's can be used.
      //
      //  insert keeps the object in an internal store, but does not 
      //  immediately insert it. The object in the internal store are only inserted
      //  when CompleteDelayedInsert is called or when one of the search functions
      //  is invoked (they call CompleteDelayedInsert). When that is called, all
      //  of the stored objects are then inserted into the list in a way designed
      //  to give a relatively more balanced tree even if the data are strongly sorted.
      //
      //=======================================================================
      template< typename InputContainer >
      void insert ( const InputContainer& o )
    {
        size_t localDepth = 0;
          typename InputContainer::const_iterator it;

        for( it=o.begin(); it!=o.end(); ++it )
        {
              m_ObjectStore    .push_back( *it );
              m_DelayedIndices .push_back( (long)m_ObjectStore.size( ) - 1 );
        }
        m_DeepestDepth = std::max( localDepth, m_DeepestDepth );
    }

      //=======================================================================
      //  void ImmediateInsert ( const T& t )
      //
      //  Function to insert some "point" as an object into a CNearTree for
      //  later searching. Data is immediately put into the tree, instead of 
      //  being delayed (as function insert does). Use insert unless there is
      //  some known need to insert immediately.
      //
      //     t is an object of the templated type which is to be inserted into a
      //     NearTree
      //
      //  Three possibilities exist: put the datum into the left
      //  position (first test),into the right position, or else
      //  into a node descending from the nearer of those positions
      //  when they are both already used.
      //
    //=======================================================================
      void ImmediateInsert ( const T& t )
    {
        size_t localDepth = 0;
          m_BaseNode.Inserter( t, localDepth, m_ObjectStore );
        m_DeepestDepth = std::max( localDepth, m_DeepestDepth );
    }

    //=======================================================================
      // ImmediateInsert ( const InputContainer& o )
      //
      // see the description of ImmediateInsert above
      //
      //=======================================================================
      template< typename InputContainer >
      void ImmediateInsert ( const InputContainer& o )
    {
        size_t localDepth = 0;
          typename InputContainer::const_iterator it;

        for( it=o.begin(); it!=o.end(); ++it )
        {
            m_BaseNode.Inserter( *it, localDepth, m_ObjectStore );
        }
        m_DeepestDepth = std::max( localDepth, m_DeepestDepth );
    }
    
//=======================================================================
//  iterator NearestNeighbor ( const DistanceType &radius, const T& t ) const
//
//  Function to search a NearTree for the object closest to some probe point, t. This function
//  is only here so that the function Nearest can be called without having the radius const.
//  This was necessary because Nearest is recursive, but needs to keep the current smallest radius.
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    t  is the probe point
//
//    the return is an iterator to the templated type and is the returned nearest point
//             to the probe point (t) that can be found in the NearTree
//             or iterator::end if no point was found
//
//=======================================================================
    iterator NearestNeighbor ( const DistanceType& radius, const T& t ) const
    {
       T closest;
       size_t index = ULONG_MAX;
        DistanceType tempRadius = radius;
        const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if( this->empty( ) || radius < DistanceType( 0 ) )
        {
              return ( iterator(end( )) );
        }
        else if ( m_BaseNode.Nearest( tempRadius, closest, t, index, m_ObjectStore ) )
        {
            return ( iterator( (long)index, this ) );
       }
       else
       {
              return ( iterator(end( )) );
       }
    }

//=======================================================================
//  bool NearestNeighbor ( const DistanceType& dRadius,  T& tClosest,   const T& t ) const
//
//  Function to search a NearTree for the object closest to some probe point, t. This function
//  is only here so that the function Nearest can be called without having the radius const.
//  This was necessary because Nearest is recursive, but needs to keep the current smallest radius.
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
//    tClosest is an object of the templated type and is the returned nearest point
//             to the probe point that can be found in the NearTree
//    t  is the probe point
//
//    the return value is true only if a point was found
//
//=======================================================================
    bool NearestNeighbor ( const DistanceType& dRadius,  T& tClosest,   const T& t ) const
    {
          const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if ( dRadius < DistanceType(0) ) 
        {
            return ( false );
        }
        else if ( this->empty( ) )
        {
            return ( false );
        }
        else
        {
            DistanceType dSearchRadius = dRadius;
         size_t index = ULONG_MAX;
              return ( this->m_BaseNode.Nearest ( dSearchRadius, tClosest, t, index, m_ObjectStore ) );
        }
    }  //  NearestNeighbor

//=======================================================================
//  iterator FarthestNeighbor ( const T& t ) const
//
//  Function to search a NearTree for the object farthest from some probe point, t. This function
//  is only here so that the function Farthest can be called without having the radius const.
//  This was necessary because Farthest is recursive, but needs to keep the current largest radius.
//
//    t  is the probe point
//
//    the return is an iterator to the templated type and is the returned farthest point
//             from the probe point (t) that can be found in the NearTree
//             or iterator::end if no point was found
//
//=======================================================================
    iterator FarthestNeighbor ( const T& t ) const
    {
       T farthest;
       size_t index = ULONG_MAX;
        DistanceType radius = DistanceType( distMinValue );
        const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if( this->empty( ) )
        {
            return ( end( ) );
        }
        else if ( m_BaseNode.Farthest( radius, farthest, t, index, m_ObjectStore ) )
        {
            return ( iterator( (long)index, this ) );
       }
       else
       {
            return ( end( ) );
       }
    }
//=======================================================================
//  bool FarthestNeighbor ( T& tFarthest, const T& t ) const
//
//  Function to search a NearTree for the object farthest from some probe point, t. This function
//  is only here so that the function FarthestNeighbor can be called without the user
//  having to input a search radius and so the search radius can be guaranteed to be
//  negative at the start.
//
//    tFarthest is an object of the templated type and is the returned farthest point
//             from the probe point that can be found in the NearTree
//    t  is the probe point
//
//    the return value is true only if a point was found (should only be false for
//             an empty tree)
//
//=======================================================================
    bool FarthestNeighbor ( T& tFarthest, const T& t ) const
    {
          const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if ( this->empty( ) )
        {
            return ( false );
        }
        else
        {
            DistanceType dSearchRadius = DistanceType( distMinValue );
         size_t index = ULONG_MAX;
              return ( this->m_BaseNode.Farthest ( dSearchRadius, tFarthest, t, index, m_ObjectStore ) );
        }
    }  //  FarthestNeighbor
    
//=======================================================================
      //  long FindInSphere ( const DistanceType& dRadius,  OutputContainerType& tClosest,   const T& t ) const
//
//  Function to search a NearTree for the set of objects closer to some probe point, t,
      //  than dRadius. This is only here so that tClosest can be cleared before starting the work.
//
//    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
      //    tClosest is returned as a container of objects of the templated type and is the 
      //             returned set of nearest points to the probe point that can be found 
      //             in the NearTree. The container can be a Standard Library container or
      //             a CNearTree
//    t  is the probe point
      //
      // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
      template<typename OutputContainerType> 
      long FindInSphere ( const DistanceType& dRadius,  OutputContainerType& tClosest,   const T& t ) const
    {
        // clear the contents of the return vector so that things don't accidentally accumulate
        tClosest.clear( );
        const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if( this->empty( ) )
        {
           return( 0L );
        }
        else
        {
              return ( m_BaseNode.InSphere( dRadius, tClosest, t, m_ObjectStore ) );
        }
    }  //  FindInSphere
    
//=======================================================================
      //  long FindOutSphere ( const DistanceType& dRadius,  OutputContainerType& tFarthest,   const T& t ) const
//
//  Function to search a NearTree for the set of objects farther from some probe point, t,
//  than dRadius. This is only here so that tFarthest can be cleared before starting the work.
//
//    dRadius is the maximum search radius - any point nearer than dRadius from the probe
//             point will be ignored
      //    tFarthest is returned as a container of objects of the templated type and is the 
      //             returned set of nearest points to the probe point that can be found 
      //             in the NearTree. The container can be a Standard Library container or
      //             a CNearTree
//    t  is the probe point
      //
      // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
      template<typename OutputContainerType> 
      long FindOutSphere ( 
                            const DistanceType& dRadius,  
                            OutputContainerType& tFarthest,   
                            const T& t 
                         ) const
    {
        // clear the contents of the return vector so that things don't accidentally accumulate
        tFarthest.clear( );
          const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if( this->empty( ) )
        {
           return( 0L );
        }
        else
        {
            return ( m_BaseNode.OutSphere( dRadius, tFarthest, t, m_ObjectStore ) );
        }
    }  //  FindOutSphere


//=======================================================================
      //  long FindInAnnulus ( const DistanceType& dRadius1, const DistanceType dRadius2, OutputContainerType& tAnnular, const T& t ) const
//
//  Function to search a NearTree for the set of objects within a "spherical" annulus
//
//    dRadius1 is the minimum search radius - any point nearer than dRadius1 from the probe
//             point will be ignored
//    dRadius2 is the maximum search radius - any point farther than dRadius2 from the probe
//             point will be ignored
      //    tAnnular is returned as a container of objects of the templated type and is the 
      //             returned set of nearest points to the probe point that can be found 
      //             in the NearTree. The container can be a Standard Library container or
      //             a CNearTree
//    t  is the probe point
      //
      // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
      template<typename OutputContainerType> 
    long FindInAnnulus ( 
                       const DistanceType& dRadius1, 
                       const DistanceType& dRadius2, 
                              OutputContainerType& tAnnular, 
                              const T& t 
                         ) const
      {
          long lReturn = 0;
        // clear the contents of the return vector so that things don't accidentally accumulate
        tAnnular.clear( );
          const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if( this->empty( ) )
        {
        }
        else if ( dRadius1 > dRadius2 )
        {
            // Make sure that r1 < r2
            return ( FindInAnnulus( dRadius2, dRadius1, tAnnular, t ) );
        }
        else
        {
              lReturn = this->m_BaseNode.InAnnulus( dRadius1, dRadius2, tAnnular, t, m_ObjectStore );
          }

        return ( lReturn );
      }  //  FindInAnnulus

//=======================================================================
      //  long FindK_NearestNeighbors(  const size_t k, const DistanceType& dRadius, OutputContainerType& tClosest, const T& t ) const
//
      //  Function to search a NearTree for the set of objects closer to some probe point, t,
      //  than dRadius. This is only here so that tClosest can be cleared before starting the work
      //   and radius can be updated while processing.
//
      //    k is the maximum number of points to return
      //    dRadius is the maximum search radius - any point farther than dRadius from the probe
//             point will be ignored
      //    tClosest is returned as a container of objects of the templated type and is the 
      //             returned set of nearest points to the probe point that can be found 
      //             in the NearTree. The container can be a Standard Library container or
      //             a CNearTree
//    t  is the probe point
      //
      // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
      template<typename OutputContainerType> 
      long FindK_NearestNeighbors ( const size_t k, const DistanceType& radius,  OutputContainerType& tClosest,   const T& t )
      {
        // clear the contents of the return vector so that things don't accidentally accumulate
          tClosest.clear( );
          const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

        if( this->empty( ) )
        {
              return( 0L );
        }
        else
        {
              std::vector<std::pair<double, T> > K_Storage;
              double dRadius = radius;
              const long lFound = m_BaseNode.K_Near( k, dRadius, K_Storage, t, this->m_ObjectStore );
              for( unsigned int i=0; i<K_Storage.size( ); ++i )
              {
                  tClosest.insert( tClosest.end( ), K_Storage[i].second );
              }
              return( lFound );
          }
      }  //  FindK_NearestNeighbors
    
//=======================================================================
      //  long FindK_FarthestNeighbors const size_t k,OutputContainerType& tFarthest, const T& t ) const
//
      //  Function to search a NearTree for the set of objects farthest from some probe point, t. 
      //  This is only here so that tClosest can be cleared before starting the work.
//
      //    k is the maximum number of points to return
      //    tClosest is returned as a container of objects of the templated type and is the 
      //             returned set of nearest points to the probe point that can be found 
      //             in the NearTree. The container can be a Standard Library container or
      //             a CNearTree
      //    t  is the probe point
//
      // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
      template<typename OutputContainerType> 
      long FindK_FarthestNeighbors ( const size_t k, OutputContainerType& tFarthest,   const T& t )
      {
          // clear the contents of the return vector so that things don't accidentally accumulate
          tFarthest.clear( );
          const_cast<CNearTree*>(this)->CompleteDelayedInsert( );

          if( this->empty( ) )
          {
              return( 0L );
          }
          else
          {
              std::vector<std::pair<double, T> > K_Storage;
              double dRadius = 0;
              const long lFound = m_BaseNode.K_Far( k, dRadius, K_Storage, t, this->m_ObjectStore );
              for( unsigned int i=0; i<K_Storage.size( ); ++i )
              {
                  tFarthest.insert( tFarthest.end( ), K_Storage[i].second );
              }
              return( lFound );
          }
      }  //  FindK_FarthestNeighbors
    
//=======================================================================
//  void CompleteDelayedInsert ( void )
//
//  When CompleteDelayedInsert is invoked, if there are any objects in the 
//  delayed store they are then inserted into the neartree. CompleteDelayedInsert
//  randomly selects enough objects to partly fill a well-balanced tree. Specifically,
//  if there are n objects in the internal store, it randomly selects and inserts
//  sqrt(n) objects. After that, the remaining objects are inserted in a linear
//  sequence as they were entered.
//
//=======================================================================
    void CompleteDelayedInsert ( void )
    {
        if ( m_DelayedIndices.empty( ) )
        {
            return;
        }
        
        // insert a random selection of the objects
        const size_t vectorSize = m_DelayedIndices.size( );
      const size_t toRandomlyInsert = (size_t)::sqrt( (double)vectorSize );
        for ( size_t i=0; i<toRandomlyInsert; ++i )
        {
            size_t n = (size_t)((double)(vectorSize-1u) * (DistanceType)(random( )%MYRAND_MAX) / (double) MYRAND_MAX);
            random( ); random( );
            // Find the next pointer that hasn't already had its object "insert"ed
            // We can do this blindly since sqrt(n)<=n for all cases. n=1 would be the only 
            // bad case here, and that will not trigger the later loop.
            while ( m_DelayedIndices[n] == -1 )
            {
                ++n;
                n = n% vectorSize;
            }
            insertDelayed( (long)m_DelayedIndices[n] );
            m_DelayedIndices[n] = -1;         
        }
        
        // finish by inserting all the remaining objects
        for ( size_t i=0; i<vectorSize; ++i )
        {
            if ( m_DelayedIndices[i] != -1 )
            {
                insertDelayed( (long)m_DelayedIndices[i] );
            }
        }
        
        // now get rid of the temporary storage that was used for delayed 
        // insertions (fast way, faster than clear() )
        std::vector<long> DelayedPointersTemp;   
        DelayedPointersTemp  .swap( m_DelayedIndices );
    };
//=======================================================================
//  void CompleteDelayedInsertRandom ( void )
//
//  When CompleteDelayedInsertRandom is invoked, if there are any objects in the 
//  delayed store they are then inserted into the neartree. CompleteDelayedInsertRandom
      //  randomly selects all objects and inserts them. m_DelayedIndices is empty after
      //  a call to CompleteDelayedInsertRandom.
//
//=======================================================================
    void CompleteDelayedInsertRandom ( void )
    {
        if ( m_DelayedIndices.empty( ) )
        {
            return;
        }

        // insert a random selection of the objects
        const size_t vectorSize = m_DelayedIndices.size( );
        for ( size_t i=0; i<m_DelayedIndices.size( ); ++i )
        {
            size_t n = (size_t)((double)(vectorSize-1u) * (DistanceType)(random( )%MYRAND_MAX) / (double) MYRAND_MAX);
            random( ); random( );
            // Find the next pointer that hasn't already had its object "insert"ed
            // We can do this blindly since sqrt(n)<=n for all cases. n=1 would be the only 
            // bad case here, and that will not trigger the later loop.
            while ( m_DelayedIndices[n] == -1 )
            {
                ++n;
                n = n% vectorSize;
            }
            insertDelayed( (long)m_DelayedIndices[n] );
            m_DelayedIndices[n] = -1;         
        }

        // now get rid of the temporary storage that was used for delayed 
        // insertions (fast way, faster than clear() )
        std::vector<long> DelayedPointersTemp;   
        DelayedPointersTemp.swap( m_DelayedIndices );
    };
    
//=======================================================================
//  size_t GetDeferredSize (  void )
//
//  The number of objects currently queued for insertion.
//
//=======================================================================
    size_t GetDeferredSize ( void ) const
    {
        return ( m_DelayedIndices.size( ) );
    };
    
//=======================================================================
//  size_t GetTotalSize (  void )
//
//  The total number of objects that have been inserted plus those
//  queued for insertion.
//
//=======================================================================
    size_t GetTotalSize ( void ) const
    {
        return ( m_ObjectStore.size( ) );
    };

//=======================================================================
//  size_t size ( void )
//
//  The total number of objects that have been inserted plus those
//  queued for insertion.
//
//=======================================================================
    size_t size ( void ) const
    {
        return ( GetTotalSize( ) );
   };

//=======================================================================
//  size_t GetDepth ( void ) const
//
//  The greatest depth of the tree (1-based) from the root.
//
//=======================================================================
    size_t GetDepth ( void ) const
    {
        return ( m_DeepestDepth );
    };

//=======================================================================
      //  operator ContainerType ( void ) const
//
      //  Utility function to copy the data object to a user's container object.
//
//=======================================================================
      template<typename ContainerType>
      operator ContainerType ( void ) const
      {
        return ( m_ObjectStore );
    }
    
public:
    iterator begin ( void ) { return ( iterator( 0, this ) ); };
    iterator end   ( void ) { return ( iterator( m_ObjectStore.empty( )? 1 :(long)m_ObjectStore.size( )  , this ) ); };
    iterator back  ( void ) { return ( iterator( m_ObjectStore.empty( )? 1 :(long)m_ObjectStore.size( )-1, this ) ); };
    const_iterator begin ( void ) const { return ( const_iterator( 0, this ) ); };
    const_iterator end   ( void ) const { return ( const_iterator( m_ObjectStore.empty( )? 1 :(long)m_ObjectStore.size( )  , this ) ); };
    const_iterator back  ( void ) const { return ( const_iterator( m_ObjectStore.empty( )? 1 :(long)m_ObjectStore.size( )-1, this ) ); };

    T at( const size_t n ) const { return ( m_ObjectStore[n] ); };
    T operator[] ( const size_t position ) const { return ( m_ObjectStore[position] ); };

    
private:

    //=======================================================================
        void insertDelayed ( const long n )
        {
            size_t localDepth = 0;
            m_BaseNode.InserterDelayed( n, localDepth, m_ObjectStore );
            m_DeepestDepth = std::max( localDepth, m_DeepestDepth );
        }

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//  end of CNearTree
//=======================================================================
    template<typename U> class KT : U
    {
    private:
        DistanceType d;
        KT( void ) : U() { };
    };



// start of nested class NearTreeNode
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

//=======================================================================
//
// NEARTREENODE - nested class to hold the actual neartree and the indices to its data objects
//
// For the description of the template parameters, see the text above in the description
// of CNearTree.
//=======================================================================
    template <typename TNode, typename DistanceTypeNode=double, int distMinValueNode=-1 > 
    class NearTreeNode
    {
        size_t            m_ptLeft;            // index of left object (of type TNode) stored in this node
        size_t            m_ptRight;           // index of right object (of type TNode) stored in this node
        DistanceTypeNode      m_dMaxLeft;      // longest distance from the left object to
                                               // anything below it in the tree
        DistanceTypeNode      m_dMaxRight;     // longest distance from the right object to
                                               // anything below it in the tree
        NearTreeNode *    m_pLeftBranch;       // tree descending from the left object
        NearTreeNode *    m_pRightBranch;      // tree descending from the right object

    public:

        NearTreeNode( void ) :  //  NearTreeNode constructor
               m_ptLeft            ( ULONG_MAX ),
               m_ptRight           ( ULONG_MAX ),
               m_dMaxLeft          ( DistanceTypeNode( distMinValueNode ) ),
               m_dMaxRight         ( DistanceTypeNode( distMinValueNode ) ),
               m_pLeftBranch       ( 0 ),
               m_pRightBranch      ( 0 )
        {
        };  //  NearTreeNode constructor

//=======================================================================
        ~NearTreeNode( void )  //  NearTreeNode destructor
        {
            clear( );
          };  //  end NearTreeNode destructor

//=======================================================================
        void clear( void )
        {
            if ( m_pLeftBranch  ) m_pLeftBranch ->clear( );
            if ( m_pRightBranch ) m_pRightBranch->clear( );
            delete m_pLeftBranch  ;  m_pLeftBranch  =0;
            delete m_pRightBranch ;  m_pRightBranch =0;

            m_ptLeft     = ULONG_MAX;
            m_ptRight    = ULONG_MAX;
            m_dMaxLeft   = DistanceTypeNode( distMinValueNode );
            m_dMaxRight  = DistanceTypeNode( distMinValueNode );
          };  //  end clear

    //=======================================================================
    //  void Inserter ( const TNode& t( const TNode& t, size_t& localDepth )
//
//  Function to insert some "point" as an object into a CNearTree for
//  later searching
//
//     t is an object of the templated type which is to be inserted into a
    //     NearTree
//
//     localDepth is the returned deepest tree level reached for the current insert
//
//  Three possibilities exist: put the datum into the left
    //  position (first test),into the right position, or else
//  into a node descending from the nearer of those positions
//  when they are both already used.
//
//=======================================================================
        void Inserter ( const TNode& t, size_t& localDepth, std::vector<TNode>& objectStore )
    {
        // do a bit of precomputing if it is possible so that we can
            // reduce the number of calls to operator 'DistanceTypeNode' as much as possible;
            // 'DistanceTypeNode' might use square roots in some cases
            DistanceTypeNode dTempRight =  DistanceTypeNode(0);
            DistanceTypeNode dTempLeft  =  DistanceTypeNode(0);
        ++localDepth;
        
      if ( m_ptRight  != ULONG_MAX )
      {
                dTempRight  = DistanceBetween( t, objectStore[m_ptRight] );
                dTempLeft   = DistanceBetween( t, objectStore[m_ptLeft]  );
        }
        
      if ( m_ptLeft == ULONG_MAX )
      {
                m_ptLeft = objectStore.size( );
                objectStore.push_back( t );
      }
      else if ( m_ptRight == ULONG_MAX )
      {
                m_ptRight = objectStore.size( );
                objectStore.push_back( t );
        }
        else if ( dTempLeft > dTempRight )
        {
                if ( m_pRightBranch == 0 ) m_pRightBranch = new NearTreeNode;
            // note that the next line assumes that m_dMaxRight is negative for a new node
            if ( m_dMaxRight < dTempRight ) m_dMaxRight = dTempRight;
                m_pRightBranch->Inserter( t, localDepth, objectStore );
            }
            else  // ((DistanceTypeNode)(t - *m_tLeft) <= (DistanceTypeNode)(t - *m_tRight) )
            {
                if ( m_pLeftBranch  == 0 ) m_pLeftBranch  = new NearTreeNode;
            // note that the next line assumes that m_dMaxLeft is negative for a new node
            if ( m_dMaxLeft < dTempLeft ) m_dMaxLeft  = dTempLeft;
                m_pLeftBranch->Inserter( t, localDepth, objectStore );
        }
          }  //  end Inserter

//=======================================================================
        void InserterDelayed ( const long n, size_t& localDepth, std::vector<TNode>& objectStore )
        {
            // do a bit of precomputing if it is possible so that we can
            // reduce the number of calls to operator 'DistanceTypeNode' as much as possible;
            // 'DistanceTypeNode' might use square roots in some cases
            DistanceTypeNode dTempRight =  DistanceTypeNode(0);
            DistanceTypeNode dTempLeft  =  DistanceTypeNode(0);
            ++localDepth;

            if ( m_ptRight  != ULONG_MAX )
            {
                dTempRight  = DistanceBetween( objectStore[n], objectStore[m_ptRight] );
                dTempLeft   = DistanceBetween( objectStore[n], objectStore[m_ptLeft]  );
            }

            if ( m_ptLeft == ULONG_MAX )
            {
                m_ptLeft = n;
            }
            else if ( m_ptRight == ULONG_MAX )
            {
                m_ptRight = n;
            }
            else if ( dTempLeft > dTempRight )
            {
                if ( m_pRightBranch == 0 ) m_pRightBranch = new NearTreeNode;
                // note that the next line assumes that m_dMaxRight is negative for a new node
                if ( m_dMaxRight < dTempRight ) m_dMaxRight = dTempRight;
                m_pRightBranch->InserterDelayed( n, localDepth, objectStore );
            }
            else  // ((DistanceTypeNode)(t - *m_tLeft) <= (DistanceTypeNode)(t - *m_tRight) )
            {
                if ( m_pLeftBranch  == 0 ) m_pLeftBranch  = new NearTreeNode;
                // note that the next line assumes that m_dMaxLeft is negative for a new node
                if ( m_dMaxLeft < dTempLeft ) m_dMaxLeft  = dTempLeft;
                m_pLeftBranch->InserterDelayed( n, localDepth, objectStore );
            }
          }  //   end InserterDelayed

    //=======================================================================
    //  bool Nearest ( DistanceTypeNode& dRadius,  TNode& tClosest,   const TNode& t ) const
//
    //  Private function to search a NearTree for the object closest to some probe point, t.
    //  This function is only called by NearestNeighbor.
//
    //    dRadius is the smallest currently known distance of an object from the probe point.
    //    tClosest is an object of the templated type and is the returned closest point
    //             to the probe point that can be found in the NearTree
//    t  is the probe point
          //
    //    the return value is true only if a point was found within dRadius
//
   //=======================================================================
        bool Nearest ( 
                             DistanceTypeNode& dRadius, 
                             TNode& tClosest, 
                             const TNode& t, 
                             size_t& pClosest, 
                             const std::vector<TNode>& objectStore 
                     ) const
        {
            std::vector <NearTreeNode* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
            NearTreeNode* pt = const_cast<NearTreeNode*>(this);
            pClosest = ULONG_MAX;
            if ( pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
                    const DistanceTypeNode dDR = DistanceBetween( t, objectStore[pt->m_ptRight] );
                if ( dDR <= dRadius )
                {
                        dRadius = dDR;
                        pClosest = pt->m_ptRight;
                    }
                    if ( pt->m_pRightBranch != 0 && TRIANG(dDR,pt->m_dMaxRight,dRadius))
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
                    const DistanceTypeNode dDL = DistanceBetween( t, objectStore[pt->m_ptLeft]  );
                if ( dDL <= dRadius )
                {
                        dRadius = dDL;
                        pClosest = pt->m_ptLeft;
                }
                if ( pt->m_ptRight != ULONG_MAX ) // only stack if there's a right object
                {
                    sStack.push_back( pt );
                }
                    if ( pt->m_pLeftBranch != 0 && TRIANG(dDL,pt->m_dMaxLeft,dRadius))
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
            if ( !sStack.empty( ) ) // for safety !!!
            {
                std::vector <NearTreeNode* > sTemp;
                sTemp.swap( sStack );
            }
            if ( pClosest != ULONG_MAX )
                tClosest = objectStore[pClosest];
            return ( pClosest != ULONG_MAX );
          };   // end Nearest

//=======================================================================
    //  bool Farthest ( DistanceTypeNode& dRadius,  TNode& tFarthest,   const TNode& t ) const
//
    //  Private function to search a NearTree for the object farthest from some probe point, t.
    //  This function is only called by FarthestNeighbor.
//
    //    dRadius is the largest currently known distance of an object from the probe point.
    //    tFarthest is an object of the templated type and is the returned farthest point
    //             from the probe point that can be found in the NearTree
//    t  is the probe point
          //
    //    the return value is true only if a point was found (should only be false for
    //             an empty tree)
//
//=======================================================================
        bool Farthest ( 
                            DistanceTypeNode& dRadius,  
                            TNode& tFarthest,   
                            const TNode& t, size_t& pFarthest, 
                            const std::vector<TNode>& objectStore 
                      ) const
        {
            std::vector <NearTreeNode* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
            NearTreeNode* pt = const_cast<NearTreeNode*>(this);
            pFarthest = ULONG_MAX;
      if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
                    const DistanceTypeNode dDR = DistanceBetween( t , objectStore[pt->m_ptRight] );
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
                    const DistanceTypeNode dDL = DistanceBetween( t , objectStore[pt->m_ptLeft] );
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
            if ( !sStack.empty( ) ) // for safety !!!
            {
                std::vector <NearTreeNode* > sTemp;
                sTemp.swap( sStack );
            }
            if ( pFarthest != ULONG_MAX )
                tFarthest = objectStore[pFarthest];
            return ( pFarthest != ULONG_MAX );
          };   //  end Farthest

          //=======================================================================
    //  long InSphere ( const DistanceTypeNode dRadius,  CNearTree<  TNode >& tClosest,   const TNode& t ) const
//
    //  Private function to search a NearTree for the objects inside of the specified radius
//     from the probe point
    //  This function is only called by FindInSphere.
//
//    dRadius is the search radius
    //    tClosest is a CNearTree of objects of the templated type found within dRadius of the
//         probe point
//    t  is the probe point
          //
          // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
          template<typename ContainerType>
        long InSphere ( 
                               const DistanceTypeNode& dRadius, 
                              ContainerType& tClosest, 
                               const TNode& t,
                               const std::vector<TNode>& objectStore
                      ) const
        {
            std::vector <NearTreeNode* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
            NearTreeNode* pt = const_cast<NearTreeNode*>(this);
        if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
                    const DistanceTypeNode dDR =  DistanceBetween( t, objectStore[pt->m_ptRight] );
                    if ( dDR <= dRadius )
                {
                          tClosest.insert( tClosest.end(), objectStore[pt->m_ptRight] );
                    }
                    if ( pt->m_pRightBranch != 0 && TRIANG(dDR,pt->m_dMaxRight,dRadius) )
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
                    const DistanceTypeNode dDL = DistanceBetween( t, objectStore[pt->m_ptLeft]  );
                    if ( dDL <= dRadius )
                {
                          tClosest.insert( tClosest.end(), objectStore[pt->m_ptLeft] );
                }
                if ( pt->m_ptRight != ULONG_MAX ) // only stack if there's a right object
                {
                    sStack.push_back( pt );
                }
                    if ( pt->m_pLeftBranch != 0 && TRIANG(dDL,pt->m_dMaxLeft,dRadius) )
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

              return ( (long)tClosest.size() );
          }  //  end InSphere

//=======================================================================
    //  long OutSphere ( const DistanceTypeNode& dRadius,  CNearTree<  TNode >& tFarthest,   const TNode& t ) const
//
    //  Private function to search a NearTree for the objects outside of the specified radius
    //     from the probe point
    //  This function is only called by FindOutSphere.
//
    //    dRadius is the search radius
    //    tFarthest is a CNearTree of objects of the templated type found within dRadius of the
    //         probe point
//    t  is the probe point
          //
          // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
          template<typename ContainerType>
        long OutSphere ( 
                            const DistanceTypeNode& dRadius, 
                              ContainerType& tFarthest, 
                            const TNode& t,
                            const std::vector<TNode> objectStore
                       ) const
        {
            std::vector <NearTreeNode* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
            NearTreeNode* pt = const_cast<NearTreeNode*>(this);
        if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
                    const DistanceTypeNode dDR = DistanceBetween( t, objectStore[pt->m_ptRight] );
                    if ( dDR >= dRadius )
                {
                          tFarthest.insert( tFarthest.end(), objectStore[pt->m_ptRight] );
                    }
                    if ( pt->m_pRightBranch != 0 && TRIANG(dRadius,dDR,pt->m_dMaxRight) )
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
                    const DistanceTypeNode dDL = DistanceBetween( t, objectStore[pt->m_ptLeft]  );
                    if ( dDL >= dRadius )
                {
                          tFarthest.insert( tFarthest.end(), objectStore[pt->m_ptLeft] );
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

              return ( (long)tFarthest.size() );
          }  //  end OutSphere

          //=======================================================================
          //  long InAnnulus ( const DistanceTypeNode& dRadius1, const DistanceTypeNode& dRadius2, CNearTree< TNode >& tAnnular,   const TNode& t ) const
//
          //  Private function to search a NearTree for the objects within a specified annulus from probe point
          //  This function is only called by FindInAnnulus.
//
          //    dRadius1, dRadius2 specifies the range of the annulus
          //    tAnnular is a NearTree of objects of the templated type found between the two radii
//    t  is the probe point
          //
          // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
//=======================================================================
          template<typename ContainerType>
          long InAnnulus ( 
                              const DistanceTypeNode& dRadius1, 
                              const DistanceTypeNode& dRadius2, 
                              ContainerType& tAnnular, 
                             const TNode& t,
                             const std::vector<TNode> objectStore
                       ) const
        {
            std::vector <NearTreeNode* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
            NearTreeNode* pt = const_cast<NearTreeNode*>(this);
        if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
                    const DistanceTypeNode dDR = DistanceBetween( t, objectStore[pt->m_ptRight] );
                      if ( dDR <= dRadius2 && dDR >= dRadius1 )
                      {
                          tAnnular.insert( tAnnular.end( ), objectStore[pt->m_ptRight] );
                      }
            if ( pt->m_pRightBranch != 0 && (TRIANG(dRadius1,dDR,pt->m_dMaxRight)) && (TRIANG(dDR,pt->m_dMaxRight,dRadius2) ) )
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
                    const DistanceTypeNode dDL = DistanceBetween( t, objectStore[pt->m_ptLeft]  );
                      if ( dDL <= dRadius2 && dDL >= dRadius1 )
                      {
                          tAnnular.insert( tAnnular.end(), objectStore[pt->m_ptLeft] );
                }
                if ( pt->m_ptRight != ULONG_MAX ) // only stack if there's a right object
                {
                    sStack.push_back( pt );
                }
            if ( pt->m_pLeftBranch != 0 && (TRIANG(dRadius1,dDL,pt->m_dMaxLeft)) && (TRIANG(dDL,pt->m_dMaxLeft,dRadius2)  ) )
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

              return ( (long)tAnnular.size() );
          }  // end InAnnulus

//=======================================================================
          //  long K_Near ( const DistanceTypeNode dRadius,  std::vector<T>& tClosest,   const TNode& t ) const
          //
          //  Private function to search a NearTree for the objects inside of the specified radius
          //     from the probe point
          //  This function is only called by FindK_Nearest.
//
          // k:           the maximum number of object to return, giving preference to the nearest
          // dRadius:     the search radius, which will be updated when the internal store is resized 
          // tClosest:    is a vector of objects of the templated type found within dRadius of the
          //                 probe point, limited by the k-near search
          // t:           is the probe point
          // objectStore: the internal vector storing the object in CNearTree
//
          // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
          /*=======================================================================*/
          long K_Near ( 
                              const size_t k,
                              DistanceTypeNode& dRadius, 
                              std::vector<std::pair<double,T> >& tClosest, 
                              const TNode& t,
                              const std::vector<TNode>& objectStore
                        )
        {
            std::vector <NearTreeNode* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
            NearTreeNode* pt = const_cast<NearTreeNode*>(this);
            if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
                      const DistanceTypeNode dDR =  DistanceBetween( t, objectStore[pt->m_ptRight] );
                      if ( dDR <= dRadius )
                      {
                          tClosest.insert( tClosest.end(), std::make_pair( dDR, objectStore[pt->m_ptRight] ) );
                          if( tClosest.size( ) > 2*k ) K_Resize( k, t, tClosest, dRadius );
                      }
                      if ( pt->m_pRightBranch != 0 && TRIANG(dDR,pt->m_dMaxRight,dRadius) )
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
                    const DistanceTypeNode dDL = DistanceBetween( t, objectStore[pt->m_ptLeft]  );
                      if ( dDL <= dRadius )
                      {
                          tClosest.insert( tClosest.end(), std::make_pair( dDL, objectStore[pt->m_ptLeft] ) );
                          if( tClosest.size( ) > 2*k ) K_Resize( k, t, tClosest, dRadius );
                }
            if ( pt->m_ptRight != ULONG_MAX ) // only stack if there's a right object
                {
                    sStack.push_back( pt );
                }
                      if ( pt->m_pLeftBranch != 0 && TRIANG(dDL,pt->m_dMaxLeft,dRadius) )
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

              if( tClosest.size( ) > k ) K_Resize( k, t, tClosest, dRadius );
              return ( (long)tClosest.size( ) );
          }  // end K_Near

//=======================================================================
          //  long K_Far ( const DistanceTypeNode dRadius, std::vector<T>& tFarthest, const TNode& t ) const
          //
          //  Private function to search a NearTree for the objects inside of the specified radius
          //     from the probe point. Distances are stored in an intermediate array as negative values
          //     so that the same logic as K_Near can be used.
          //  This function is only called by FindK_Farthest.
//
          // k:           the maximum number of object to return, giving preference to the nearest
          // dRadius:     the search radius, which will be updated when the internal store is resized 
          // tClosest:    is a vector of objects of the templated type found within dRadius of the
          //                 probe point, limited by the k-near search
          // t:           is the probe point
          // objectStore: the internal vector storing the object in CNearTree
//
          // returns the number of objects returned in the container (for sets, that may not equal the number found)
//
          /*=======================================================================*/
          long K_Far ( 
                              const size_t k,
                              DistanceTypeNode& dRadius, 
                              std::vector<std::pair<double,T> >& tFarthest, 
                              const TNode& t,
                              const std::vector<TNode>& objectStore
                        )
        {
            std::vector <NearTreeNode* > sStack;
        enum  { left, right, end } eDir;
        eDir = left; // examine the left nodes first
            NearTreeNode* pt = const_cast<NearTreeNode*>(this);
      if (pt->m_ptLeft == ULONG_MAX) return false; // test for empty
        while ( ! ( eDir == end && sStack.empty( ) ) )
        {
            if ( eDir == right )
            {
                      const DistanceTypeNode dDR =  DistanceBetween( t, objectStore[pt->m_ptRight] );
                      if ( dDR >= dRadius )
                      {
                          tFarthest.insert( tFarthest.end(), std::make_pair( -dDR, objectStore[pt->m_ptRight] ) );
                          if( tFarthest.size( ) > 2*k ) K_Resize( k, t, tFarthest, dRadius );
                      }
                      if ( pt->m_pRightBranch != 0 && TRIANG(dRadius,dDR,pt->m_dMaxRight) )
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
                    const DistanceTypeNode dDL = DistanceBetween( t, objectStore[pt->m_ptLeft]  );
                      if ( dDL >= dRadius )
                      {
                          tFarthest.insert( tFarthest.end(), std::make_pair( -dDL, objectStore[pt->m_ptLeft] ) );
                          if( tFarthest.size( ) > 2*k ) K_Resize( k, t, tFarthest, dRadius );
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

              if( tFarthest.size( ) > k ) K_Resize( k, t, tFarthest, dRadius );
              return ( (long)tFarthest.size( ) );
          }  //  end K_Far


          //=======================================================================
          //  void K_Resize ( const size_t k, const TNode& t, std::vector<T>& tClosest, double& dRadius ) const
          //
          //  Private function to limit the size of internally stored data for K-nearest/farthest-neighbor searches
          //  This function is only called by K_Near and K_Far.
          //
          //    dRadius is the search radius, updated to the best-known value
          //    tClosest is a CNearTree of objects of the templated type found within dRadius of the
          //         probe point
          //    t  is the probe point
          //
          //=======================================================================
          void K_Resize( const size_t k, const TNode& t, std::vector<std::pair<double, T> >& tClosest, double& dRadius )
          {
              std::sort( tClosest.begin(), tClosest.end() );
              tClosest.resize( k );
              dRadius = DistanceBetween( t, tClosest[tClosest.size()-1].second );
          }  // end K_Resize



    }; // end NearTreeNode
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
// end NearTreeNode nested class in CNearTree
//=======================================================================
// start of iterator, a nested class in CNearTree
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

      //====================================================================================
    // class iterator
    //
    // nested class within CNearTree
    //====================================================================================
public:
      class iterator
      {
    public:
        friend class CNearTree< T, DistanceType, distMinValue >;
    private:
         long position;
        const CNearTree< T, DistanceType, distMinValue >* parent;

      public:
         iterator( void ) { }; // constructor
        iterator( const const_iterator& s ) { position = ((const_iterator)s).get_position(); parent = ((const_iterator)s).get_parent(); };// constructor

        iterator& operator=  ( const iterator& s )       { position = s.position; parent = s.parent; return ( *this ); };
        iterator& operator=  ( const const_iterator& s ) { position = s.position; parent = s.parent; return ( *this ); };
        iterator  operator++ ( const int n )             { iterator it(*this); position+=1+n; return ( it ); };
        iterator  operator-- ( const int n )             { iterator it(*this); position-=1+n; return ( it ); };
        iterator& operator++ ( void )                    { ++position; return ( *this ); };
        iterator& operator-- ( void )                    { --position; return ( *this ); };
        iterator  operator+  ( const long n ) const      { iterator it( position+n, parent); return ( it ); };
        iterator  operator-  ( const long n ) const      { iterator it( position-n, parent); return ( it ); };
        iterator& operator+= ( const long n )            { position += n; return ( *this ); };
        iterator& operator-= ( const long n )            { position -= n; return ( *this ); };
        T         operator*  ( void )         const      { return ( parent->m_ObjectStore[position] ); };

        bool      operator== ( const iterator& t ) const { return ( t.position==(parent->m_ObjectStore.empty( )?1:position) && t.parent==parent ); };
        bool      operator!= ( const iterator& t ) const { return ( ! (*this==t )); };
        bool      operator== ( const const_iterator& t ) const { return ( t.position==(parent->m_ObjectStore.empty( )?1:position) && t.parent==parent ); };
        bool      operator!= ( const const_iterator& t ) const { return ( ! (*this==t )); };

        const T * const operator-> ( void )   const      { return ( &(const_cast<CNearTree*>(parent)->m_ObjectStore[position]) ); };
        long get_position( void ) {return position;};
        const CNearTree< T, DistanceType, distMinValue >* get_parent( void ) {return parent;};

      private:
        iterator ( const long s, const CNearTree* const nt ) { position = s; parent = nt; }; // constructor

      }; // class iterator
    //====================================================================================
    // end of nested class "iterator"
    //====================================================================================

    class const_iterator
    {
    public:
        friend class CNearTree< T, DistanceType, distMinValue >;
    private:
        long position;
        const CNearTree< T, DistanceType, distMinValue >* parent;

    public:
        const_iterator( void ) { }; // constructor
        const_iterator( const iterator& s ) { position = ((iterator)s).get_position(); parent = ((iterator)s).get_parent(); }; // constructor

        const_iterator& operator=  ( const const_iterator& s ) { position = s.position; parent = s.parent; return ( *this ); };
        const_iterator& operator=  ( const       iterator& s ) { position = s.position; parent = s.parent; return ( *this ); };
        const_iterator  operator++ ( const int n )             { const_iterator it(*this); position+=1+n; return ( it ); };
        const_iterator  operator-- ( const int n )             { const_iterator it(*this); position-=1+n; return ( it ); };
        const_iterator& operator++ ( void )                    { ++position; return ( *this ); };
        const_iterator& operator-- ( void )                    { --position; return ( *this ); };
        const_iterator  operator+  ( const long n ) const      { const_iterator it( position+n, parent); return ( it ); };
        const_iterator  operator-  ( const long n ) const      { const_iterator it( position-n, parent); return ( it ); };
        const_iterator& operator+= ( const long n )            { position += n; return ( *this ); };
        const_iterator& operator-= ( const long n )            { position -= n; return ( *this ); };
        T               operator*  ( void )         const      { return ( parent->m_ObjectStore[position] ); };
                 
        bool            operator== ( const const_iterator& t ) const { return ( t.position==(parent->m_ObjectStore.empty( )?1:position) && t.parent==parent ); };
        bool            operator!= ( const const_iterator& t ) const { return ( ! (*this==t )); };
        bool            operator== ( const iterator& t ) const { return ( t.position==(parent->m_ObjectStore.empty( )?1:position) && t.parent==parent ); };
        bool            operator!= ( const iterator& t ) const { return ( ! (*this==t )); };

        const T * const operator-> ( void )   const      { return ( &(const_cast<CNearTree*>(parent)->m_ObjectStore[position]) ); };
        long get_position          ( void ) {return position;};
        const CNearTree< T, DistanceType, distMinValue >* get_parent( void ) {return parent;};

    private:
        const_iterator ( const long s, const CNearTree* nt ) { position = s; parent = nt; }; // constructor

    }; // end class const_iterator
    //====================================================================================
    // end of nested class "const_iterator"
    //====================================================================================


}; // template class CNearTree

#endif // !defined(TNEAR_H_INCLUDED)
