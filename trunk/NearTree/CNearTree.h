/*
 *  CNearTree.h
 *  NearTree
 *
 *  Based on TNear.h C++ Template
 *  Copyright 2001 Larry Andrews.  All rights reserved
 *
 *  C Version created by Herbert J. Bernstein on 11/29/08
 *  with permission from Larry Andrews.
 *  Copyright 2008 Larry Andrews and Herbert J. Bernstein. 
 *  All rights reserved.
 *
 */

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE CNearTree API UNDER THE TERMS OF THE LGPL *
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

/* Notices from original C++ template:
 
 Nearest Neighbor algorithm after Kalantari and McDonald,
 (IEEE Transactions on Software Engineering, v. SE-9, pp.
 631-634,1983)
 modified to use recursion instead of a double-linked tree
 and simplified so that it does a bit less checking for
 things like is the distance to the right less than the
 distance to the left; it was found that these checks little
 to no difference.
 
 copyright by Larry Andrews, 2001
 may be freely distributed or used as long as this copyright notice
 is included
 */

#ifndef CNEARTREE_H_INCLUDED
#define CNEARTREE_H_INCLUDED

#ifdef __cplusplus

extern "C" {
    
#endif
    
#ifdef CNEARTREE_USE_FAR
#include <malloc.h>
#define FAR __far
#define MALLOC _fmalloc
#define FREE _ffree
#define MEMSET _fmemset
#define MEMMOVE _fmemmove
#else
#include <stdlib.h>
#define FAR
#define MALLOC malloc
#define FREE free
#define MEMSET memset
#define MEMMOVE memmove
#endif
    
#include <limits.h>
#include <float.h>
#include <math.h>
#ifndef USE_LOCAL_HEADERS
#include <CVector.h>
#else
#include "CVector.h"
#endif    
    
#ifdef USE_MINGW_RAND
#define random(x) rand(x)
#define srandom(x) srand(x)
#endif
    
    
    /* function returns */
#define CNEARTREE_SUCCESS       0
#define CNEARTREE_MALLOC_FAILED 1
#define CNEARTREE_BAD_ARGUMENT  2
#define CNEARTREE_NOT_FOUND     4
#define CNEARTREE_FREE_FAILED   8
#define CNEARTREE_CVECTOR_FAILED   16

    
    /* flags 
     0 for n data, n children
     */
    
#define CNEARTREE_FLAG_LEFT_DATA   1          /* 0x0001 */
#define CNEARTREE_FLAG_RIGHT_DATA  2          /* 0x0002 */
#define CNEARTREE_FLAG_LEFT_CHILD  4          /* 0x0004 */
#define CNEARTREE_FLAG_RIGHT_CHILD 8          /* 0x0008 */
    
#define CNEARTREE_DATA_OR_CHILDREN 15         /* 0x000F */
    
#define CNEARTREE_TYPE_DOUBLE      16         /* 0x0010 */
#define CNEARTREE_TYPE_INTEGER     32         /* 0x0020 */
    
#define CNEARTREE_TYPE             48         /* 0x0030 */
    
#define CNEARTREE_NORM_UNKNOWN     64         /* 0x0040 */
#define CNEARTREE_NORM_L1         128         /* 0x0080 */
#define CNEARTREE_NORM_L2         256         /* 0x0100 */
#define CNEARTREE_NORM_LINF       512         /* 0x0200 */
    
#define CNEARTREE_NORM            960         /* 0x03C0 */
    
#define CNEARTREE_FLIP           1024         /* 0x0400 */
#define CNEARTREE_DEFER_ALL      2048         /* 0x0800 */

    
    
    typedef struct _CNearTreeNode {
        size_t           m_indexLeft;    /* index of left coords in m_CoordStore  
                                            and of left object in m_ObjectStore     */
        size_t           m_indexRight;   /* index of right coords in m_CoordStore 
                                            and of right object in m_ObjectStore     */
        double           m_dMaxLeft;     /* longest distance from the left object
                                            to anything below it in the tree            */
        double           m_dMaxRight;    /* longest distance from the right object 
                                            to anything below it in the tree            */
        struct _CNearTreeNode FAR * m_pLeftBranch;  
                                         /* tree descending from the left object        */
        struct _CNearTreeNode FAR * m_pRightBranch; 
                                         /* tree descending from the right object       */
        int              m_iflags;       /* flags                                       */
    } CNearTreeNode;
    
    
    typedef CNearTreeNode FAR * CNearTreeNodeHandle;   
    
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
    
    /*
     =======================================================================
     double CNearTreeDistsq(void FAR * coord1, 
                            void FAR * coord2,  
                            size_t treedim, 
                            int treetype)
     
     function to return the square of the Euclidean distance between two 
     coordinate vectors.  
     
     THIS FUNCTION IS DEPRECATED
     
     treedim -- the dimension of the vectors

     treetype -- and integer flag for type of the vectors
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer
     
     =======================================================================
     */
    
     double CNearTreeDistsq(void FAR * coord1,
                            void FAR * coord2, 
                            size_t treedim, 
                            int treetype);
    /*
     =======================================================================
     double CNearTreeDist(CNearTreeHandle treehandle, void FAR * coord1, 
                                      void FAR * coord2)
     
     function to return the distance (L1, L2 or L-infinity) between two 
     coordinate vectors according to the parameters of the given tree  
     
     =======================================================================
     */
     double CNearTreeDist(CNearTreeHandle treehandle, void FAR * coord1, 
                         void FAR * coord2);
     
     /*
     =======================================================================
     int CNearTreeSetNorm(CNearTreeHandle treehandle, int treenorm);
     
     function to set the norm to use used for this tree
     
     treenorm should be one of CNEARTREE_NORM_L1 for an L-1 norm
                               CNEARTREE_NORM_L2 for an L-2 norm
                               CNEARTREE_NORM_LINF for an L-infinity norm
     
     the function returns CNEARTREE_BAD_ARGUMENT for an invalid argument
     CNEARTREE_SUCCESS (0) otherwise
     =======================================================================
     */
    
     int CNearTreeSetNorm(CNearTreeHandle treehandle, int treenorm);
    

    /*
     =======================================================================
     int CNearTreeNodeCreate ( CNearTreeHandle treehandle,  
                               CNearTreeNodeHandle FAR * treenodehandle) 
     
     Create a CNearTreeNode
     
     returns a pointer to the newly allocated block of memory as a 
     CNearTreeNodeHandle in *treenodehandle
     
     flags are inherited from the treehandle  
     
     creates an empty tree with no right or left node and with the dMax-below
     set to negative values so that any match found will be stored since it will
     greater than the negative value
     
     =======================================================================
     */
    
    int CNearTreeNodeCreate ( CNearTreeHandle treehandle,  
                              CNearTreeNodeHandle FAR * treenodehandle);
    
    
    /*
     =======================================================================
     int CNearTreeCreate ( CNearTreeHandle FAR * treehandle, 
                          size_t treedim, 
                          int treetype )
     
     Create a CNearTree
     
     treedim -- the dimension of the vectors
     treetype -- and integer flag for type of the vectors
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer     
     
     returns a pointer to the newly allocated block of memory as a 
     CNearTreeHandle in *treehandle, assuming coordinates of dimension treedim
     
     creates an empty tree with no right or left node and with the dMax-below
     set to negative values so that any match found will be stored since it will
     greater than the negative value
     
     =======================================================================
     */
    
    int CNearTreeCreate(CNearTreeHandle FAR * treehandle, 
                        size_t treedim, 
                        int treetype);
    
    /*
     =======================================================================
     int CNearTreeFree ( CNearTreeHandle FAR * treehandle )
     
     Free a CNearTree
     
     recursively frees the NearTree with the handle *treehandle
     and nulls the treehandle.
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeFree(CNearTreeHandle FAR * treehandle);
    
    /*
     =======================================================================
     int CNearTreeGetSize (CNearTreeHandle treehandle, size_t FAR * size)
     
     Return the number of objects in the tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetSize (CNearTreeHandle treehandle, size_t FAR * size);
    
    /*
     =======================================================================
     int CNearTreeGetDelayedSize (CNearTreeHandle treehandle, size_t FAR * size)
     
     Return the number of objects in the delay queue tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetDelayedSize (CNearTreeHandle treehandle, size_t FAR * size);
    
    /*
     =======================================================================
     int CNearTreeGetTotalSize (CNearTreeHandle treehandle, size_t FAR * size)
     
     Return the number of objects in both in the tree and in the delay queue tree in size
     Identical to CNearTreeGetSize, retained to support older code
     
     =======================================================================
     */
    
#define CNearTreeGetTotalSize(treehandle,size) CNearTreeGetSize(treehandle,size)
    /*
     =======================================================================
     int CNearTreeGetDepth (CNearTreeHandle treehandle, size_t FAR * depth)
     
     Return the depth of the tree in depth
     
     =======================================================================
     */
    
    int CNearTreeGetDepth (CNearTreeHandle treehandle, size_t FAR * depth);
    
    /*
     =======================================================================
     int CNearTreeInsert ( CNearTreeHandle treehandle, 
                          const void FAR * coord, 
                          const void * obj )
     
     Function to insert some "point" as an object into a CNearTree for
     later searching
     
     coord is a coordinate vector for an object, obj, to be inserted into a
     Neartree.  A static copy of the coordinates and a pointer to the
     object are inserted
     
     Three possibilities exist: put the datum into the left
     position (first test),into the right position, or else
     into a node descending from the nearer of those positions
     when they are both already used.
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
 
     int CNearTreeInsert( CNearTreeHandle treehandle,
                        const void FAR * coord, 
                        const void * obj );
    
 
    /*
     =======================================================================
     int CNearTreeNodeInsert( CNearTreeHandle treehandle,
                              CNearTreeNodeHandle treenodehandle,
                              const size_t index,
                              size_t FAR * depth );
     
     Function to insert some "point" as an object into a CNearTree for
     later searching, starting at a given treenode
     
     coord is a coordinate vector for an object, obj, to be inserted into a
     Neartree.  A static copy of the coordinates and a pointer to the
     object are inserted
     
     Three possibilities exist: put the datum into the left
     position (first test),into the right position, or else
     into a node descending from the nearer of those positions
     when they are both already used.
 
     depth is used to keep track of the depth at which the insertion is done

     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
     int CNearTreeNodeInsert( CNearTreeHandle treehandle,
                             CNearTreeNodeHandle treenodehandle,
                             const size_t index, 
                             size_t FAR * depth );
    
    
    /*
     =======================================================================
     int CNearTreeDelayedInsert ( CNearTreeHandle treehandle, 
     const void FAR * coord, 
     const void * obj )
     
     Function to queue some "point" as an object for future insertion
     into a CNearTree for later searching
     
     coord is a coordinate vector for an object, obj, to be inserted into a
     Neartree.  A static copy of the coordinates and a pointer to the
     object are queued for insertion.  The exact order of insertion
     is not predetermined.  It will be partially randomized to help to
     balance the tree.
     
     The insertions will be completed by a call to 
     CNearTreeCompleteDelayedInsert(CNearTreeHandle treehandle) 
     or by execution of any search.
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
     int CNearTreeDelayedInsert( CNearTreeHandle treehandle,
                        const void FAR * coord, 
                        const void * obj );
    
    /*
     =======================================================================
     int CNearTreeCompleteDelayedInsert ( CNearTreeHandle treehandle )
     
     Function to dequeue the "points" queued as an objects for future insertion
     into a CNearTree for later searching
      
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
     int CNearTreeCompleteDelayedInsert( CNearTreeHandle treehandle );
    
    /*
     =======================================================================
     int CNearTreeZeroIfEmpty (CNearTreeHandle treehandle)
     
     Test for an empty CNearTree, returning 0 in that case
     
     =======================================================================
     */
    
    int CNearTreeZeroIfEmpty (CNearTreeHandle treehandle);
    
    /*
     =======================================================================
     int CNearTreeNearestNeighbor ( CNearTreeHandle treehandle, 
                                    const double dRadius,  
                                    void FAR *  FAR * coordClosest,
                                    void FAR * FAR * objClosest,   
                                    const void FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord. This function
     is only here so that the function CNearTreeNearest can be called without having dRadius const
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
    
     coordClosest is the coordinate vector into which the coordinates of the nearest point
     will be stored
    
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
  
     coord  is the probe point
  
     the return value is true only if a point was found
     
     =======================================================================
     */
    
    int CNearTreeNearestNeighbor (CNearTreeHandle treehandle, 
                                  const double dRadius,  
                                  void FAR *  FAR * coordClosest,
                                  void FAR * FAR * objClosest, 
                                  const void FAR * coord );
    /*
     =======================================================================
     int CNearTreeFarthestNeighbor ( CNearTreeHandle treehandle, 
                                    void FAR *  FAR * coordFarthest,
                                    void FAR * FAR * objFarthest,   
                                     const void FAR * coord )
     
     Function to search a Neartree for the object farthest some probe point, coord.
     
     coordClosest is the coordinate vector into which the coordinates of the nearest point
     will be stored

     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored

     coord  is the probe point
  
     the return value is true only if a point was found
     
     =======================================================================
     */
    
    int CNearTreeFarthestNeighbor (CNearTreeHandle treehandle, 
                                   void FAR * FAR * coordFarthest,
                                   void FAR * FAR * objFarthest,   
                                   const void FAR * coord );
    /*
     =======================================================================
     int CNearTreeFindInSphere ( CNearTreeHandle treehandle, 
                                const double FAR * dRadius,
                                CVectorHandle coordClosest,
                                CVectorHandle objClosest,
                                const void FAR * coord,
                                int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius. This is only here so that objClosest can be cleared before starting the work.
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored

     coordClosest is a vector of pointers to coordinate tuples and is the 
     returned set of nearest points to the probe point that can be found in the Neartree

     objClosest is a vector of objects and is the returned set of nearest points
     to the probe point that can be found in the Neartree

     coord  is the probe point

     resetcount should be non-zero to clear objClosest on entry
 
     return value is 0 if points were found
     
     =======================================================================
     */
    
    int CNearTreeFindInSphere ( CNearTreeHandle treehandle,
                               const double dRadius,
                               CVectorHandle coordClosest,
                               CVectorHandle objClosest,
                               const void * coord,
                               int resetcount);
    /*
     =======================================================================
     int CNearTreeNearest ( CNearTreeHandle treehandle, 
                           double FAR *dRadius,  
                           void FAR * FAR * coordClosest,
                           void FAR * FAR * objClosest,
                           const void FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord.
     This function is called by CNearTreeNearestNeighbor.
     
     *dRadius is the smallest currently known distance of an object from the probe point
     and is modified as the search progresses.

     coordClosest is the returned closest point
     to the probe point that can be found in the Neartree

     objClosest is a pointer to a pointer to hold the corresponding object or is NULL

     coord  is the probe point

     the return value is 0 only if a point was found within dRadius
     
     =======================================================================
     */
    int CNearTreeNearest ( CNearTreeHandle treehandle, 
                          double FAR * dRadius,  
                          void FAR * FAR * coordClosest,
                          void FAR * FAR * objClosest,
                          const void FAR * coord );
    /*
     =======================================================================
     int CNearTreeFindFarthest ( CNearTreeHandle treehandle,
                                double FAR * dRadius,  
                                void FAR * FAR * coordFarthest,
                                void FAR * FAR * objFarthest,
                                const void FAR * coord )
     
     Function to search a Neartree for the object farthest from some probe point, coord.
     This function is called by CNearTreeFarthestNeighbor.
     
     *dRadius is the largest currently known distance of an object from the probe point
     and is modified as the search progresses.

     coordFarthest is the returned farthest point
     from the probe point that can be found in the Neartree

     objFarthest is a pointer to a pointer to hold the corresponding object or is NULL
     coord  is the probe point

     the return value is 0 only if a point was found.
     
     =======================================================================
     */
    
    int CNearTreeFindFarthest ( CNearTreeHandle treehandle,
                               double FAR * dRadius,  
                               void FAR * FAR * coordFarthest,
                               void FAR * FAR * objFarthest,
                               const void FAR * coord );
    
#ifdef __cplusplus
    
}

#endif


#endif