/*
 *  CNearTree.c
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


#ifdef __cplusplus

extern "C" {
    
#endif
    
#ifndef USE_LOCAL_HEADERS
#include <CNearTree.h>
#else
#include "CNearTree.h"
#endif
#include <math.h>
    
#ifdef CNEARTREE_SAFE_TRIANG
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0) \
                      || ((b)-((a)-(c)) >= 0) \
                      || ((c)-((a)-(b)) >= 0))    
#else
#define TRIANG(a,b,c) (  (((b)+(c))-(a) >= 0))
#endif
    
#define max2(x,y) ((x)>=(y)?(x):(y))    
    
    /*
     =======================================================================
     double CNearTreeDistsq(void FAR * coord1, 
     void FAR * coord2,
     size_t treedim, int treetype)
     
     function to return the square of the Euclidean distance between two 
     coordinate vectors.
     
     treedim -- the dimension of the vectors
     treetype -- and integer flag for type of the vectors
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer
     
     =======================================================================
     */
    
    double CNearTreeDistsq(void FAR * coord1, 
                           void FAR * coord2,
                           size_t treedim, int treetype) {
        size_t index;
        double distsq;
        
        if (treetype == CNEARTREE_TYPE_DOUBLE) {
            
            double FAR * dcoord1;
            double FAR * dcoord2;
            
            dcoord1 = (double FAR *)coord1;
            dcoord2 = (double FAR *)coord2;
            
            distsq = (dcoord1[0]-dcoord2[0])*(dcoord1[0]-dcoord2[0]);
            
            for (index=1; index < treedim; index++) {
                distsq += (dcoord1[index]-dcoord2[index])
                *(dcoord1[index]-dcoord2[index]);
            }
            
            return distsq;
            
        } else if (treetype == CNEARTREE_TYPE_INTEGER) {
            
            int FAR * icoord1;
            int FAR * icoord2;
            
            icoord1 = (int FAR *)coord1;
            icoord2 = (int FAR *)coord2;
            
            distsq = ((double)(icoord1[0]-icoord2[0]))*((double)(icoord1[0]-icoord2[0]));
            
            for (index=1; index < treedim; index++) {
                distsq += ((double)(icoord1[index]-icoord2[index]))
                *((double)(icoord1[index]-icoord2[index]));
            }
            
            return distsq;
            
            
        } else return -1.0;
        
    }
    
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
    
    int CNearTreeSetNorm(CNearTreeHandle treehandle, int treenorm) {
        
        if (!treehandle ||
            (treenorm != CNEARTREE_NORM_L1
             && treenorm != CNEARTREE_NORM_L2
             && treenorm != CNEARTREE_NORM_LINF)
            || (treehandle->m_iflags & CNEARTREE_NORM)!=CNEARTREE_NORM_UNKNOWN ) return CNEARTREE_BAD_ARGUMENT;
        treehandle->m_iflags &= ~CNEARTREE_NORM_UNKNOWN;
        treehandle->m_iflags |= treenorm;
        return CNEARTREE_SUCCESS;
    }
    /*
     =======================================================================
     double CNearTreeDist(CNearTreeHandle treehandle, 
                          void FAR * coord1, 
                          void FAR * coord2))
     
     function to return the distance (L1, L2 or L-infinity) between two 
     coordinate vectors according to the parameters of the given tree  
     
     =======================================================================
     */
    
    double CNearTreeDist(CNearTreeHandle treehandle, 
                         void FAR * coord1,
                         void FAR * coord2)  {
        size_t index;
        double dist, distsq;
        
        size_t treedim;
        int treetype;
        int treenorm;
        
        double FAR * dcoord1 = NULL;
        double FAR * dcoord2 = NULL;
        int FAR * icoord1 = NULL;
        int FAR * icoord2 = NULL;
        
        
        treedim = treehandle->m_szdimension;
        treetype = treehandle->m_iflags&CNEARTREE_TYPE;
        treenorm = treehandle->m_iflags&CNEARTREE_NORM;
        
        if (treenorm == CNEARTREE_NORM_UNKNOWN || treenorm == 0) {
            treehandle->m_iflags &= ~CNEARTREE_NORM;
            treehandle->m_iflags |= CNEARTREE_NORM_L2;
            treenorm = CNEARTREE_NORM_L2;
        }
        
        if (!treehandle || !coord1 || !coord2) return -1.;
        
        
        if (treetype == CNEARTREE_TYPE_DOUBLE) {
            
            dcoord1 = (double FAR *)coord1;
            dcoord2 = (double FAR *)coord2;
        } else if (treetype == CNEARTREE_TYPE_INTEGER) {
            icoord1 = (int FAR *)coord1;
            icoord2 = (int FAR *)coord2;
        } else return -1.0;
        
        if (treedim == 1) {
            if (treetype == CNEARTREE_TYPE_DOUBLE) {
                return fabs(dcoord1[0]-dcoord2[0]);
            } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                return fabs((double)(icoord1[0]-icoord2[0]));
            } else {
                return -1.0;
            }
        }
        
        switch (treenorm) {
            case CNEARTREE_NORM_L1:
                if (treetype == CNEARTREE_TYPE_DOUBLE) {
                    dist= fabs(dcoord1[0]-dcoord2[0]);
                    for (index=1; index < treedim; index++) {
                        dist += fabs(dcoord1[index]-dcoord2[index]);
                    }
                } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                    dist = fabs((double)(icoord1[0]-icoord2[0]));
                    for (index=1; index < treedim; index++) {
                        dist += fabs((double)(dcoord1[index]-dcoord2[index]));
                    }
                } else return -1.0;
                return dist;
            case CNEARTREE_NORM_LINF:
                if (treetype == CNEARTREE_TYPE_DOUBLE) {
                    dist= fabs(dcoord1[0]-dcoord2[0]);
                    for (index=1; index < treedim; index++) {
                        dist = max2(dist,fabs(dcoord1[index]-dcoord2[index]));
                    }
                } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                    dist = fabs((double)(icoord1[0]-icoord2[0]));
                    for (index=1; index < treedim; index++) {
                        dist = max2(dist,fabs((double)(dcoord1[index]-dcoord2[index])));
                    }
                } else return -1.0;
                return dist;
                
            case CNEARTREE_NORM_L2:
            default:
                if (treetype == CNEARTREE_TYPE_DOUBLE) {
                    distsq = (dcoord1[0]-dcoord2[0])*(dcoord1[0]-dcoord2[0]);
                    
                    for (index=1; index < treedim; index++) {
                        distsq += (dcoord1[index]-dcoord2[index])
                        *(dcoord1[index]-dcoord2[index]);
                    }
                } else if (treetype == CNEARTREE_TYPE_INTEGER) {
                    distsq = ((double)(icoord1[0]-icoord2[0]))*((double)(icoord1[0]-icoord2[0]));
                    
                    for (index=1; index < treedim; index++) {
                        distsq += ((double)(icoord1[index]-icoord2[index]))
                        *((double)(icoord1[index]-icoord2[index]));
                    }
                } else return -1;
                return sqrt(distsq);
        }
        
    }
    
    
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
                             CNearTreeNodeHandle FAR * treenodehandle) {
        
        size_t index;
        
        size_t ntsize, nvsize;
        
        int treenorm, treetype;
        
        if (!treehandle || !treenodehandle) return CNEARTREE_BAD_ARGUMENT;
        
        treenorm = (treehandle->m_iflags) & CNEARTREE_NORM;
        
        if (!treenorm) treenorm = CNEARTREE_NORM_UNKNOWN;
        
        treetype = (treehandle->m_iflags) & CNEARTREE_TYPE;
        
        if ( (treetype != CNEARTREE_TYPE_DOUBLE 
              && treetype != CNEARTREE_TYPE_INTEGER)) return CNEARTREE_BAD_ARGUMENT;
        
        
        if (treetype == CNEARTREE_TYPE_INTEGER) {
            ntsize = (sizeof(CNearTreeNode)+(sizeof(int)-1))/sizeof(int);
            ntsize *= sizeof(int);
            nvsize = sizeof(int)*treehandle->m_szdimension;
            *treenodehandle = (CNearTreeNodeHandle)MALLOC(ntsize+2*nvsize);
        } else {
            ntsize = (sizeof(CNearTreeNode)+(sizeof(double)-1))/sizeof(double);
            ntsize *= sizeof(double);
            nvsize = sizeof(double)*treehandle->m_szdimension;
            *treenodehandle = (CNearTreeNodeHandle)MALLOC(ntsize+2*nvsize);
        }
        
        if (!(*treenodehandle)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        (*treenodehandle)->m_coordLeft = ((char *)(*treenodehandle))+ntsize;
        (*treenodehandle)->m_coordRight = ((char *)(*treenodehandle))+ntsize+nvsize;
        if (treetype == CNEARTREE_TYPE_DOUBLE) {
            for (index=0; index < treehandle->m_szdimension; index++) {
                ((double *)(*treenodehandle)->m_coordLeft)[index]    = 0;
                ((double *)(*treenodehandle)->m_coordRight)[index]   = 0;
            }
        } else {
            for (index=0; index < treehandle->m_szdimension; index++) {
                ((int *)(*treenodehandle)->m_coordLeft)[index]    = 0;
                ((int *)(*treenodehandle)->m_coordRight)[index]   = 0;
            }
        }
        (*treenodehandle)->m_pLeftBranch  = 0;
        (*treenodehandle)->m_pRightBranch = 0;
        (*treenodehandle)->m_dMaxLeft     = -1.;        /* negative for an empty branch */
        (*treenodehandle)->m_dMaxRight    = -1.;        /* negative for an empty branch */
        (*treenodehandle)->m_ptobjLeft    = NULL;
        (*treenodehandle)->m_ptobjRight   = NULL;
        (*treenodehandle)->m_iflags       = treetype;  /* no data, no children */
        (*treenodehandle)->m_iflags      |= treenorm;
        return CNEARTREE_SUCCESS;
    }
    
    
    /*
     =======================================================================
     int CNearTreeCreate ( CNearTreeHandle FAR * treehandle, 
     size_t treedim, int treetype)
     
     Create a CNearTree
     
     returns a pointer to the newly allocated block of memory as a 
     CNearTreeHandle in *treehandle
     
     
     treedim -- the dimension of the vectors
     treetype -- double or integer flag for type of the vectors ored with norm
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer
       ored with
     CNEARTREE_NORM_L1        for the sum of the absolute values
     CNEARTREE_NORM_L2        for the square root of the sum of the squares
     CNEARTREE_NORM_LINF      for the max
     
     
     creates an empty tree with no right or left node and with the dMax-below
     set to negative values so that any match found will be stored since it will
     greater than the negative value
     
     =======================================================================
     */
    
    int CNearTreeCreate ( CNearTreeHandle FAR * treehandle, 
                         size_t treedim, int treetype) {
        
        int treenorm;
        
        int treeflip;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        treenorm = treetype & CNEARTREE_NORM;
        
        treeflip = treetype & CNEARTREE_FLIP;
        
        if (!treenorm) treenorm = CNEARTREE_NORM_UNKNOWN;
        
        treetype &= CNEARTREE_TYPE;
        
        if (treedim < 0 
            || (treetype != CNEARTREE_TYPE_DOUBLE 
                && treetype != CNEARTREE_TYPE_INTEGER)) return CNEARTREE_BAD_ARGUMENT;
        
        *treehandle = (CNearTreeHandle)MALLOC(sizeof(CNearTree));
        if (!(*treehandle)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        (*treehandle)->m_iflags        = treetype; /* no data, no children         */
        (*treehandle)->m_iflags       |= treenorm; /* record the chosen norm       */
        (*treehandle)->m_iflags       |= treeflip; /* record whether to flip       */
        (*treehandle)->m_szdimension   = treedim;  /* number of ints or doubles    */
        (*treehandle)->m_szsize        = 0;        /* number of nodes in the tree  */
        (*treehandle)->m_szdepth       = 0;        /* depth of in the tree         */
        if( CNearTreeNodeCreate(*treehandle,&((*treehandle)->m_ptTree))) {
            FREE(*treehandle);
            return CNEARTREE_MALLOC_FAILED;
        }
        (*treehandle)->m_ptDelayCoords = NULL;
        (*treehandle)->m_ptDelayObjs = NULL;
        return CNEARTREE_SUCCESS;
        
    }
    
    /*
     =======================================================================
     int CNearTreeNodeFree ( CNearTreeNodeHandle FAR * treenodehandle )
     
     Free a CNearTreeNode
     
     recursively frees the NearTreeNode tree with the handle *trenodeehandle
     and nulls the treenodehandle.
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeNodeFree ( CNearTreeNodeHandle FAR * treenodehandle ) {
        
        int errorcode;
        
        if (!treenodehandle) return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = CNEARTREE_SUCCESS;
        
        if ((*treenodehandle)->m_pLeftBranch) {
            errorcode |= CNearTreeNodeFree(&((*treenodehandle)->m_pLeftBranch));
        }
        
        if ((*treenodehandle)->m_pRightBranch) {
            errorcode |= CNearTreeNodeFree(&((*treenodehandle)->m_pRightBranch));
        }
        
        FREE(*treenodehandle);
        
        *treenodehandle = NULL;
        
        return errorcode;
    }
    /*
     =======================================================================
     int CNearTreeFree ( CNearTreeHandle FAR * treehandle )
     
     Free a CNearTree
     
     Frees the NearTree with the handle *treehandle
     and nulls the treehandle.
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeFree ( CNearTreeHandle FAR * treehandle ) {
        
        int errorcode, errorcodev;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = CNEARTREE_SUCCESS;
        
        if ((*treehandle)->m_ptTree) {
            errorcode |= CNearTreeNodeFree(&((*treehandle)->m_ptTree));
        }
        
        if ((*treehandle)->m_ptDelayCoords) {
            errorcodev = CVectorFree(&((*treehandle)->m_ptDelayCoords));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
        }
        
        if ((*treehandle)->m_ptDelayObjs) {
            errorcodev = CVectorFree(&((*treehandle)->m_ptDelayObjs));
            if (errorcodev) errorcode |= CNEARTREE_FREE_FAILED ;
        }
        
        FREE(*treehandle);
        
        *treehandle = NULL;
        
        return errorcode;
    }
    
    /*
     =======================================================================
     int CNearTreeZeroIfEmpty (CNearTreeHandle treehandle)
     
     Test for an empty CNearTree, returning 0 in that case
     
     =======================================================================
     */
    
    int CNearTreeZeroIfEmpty (CNearTreeHandle treehandle)
    {
        return (treehandle==NULL||(((treehandle->m_ptTree->m_iflags)&CNEARTREE_DATA_OR_CHILDREN)  == 0))?0:1;
    }
    
    /*
     =======================================================================
     int CNearTreeGetSize (CNearTreeHandle treehandle, size_t FAR * size)
     
     Return the number of objects in the tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetSize (CNearTreeHandle treehandle, size_t FAR * size)
    {
        if (!treehandle || !size ) return CNEARTREE_BAD_ARGUMENT;
        
        *size = treehandle->m_szsize;
        
        return CNEARTREE_SUCCESS;

    }

    /*
     =======================================================================
     int CNearTreeGetDelayedSize (CNearTreeHandle treehandle, size_t FAR * size)
     
     Return the number of objects in the delay queue tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetDelayedSize (CNearTreeHandle treehandle, size_t FAR * size)
    {
        if (!treehandle || !size ) return CNEARTREE_BAD_ARGUMENT;
        
        *size = 0;
        
        if (treehandle->m_ptDelayCoords) *size = CVectorSize(treehandle->m_ptDelayCoords);
        
        return CNEARTREE_SUCCESS;

    }
    
       /*
     =======================================================================
     int CNearTreeGetTotalSize (CNearTreeHandle treehandle, size_t FAR * size)
     
     Return the number of objects in both in the tree and in the delay queue tree in size
     
     =======================================================================
     */
    
    int CNearTreeGetTotalSize (CNearTreeHandle treehandle, size_t FAR * size)
    {
        if (!treehandle || !size ) return CNEARTREE_BAD_ARGUMENT;
        
        *size = treehandle->m_szsize;;
        
        if (treehandle->m_ptDelayCoords) *size += CVectorSize(treehandle->m_ptDelayCoords);
        
        return CNEARTREE_SUCCESS;

    }


    /*
     =======================================================================
     int CNearTreeGetDepth (CNearTreeHandle treehandle, size_t FAR * depth)
     
     Return the depth of the tree in depth
     
     =======================================================================
     */
    
    int CNearTreeGetDepth (CNearTreeHandle treehandle, size_t FAR * depth)
    {
        if (!treehandle || !depth ) return CNEARTREE_BAD_ARGUMENT;
        
        *depth = treehandle->m_szdepth;
        
        return CNEARTREE_SUCCESS;
        
    }
    
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
                        const void * obj ) {
        
        /* do a bit of precomputing if it is possible so that we can
         reduce the number of calls to sqrt */
        
        size_t depth;
        
        int errorcode;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( !(treehandle->m_ptTree) ) return CNEARTREE_BAD_ARGUMENT;
        
        depth = 1;
        
        if ((errorcode = 
             CNearTreeNodeInsert( treehandle, treehandle->m_ptTree, coord, obj, &depth)) == 0) {
            
            if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
            
            (treehandle->m_szsize)++;
            
        }
        
        return errorcode;
        
    }
    
    
    /*
     =======================================================================
     int CNearTreeNodeInsert ( CNearTreeHandle treehandle,
     CNearTreeNodeHandle treenodehandle, 
     const void FAR * coord, 
     const void * obj,
     size_t FAR * depth)
     
     Function to insert some "point" as an object into a a node if a CNearTree for
     later searching
     
     treenodehandle is the handle of the node at which to start the insertion
     
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
                            const void FAR * coord, 
                            const void * obj,
                            size_t FAR * depth) {
        
        /* do a bit of precomputing if it is possible so that we can
         reduce the number of calls to sqrt */
        
        double dTempRight =  0.;
        double dTempLeft  =  0.;
        int errorcode = 0;
        size_t index;
        
        if ( !treehandle || !treenodehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( (treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA ) {
            dTempRight = CNearTreeDist(treehandle, (void FAR *)coord,treenodehandle->m_coordRight);
        }
        
        if ( (treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_DATA ) {
            dTempLeft = CNearTreeDist(treehandle, (void FAR *)coord,treenodehandle->m_coordLeft);
        }
        
        if ( !((treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) ) {
            if (treehandle->m_iflags&CNEARTREE_TYPE_DOUBLE) {
                for (index=0; index<treehandle->m_szdimension; index++)  {
                    ((double *)(treenodehandle->m_coordLeft))[index]=((double *)coord)[index];    	
                }
            } else {
                for (index=0; index<treehandle->m_szdimension; index++)  {
                    ((int *)(treenodehandle->m_coordLeft))[index]=((int *)coord)[index];    	
                }
            }
            treenodehandle->m_dMaxLeft = -1.;
            treenodehandle->m_ptobjLeft = (void *)obj;
            treenodehandle->m_iflags |= CNEARTREE_FLAG_LEFT_DATA;
            return CNEARTREE_SUCCESS;
        } else if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA)  ){ 
            size_t index;
            if (treehandle->m_iflags&CNEARTREE_TYPE_DOUBLE) {
                for (index=0; index<treehandle->m_szdimension; index++)  {
                    ((double *)treenodehandle->m_coordRight)[index]=((double *)coord)[index];    	
                }
            } else {
                for (index=0; index<treehandle->m_szdimension; index++)  {
                    ((int *)treenodehandle->m_coordRight)[index]=((int *)coord)[index];    	
                }
            }
            treenodehandle->m_dMaxRight = -1.;
            treenodehandle->m_ptobjRight = (void *)obj;
            treenodehandle->m_iflags |= CNEARTREE_FLAG_RIGHT_DATA;
            return CNEARTREE_SUCCESS;
        } else if ( dTempLeft > dTempRight ) {
            if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_CHILD) ) {
                if ( (errorcode = CNearTreeNodeCreate(treehandle, &(treenodehandle->m_pRightBranch)))) return errorcode;
                treenodehandle->m_iflags |= CNEARTREE_FLAG_RIGHT_CHILD;
                treenodehandle->m_dMaxRight = dTempRight;
                if (((treehandle->m_iflags)&CNEARTREE_FLIP)&&!((treenodehandle->m_iflags)&CNEARTREE_FLAG_RIGHT_CHILD)
                    && (CNearTreeDist(treehandle,treenodehandle->m_coordLeft,treenodehandle->m_coordRight))< dTempRight) {
                        if ( treenodehandle->m_dMaxRight < dTempRight ) 
                            treenodehandle->m_dMaxRight = dTempRight;
                        depth++;
                        errorcode = CNearTreeNodeInsert(treehandle, treenodehandle->m_pRightBranch, 
                                                        treenodehandle->m_coordRight,  
                                                        treenodehandle->m_ptobjRight, depth);
                        treenodehandle->m_ptobjRight = (void *)obj;    
                        if (treehandle->m_iflags&CNEARTREE_TYPE_DOUBLE) {
                            for (index=0; index<treehandle->m_szdimension; index++)  {
                                ((double *)treenodehandle->m_coordRight)[index]=((double *)coord)[index];    	
                            }
                        } else {
                            for (index=0; index<treehandle->m_szdimension; index++)  {
                                ((int *)treenodehandle->m_coordRight)[index]=((int *)coord)[index];    	
                            }
                        }
                        return errorcode;
                    }
            }
            /* note that the next line assumes that m_dMaxRight is negative for a new node */
            if ( treenodehandle->m_dMaxRight < dTempRight ) 
                treenodehandle->m_dMaxRight = dTempRight;
            (*depth)++;
            return CNearTreeNodeInsert(treehandle, treenodehandle->m_pRightBranch, coord, obj, depth);
        } else { /* ((double)(t - *m_tLeft) <= (double)(t - *m_tRight) ) */
            if (  !((treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_CHILD) ) {
                if ( (errorcode = CNearTreeNodeCreate(treehandle, &(treenodehandle->m_pLeftBranch)))) return errorcode;
                treenodehandle->m_iflags |= CNEARTREE_FLAG_LEFT_CHILD;
                treenodehandle->m_dMaxLeft = dTempLeft;
                if (((treehandle->m_iflags)&CNEARTREE_FLIP)&&!((treenodehandle->m_iflags)&CNEARTREE_FLAG_LEFT_CHILD)
                    && (CNearTreeDist(treehandle,treenodehandle->m_coordLeft,treenodehandle->m_coordRight))< dTempLeft) {
                        if ( treenodehandle->m_dMaxLeft < dTempLeft ) 
                            treenodehandle->m_dMaxLeft = dTempLeft;
                        depth++;
                        errorcode = CNearTreeNodeInsert(treehandle, treenodehandle->m_pLeftBranch, 
                                                        treenodehandle->m_coordLeft,  
                                                        treenodehandle->m_ptobjLeft, depth);
                        treenodehandle->m_ptobjLeft = (void *)obj;    
                        if (treehandle->m_iflags&CNEARTREE_TYPE_DOUBLE) {
                            for (index=0; index<treehandle->m_szdimension; index++)  {
                                ((double *)treenodehandle->m_coordLeft)[index]=((double *)coord)[index];    	
                            }
                        } else {
                            for (index=0; index<treehandle->m_szdimension; index++)  {
                                ((int *)treenodehandle->m_coordLeft)[index]=((int *)coord)[index];    	
                            }
                        }
                        return errorcode;
                    }
            }
            /* note that the next line assumes that m_dMaxLeft is negative for a new node */
            if ( treenodehandle->m_dMaxLeft < dTempLeft ) 
                treenodehandle->m_dMaxLeft  = dTempLeft;
            (*depth)++;
            return CNearTreeNodeInsert(treehandle,treenodehandle->m_pLeftBranch, coord, obj, depth);
        }
        
    }
    
    
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
                               const void * obj ) {
        
        int treetype;
        
        size_t nvsize;
        
        size_t treedim;
        
        if (!treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        treetype = (treehandle->m_iflags) & CNEARTREE_TYPE;
        
        treedim = treehandle->m_szdimension;
        
        if (treetype == CNEARTREE_TYPE_INTEGER) {
            nvsize = sizeof(int)*treedim;
        } else {
            nvsize = sizeof(double)*treedim;
        }
        
        if (!(treehandle->m_ptDelayCoords)) {
            if (CVectorCreate(&(treehandle->m_ptDelayCoords),nvsize,1)) {
                return CNEARTREE_MALLOC_FAILED;
            }
        }
        
        if (!(treehandle->m_ptDelayObjs)) {
            if (CVectorCreate(&(treehandle->m_ptDelayObjs),sizeof(void FAR *),1)) {
                return CNEARTREE_MALLOC_FAILED;
            }
        }
        
        if (CVectorAddElement( treehandle->m_ptDelayObjs, &obj)) {
            return CNEARTREE_BAD_ARGUMENT;
        }
        
        if (CVectorAddElement( treehandle->m_ptDelayCoords, coord) ) {
            return CNEARTREE_BAD_ARGUMENT;
        }
        
        return CNEARTREE_SUCCESS;
        
    }
    
    /*
     =======================================================================
     int CNearTreeCompleteDelayedInsert ( CNearTreeHandle treehandle )
     
     Function to dequeue the "points" queued as an objects for future insertion
     into a CNearTree for later searching
     
     return 0 for success, nonzero for an error
     
     =======================================================================
     */
    
    int CNearTreeCompleteDelayedInsert ( CNearTreeHandle treehandle ) {
        
        size_t nqueued;
        
        size_t nrandom;
        
        size_t ielement, kelement=0;
        
        void FAR * obj;
        
        void FAR * dummyobj;
        
        void FAR * coord;
        
        int errorcode;
        
        size_t depth;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        if (treehandle->m_ptDelayCoords == NULL && treehandle->m_ptDelayObjs == NULL ) 
            return CNEARTREE_SUCCESS;
        
        if (treehandle->m_ptDelayCoords == NULL || treehandle->m_ptDelayObjs == NULL ||
            CVectorSize(treehandle->m_ptDelayCoords) != CVectorSize(treehandle->m_ptDelayObjs))
            return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = 0;
        
        nqueued = CVectorSize(treehandle->m_ptDelayCoords);
        
        nrandom = (size_t)sqrt((double)nqueued);
        
        dummyobj = (void FAR *)&nrandom;
        
        for (ielement = 0; ielement < nrandom; ielement++) {
            kelement = random()%nqueued;
            do {
                if (kelement >= nqueued) kelement = 0;
                if (CVectorGetElement(treehandle->m_ptDelayObjs,&obj,kelement)) return  CNEARTREE_BAD_ARGUMENT;
                kelement ++;
            } while (obj == dummyobj);
            kelement--;
            if (CVectorSetElement(treehandle->m_ptDelayObjs,&dummyobj,kelement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
            if (CVectorGetElementptr(treehandle->m_ptDelayCoords,&coord,kelement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
            depth = 1;
            errorcode |= CNearTreeNodeInsert(treehandle,treehandle->m_ptTree,coord,obj,&depth);
            if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
            (treehandle->m_szsize)++;
        }
        
        for (ielement = 0; ielement < nqueued; ielement++) {
            if (CVectorGetElement(treehandle->m_ptDelayObjs,&obj,ielement)) return CNEARTREE_CVECTOR_FAILED;
            if (obj == dummyobj) continue;
            if (CVectorGetElementptr(treehandle->m_ptDelayCoords,&coord,ielement)) errorcode |= CNEARTREE_CVECTOR_FAILED;
            depth = 1;
            errorcode |= CNearTreeNodeInsert(treehandle,treehandle->m_ptTree,coord,obj,&depth);
            if (depth > treehandle->m_szdepth) treehandle->m_szdepth = depth;
            (treehandle->m_szsize)++;
        }
        
        if (CVectorFree(&(treehandle->m_ptDelayObjs))) errorcode |= CNEARTREE_CVECTOR_FAILED;
        if (CVectorSetFlags(treehandle->m_ptDelayCoords,0)) errorcode |= CNEARTREE_CVECTOR_FAILED;
        if (CVectorFree(&(treehandle->m_ptDelayCoords))) errorcode |= CNEARTREE_CVECTOR_FAILED;
        
        return (errorcode != 0)?errorcode:CNEARTREE_SUCCESS;
    }
    
    /*
     =======================================================================
     int CNearTreeNearestNeighbor ( CNearTreeHandle treehandle, 
     const double dRadius,  
     void FAR * FAR * coordClosest,
     void FAR * FAR * objClosest,   
     const void FAR * coord )
     
     Function to search a Neartree for the object closest to some probe point, coord. This function
     is only here so that the function CNearTreeNearest can be called without having dRadius const
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     coord  is the probe point
     the return value is true only if a point was found
     
     =======================================================================
     */
    int CNearTreeNearestNeighbor (CNearTreeHandle treehandle, 
                                  const double dRadius,  
                                  void FAR * FAR * coordClosest,
                                  void FAR * FAR * objClosest, 
                                  const void FAR * coord ) {
        
        double dSearchRadius = dRadius;
        if (!treehandle || ! coord ) return CNEARTREE_BAD_ARGUMENT;
        if (!(treehandle->m_ptTree->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) 
            return CNEARTREE_NOT_FOUND;
        return ( CNearTreeNearest ( treehandle, &dSearchRadius, coordClosest, objClosest, coord ) );
    }
    
    
    /*
     =======================================================================
     int CNearTreeFarthestNeighbor ( CNearTreeHandle treehandle, 
     void FAR * FAR * coordFarthest,
     void FAR * FAR * objFarthest,   
     const void FAR * coord )
     
     Function to search a Neartree for the object farthest some probe point, coord.
     
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is the address into which a pointer to the object associated with coordClosest
     will be stored
     coord  is the probe point
     the return value is 0 only if a point was found
     
     =======================================================================
     */
    int CNearTreeFarthestNeighbor (CNearTreeHandle treehandle, 
                                   void FAR *  FAR * coordFarthest,
                                   void FAR * FAR * objFarthest,   
                                   const void FAR * coord ) {
        double dSearchRadius = DBL_MIN;
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        if (!(treehandle->m_ptTree->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) 
            return CNEARTREE_NOT_FOUND;
        return ( CNearTreeFindFarthest ( treehandle, &dSearchRadius, coordFarthest, objFarthest, coord ) );
    }
    
    /*
     =======================================================================
     int CNearTreeFindInSphere ( CNearTreeHandle treehandle, 
     const double FAR * dRadius,
     CVectorHandle coordclosest,
     CVectorHandle objClosest,
     const void FAR * coord,
     int resetcount)
     
     Function to search a Neartree for the set of objects closer to some probe point, coord,
     than dRadius. This is only here so that objClosest can be cleared before starting the work.
     
     dRadius is the maximum search radius - any point farther than dRadius from the probe
     point will be ignored
     coordClosest is a vector of pointers to coordinate tuples of nearest points
     objClosest is a vector of objects and is the returned set of nearest points
     to the probe point that can be found in the Neartree
     coord  is the probe point
     resetcount should be non-zero to clear coordclosest and objClosest on entry
     return value is 0 if points were found
     
     =======================================================================
     */
    int CNearTreeFindInSphere ( CNearTreeHandle treehandle,
                               const double dRadius,
                               CVectorHandle coordClosest,
                               CVectorHandle objClosest,
                               const void FAR * coord,
                               int resetcount) {
        double dDR, dDL;
        int nopoints;
        CVectorHandle sStack;
        CNearTreeNodeHandle pt;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        nopoints = 1;
        
        dDR = dDL = DBL_MAX;
        
        if (dRadius < 0.) return 1;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_ptDelayObjs ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordClosest) CVectorClear( coordClosest );
            if (objClosest) CVectorClear( objClosest );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = DBL_MAX;
                if ((pt->m_iflags)&CNEARTREE_FLAG_RIGHT_DATA) {
                    dDR = CNearTreeDist(treehandle, (void FAR *)coord, pt->m_coordRight);
                    if (dDR <= dRadius ) {
                        nopoints = 0;
                        if (coordClosest) CVectorAddElement(coordClosest,&(pt->m_coordRight));
                        if (objClosest && pt->m_ptobjRight) CVectorAddElement(objClosest,&(pt->m_ptobjRight));
                    }
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dDR,pt->m_dMaxRight,dRadius))){
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_iflags)&CNEARTREE_FLAG_LEFT_DATA) {
                    dDL = CNearTreeDist(treehandle, (void FAR *)coord, pt->m_coordLeft);
                    if (dDL <= dRadius ) {
                        nopoints = 0;
                        if (coordClosest) CVectorAddElement(coordClosest,&(pt->m_coordLeft));
                        if (objClosest &&  pt->m_ptobjLeft) CVectorAddElement(objClosest,&(pt->m_ptobjLeft));
                    }
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,dRadius))){
                    pt = pt->m_pLeftBranch;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        return  nopoints?CNEARTREE_NOT_FOUND:CNEARTREE_SUCCESS;
    }     
    
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
     coordClosest is a pointer to the coordinate vector of the nearest point
     objClosest is a pointer to a pointer to hold the corresponding object or is NULL
     coord  is the probe point
     the return value is 0 only if a point was found within dRadius
     
     =======================================================================
     */
    int CNearTreeNearest ( CNearTreeHandle treehandle, 
                          double FAR * dRadius,  
                          void FAR * FAR * coordClosest,
                          void FAR * FAR * objClosest,
                          const void FAR * coord ) {
        double   dDR, dDL;
        CVectorHandle sStack;
        CNearTreeNodeHandle pt;
        void FAR * pobjClosest;
        void FAR * pcoordClosest;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        pobjClosest = NULL;
        pcoordClosest = NULL;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_ptDelayObjs ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = CNearTreeDist(treehandle, (void FAR *)coord,  pt->m_coordRight);
                if (dDR <= *dRadius ) {
                    *dRadius = dDR;
                    pobjClosest = pt->m_ptobjRight;
                    pcoordClosest = pt->m_coordRight;
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(dDR,pt->m_dMaxRight,*dRadius))) {
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                dDL = CNearTreeDist(treehandle, (void FAR *)coord, pt->m_coordLeft);
                if (dDL <= *dRadius ) {
                    *dRadius = dDL;
                    pobjClosest = pt->m_ptobjLeft;
                    pcoordClosest = pt->m_coordLeft;
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(dDL,pt->m_dMaxLeft,*dRadius))) {
                    pt = pt->m_pLeftBranch;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        if (coordClosest) *coordClosest = pcoordClosest;
        if (objClosest) *objClosest = pobjClosest;
        return  pcoordClosest?CNEARTREE_SUCCESS:CNEARTREE_NOT_FOUND;
    }   /* Nearest */
    
    
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
     coordFarthest is a pointer to the returned farthest point
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
                               const void FAR * coord ) {
        double   dDR, dDL;
        CVectorHandle sStack;
        
        CNearTreeNodeHandle pt;
        void FAR * pobjFarthest;
        void FAR * pcoordFarthest;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        pobjFarthest = NULL;
        pcoordFarthest = NULL;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( treehandle->m_ptDelayObjs ) {
            if (CNearTreeCompleteDelayedInsert(treehandle)!=CNEARTREE_SUCCESS) return CNEARTREE_BAD_ARGUMENT;
        }
        
        pt = treehandle->m_ptTree;
        
        if (!pt) return CNEARTREE_BAD_ARGUMENT;
        
        if (!(pt->m_iflags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        CVectorCreate(&sStack,sizeof(CNearTreeNodeHandle),10);
        
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = CNearTreeDist(treehandle, (void FAR *)coord,  pt->m_coordRight);
                if (dDR >= *dRadius ) {
                    *dRadius = dDR;
                    pobjFarthest = pt->m_ptobjRight;
                    pcoordFarthest = pt->m_coordRight;
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (TRIANG(*dRadius,dDR,pt->m_dMaxRight))) {
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                dDL = CNearTreeDist(treehandle, (void FAR *)coord, pt->m_coordLeft);
                if (dDL >= *dRadius ) {
                    *dRadius = dDL;
                    pobjFarthest = pt->m_ptobjLeft;
                    pcoordFarthest = pt->m_coordLeft;
                }
                if (pt->m_iflags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_iflags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (TRIANG(*dRadius,dDL,pt->m_dMaxLeft))) {
                    pt = pt->m_pLeftBranch;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == end && CVectorSize(sStack) != 0 ) {
                CVectorGetElement(sStack,&pt,CVectorSize(sStack)-1);
                CVectorRemoveElement(sStack,CVectorSize(sStack)-1);
                eDir = right;
            }
            
        }    
        
        CVectorFree(&sStack);
        if (coordFarthest) *coordFarthest = pcoordFarthest;
        if (objFarthest) *objFarthest = pobjFarthest;
        return  pcoordFarthest?CNEARTREE_SUCCESS:CNEARTREE_NOT_FOUND;
    }
    
#ifdef __cplusplus
    
}

#endif

