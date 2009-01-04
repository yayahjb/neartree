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
    
#include "CNearTree.h"
#include <math.h>
    
    
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
     int CNearTreeCreate ( CNearTreeHandle FAR * treehandle, 
     size_t treedim, int treetype)
     
     Create a CNearTree
     
     returns a pointer to the newly allocated block of memory as a 
     CNearTreeHandle in *treehandle
     
     
     treedim -- the dimension of the vectors
     treetype -- double or integer flag for type of the vectors
     CNEARTREE_TYPE_DOUBLE for double
     CNEARTREE_TYPE_INTEGER for integer
     
     creates an empty tree with no right or left node and with the dMax-below
     set to negative values so that any match found will be stored since it will
     greater than the negative value
     
     =======================================================================
     */
    
    int CNearTreeCreate ( CNearTreeHandle FAR * treehandle, 
                         size_t treedim, int treetype) {
        
        size_t index;
        
        size_t ntsize, nvsize;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        if (treedim < 0 
            || (treetype != CNEARTREE_TYPE_DOUBLE 
                && treetype != CNEARTREE_TYPE_INTEGER)) return CNEARTREE_BAD_ARGUMENT;
        
        
        if (treetype == CNEARTREE_TYPE_INTEGER) {
            ntsize = (sizeof(CNearTree)+(sizeof(int)-1))/sizeof(int);
            ntsize *= sizeof(int);
            nvsize = sizeof(int)*treedim;
            *treehandle = (CNearTreeHandle)MALLOC(ntsize+2*nvsize);
        } else {
            ntsize = (sizeof(CNearTree)+(sizeof(double)-1))/sizeof(double);
            ntsize *= sizeof(double);
            nvsize = sizeof(double)*treedim;
            *treehandle = (CNearTreeHandle)MALLOC(ntsize+2*nvsize);
        }
        
        if (!(*treehandle)) {
            return CNEARTREE_MALLOC_FAILED;
        }
        
        (*treehandle)->m_coordLeft = ((char *)(*treehandle))+ntsize;
        (*treehandle)->m_coordRight = ((char *)(*treehandle))+ntsize+nvsize;
        if (treetype == CNEARTREE_TYPE_DOUBLE) {
            for (index=0; index < treedim; index++) {
                ((double *)(*treehandle)->m_coordLeft)[index]    = 0;
                ((double *)(*treehandle)->m_coordRight)[index]   = 0;
            }
        } else {
            for (index=0; index < treedim; index++) {
                ((int *)(*treehandle)->m_coordLeft)[index]    = 0;
                ((int *)(*treehandle)->m_coordRight)[index]   = 0;
            }
        }
        (*treehandle)->m_pLeftBranch  = 0;
        (*treehandle)->m_pRightBranch = 0;
        (*treehandle)->m_dMaxLeft     = -1.;        /* negative for an empty branch */
        (*treehandle)->m_dMaxRight    = -1.;        /* negative for an empty branch */
        (*treehandle)->m_ptobjLeft    = NULL;
        (*treehandle)->m_ptobjRight   = NULL;
        (*treehandle)->m_flags        = treetype;   /* no data, no children         */
        (*treehandle)->m_dimension    = treedim;    /* number of ints or doubles    */
        return CNEARTREE_SUCCESS;
    }
    
    /*
     =======================================================================
     int CNearTreeFree ( CNearTreeHandle FAR * treehandle )
     
     Free a CNearTree
     
     recursively frees the NearTree with the handle *treehandle
     and nulls the treehandle.
     
     note that the objects referenced are not freed.
     =======================================================================
     */
    
    int CNearTreeFree ( CNearTreeHandle FAR * treehandle ) {
        
        int errorcode;
        
        if (!treehandle) return CNEARTREE_BAD_ARGUMENT;
        
        errorcode = CNEARTREE_SUCCESS;
        
        if ((*treehandle)->m_pLeftBranch) {
            errorcode |= CNearTreeFree(&((*treehandle)->m_pLeftBranch));
        }
        
        if ((*treehandle)->m_pRightBranch) {
            errorcode |= CNearTreeFree(&((*treehandle)->m_pRightBranch));
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
        return (((treehandle->m_flags)&CNEARTREE_DATA_OR_CHILDREN)  == 0)?0:1;
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
        
        double dTempRight =  0.;
        double dTempLeft  =  0.;
        int errorcode = 0;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        
        if ( (treehandle->m_flags)&CNEARTREE_FLAG_RIGHT_DATA ) {
            dTempRight = CNearTreeDistsq((void FAR *)coord,treehandle->m_coordRight,
                                         treehandle->m_dimension, treehandle->m_flags&CNEARTREE_TYPE);
        }
        
        if ( (treehandle->m_flags)&CNEARTREE_FLAG_LEFT_DATA ) {
            dTempLeft = CNearTreeDistsq((void FAR *)coord,treehandle->m_coordLeft,
                                        treehandle->m_dimension, treehandle->m_flags&CNEARTREE_TYPE);
        }
        
        if ( !((treehandle->m_flags)&CNEARTREE_FLAG_LEFT_DATA) ) {
            size_t index;
            if (treehandle->m_flags&CNEARTREE_TYPE_DOUBLE) {
                for (index=0; index<treehandle->m_dimension; index++)  {
                    ((double *)(treehandle->m_coordLeft))[index]=((double *)coord)[index];    	
                }
            } else {
                for (index=0; index<treehandle->m_dimension; index++)  {
                    ((int *)(treehandle->m_coordLeft))[index]=((int *)coord)[index];    	
                }
            }
            treehandle->m_dMaxLeft = -1.;
            treehandle->m_ptobjLeft = (void *)obj;
            treehandle->m_flags |= CNEARTREE_FLAG_LEFT_DATA;
            return CNEARTREE_SUCCESS;
        } else if (  !((treehandle->m_flags)&CNEARTREE_FLAG_RIGHT_DATA)  ){ 
            size_t index;
            if (treehandle->m_flags&CNEARTREE_TYPE_DOUBLE) {
                for (index=0; index<treehandle->m_dimension; index++)  {
                    ((double *)treehandle->m_coordRight)[index]=((double *)coord)[index];    	
                }
            } else {
                for (index=0; index<treehandle->m_dimension; index++)  {
                    ((int *)treehandle->m_coordRight)[index]=((int *)coord)[index];    	
                }
            }
            treehandle->m_dMaxRight = -1.;
            treehandle->m_ptobjRight = (void *)obj;
            treehandle->m_flags |= CNEARTREE_FLAG_RIGHT_DATA;
            return CNEARTREE_SUCCESS;
        } else if ( dTempLeft > dTempRight ) {
            dTempRight = sqrt(dTempRight);
            if (  !((treehandle->m_flags)&CNEARTREE_FLAG_RIGHT_CHILD) ) {
                if ( (errorcode = CNearTreeCreate(&(treehandle->m_pRightBranch), 
                                                    treehandle->m_dimension, 
                                                    treehandle->m_flags&CNEARTREE_TYPE))) return errorcode;
                treehandle->m_flags |= CNEARTREE_FLAG_RIGHT_CHILD;
                treehandle->m_dMaxRight = dTempRight;
            }
            /* note that the next line assumes that m_dMaxRight is negative for a new node */
            if ( treehandle->m_dMaxRight < dTempRight ) 
                treehandle->m_dMaxRight = dTempRight;
            return CNearTreeInsert(treehandle->m_pRightBranch, coord, obj);
        } else { /* ((double)(t - *m_tLeft) <= (double)(t - *m_tRight) ) */
            dTempLeft = sqrt(dTempLeft);
            if (  !((treehandle->m_flags)&CNEARTREE_FLAG_LEFT_CHILD) ) {
                if ( (errorcode = CNearTreeCreate(&(treehandle->m_pLeftBranch),
                                                    treehandle->m_dimension, 
                                                    treehandle->m_flags&CNEARTREE_TYPE))) return errorcode;
                treehandle->m_flags |= CNEARTREE_FLAG_LEFT_CHILD;
                treehandle->m_dMaxLeft = dTempLeft;
            }
            /* note that the next line assumes that m_dMaxLeft is negative for a new node */
            if ( treehandle->m_dMaxLeft < dTempLeft ) 
                treehandle->m_dMaxLeft  = dTempLeft;
            return CNearTreeInsert(treehandle->m_pLeftBranch, coord, obj );
        }
        
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
        if (!(treehandle->m_flags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
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
        if (!(treehandle->m_flags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
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
        double dradiussq;
        double dDRsq, dDLsq;
        int nopoints;
        CVectorHandle sStack;
        CNearTreeHandle pt;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        pt = treehandle;
        
        nopoints = 1;
        
        if (dRadius < 0.) return 1;
        
        dradiussq = dRadius*dRadius;

        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;

        if (!(treehandle->m_flags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        CVectorCreate(&sStack,sizeof(CNearTreeHandle),10);
        
        /* clear the contents of the return vector so that things don't accidentally accumulate */
        if (resetcount) {
            if (coordClosest) CVectorClear( coordClosest );
            if (objClosest) CVectorClear( objClosest );
        }
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
               dDRsq = DBL_MAX;
                if ((pt->m_flags)&CNEARTREE_FLAG_RIGHT_DATA) {
                  dDRsq = CNearTreeDistsq((void FAR *)coord,
                                        pt->m_coordRight,
                                        pt->m_dimension,
                                        pt->m_flags&CNEARTREE_TYPE);
                  if (dDRsq <= dradiussq ) {
                    nopoints = 0;
                    if (coordClosest) CVectorAddElement(coordClosest,&(pt->m_coordRight));
                    if (objClosest && pt->m_ptobjRight) CVectorAddElement(objClosest,&(pt->m_ptobjRight));
                  }
                }
                if ((pt->m_flags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (pt->m_dMaxRight+dRadius)*(pt->m_dMaxRight+dRadius) >= dDRsq) {
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                if ((pt->m_flags)&CNEARTREE_FLAG_LEFT_DATA) {
                  dDLsq = CNearTreeDistsq((void FAR *)coord,
                                        pt->m_coordLeft,
                                        pt->m_dimension,
                                        pt->m_flags&CNEARTREE_TYPE);
                  if (dDLsq <= dradiussq ) {
                    nopoints = 0;
                    if (coordClosest) CVectorAddElement(coordClosest,&(pt->m_coordLeft));
                    if (objClosest &&  pt->m_ptobjLeft) CVectorAddElement(objClosest,&(pt->m_ptobjLeft));
                  }
                }
                if (pt->m_flags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_flags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (pt->m_dMaxLeft+dRadius)*(pt->m_dMaxLeft+dRadius) >= dDLsq){
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
        CNearTreeHandle pt;
        void FAR * pobjClosest;
        void FAR * pcoordClosest;
        enum  { left, right, end } eDir;
        
        eDir = left; /* examine the left nodes first */
        
        pt = treehandle;
        pobjClosest = NULL;
        pcoordClosest = NULL;
        
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        if (!(treehandle->m_flags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        CVectorCreate(&sStack,sizeof(CNearTreeHandle),10);
        
        
        while (!(eDir == end && CVectorSize(sStack) == 0)) {
            
            if ( eDir == right ) {
                dDR = sqrt(CNearTreeDistsq((void FAR *)coord,
                                        pt->m_coordRight,
                                        pt->m_dimension,
                                        pt->m_flags&CNEARTREE_TYPE));
                if (dDR < *dRadius ) {
                    *dRadius = dDR;
                    pobjClosest = pt->m_ptobjRight;
                    pcoordClosest = pt->m_coordRight;
                }
                if ((pt->m_flags&CNEARTREE_FLAG_RIGHT_CHILD)&& 
                    (pt->m_dMaxRight+*dRadius) >= dDR) {
                    /* we did the left and now we finished the right, go down */
                    pt = pt->m_pRightBranch;
                    eDir = left;
                } else {
                    eDir = end;
                }
            }
            if ( eDir == left ) {
                dDL = sqrt(CNearTreeDistsq((void FAR *)coord,
                                        pt->m_coordLeft,
                                        pt->m_dimension,
                                        pt->m_flags&CNEARTREE_TYPE));
                if (dDL < *dRadius ) {
                    *dRadius = dDL;
                    pobjClosest = pt->m_ptobjLeft;
                    pcoordClosest = pt->m_coordLeft;
                }
                if (pt->m_flags&CNEARTREE_FLAG_RIGHT_DATA) {
                    CVectorAddElement(sStack,&pt);
                }
                if ((pt->m_flags&CNEARTREE_FLAG_LEFT_CHILD)&&
                    (pt->m_dMaxLeft+*dRadius) >= dDL){
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
        double   dTempRadius = DBL_MAX;
        double   dradiussq;
        int  bRet = 0;
                
        if ( !treehandle || !coord ) return CNEARTREE_BAD_ARGUMENT;
        if (!(treehandle->m_flags&CNEARTREE_FLAG_LEFT_DATA)) return CNEARTREE_NOT_FOUND;
        
        dradiussq = (*dRadius)*(*dRadius);
        
        /* first test each of the left and right positions to see if
         one holds a point farther than the farthest so far discovered.
         the calling function is presumed initially to have set dRadius to a
         negative value before the recursive calls to  CNearTreeFindFarthest */
        
        if ((treehandle->m_flags&CNEARTREE_FLAG_LEFT_DATA ) && 
            (( dTempRadius = CNearTreeDistsq((void FAR *)coord,
                                             treehandle->m_coordLeft,
                                             treehandle->m_dimension,
                                             treehandle->m_flags&CNEARTREE_TYPE))>= dradiussq )) {
            *dRadius  = sqrt(dTempRadius);
            dradiussq = dTempRadius;
            *coordFarthest = treehandle->m_coordLeft;
            if ( objFarthest) *objFarthest = treehandle->m_ptobjLeft;
            bRet     =  1;
        }
        if ((treehandle->m_flags&CNEARTREE_FLAG_RIGHT_DATA ) && 
            (( dTempRadius = CNearTreeDistsq((void FAR *)coord,
                                             treehandle->m_coordRight,
                                             treehandle->m_dimension,
                                             treehandle->m_flags&CNEARTREE_TYPE))>= dradiussq )) {            
            *dRadius   = sqrt(dTempRadius);
            *coordFarthest = treehandle->m_coordRight;
            if (objFarthest) *objFarthest = treehandle->m_ptobjRight;
            bRet      = 1;
        }
        /*
         Now we test to see if the branches below might hold an object
         farther than the best so far found. The triangle rule is used
         to test whether it's even necessary to descend.
         */
        if (( treehandle->m_pLeftBranch  != 0 )  && ( treehandle->m_dMaxLeft>=0. )
            && (((*dRadius) - treehandle->m_dMaxLeft< 0.)||( (*dRadius) - treehandle->m_dMaxLeft )*( (*dRadius) - treehandle->m_dMaxLeft ) <= 
                CNearTreeDistsq ( (void FAR *)coord, 
                                 treehandle->m_coordLeft,
                                 treehandle->m_dimension,
                                 treehandle->m_flags&CNEARTREE_TYPE)))
        {
            if (!CNearTreeFindFarthest( treehandle->m_pLeftBranch,
                                       dRadius, coordFarthest, objFarthest, coord)) bRet = 1;
        }
        
        if (( treehandle->m_pRightBranch  != 0 ) && ( treehandle->m_dMaxRight>=0. )
            && (((*dRadius) - treehandle->m_dMaxRight< 0.)||( (*dRadius) - treehandle->m_dMaxRight )*( (*dRadius) - treehandle->m_dMaxRight ) <= 
                CNearTreeDistsq ( (void FAR *)coord,
                                 treehandle->m_coordRight,
                                 treehandle->m_dimension,
                                 treehandle->m_flags&CNEARTREE_TYPE)))        {
            if (!CNearTreeFindFarthest( treehandle->m_pRightBranch,
                                       dRadius, coordFarthest, objFarthest, coord)) bRet = 1;
        }
        
        return ( bRet?0:1 );
    }
    
#ifdef __cplusplus
    
}

#endif

