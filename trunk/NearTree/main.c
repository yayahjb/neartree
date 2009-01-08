/*
 *  main.c
 *  NearTree
 *
 *  Copyright 2001, 2008 Larry Andrews.  All rights reserved
 *  Revised 12 Dec 2008 for sourceforge release -- H. J. Bernstein
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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <CVector.h>


#include "CNearTree.h"

int main ( int argc, char** argv )
{
    
    CNearTreeHandle treehandle;
    double v[3],vSearch[3];
    double * vBest;
    void * vvBest;
    CVectorHandle vReturn;
    CVectorHandle oReturn;
    long i,j,k;
    const long lMaxRow = 10;
    double   dRad = 0.6;
    
    
    CNearTreeCreate(&treehandle,3,CNEARTREE_TYPE_DOUBLE);
    CVectorCreate(&vReturn,sizeof(void *),10);
    CVectorCreate(&oReturn,sizeof(void *),10);
    
    if (argc <= 1) {
        srand( (unsigned int)time( NULL ) );  /* use the current time to seed the
                                                 random number generator */
    } else {
        srand((unsigned int)atoi(argv[1]));
    }
    /*---------------------------------------
     build up a library of points to search among
     ---------------------------------------*/
    for ( k=-1; k<=lMaxRow; k++ )
    {
        for ( j=-1; j<=lMaxRow; j++) 
        {
            for ( i= lMaxRow ; i>=-1;  i-- ) 
            {  v[0] = (double)i; v[1] = (double)j; v[2] = (double)k;
                CNearTreeInsert(treehandle,&v[0],NULL);
            }  /* for i */
        }     /* for j */
    }        /* for k */
    fprintf(stdout,"\n");
    
    /*---------------------------------------
     Done building the tree; now try a retrieval
     ---------------------------------------*/
    
    for ( i=0;  i<10; i++ )
    {  double x, y, z;
        dRad += 0.05;
        x = rand( ) * ((double) lMaxRow ) / RAND_MAX;
        y = x;
        z = ( 1.25 * ((double) lMaxRow) - 1.5 * x );
        vSearch[0] = x; vSearch[1] = 0.5*(x+y); vSearch[2] = z;
        fprintf(stdout,"Trial %ld from probe point [%g, %g, %g]\n",
                i, vSearch[0], vSearch[1], vSearch[2]);
        
        
        /* find the nearest point to vSearch */
        if ( !CNearTreeNearestNeighbor(treehandle,dRad,&vvBest,NULL,vSearch))
        {   vBest = (double *)vvBest;
            fprintf(stdout," Closest distance %g to [%g, %g, %g]\n",
                   sqrt(CNearTreeDistsq((void *)vSearch,vvBest,3,CNEARTREE_TYPE_DOUBLE)),
                   vBest[0],vBest[1],vBest[2]);
        }
        else
        {  fprintf(stdout," ***** nothing within %g of [%g, %g, %g]\n",
                   dRad,  vSearch[0], vSearch[1], vSearch[2]);
        }
        
        /* find the farthest point from vSearch */
        if ( !CNearTreeFarthestNeighbor(treehandle,&vvBest,NULL,vSearch))
        {   vBest = (double *)vvBest; 
            fprintf(stdout," Farthest distance %g to [%g, %g, %g]\n",
                   sqrt(CNearTreeDistsq(vSearch,vBest,3,CNEARTREE_TYPE_DOUBLE)),
                   vBest[0],vBest[1],vBest[2]);
        }
        else
        {  fprintf(stdout," No Farthest object found\n");
        }
        
        /* search for all points within a "sphere" out to radius dRad */
        
        if ( !CNearTreeFindInSphere( treehandle, dRad, vReturn, oReturn, vSearch,1 ) ) 
        {
            size_t index;
            void * foundpoint;
            
            fprintf(stdout," Returned %lu items within %g of [%g,%g,%g]\n",
                    (long unsigned int)vReturn->size, dRad, vSearch[0], vSearch[1], vSearch[2]);
            for (index=0; index < CVectorSize(vReturn); index++) {
                CVectorGetElement(vReturn,&foundpoint,index);
                fprintf (stdout,"\t [%g,%g,%g]\n",((double *)foundpoint)[0],
                         ((double *)foundpoint)[1],
                         ((double *)foundpoint)[2]);
            }
        }
        
        fprintf(stdout," -------------------------------------------------------------\n");
    }  /* for i */
    
    CNearTreeFree(&treehandle);
    
    return ( EXIT_SUCCESS );
}
