//*
//*  main.cpp
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

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <vector>

#include "TNear.h"
#include "v.h"

int main ( )
{
   CNearTree <v> vTree;
   v vBest = DBL_MAX;
   long i,j,k;
   const long lMaxRow = 10;
   std::vector <v> vReturn;

   srand( (unsigned)time( NULL ) );  // use the current time to seed the
                              // random number generator
//---------------------------------------
// build up a library of points to search among
//---------------------------------------
   for ( k=-1; k<=lMaxRow; k++ )
   {
      for ( j=-1; j<=lMaxRow; j++) 
      {
         for ( i= lMaxRow ; i>=-1;  i-- ) 
         {
            vTree.Insert( v((double)i, (double)j, (double)k) );
         }  // for i
      }     // for j
   }        // for k
    cout
      << endl;

//---------------------------------------
//   Done building the tree; now try a retrieval
//---------------------------------------

   double   dRad = 0.6;
   long lReturn;
   for ( i=0;  i<10; i++ )
   {
      dRad += 0.05;
      double x = rand( ) * double( lMaxRow ) / RAND_MAX;
      double y = x;
      double z = ( 1.25 * double(lMaxRow) - 1.5 * x );
      v vSearch( x, 0.5*(x+y), z );
      cout
         << "Trial " << i << " from probe point " << vSearch << endl;




// find the nearest point to vSearch
      if ( vTree.NearestNeighbor( dRad, vBest, vSearch ) )
      {
         cout
            << " Closest distance " << (double) ( vSearch - vBest ) << " to " << vBest << endl;
      }
      else
      {
         cout
            << " ***** nothing within " << dRad << " of " << vSearch << endl;
      }

// find the farthest point from vSearch
      if ( vTree.FarthestNeighbor ( vBest, vSearch ) )
      {
         cout 
            << " Farthest distance " << (double) ( vSearch - vBest ) << " to " << vBest << endl;
      }
      else
      {
         cout << " No Farthest object found" << endl;
      }

// search for all points within a "sphere" out to radius dRad
      if ( ( lReturn = vTree.FindInSphere( dRad, vReturn, vSearch )) > 0 ) 
      {
         cout << " Returned " << lReturn << " items within " << dRad << " of "<< vSearch << endl;
         std::vector <v>::iterator iv;
         for ( iv=vReturn.begin( ); iv<vReturn.end( ); iv++ )
         {
            cout << "\t" << *iv << endl;
         }
      }

      cout << " -------------------------------------------------------------"<<endl;
   }  // for i

   vTree.CNearTree<v>::~CNearTree( );

   return ( EXIT_SUCCESS );
}