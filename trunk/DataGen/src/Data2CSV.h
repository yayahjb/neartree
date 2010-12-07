// Data2CSV.h
//
// utility code for the various small programs intended to be
// used with other programs for generating test data
// sets for NearTree. Also included is a limited n-dimensional
// vector class
//

#ifndef DATA2CSV_H
#define DATA2CSV_H

#include <string>
#include <vector>
#include <iostream>

#include "vector_3d.h"

template <typename T>
T MAX( const T& t1, const T& t2 )
{
    if ( t1 > t2 )
    {
        return( t1 );
    }
    else
    {
        return( t2 );
    }
}

template <typename T>
T MIN( const T& t1, const T& t2 )
{
    if ( t1 < t2 )
    {
        return( t1 );
    }
    else
    {
        return( t2 );
    }
}

template <typename T>
T sign_t( const T& t )
{
    return( (t<T(0)) ? T(-1) : T(1) );
}

class vecN;

std::vector<std::string> GetNextLine( std::string& line );

std::vector<Vector_3> ReadVector_3File( std::istream& str );

void WriteVector_3File( const std::vector<Vector_3>& v );

void WriteVectorFile( const std::vector<vecN>& v );

std::pair<vecN,vecN> FindBox( const std::vector<vecN>& v );
void BoxIt( std::vector<vecN>& v );

class vecN;
std::vector<vecN> ReadGeneralFile( std::istream& str );


/*=======================================================================*/
/* make an N-dimension vector class for testing */
class vecN
{
public:
    friend vecN operator* ( const double&, const vecN& );

private:
    std::vector<double> pd;
public:
    int dim;
    double length;  // just for a signature for debugging

    explicit vecN( const long n );
    explicit vecN( const int n );

    vecN( void )
        : dim( 0 )
        , length( 0.0 )
    {
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    explicit vecN( const long n, const double d ):
    dim(n), length(0.0)
    //-------------------------------------------------------------------------------------
    {
        pd.resize( n );
        for( int i=0; i<dim; ++i )
        {
            pd[i] = d;
        }
        length = Norm( );
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    explicit vecN( const std::vector<double>& v )
    //-------------------------------------------------------------------------------------
    {
        pd = v;
        dim = (int)v.size( );
        length = Norm( );
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    ~vecN( void )
    //-------------------------------------------------------------------------------------
    {
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    vecN( const vecN& v )
    //-------------------------------------------------------------------------------------
    {
        length = v.length;
        pd = v.pd;
        dim = v.dim;
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    vecN operator+( const vecN& v ) const /* USERS: be sure to make both const */
    //-------------------------------------------------------------------------------------
    {
        vecN vtemp = v;
        for( int i=0; i<dim; ++i )
        {
            vtemp.pd[i] = this->pd[i] + v.pd[i];
        }
        vtemp.length = Norm(  );
        return( vtemp );
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    vecN operator-( const vecN& v ) const /* USERS: be sure to make both const */
    //-------------------------------------------------------------------------------------
    {
        vecN vtemp = v;
        for( int i=0; i<dim; ++i )
        {
            vtemp.pd[i] = this->pd[i] - v.pd[i];
        }
        vtemp.length = Norm(  );
        return( vtemp );
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    double Norm( void ) const
    //-------------------------------------------------------------------------------------
    {
        double dtemp = 0.0;
        for( int i=0; i<dim; ++i )
        {
            dtemp += pd[i]*pd[i];  //  L2 measure here
        }
        return( double(sqrt( dtemp )) );
    }

    //-----------------------------------------------------------------------------
    // Name: operator[]()
    // Description: access function for the components of a vector
    //
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    inline double operator[] ( const int& i ) const
    //-------------------------------------------------------------------------------------
    { // used for access (e.g. b = a[i];)
        unsigned int n = ( i<0 ) ? 0 : i;
        if( (unsigned)i > pd.size( )-1 ) n = (unsigned int)pd.size( )-1;
        return (pd[n]);
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    inline double& operator[] ( const int& i )
    //-------------------------------------------------------------------------------------
    { // for use as an assigment (e.g. a[i] = 3.0;)
        unsigned int n = ( i<0 ) ? 0 : i;
        if( (size_t)i > pd.size( )-1 ) n = (unsigned int)pd.size( )-1;
        return (pd[n]);
    }

     /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
   size_t size( void ) const
    //-------------------------------------------------------------------------------------
    {
        return( pd.size( ) );
    }

////-----------------------------------------------------------------------------
//// Name: operator+=()
//// Description: add two vectors
////
///*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//    inline vecN operator+ ( const vecN& vv )
//        //-----------------------------------------------------------------------------
//    {
//        vecN v( vv.pd.size( ), 0.0 );
//        for ( unsigned int i=0; i<vv.pd.size(); ++i )
//        {
//            v.pd[i] = (*this)[i] + vv[i];
//        }
//
//        return( v );
//    }

};  // end vecN



#endif  // DATA2CSV_H
