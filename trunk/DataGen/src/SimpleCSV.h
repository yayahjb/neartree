//=============================================================================
//      Micro Encoder Inc.
//      Copyright © 1994 - 2010. All Rights Reserved.
//
//      File    :       SimpleCSV.h
//      Name    :       Larry Andrews
//      Date    :       01/12/2010
//      Desc    :       utility functions for creating and parsing comma 
//                      separated variable text
//=============================================================================
#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace FXTK {

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static inline bool comma ( const wchar_t c )
//-----------------------------------------------------------------------------
{
   return( c == wchar_t(',') );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static inline bool notblank( const wchar_t c )
//-----------------------------------------------------------------------------
{
   return( c != wchar_t(' ') );
}

static std::vector<std::string> VectorFromCSV( const std::string& s );
static std::string              CSVFromVector( const std::vector<std::string>& v );
static std::string              GetNthField  ( const std::string& s, const size_t n );

static std::vector<std::wstring> VectorFromCSV( const std::wstring& s );
static std::wstring              CSVFromVector( const std::vector<std::wstring>& v );
static std::wstring              GetNthField  ( const std::wstring& s, const size_t n );

// THE FOLLOWING ARE FOR CString ONLY
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template< typename T>
static std::vector<T> VectorFromCSV( const T& s )
//--------------------------------------------------------------------------
{
   std::vector<T> v;
   if ( s.IsEmpty( ) )
      return( v );

   int i = 0;
   while ( i < s.GetLength( ) )
   {
      // find the end of the next token
      int j = s.Find(T(','), i);
      if ( j <= -1 )
      {
         // no comma found
         j = s.GetLength();
      }

      // copy the characters in [i,j]
      v.push_back( s.Mid(i, j-i) );

      // no terminating comma !!!
      if ( j >= s.GetLength() )
      {
         break;
      }

      ++j;

      i = j;
   }
   return( v );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template< typename T>
static T CSVFromVector( const std::vector<T>& v )
//--------------------------------------------------------------------------
{
   T s;
   T sTemp;
   if ( v.size( ) >  0 )
   {
      s = v[0];
      for ( size_t i=1U; i<v.size( ); ++i )
      {
         sTemp = v[i];
         s += T(",");
         s += sTemp;
      }
   }
   return( s );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template< typename T>
static T GetNthField  ( const T& s, const size_t n )
//--------------------------------------------------------------------------
{
   const std::vector<T> sv = VectorFromCSV( s );
   return( (n<sv.size()-1) ? sv[n] : T( ) ); // construct an empty string
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Name: CSVFromVector()
// Original Defect: CR08883(MRW)
// Description:
//   Composes a csv wstring using the elements of an input vector of wstrings
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::wstring CSVFromVector( const std::vector<std::wstring>& v )
//-----------------------------------------------------------------------------
{
   std::wstring s;
   std::wstring sTemp;
   if ( v.size( ) >  0 )
   {
      s = v[0];
      for ( size_t i=1U; i<v.size( ); ++i )
      {
         sTemp = v[i];
         s += wchar_t(',');
         s += sTemp;
      }
   }
   return( s );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Name: VectorFromCSV()
// Original Defect: CR08883(MRW)
// Description:
//   Takes a text string and returns a vector of strings containing the tokens 
//   found usings only commas as separators. Leading and trailing blanks are 
//   deleted.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::vector<std::string> VectorFromCSV( const std::string& s )
//-----------------------------------------------------------------------------
{
   std::vector<std::string> v;
   if ( s.empty( ) )
      return( v );

   std::string::const_iterator i;
   i = s.begin( );
   while ( i != s.end( ) )
   {
      // find the end of the next token
      std::string::const_iterator j = std::find_if( i, s.end( ), comma );

      // copy the characters in [i,j]
      v.push_back( std::string( i, j ) );
      // no terminating comma !!!
      if ( j == s.end( ) )
      {
         break;
      }
      ++j;
      i = std::find_if( j, s.end( ), notblank );
   }
   return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Name: CSVFromVector()
// Original Defect: CR08883(MRW)
// Description:
//   Composes a csv string using the elements of an input vector of strings
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::string CSVFromVector( const std::vector<std::string>& v )
//-----------------------------------------------------------------------------
{
   std::string s;
   std::string sTemp;
   if ( v.size( ) >  0 )
   {
      s = v[0];
      for ( size_t i=1U; i<v.size( ); ++i )
      {
         sTemp = v[i];
         s += ',';
         s += sTemp;
      }
   }
   return( s );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Name: GetNthField()
// Original Defect: CR08883(MRW)
// Description:
//   Takes as input a CSV string and a field number, parses the CSV string, and
//   return the n-th (zero-based) field. If that field does not exist, returns 
//   an empty string.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::string GetNthField( const std::string& s, const size_t n )
//-----------------------------------------------------------------------------
{
   const std::vector<std::string> sv = VectorFromCSV( s );
   return( (n<sv.size()-1) ? sv[n] : std::string( ) );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Name: VectorFromCSV()
// Original Defect: CR08883(MRW)
// Description:
//   Takes a text wstring and returns a vector of wstrings containing the tokens 
//   found usings only commas as separators. Leading and trailing blanks are 
//   deleted.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::vector<std::wstring> VectorFromCSV( const std::wstring& s )
//-----------------------------------------------------------------------------
{
   std::vector<std::wstring> v;
   if ( s.empty( ) )
      return( v );

   std::wstring::const_iterator i;
   i = s.begin( );
   while ( i != s.end( ) )
   {
      // find the end of the next token
      std::wstring::const_iterator j = std::find_if( i, s.end( ), comma );

      // copy the characters in [i,j]
      v.push_back( std::wstring( i, j ) );
      // no terminating comma !!!
      if ( j == s.end( ) )
      {
         break;
      }
      ++j;
      i = std::find_if( j, s.end( ), notblank );
   }
   return( v );
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Name: GetNthField()
// Original Defect: CR08883(MRW)
// Description:
//   Takes as input a CSV wstring and a field number, parses the CSV wstring,
//   and return the n-th (zero-based) field. If that field does not exist, 
//   returns an empty wstring.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
static std::wstring GetNthField( const std::wstring& s, const size_t n )
//-----------------------------------------------------------------------------
{
   const std::vector<std::wstring> sv = VectorFromCSV( s );
   return( (n<sv.size()-1) ? sv[n] : std::wstring( ) );
}

}//namespace FXTK;
