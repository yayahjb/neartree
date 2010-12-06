// svg.cpp
//
// svg reads from standard input and writes to standard output. It writes
// an svg file for displaying the data set (input as a n-space vector in
// csv). The command line is read to get a scale factor for the display.
// Only the first two coordinates are plotted.
//
// This is intended to be used with other programs for generating test data
// sets for NearTree.
//

#include <vector>
#include <iostream>

#include "Data2CSV.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SVG_Output(
    const long picSize,
    const std::vector<vecN>& points )
//---------------------------------------------------------------------
{
    const long sqSide = (long)(2.1 * picSize);
    fprintf( stdout, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
    fprintf( stdout, "<svg width=\"%ld\" height=\"%ld\" version=\"1.1\"\n", sqSide, sqSide);
    fprintf( stdout, "xmlns=\"http://www.w3.org/2000/svg\">\n");
    fprintf( stdout, "<rect width=\"%ld\" height=\"%ld\" style=\"fill:rgb(255,255,255); stroke-width:10;  stroke:rgb(255,255,255)\"/>\n",
        sqSide, sqSide);


    double maxx = -DBL_MAX;
    double minx = DBL_MAX;
    double maxy = -DBL_MAX;
    double miny = DBL_MAX;
    for( unsigned int i=0; i<points.size( ); ++i )
    {
       maxx = MAX(maxx,points[i][0]);
       maxy = MAX(maxy,points[i][1]);
       minx = MIN(minx,points[i][0]);
       miny = MIN(miny,points[i][1]);
    }

    const double scale = MAX( (maxx-minx), (maxy-miny) ) / picSize;

    for( unsigned int i=0; i<points.size( ); ++i )
    {
        const int r = 2;
        const double xd = (points[i][0]-minx)/scale;
        const double yd = (points[i][1]-miny)/scale;
        const int x = int( sign_t(xd)*xd );
        const int y = int( sign_t(yd)*yd );
        fprintf( stdout, "<circle r=\"%ld\" cx=\"%ld\" cy=\"%ld\" stroke=\"black\" stroke-width=\"2\" fill=\"red\"/>\n",
            r, x, y );
    }

    fprintf( stdout, "</svg>\n");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char* argv[])
//---------------------------------------------------------------------
{
    double scale = 0.0;
    // see if there is a scale factor on the command line
    if ( argc>1 )
        scale = atof( argv[1] );

    if ( scale <= 0.0 ) scale = 500.0;

    const std::vector<vecN> v = ReadGeneralFile( std::cin );

    SVG_Output( (long)scale, v );

    return 0;}

