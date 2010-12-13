// testmaxdist.cpp
//
// testMaxDist is designed as a test harness for NearTree. It reads a csv file
// of general vectors. The number of data elements in each vector must be
// consistent throughout the entire input file, but the dimension is not
// constrained. Various tests are/can be executed and reported. Originally
// designed to connect to DataGen and a few related programs. However, it
// still contains functions such as testHamm that generate their own data.

#include <cfloat>
#include <climits>
#include <ctime>

#include "hammersley.h"
#include "vector_3d.h"
#include "Data2CSV.h"

#include "TNear.h"
#include "rhrand.h"

#ifdef _MSC_VER
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#endif

RHrand rhr;
void testBigVector(  );
void testDouble(  );
int g_errorCount;

/*=======================================================================*/
// testGeneral
//
// Designed to test input data sets of arbitrary size and dimension.
// The input can come from any source, but the main in testMaxDist
// gets its input from standard input.
/*=======================================================================*/
void testGeneral( const bool Do_Random_Insertion, const std::vector<vecN>& v )
//---------------------------------------------------------------------
{
    CNearTree<vecN> nt( v );

    if ( Do_Random_Insertion )
    {
        nt.CompleteDelayedInsertRandom( );
    }
    else
    {
        nt.CompleteDelayedInsert( );
    }

    if ( nt.empty( ) )
    {
        fprintf( stdout, "No input found\n" );
        return;
    }

    const int nTests = 10000;
    nt.SetFlags( CNearTree<vecN>::NTF_NoPrePrune );

    long nodevisits1 = (long)nt.GetNodeVisits( );
    const clock_t tc1 = std::clock();
    for ( int i=0; i<nTests; ++i )
    {
        vecN newPoint(v[0].dim);
        const vecN probe = vecN(v[0].dim);
        nt.NearestNeighbor( DBL_MAX, newPoint, probe );
    }
    {

        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const clock_t tc2 = std::clock();
        long nodevisits2 = (long)nt.GetNodeVisits( );

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.NearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( 0 );
        long nodevisits3 = (long)nt.GetNodeVisits( );
        //const clock_t tc3 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe(v[0].dim);
            nt.NearestNeighbor( DBL_MAX, newPoint, probe );
        }

        long nodevisits4 = (long)nt.GetNodeVisits( );
        //const clock_t tc4 = std::clock();

        //fprintf( stdout, "testGeneral treesize:%ld  tests:%d old: %f, new: %f, default: %f, ratio old:new %f\n", (long)nt.size( ), nTests,
        //        ((double)(tc2-tc1))/CLOCKS_PER_SEC,
        //        ((double)(tc3-tc2))/CLOCKS_PER_SEC,
        //        ((double)(tc4-tc3))/CLOCKS_PER_SEC,
        //        ((double)(tc2-tc1))/((double)(tc3-tc2)) );
        //fprintf( stdout, "testHGeneral tree size: %ld, dimension estimate: %ld, diameter est,: %g, mean spacing est.: %g, variance spacing est.: %g\n",
        //        (long)nt.size(), (long)nt.GetDimEstimate(), (double)nt.GetDiamEstimate(),
        //        (double)nt.GetMeanSpacing(), (double)nt.GetVarSpacing());

        fprintf( stdout, "CSV-balanced,%ld,%d,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,%.4f,%.4f\n",
            (long)nt.size( ),
            nTests,
            ((double)(tc2-tc1))/CLOCKS_PER_SEC,
            (double)(nodevisits2-nodevisits1)/(double)nTests,
            (double)(nodevisits3-nodevisits2)/(double)nTests,
            ((double)nodevisits4-nodevisits3)/(double)nTests,
            ((double)(nodevisits2-nodevisits1))/((double)(nodevisits4-nodevisits3)),
            nt.GetDimEstimate(),
            (double)nt.GetDiamEstimate(),
            (double)nt.GetMeanSpacing(),
            (double)nt.GetVarSpacing());
    }
    /*----------------------------end balanced tree test--------------------------------------------*/
    /*----------------------------end left branch first test--------------------------------------------*/
    {
        nt.SetFlags( CNearTree<vecN>::NTF_NoPrePrune );


        const long nodevisits1 = (long)nt.GetNodeVisits( );
        const clock_t tc1 = std::clock();
        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.LeftNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const long nodevisits2 = (long)nt.GetNodeVisits( );
        const clock_t tc2 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.LeftNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( 0 );
        const long nodevisits3 = (long)nt.GetNodeVisits( );
        //const clock_t tc3 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe(v[0].dim);
            nt.LeftNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        const long nodevisits4 = (long)nt.GetNodeVisits( );
        //const clock_t tc4 = std::clock();

        //fprintf( stdout, "testGeneral treesize:%ld  tests:%d old: %f, new: %f, default: %f, ratio old:new %f\n", (long)nt.size( ), nTests,
        //        ((double)(tc2-tc1))/CLOCKS_PER_SEC,
        //        ((double)(tc3-tc2))/CLOCKS_PER_SEC,
        //        ((double)(tc4-tc3))/CLOCKS_PER_SEC,
        //        ((double)(tc2-tc1))/((double)(tc3-tc2)) );
        //fprintf( stdout, "testHGeneral tree size: %ld, dimension estimate: %ld, diameter est,: %g, mean spacing est.: %g, variance spacing est.: %g\n",
        //        (long)nt.size(), (long)nt.GetDimEstimate(), (double)nt.GetDiamEstimate(),
        //        (double)nt.GetMeanSpacing(), (double)nt.GetVarSpacing());

        fprintf( stdout, "CSV-LEFT,%ld,%d,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,%.4f,%.4f\n",
            (long)nt.size( ),
            nTests,
            ((double)(tc2-tc1))/CLOCKS_PER_SEC,
            (double)(nodevisits2-nodevisits1)/(double)nTests,
            (double)(nodevisits3-nodevisits2)/(double)nTests,
            (double)(nodevisits4-nodevisits3)/(double)nTests,
            ((double)(nodevisits2-nodevisits1))/((double)(nodevisits4-nodevisits3)),
            nt.GetDimEstimate(),
            (double)nt.GetDiamEstimate(),
            (double)nt.GetMeanSpacing(),
            (double)nt.GetVarSpacing());
    }
    /*----------------------------end left branch first test--------------------------------------------*/
    /*----------------------------start short radius test--------------------------------------------*/
    {
        nt.SetFlags( CNearTree<vecN>::NTF_NoPrePrune );

        const clock_t tc1 = std::clock();
        const long nodevisits1 = (long)nt.GetNodeVisits( );
        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.ShortNearestNeighbor( DBL_MAX, newPoint, probe );
        }


        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const long nodevisits2 = (long)nt.GetNodeVisits( );
        const clock_t tc2 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.ShortNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( 0 );
        const long nodevisits3 = (long)nt.GetNodeVisits( );
        //const clock_t tc3 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe(v[0].dim);
            nt.ShortNearestNeighbor( DBL_MAX, newPoint, probe );
        }
        const long nodevisits4 = (long)nt.GetNodeVisits( );

        //const clock_t tc4 = std::clock();

        //fprintf( stdout, "testGeneral treesize:%ld  tests:%d old: %f, new: %f, default: %f, ratio old:new %f\n", (long)nt.size( ), nTests,
        //        ((double)(tc2-tc1))/CLOCKS_PER_SEC,
        //        ((double)(tc3-tc2))/CLOCKS_PER_SEC,
        //        ((double)(tc4-tc3))/CLOCKS_PER_SEC,
        //        ((double)(tc2-tc1))/((double)(tc3-tc2)) );
        //fprintf( stdout, "testHGeneral tree size: %ld, dimension estimate: %ld, diameter est,: %g, mean spacing est.: %g, variance spacing est.: %g\n",
        //        (long)nt.size(), (long)nt.GetDimEstimate(), (double)nt.GetDiamEstimate(),
        //        (double)nt.GetMeanSpacing(), (double)nt.GetVarSpacing());
        fprintf( stdout, "CSV-SHORT,%ld,%d,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,%.4f,%.4f\n",
            (long)nt.size( ),
            nTests,
            ((double)(tc2-tc1))/CLOCKS_PER_SEC,
            (double)(nodevisits2-nodevisits1)/(double)nTests,
            (double)(nodevisits3-nodevisits2)/(double)nTests,
            (double)(nodevisits4-nodevisits3)/(double)nTests,
            ((double)(nodevisits2-nodevisits1))/((double)(nodevisits4-nodevisits3)),
            nt.GetDimEstimate(),
            (double)nt.GetDiamEstimate(),
            (double)nt.GetMeanSpacing(),
            (double)nt.GetVarSpacing());
    }
    /*----------------------------end short radius test--------------------------------------------*/

}

/*=======================================================================*/
// testInputVector
//
// confirm that the input data vector elements have consistent dimension
/*=======================================================================*/
bool testInputVector( const std::vector<vecN>& v )
//---------------------------------------------------------------------
{
    int maxDim = INT_MIN;
    int minDim = INT_MAX;

    for ( unsigned int i=0; i<v.size( ); ++i )
    {
        const vecN& vi = v[i];
        if ( maxDim < vi.dim ) maxDim = vi.dim;
        if ( minDim > vi.dim ) minDim = vi.dim;
    }

    return ( (maxDim==minDim) && (maxDim>0) );
}

/*=======================================================================*/
int main( int argc, char* argv[] )
//---------------------------------------------------------------------
{
    bool Do_Random_Insertion = true;

    if ( argc > 1 )
    {
        if ( argv[1][0] == '?' )
        {
            fprintf( stdout, "CSV,size,test,oldTime,old,new,default,old/default,dim. Est,diam. Est,meanSpacing est,var spacing est\n" );
        }
        else if ( argv[1][0] == 'i' || argv[1][0] == 'I'  )
        {
            Do_Random_Insertion = false;
        }

        exit(0);
    }

    const std::vector<vecN> v = ReadGeneralFile( std::cin );

    if ( ! testInputVector( v ) )
    {
        fprintf( stdout, "Inconsistent csv input\n" );
    }
    else
    {
        testGeneral( Do_Random_Insertion, v );
    }

    return 0;
}

