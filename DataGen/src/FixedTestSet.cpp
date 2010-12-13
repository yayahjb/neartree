// FixedTestSet.cpp
//
// Test code for TNear (neartree). The input is expected to contain two parts.
// The first part (order not required) is the "probe" section. The points in
// the probe data will be used to search in the "target" section for their
// nearest neighbor. The probe section will be headed by a line with the
// the text "probe_data" and the target section will have a header line
// that contains "target_data". The text is case sensitive. Various
// diagnostic information is printed out. The command line has one
// possible argument. If it is "i" or "I", then the neartree is constructed
// by immediate insertion, instead of delaying and attempting a more
// evenly balanced tree. If the argument is ?, then the header
// descriptions are output. They can be used as a title line
// in Excel.

#include "TNear.h"

#include <cfloat>
#include <climits>

#include <ctime>

#include "vector_3d.h"
#include <utility>
#include <vector>

#include "Data2CSV.h"

#ifdef _MSC_VER
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#endif

int g_errorCount;


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void testProbes( const bool Do_Random_Insertion, const std::pair<std::vector<vecN>, std::vector<vecN> >& v )
//---------------------------------------------------------------------
{
    const long requiredTests = 1;
    long totalTests = 0;
    const std::vector<vecN>& vProbe = v.first;
    const std::vector<vecN>& vTarget = v.second;

    CNearTree<vecN> nt( vTarget );

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

    /*----------------------------begin balanced tree test--------------------------------------------*/

    nt.SetFlags( CNearTree<vecN>::NTF_NoPrePrune );

    long nodevisits1 = (long)nt.GetNodeVisits( );
    const clock_t tc1 = std::clock();
    vecN newPoint(vTarget[0]);

    totalTests = 0;
    while ( totalTests < requiredTests )
    {
        totalTests += (long)vProbe.size( );
        for ( unsigned int i=0; i<vProbe.size( ); ++i )
        {
            const vecN& probe = vProbe[i];
            nt.NearestNeighbor( DBL_MAX, newPoint, probe );
        }
    }

    {

        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const clock_t tc2 = std::clock();
        long nodevisits2 = (long)nt.GetNodeVisits( );

        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                const vecN& probe = vProbe[i];
                nt.NearestNeighbor( DBL_MAX, newPoint, probe );
            }
        }

        nt.SetFlags( 0 );
        long nodevisits3 = (long)nt.GetNodeVisits( );
        //const clock_t tc3 = std::clock();

        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                const vecN& probe = vProbe[i];
                nt.NearestNeighbor( DBL_MAX, newPoint, probe );
            }
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

        fprintf( stdout, "CSV-balanced,%ld,%ld,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,%.4f,%.4f\n",
            (long)nt.size( ),
            totalTests,
            ((double)(tc2-tc1))/CLOCKS_PER_SEC,
            (double)(nodevisits2-nodevisits1)/(double)totalTests,
            (double)(nodevisits3-nodevisits2)/(double)totalTests,
            ((double)nodevisits4-nodevisits3)/(double)totalTests,
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
        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                vecN newPoint(vTarget[0]);
                nt.LeftNearestNeighbor( DBL_MAX, newPoint, vProbe[i] );
            }
        }

        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const long nodevisits2 = (long)nt.GetNodeVisits( );
        const clock_t tc2 = std::clock();

        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                vecN newPoint(vTarget[0]);
                nt.LeftNearestNeighbor( DBL_MAX, newPoint, vProbe[i] );
            }
        }

        nt.SetFlags( 0 );
        const long nodevisits3 = (long)nt.GetNodeVisits( );
        //const clock_t tc3 = std::clock();

        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                vecN newPoint(vTarget[0]);
                nt.LeftNearestNeighbor( DBL_MAX, newPoint, vProbe[i] );
            }
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

        fprintf( stdout, "CSV-LEFT,%ld,%ld,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,%.4f,%.4f\n",
            (long)nt.size( ),
            totalTests,
            ((double)(tc2-tc1))/CLOCKS_PER_SEC,
            (double)(nodevisits2-nodevisits1)/(double)totalTests,
            (double)(nodevisits3-nodevisits2)/(double)totalTests,
            (double)(nodevisits4-nodevisits3)/(double)totalTests,
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
        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                vecN newPoint(vTarget[0]);
                nt.ShortNearestNeighbor( DBL_MAX, newPoint, vProbe[i] );
            }
        }


        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const long nodevisits2 = (long)nt.GetNodeVisits( );
        const clock_t tc2 = std::clock();

        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                vecN newPoint(vTarget[0]);
                nt.ShortNearestNeighbor( DBL_MAX, newPoint, vProbe[i] );
            }
        }

        nt.SetFlags( 0 );
        const long nodevisits3 = (long)nt.GetNodeVisits( );
        //const clock_t tc3 = std::clock();

        totalTests = 0;
        while ( totalTests < requiredTests )
        {
            totalTests += (long)vProbe.size( );
            for ( unsigned int i=0; i<vProbe.size( ); ++i )
            {
                vecN newPoint(vTarget[0]);
                nt.ShortNearestNeighbor( DBL_MAX, newPoint, vProbe[i] );
            }
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
        fprintf( stdout, "CSV-SHORT,%ld,%ld,%.3f,%.3f,%.3f,%.3f,%.2f,%.4f,%.4f,%.4f,%.4f\n",
            (long)nt.size( ),
            totalTests,
            ((double)(tc2-tc1))/CLOCKS_PER_SEC,
            (double)(nodevisits2-nodevisits1)/(double)totalTests,
            (double)(nodevisits3-nodevisits2)/(double)totalTests,
            (double)(nodevisits4-nodevisits3)/(double)totalTests,
            ((double)(nodevisits2-nodevisits1))/((double)(nodevisits4-nodevisits3)),
            nt.GetDimEstimate(),
            (double)nt.GetDiamEstimate(),
            (double)nt.GetMeanSpacing(),
            (double)nt.GetVarSpacing());
    }
    /*----------------------------end short radius test--------------------------------------------*/

}




/*=======================================================================*/
// read the probe and target sections of a test data file
//
// everything will be ignored until the first of the data indicators ("probe_data"
// and "target_data") is found. Then that section of data is read until the
// other indicator is found. The rest of the data to eof is then read for
// that section.
/*=======================================================================*/
std::pair<std::vector<vecN>,std::vector<vecN> > ReadProbeTarget( std::istream& str )
{
    std::vector<vecN> vProbe;
    std::vector<vecN> vTarget;
    bool isProbe = false;
    bool isTarget = false;

    while ( ! str.fail() && ! str.eof( ) )
    {
        std::string line;
        getline( str, line, '\n' );
        if ( ! str.fail( ) && ! str.eof( ) )
        {
            if ( line.empty( ) )
            {
            }
            else if ( std::string::npos != line.find( std::string("probe_data") ) )
            {
                isProbe  = true;
                isTarget = false;
            }
            else if ( std::string::npos != line.find( std::string("target_data") ) )
            {
                isProbe  = false;
                isTarget = true;
            }
            else
            {
                const std::vector<std::string> vs = GetNextLine( line );
                std::vector<double> vd(vs.size( ));

                for ( unsigned int i=0; i<vs.size( ); ++i )
                {
                    vd[i] = atof(vs[i].c_str( ));
                }

                const vecN v1( vd );
                if ( isProbe )
                { // it's for the probe
                    vProbe.push_back( v1 );
                }
                else if ( isTarget )
                { // it's for the target
                    vTarget.push_back( v1 );
                }
            }
        }
    }
    return( std::make_pair( vProbe, vTarget ) );
}

/*=======================================================================*/
// testInputVector
//
// confirm that the input data vector elements have consistent dimension
/*=======================================================================*/
bool testInputVector( const std::vector<vecN>& v )
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
{
    bool Do_Random_Insertion = true;

    if ( argc > 1 )
    {
        if ( argv[1][0] == '?' )
        {
            fprintf( stdout, "CSV,targets,probes,oldTime,old,new,default,old/default,dim. Est,diam. Est,meanSpacing est,var spacing est\n" );
            exit(0);
        }
        else if ( argv[1][0] == 'i' || argv[1][0] == 'I'  )
        {
            Do_Random_Insertion = false;
        }

    }

    const std::pair<std::vector<vecN>,std::vector<vecN> > v = ReadProbeTarget( std::cin );

    testProbes( Do_Random_Insertion, v );

    return 0;
}


