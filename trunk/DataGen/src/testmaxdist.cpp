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

#include "tnear.h"
#include "rhrand.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

RHrand rhr;
void testBigVector(  );
void testDouble(  );
int g_errorCount;

/*=======================================================================*/
/* make a 17-dimension vector class for testing */
class vec17
{
public:
    double pd[17];
    int dim;
    double length;  // just for a signature for debugging
    static RHrand srhr;

    vec17( ) :
    dim(17), length(0)
    {
        for( int i=0; i<dim; ++i )
        {
            pd[i] = srhr.urand()*((double)RHrand::RHRAND_MAX);
        }
        length = Norm( );
    }
    explicit vec17( const double d ):
    dim(17), length(0)
    {
        for( int i=0; i<dim; ++i )
        {
            pd[i] = d;
        }
        length = Norm( );
    }
    explicit vec17( const int n ):
    dim(17), length(0)
    {
        for( int i=0; i<dim; ++i )
        {
            pd[i] = double(n);
        }
        length = Norm( );
    }
    ~vec17( void )
    {
    }

    vec17 operator-( const vec17& v ) const /* USERS: be sure to make both const */
    {
        vec17 vtemp;
        for( int i=0; i<dim; ++i )
        {
            vtemp.pd[i] = this->pd[i] - v.pd[i];
        }
        vtemp.length = Norm(  );
        return( vtemp );
    }
    double Norm( void ) const
    {
        double dtemp = 0.0;
        for( int i=0; i<dim; ++i )
        {
            dtemp += pd[i]*pd[i];  //  L2 measure here
        }
        return( double(sqrt( dtemp )) );
    }

    vec17& operator-= ( const vec17 ) { return ( *this ); }; // just to keep LINT happy
};  // end vec17

RHrand vec17::srhr;

/*=======================================================================*/
// HAM_Gen
//
// Create nPoints vectors of dimension dim with that obey the
// Hammersley quasirandom distribution.
/*=======================================================================*/
std::vector<double> HAM_Gen( const long dim, const long nPoints )
//---------------------------------------------------------------------
{
    const int step = 1;
    std::vector<int> seed( dim, 0 );
    std::vector<int> leap( dim, 1 );
    std::vector<int> base;
    base.push_back( -1000 );
    for ( int i=1; i<dim; ++i )
    {
        base.push_back( prime(i) );
    }
    std::vector<double> r(nPoints*dim);
    i4_to_hammersley_sequence ( dim, nPoints, step, &seed[0], &leap[0], &base[0], &r[0] );
    return( r );
}


/*=======================================================================*/
// testDouble
//
// Generate a comb of double values
/*=======================================================================*/
void testDouble( )
//---------------------------------------------------------------------
{

    CNearTree<double> nt;

    const int maxData = 50000;
    //for ( int i=0; i<maxData; ++i )
    //{
    //    nt.insert( rand( ) );
    //}

    for ( int i=0; i<maxData; ++i )
    {
        nt.insert( double(i) );
    }

    nt.CompleteDelayedInsertRandom( );
    fprintf( stdout, "testDouble tree size: %ld, dimension estimate: %ld, diameter est,: %g, mean spacing est.: %g, variance spacing est.: %g\n",
        (long)nt.size(), (long)nt.GetDimEstimate(), (double)nt.GetDiamEstimate(),
        (double)nt.GetMeanSpacing(), (double)nt.GetVarSpacing());

    clock_t tc1 = std::clock();
    const int nTests = 10000;
    nt.SetFlags( CNearTree<double>::NTF_NoPrePrune );

    RHrand rhr;
    for ( int i=0; i<nTests; ++i )
    {
        double newPoint;
        const double probe = rhr.urand( ) * double(maxData);
        nt.NearestNeighbor( DBL_MAX, newPoint, probe );
    }

    clock_t tc2 = std::clock();
    nt.SetFlags( CNearTree<double>::NTF_ForcePrePrune );

    for ( int i=0; i<nTests; ++i )
    {
        double newPoint;
        const double probe = rhr.urand( ) * double(maxData);
        nt.NearestNeighbor( DBL_MAX, newPoint, probe );
    }

    clock_t tc3 = std::clock();
    nt.SetFlags( 0 );

    for ( int i=0; i<nTests; ++i )
    {
        double newPoint;
        const double probe = rhr.urand( ) * double(maxData);
        nt.NearestNeighbor( DBL_MAX, newPoint, probe );
    }

    clock_t tc4 = std::clock();

    fprintf( stdout, "testDouble treesize:%ld  tests:%d old: %f, new: %f, default: %f, ratio old:new %f\n", (long)nt.size( ), nTests,
        ((double)(tc2-tc1))/CLOCKS_PER_SEC,
        ((double)(tc3-tc2))/CLOCKS_PER_SEC,
        ((double)(tc4-tc3))/CLOCKS_PER_SEC,
        ((double)(tc2-tc1))/((double)(tc3-tc2)) );

}

/*=======================================================================*/
// testBigVector
//
// generate a lot of 17-D vectors for testing
/*=======================================================================*/
void testBigVector(  )
//---------------------------------------------------------------------
{
    const int vectorsize = 100000;

    CNearTree<vec17> nt;
    std::vector<vec17>vAll(vectorsize); /* keep a list of all of the input so we can find particular entries */
    double rmax = -DBL_MAX;
    double rmin =  DBL_MAX;
    vec17 v17min; /* to be the point nearest to the origin */
    vec17 v17max; /* to be the point farthest from the origin */

    /* All of the coordinate values will be in the range 0-RHrand::RHRAND_MAX. In other words,
    all of the data points will be within a 17-D cube that has a corner at
    the origin of the space.
    */
    for( int i=0; i<vectorsize; ++i )
    {
        vec17 v;
        if( v.Norm( ) < rmin )
        {
            rmin = v.Norm( );
            v17min = v;
        }

        if( v.Norm( ) > rmax )
        {
            rmax = v.Norm( );
            v17max = v;
        }

        nt.insert( v );
        vAll[i] = v;
    }

    nt.CompleteDelayedInsertRandom( );
    fprintf( stdout, "testBigVector tree size: %ld, dimension estimate: %ld, diameter est,: %g, mean spacing est.: %g, variance spacing est.: %g\n",
        (long)nt.size(), (long)nt.GetDimEstimate(), (double)nt.GetDiamEstimate(),
        (double)nt.GetMeanSpacing(), (double)nt.GetVarSpacing());

    clock_t tc1 = std::clock();
    const int nTests = 100;
    {
        nt.SetFlags( CNearTree<vec17>::NTF_NoPrePrune );
        for ( int i=0; i<nTests; ++i )
        {
            vec17 newPoint;
            vec17 returnPoint;
            //nt.NearestNeighbor( double(RHrand::RHRAND_MAX/2)*sqrt(17.), returnPoint, newPoint );
            nt.NearestNeighbor( double(RHrand::RHRAND_MAX/2)*sqrt(17.), newPoint );
        }
    }
    clock_t tc2 = std::clock();
    {
        nt.SetFlags( CNearTree<vec17>::NTF_ForcePrePrune );
        for ( int i=0; i<nTests; ++i )
        {
            vec17 newPoint;
            vec17 returnPoint;
            //nt.NearestNeighbor( double(RHrand::RHRAND_MAX/2)*sqrt(17.), returnPoint, newPoint );
            nt.NearestNeighbor( double(RHrand::RHRAND_MAX/2)*sqrt(17.), newPoint );
        }
    }
    clock_t tc3 = std::clock();
    {
        nt.SetFlags( 0 );
        for ( int i=0; i<nTests; ++i )
        {
            vec17 newPoint;
            vec17 returnPoint;
            //nt.NearestNeighbor( double(RHrand::RHRAND_MAX/2)*sqrt(17.), returnPoint, newPoint );
            nt.NearestNeighbor( double(RHrand::RHRAND_MAX/2)*sqrt(17.), newPoint );
        }
    }
    clock_t tc4 = std::clock();
    fprintf( stdout, "testBigVector treesize:%ld  tests:%d  old: %f, new: %f, default: %f, ratio old:new %f\n", (long)nt.size( ), nTests,
        ((double)(tc2-tc1))/CLOCKS_PER_SEC,
        ((double)(tc3-tc2))/CLOCKS_PER_SEC,
        ((double)(tc4-tc3))/CLOCKS_PER_SEC,
        ((double)(tc2-tc1))/((double)(tc3-tc2)));
}  //

/*=======================================================================*/
// testHamm
//
// perform tests using the Hammersley quasirandom distribution
/*=======================================================================*/
void testHamm( void )
//---------------------------------------------------------------------
{
    std::vector<double> vHAM = HAM_Gen( 3, 100000 );
    std::vector<Vector_3> pork;

    CNearTree<Vector_3> nt;

    for ( unsigned int i=0; i<vHAM.size( ); i+=3 )
    {
        nt.insert( Vector_3( vHAM[i], vHAM[i+1], vHAM[i+2] ) );
    }

    nt.CompleteDelayedInsertRandom( );
    fprintf( stdout, "testHamm tree size: %ld, dimension estimate: %ld, diameter est,: %g, mean spacing est.: %g, variance spacing est.: %g\n",
        (long)nt.size(), (long)nt.GetDimEstimate(), (double)nt.GetDiamEstimate(),
        (double)nt.GetMeanSpacing(), (double)nt.GetVarSpacing());

    clock_t tc1 = std::clock();
    const int nTests = 10000;
    {
        nt.SetFlags( CNearTree<Vector_3>::NTF_NoPrePrune );
        int iseed = 19191;
        for ( int i=0; i<nTests; ++i )
        {
            Vector_3 returnPoint;
            nt.NearestNeighbor( DBL_MAX, randVector(iseed) );
        }
    }
    clock_t tc2 = std::clock();
    {
        nt.SetFlags( CNearTree<Vector_3>::NTF_ForcePrePrune );
        int iseed = 19191;
        for ( int i=0; i<nTests; ++i )
        {
            Vector_3 returnPoint;
            nt.NearestNeighbor( DBL_MAX, randVector(iseed) );
        }
    }
    clock_t tc3 = std::clock();
    {
        nt.SetFlags( 0 );
        int iseed = 19191;
        for ( int i=0; i<nTests; ++i )
        {
            Vector_3 returnPoint;
            nt.NearestNeighbor( DBL_MAX, randVector(iseed) );
        }
    }
    clock_t tc4 = std::clock();
    fprintf( stdout, "testHamm treesize:%ld  tests:%d  old: %f, new: %f, default: %f, ratio old:new %f\n",
        (long)nt.size( ), nTests,
        ((double)(tc2-tc1))/CLOCKS_PER_SEC,
        ((double)(tc3-tc2))/CLOCKS_PER_SEC,
        ((double)(tc4-tc3))/CLOCKS_PER_SEC,
        ((double)(tc2-tc1))/((double)(tc3-tc2)) );
}  //

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

    int iseed  = 19191;
    long nodevisits1 = nt.GetNodeVisits( );
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
        long nodevisits2 = nt.GetNodeVisits( );

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.NearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( 0 );
        long nodevisits3 = nt.GetNodeVisits( );
        const clock_t tc3 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe(v[0].dim);
            nt.NearestNeighbor( DBL_MAX, newPoint, probe );
        }

        long nodevisits4 = nt.GetNodeVisits( );
        const clock_t tc4 = std::clock();

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


        const long nodevisits1 = nt.GetNodeVisits( );
        const clock_t tc1 = std::clock();
        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.LeftNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const long nodevisits2 = nt.GetNodeVisits( );
        const clock_t tc2 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.LeftNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( 0 );
        const long nodevisits3 = nt.GetNodeVisits( );
        const clock_t tc3 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe(v[0].dim);
            nt.LeftNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        const long nodevisits4 = nt.GetNodeVisits( );
        const clock_t tc4 = std::clock();

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
        const long nodevisits1 = nt.GetNodeVisits( );
        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.ShortNearestNeighbor( DBL_MAX, newPoint, probe );
        }


        nt.SetFlags( CNearTree<vecN>::NTF_ForcePrePrune );
        const long nodevisits2 = nt.GetNodeVisits( );
        const clock_t tc2 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe = vecN(v[0].dim);
            nt.ShortNearestNeighbor( DBL_MAX, newPoint, probe );
        }

        nt.SetFlags( 0 );
        const long nodevisits3 = nt.GetNodeVisits( );
        const clock_t tc3 = std::clock();

        for ( int i=0; i<nTests; ++i )
        {
            vecN newPoint(v[0].dim);
            const vecN probe(v[0].dim);
            nt.ShortNearestNeighbor( DBL_MAX, newPoint, probe );
        }
        const long nodevisits4 = nt.GetNodeVisits( );

        const clock_t tc4 = std::clock();

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

    const bool bGeneralTest = true;

    const std::vector<vecN> v = ReadGeneralFile( std::cin );

    if ( ! testInputVector( v ) )
    {
        fprintf( stdout, "Inconsistent csv input\n" );
    }
    else  if ( bGeneralTest )
    {
        testGeneral( Do_Random_Insertion, v );
    }
    else
    {
        testDouble( );
        //testBigVector( );
        //testHamm( );
    }



    return 0;
}

