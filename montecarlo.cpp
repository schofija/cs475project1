#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// print debugging messages?
#ifndef DEBUG
#define DEBUG	false
#endif

// setting the number of threads:
#ifndef NUMT
#define NUMT		    2
#endif

// setting the number of trials in the monte carlo simulation:
#ifndef NUMTRIALS
#define NUMTRIALS	50000
#endif

// how many tries to discover the maximum performance:
#ifndef NUMTIMES
#define NUMTIMES	20
#endif

/* Function Declarations */
inline float Radians(float);
float Ranf(float, float);
int Ranf(int, int);
void TimeOfDaySeed();

/* Random number ranges */
const float txmin = -10.; /* Truck x starting location in feet */
const float txmax = 10.;

const float txvmin = 10.; /* Truck x velocity in ft/sec */
const float txvmax = 30.;

const float tymin = 45.; /* Truck Y location in feet */
const float tymax = 55.;

const float svmin = 10.; /* Snowball overall velocity in feet/sec */
const float svmax = 30.;

const float sthmin = 10.; /* Snowball horizontal launch angle in degrees */
const float sthmax = 90.;

const float halflen = 20.; /* Half-length of truck in feet */

int
main( int argc, char *argv[ ] )
{
#ifndef _OPENMP
        fprintf( stderr, "No OpenMP support!\n" );
        return 1;
#endif

        TimeOfDaySeed( );               // seed the random number generator

        omp_set_num_threads( NUMT );    // set the number of threads to use in parallelizing the for-loop:`

        // better to define these here so that the rand() calls don't get into the thread timing:
        float *txs  = new float [NUMTRIALS];
        float *tys  = new float [NUMTRIALS];
        float *txvs = new float [NUMTRIALS];
        float *svs  = new float [NUMTRIALS];
        float *sths = new float [NUMTRIALS];

        // fill the random-value arrays:
        for( int n = 0; n < NUMTRIALS; n++ )
        {
                txs[n]  = Ranf(  txmin,  txmax );
                tys[n]  = Ranf(  tymin,  tymax );
                txvs[n] = Ranf(  txvmin, txvmax );
                svs[n]  = Ranf(  svmin,  svmax );
                sths[n] = Ranf(  sthmin, sthmax );
        }

        // get ready to record the maximum performance and the probability:
        double  maxPerformance = 0.;    // must be declared outside the NUMTIMES loop
        int     numHits;                // must be declared outside the NUMTIMES loop

        // looking for the maximum performance:
        for( int times = 0; times < NUMTIMES; times++ )
        {
                double time0 = omp_get_wtime( );

                numHits = 0;

                #pragma omp parallel for default(none) shared(txs,tys,txvs,svs,sths,stderr) reduction(+:numHits)
                for( int n = 0; n < NUMTRIALS; n++ )
                {
                        // randomize everything:
                        float tx   = txs[n];
                        float ty   = tys[n];
                        float txv  = txvs[n];
                        float sv   = svs[n];
                        float sthd = sths[n];	/* Snowball theta in degrees */
                        float sthr = Radians(sthd); /* Converted to radians */
                        float svx  = sv * cos(sthr);
                        float svy  = sv * sin(sthr);

                        // how long until the snowball reaches the y depth:
                        float t = ty/svy;

			// how far the truck has moved in x in that amount of time:
                        float truckx = tx + txv * t; /* (starting xpos + xvel*t) */

			// how far the snowball has moved in x in that amount of time:
                        float sbx = svx * t;

			// does the snowball hit the truck (just check x distances, not height):
                        if( fabs(truckx - sbx) < halflen )
                        {
                                numHits++;
                                if( DEBUG )  fprintf( stderr, "Hits the truck at time = %8.3f\n", t );
                        }
                } // for( # of  monte carlo trials )

                double time1 = omp_get_wtime( );
                double megaTrialsPerSecond = (double)NUMTRIALS / ( time1 - time0 ) / 1000000.;
                if( megaTrialsPerSecond > maxPerformance )
                        maxPerformance = megaTrialsPerSecond;

        } // for ( # of timing tries )

        float probability = (float)numHits/(float)( NUMTRIALS );        // just get for last NUMTIMES run
		
		/* Printed to stderr */
		#if 0
        fprintf(stderr, "%2d threads : %8d trials ; probability = %6.2f% ; megatrials/sec = %6.2lf\n",
                NUMT, NUMTRIALS, 100.*probability, maxPerformance);
		#endif
				
		/* This gets printed to the file when we use output redirection */		
		fprintf(stdout,"%d,%d,%f,%f\n",NUMT,NUMTRIALS,probability,maxPerformance);
}

// degrees-to-radians:
inline float Radians( float d )
{
	return (M_PI/180.f) * d;
}

float Ranf(float low, float high)
{
        float r = (float) rand();               // 0 - RAND_MAX
        float t = r  /  (float) RAND_MAX;       // 0. - 1.

        return   low  +  t * ( high - low );
}

int Ranf(int ilow, int ihigh)
{
        float low = (float)ilow;
        float high = ceil( (float)ihigh );

        return (int) Ranf(low,high);
}

// call this if you want to force your program to use
// a different random number sequence every time you run it:
void TimeOfDaySeed( )
{
	struct tm y2k = { 0 };
	y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

	time_t  timer;
	time( &timer );
	double seconds = difftime( timer, mktime(&y2k) );
	unsigned int seed = (unsigned int)( 1000.*seconds );    // milliseconds
	srand( seed );
}