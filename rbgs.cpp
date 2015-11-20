#include<fstream>
#include <cassert>
#include<iostream>
#include "RBGS.h"

#include "Timer.h"

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}

#endif
using namespace std;


int main(int argc, char *argv[]){
    assert(argc>2);
    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int c  = atoi(argv[3]);
    int nrthreads = atoi(argv[3]);
    
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "RBGS" );
#endif

    siwir::Timer timer;
    
    RBGS fast(nx,ny,c,nrthreads);
    double  r = fast.solve();
    
    double time = timer.elapsed();

#ifdef USE_LIKWID
   likwid_markerStopRegion( "RBGS" );
   likwid_markerClose();
#endif

  
   cout<<"time:"<<time<<'\n';
   cout<<"R:"<<r<<'\n';


}
