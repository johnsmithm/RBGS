#include<fstream>
#include <cassert>
#include<iostream>
#include <time.h>
#include "RBGS.h"


using namespace std;


int main(int argc, char *argv[]){
    assert(argc>2);
    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int c  = atoi(argv[3]);
    //cout<<nx<<" "<<ny<<" "<<c<<endl;

    clock_t t1,t2;
    t1=clock();
    
    RBGS fast(nx,ny,c,1);
    //double  r1 = fast.solve_naive();
    double  r = fast.solve();
    
    t2=clock();
    float diff ((float)t2-(float)t1);
    float time = diff / CLOCKS_PER_SEC;

   fast.print_gnuplot();

  
    cout<<"time:"<<time<<'\n';
    cout<<"R:"<<r<<'\n';


}
