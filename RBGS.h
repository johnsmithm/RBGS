#include<math.h>
#include<iostream>
#include<fstream>

class RBGS{
    public:
        RBGS(int nx,int ny,int c):nx_(nx),ny_(ny),c_(c),pi(3.141592653589793){
            //distance between grid points    
            hx_ = 2.0/nx;
            hy_ = 1.0/ny;
            //stencil    
            xst =  1.0/(hx_*hx_);
            yst =  1.0/(hy_*hy_);
            mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_)+4*pi*pi;          
            pg = nx*ny;
            //add a padding if nx_ is odd    //not for naive
            even = (nx_%2==0?true:false);
            nx_ = (even==true ? nx_ : nx_ + 1);               
            //# of grid's points    
            pg = nx_*ny_; 

            initialization();
        }

        ~RBGS(){}   
    public:

        double solve(){
            size_t nx__ = nx_ ;
            if(!even)nx__--;
            //bool  odd;
            size_t x,inx_,inx_2,ny_1 = ny_-1;
            --nx__;
	       double mst_1 = 1.0/mst;
	       size_t nx_2 = nx_/2; 
               size_t c__ = c_;      
           for(size_t it=0;it<c__;++it){
                 #pragma omp parallel for private(x) schedule( static )
                for(size_t i=1;i<ny_1;++i){//every black grid point  
                    inx_ = i*nx_;
                    inx_2 = i*(nx_2);
                    for(size_t j=((i&1)?1:2);j<nx__;j=j+2){
                        x = j/2+inx_2;
                        ub[x] = mst_1*(f[j+inx_] + yst*(ur[x+nx_2]+ur[x-nx_2])+
                                                       xst*(ur[x]+ur[((i&1)?x+1:x-1)]) );
                    }

                }           
                 //red
                 #pragma omp parallel for private(x) schedule( static )
                for(size_t i=1;i<ny_1;++i){//every red grid point
                    inx_ = i*nx_;
                    inx_2 = i*(nx_2);
                    for(size_t j=((i&1)?2:1);j<nx__;j=j+2){
                        x = j/2+inx_2;
                        ur[x] =mst_1*(f[j+inx_] + yst*(ub[x+nx_2]+ub[x-nx_2]) +
                                                       xst*(ub[x]+ub[((i&1)?x-1:x+1)]) );

                    }
            }
	      
            }
            
            return  residual_norm();
        }
        void initialization(){
            //we begin with back
            f = new double[pg];    
            ur  = new double[1+pg/2]; 
            ub  = new double[1+pg/2]; 
            double C = 4*pi*pi;   
            double freqx = 2*pi*hx_;   
            double freqy = 2*pi*hy_;     

            //TODO  threads

	        //first touch
            int nx__ = nx_ ;
            if(!even)nx__--;        
            int nx_2 = nx_/2;
            int x;
            --nx__;

           #pragma omp parallel for private(x) schedule( static )
                for(int i=1;i<ny_-1;++i){//every black grid point                
                    for(int j=(i%2==0?2:1);j<nx__;j+=2){
                        x = j/2+i*(nx_2);
                        ub[x] = 0;
                        f[j+i*nx_]=0; 
                        ur[x+nx_2]=0;
                        ur[x-nx_2]=0;
                        ur[x]=0;
                        ur[(i%2==0?x-1:x+1)]=0;
                    }             
            }     
	    
            for(int i=0;i<ny_;++i)
                for(int j=0;j<nx_;++j){
                    f[j+i*nx_]=C*sin(freqx*j)*sinh(freqy*i);//4π^2 sin(2πx) sinh(2πy)
                }  

            for(int i=0;i<pg/2+1;++i)ur[i]=ub[i]=0;
 
            double SINH = sinh(2*pi);  
            double CSINH = C*SINH;   
            //set the nonnull boundary-. check what color is the first point of the last row
            double *ptrfirst = ur, *ptrsecond = ub;
            if(ny_%2==1){
                ptrfirst  = ub;
                ptrsecond = ur;   
            }    
            for(int i=0;i<nx_;i+=2){   
                int g =(ny_-1)*(nx_/2)+i/2;
                ptrfirst[g] = sin(i*freqx) * CSINH;//sin(2πx) sinh(2πy)
                ptrsecond[g] = sin((i+1)*freqx) * CSINH;//sin(2πx) sinh(2πy)
            }
        }
        double residual_norm(){
            double r1=0;
            double f_au;
            //check if we have a padding
            int nx__ = nx_;
            if(!even)--nx__;
            int x;
            int nx_2=nx_/2;
            #pragma omp parallel for private(x,f_au) schedule( static ) reduction(+:r1)
            for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?2:1);j<nx__-1;j+=2){                   

                    x = i*(nx_/2)+j/2;
                    f_au =f[i*nx_+j] - mst * ub[x] +
                        yst*(ur[x+nx_2]+ur[x-nx_2])+
                        xst*(ur[x]+ur[((i)%2==0?x-1:x+1)]) ;
                    r1 = r1  + f_au*f_au;
                }
            }
            double r2=0;
            #pragma omp parallel for private(x,f_au) schedule( static ) reduction(+:r2)
            for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?1:2);j<nx__-1;j+=2){                 

                     x = i*(nx_/2)+j/2;
                    f_au =f[i*nx_+j] - mst * ur[x] +
                         yst*(ub[x+nx_2]+ub[x-nx_2])+
                        xst*(ub[x]+ub[((i)%2==1?x-1:x+1)]);
                    r2 = r2  + f_au*f_au;
                }
            }
            return sqrt(r1+r2);
        }


    public: void print_gnuplot(){
                std::ofstream out("solution.txt");
                int nx__ = nx_;
                if(!even)--nx__;
                for(int i=0;i<ny_;++i){//every black grid point
                    for(int j=(i%2==0?0:1);j<nx__;j+=2){               
                        
                        int x = i*(nx_/2)+j/2;
                        out<<(j*hy_)<<" "<<(i*hx_)<<" "<<ub[x]<<"\n";
                    }
                }
                for(int i=0;i<ny_;++i){//every red grid point
                    for(int j=(i%2==0?1:0);j<nx__;j+=2){                   
                        int x = i*(nx_/2)+j/2;
                        out<<(j*hy_)<<" "<<(i*hx_)<<" "<<ur[x]<<"\n";
                         
                    }
                }
            }

  

    private:
            //we make nx_ even
            bool even;
            //optimize implemantation: red point, black points
            double *ur, * ub;
            //f points,
            double *f;
            //delta x, delta y, left and right points stencil, top and bottom point stencil, middle point stencil
            double hx_,hy_, xst, yst, mst;
            //nr of iteration, nr of grid points x ax, nr of grid points y ax, #points grid, # threads
            int nx_,ny_, c_, pg;
            //pi from maths
            const double pi;
};
