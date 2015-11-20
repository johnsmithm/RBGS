#include<math.h>
#include<iostream>

class RBGS{
    public:
    RBGS(int nx,int ny,int c,int nrthreads):nx_(nx),ny_(ny),c_(c),nrthreads_(nrthreads){
       hx_ = 2.0/nx;
       hy_ = 1.0/ny;
       pg = nx*ny;
       xst =  1.0/(hx_*hx_);
       yst =  1.0/(hy_*hy_);
       mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_)+4*pi*pi;    
    }
    
    ~RBGS(){}   
    public:
    
    double solve(){
    /*TODO use two array for grid points
      ur - red points
      ub - black points
    */
    }
    void initialization(){
    /*TODO find a way to split the points in red and black (ur,ub)*/
    //split and f in fr and fb   -> fr and ur willhave same index
    //how to find neighbors ????    
    }
    double residual_norm(){
    /*TODO compute rezidual norm using ub and ur */
    }
    
    
    double solve_naive(){
       initialization_naive();
       // view(f);   
       //TODO OpenMP.    
       for(int it=0;it<c_;++it){//nr iterations
           
            for(int i=1;i<ny_-1;++i){//every black grid point
                for(int j=(i&1)+1;j<nx_-1;j+=2){
                    /*TODO non-temporal writes*/
                    u[i*nx_+j] =(1.0/mst)*(f[i*nx_+j] + xst*(u[i*nx_+j+1]+u[i*nx_+j-1])+yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]));
                }
            } 
           
           
            for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=((i+1)&1)+1;j<nx_-1;j+=2){
                    /*TODO non-temporal writes*/
                    u[i*nx_+j] =(1.0/mst)*(f[i*nx_+j] + (xst*(u[i*nx_+j+1]+u[i*nx_+j-1])+yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_])));
                }
            }
           if(it%100==0)std::cout<<rezid_norm()<<"\n";
           
       }
        //view(u);
       return  residual_norm_naive();
    }
    private:
    //grid points for right to left, bottom to top
    void initialization_naive(){    
       //using one vector for grid points    
       f = new double[pg];    
       u = new double[pg]; 
        //TODO OpenMP.
       double C = 4*pi*pi;   
       double freqx = 2*pi*hx_;   
       double freqy = 2*pi*hy_;     
        
       //TODO  threads
       for(int i=0;i<ny_;++i)
          for(int j=0;j<nx_;++j){
           f[j+i*nx_]=C*sin(freqx*j)*sinh(freqy*i);//4π^2 sin(2πx) sinh(2πy)
           u[j+i*nx_]=0;
        }       
        
       int last_row = pg-nx_-1;  
       double SINH = sinh(2*pi);  
       double CSINH = C*SINH;   
       for(int i=0;i<nx_;++i){           
           u[last_row+i] = sin(i*freqx) * CSINH;//sin(2πx) sinh(2πy)
       }
    }
    
    void view(double *a){
    for(int i=0;i<ny_;i++){
        for(int j=0;j<nx_;j++)std::cout<<a[i*nx_+j]<<" ";
        std::cout<<"\n";
    }
    }
    
    double residual_norm_naive(){
        double r=0;
        double f_au;
        //TODO OpenMP.
        for(int i=1;i<ny_-1;++i){//every  grid point
                for(int j=1;j<nx_-1;++j){                    
                    
                    f_au = f[j+i*nx_] - mst*u[j+i*nx_] + xst*(u[i*nx_+j+1]+u[i*nx_+j-1])+yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]);
                    
                    r = r  + f_au*f_au;
                }
        }
        return sqrt(r);
    }
    
    void print_gnuplot(){
    //TODO print grid point in gnuplot format
    }
   
    private:
    //optimize implemantation: red point, black points
    double *ur, * ub, * fr, * fb;
    //naive implementation:f points, u points
    double *f, *u;
    //delta x, delta y, left and right points stencil, top and bottom point stencil, middle point stencil
    double hx_,hy_, xst, yst, mst;
    //nr of iteration, nr of grid points x ax, nr of grid points y ax, #points grid, # threads
    int c_,nx_,ny_,pg,nrthreads_;
    //pi from maths
    const double pi=3.141592653589793;
};