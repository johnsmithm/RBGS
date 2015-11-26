#include<math.h>
#include<iostream>
#include<fstream>
using namespace std;
class RBGS{
    public:
    RBGS(int nx,int ny,int c,int nrthreads):nx_(nx),ny_(ny),c_(c),nrthreads_(nrthreads),pi(3.141592653589793){
       //distance between grid points    
       hx_ = 2.0/nx;
       hy_ = 1.0/ny;
       //stencil    
       xst =  1.0/(hx_*hx_);
       yst =  1.0/(hy_*hy_);
       mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_)+4*pi*pi;          
       pg = nx*ny;
    }
    
    ~RBGS(){}   
    public:
    
    double solve(){
       //add a padding if nx_ is odd    //not for naive
       even = (nx_%2==0?true:false);
       nx_ = (even==true ? nx_ : nx_ + 1);               
       //# of grid's points    
       pg = nx_*ny_; 
        
       initialization();
        //view(ub,ur);   
       //TODO OpenMP.  
        //check if we have padding
        int nx__ = nx_ ;
        if(!even)nx__--;
        
        int nx_2 = nx_/2;
        int pg2 = pg/2;
        int pg2_ = pg2 - nx_2;
        double *utop, *ubottom, *uleft , * uright, *ff;
        
       // std::cout<<nx__<<" "<<nx_<<'\n';
        int x;
        --nx__;
       for(int it=0;it<c_;++it){//nr iterations
          /* //black
            for(int i=1;i<ny_-1;++i){//every black grid point                
                for(int j=(i%2==0?2:1);j<nx__-1;j+=2){
                    
                    ub[j/2+i*(nx_2)] = (1.0/mst)*(f[j+i*nx_] + yst*(ur[(j)/2+(i-1)*(nx_2)]+ur[(j)/2+(1+i)*(nx_2)])+
                                                   xst*(ur[(j-1)/2+i*(nx_2)]+ur[(j+1)/2+i*(nx_2)]) );
                }
             
            }           
           //red
            for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?1:2);j<nx__-1;j+=2){
                    
                    ur[j/2+i*(nx_2)] = (1.0/mst)*(f[j+i*nx_] + yst*(ub[(j)/2+(i-1)*(nx_2)]+ub[(j)/2+(1+i)*(nx_2)]) +
                                                   xst*(ub[(j-1)/2+i*(nx_2)]+ub[(j+1)/2+i*(nx_2)]) );
                                                   
                }
            }
           */
           
            for(int i=1;i<ny_-1;++i){//every black grid point                
                for(int j=(i%2==0?2:1);j<nx__;j+=2){
                    x = j/2+i*(nx_2);
                    ub[x] = (1.0/mst)*(f[j+i*nx_] + yst*(ur[x+nx_2]+ur[x-nx_2])+
                                                   xst*(ur[x]+ur[(i%2==0?x-1:x+1)]) );
                }
             
            }           
           //red
            for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?1:2);j<nx__;j+=2){
                    x = j/2+i*(nx_2);
                    ur[x] = (1.0/mst)*(f[j+i*nx_] + yst*(ub[x+nx_2]+ub[x-nx_2]) +
                                                   xst*(ub[x]+ub[(i%2==1?x-1:x+1)]) );
                                                   
                }
            }
           
           /*
           //black
            double midle;
            for(int i=1;i<ny_-1;++i){//every black grid point                
                for(int j=(i%2==0?2:1);j<nx__-1;j+=4){
                    midle = ur[i*(nx_2)+(j+1)/2];
                    ub[i*(nx_/2)+j/2] = (1.0/mst)*(f[i*nx_+j] + yst*(ur[(i-1)*(nx_2)+(j)/2]+ur[(1+i)*(nx_2)+(j)/2])+
                                                   xst*(ur[i*(nx_2)+(j-1)/2]+midle) );
                    ub[i*(nx_/2)+j/2+1] = (1.0/mst)*(f[i*nx_+j+1] + yst*(ur[(i-1)*(nx_2)+(j)/2+1]+ur[(1+i)*(nx_2)+(j)/2+1])+
                                                   xst*(ur[i*(nx_2)+(j+1)/2+1]+midle) );
                    cout<<(i*(nx_/2)+j/2)<<"|"<<(i*(nx_/2)+j/2+1)<<" ";
                }               
             cout<<endl;
            }           
           //red
           cout<<"red \n";
            for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?1:2);j<nx__-1;j+=4){
                    midle = ub[i*(nx_2)+(j+1)/2];
                    ur[i*(nx_/2)+j/2] = (1.0/mst)*(f[i*nx_+j] + yst*(ub[(i-1)*(nx_2)+(j)/2]+ub[(1+i)*(nx_2)+(j)/2]) +
                                                   xst*(ub[i*(nx_2)+(j-1)/2]+midle) );
                    ur[i*(nx_/2)+j/2+1] = (1.0/mst)*(f[i*nx_+j+1] + yst*(ub[(i-1)*(nx_2)+(j)/2+1]+ub[(1+i)*(nx_2)+(j)/2+1]) +
                                                   xst*(ub[i*(nx_2)+(j-1)/2+1]+midle) );
                         cout<<(i*(nx_/2)+j/2)<<"|"<<(i*(nx_/2)+j/2+1)<<" ";                          
                }
                cout<<endl;
            }
           */
          /*
          for(int i=nx_2;i<pg2_;++i){//black  
             if((i/nx_2)%2==0){//black begin
                if(i%nx_2==0 || (!even && i%nx_2==nx_2-1))continue; //check boundary 
                uleft  = (ur+i-1);
                ff = (f+i+i); 
             }else{
                if(i%nx_==nx_-1)continue; //check boundary 
                uleft  = (ur+i+1);
                ff = (f+i+i+1); 
             }
             utop    = (ur+i+nx_2);
             ubottom = (ur+i-nx_2);
             uright  = (ur+i); 
             ub[i] = (1.0/mst)*(*ff + xst*( *uright+ *uleft )+yst*( *utop+ *ubottom ));  
          }   
           
          for(int i=nx_2;i<pg2_;++i){//red
              if((i/nx_2)%2==0){//black begin
                if(i%nx_2==nx_2-1 )continue;//check boundary  
                uleft  = (ub+i+1);
                ff = (f+i+i+1); 
             }else{//red begin
                if(i%nx_2==0 || (!even && i%nx_2==nx_2-1))continue;//check boundary 
                uleft  = (ub+i-1);
                ff = (f+i+i); 
             }
             utop    = (ub+i+nx_2);
             ubottom = (ub+i-nx_2);
             uright  = (ub+i);               
             ur[i] = (1.0/mst)*(*ff + xst*( *uright+ *uleft )+yst*( *utop+ *ubottom ));  
          }
           
           */
           
           if(it%100==0)std::cout<<residual_norm()<<"\n";
           
       }
         std::cout<<nx_2<<" "<<pg2_<<"\n";
        //view(ub,ur);
       return  residual_norm();
    }
    void initialization(){
    //we begin with back
       f = new double[pg];    
       ur  = new double[1+pg/2]; 
       ub  = new double[1+pg/2]; 
        //TODO OpenMP.
       double C = 4*pi*pi;   
       double freqx = 2*pi*hx_;   
       double freqy = 2*pi*hy_;     
        
       //TODO  threads
          
       for(int i=0;i<ny_;++i)
          for(int j=0;j<nx_;++j){
           f[j+i*nx_]=C*sin(freqx*j)*sinh(freqy*i);//4π^2 sin(2πx) sinh(2πy)
        }  
        
        for(int i=0;i<pg/2+1;++i)ur[i]=ub[i]=0;
        
       int last_row = (nx_/2)*(ny_-2);  
       double SINH = sinh(2*pi);  
       double CSINH = C*SINH;   
        //set the nonnull boundary-. check what color is the first point of the last row
       double *ptrfirst = ur, *ptrsecond = ub;
       if(ny_&1==1){
         ptrfirst  = ub;
         ptrsecond = ur;   
       }    
       // if(ptrfirst==ur)std::cout<<"red first \n";
       for(int i=0;i<nx_;i+=2){   
          // std::cout<<((ny_-1)*(nx_/2)+i/2)<<'\n';
           int g =(ny_-1)*(nx_/2)+i/2;
           //std::cout<<g<<'\n';
           ptrfirst[g] = sin(i*freqx) * CSINH;//sin(2πx) sinh(2πy)
           ptrsecond[g] = sin((i+1)*freqx) * CSINH;//sin(2πx) sinh(2πy)
       }
    }
    double residual_norm(){
        double r=0;
        double f_au;
        //TODO OpenMP.
        //check if we have a padding
        int nx__ = nx_;
        if(!even)--nx__;
        for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?2:1);j<nx__-1;j+=2){                   
                    
                    int x = i*(nx_/2)+j/2;
                    f_au =f[i*nx_+j] - mst * ub[x] +
                                                   //xst*(ur[x]+ur[x+1]) +
                                                   //yst*(ur[x+nx_/2]+ur[x-nx_/2]) ;
                                                    xst*(ur[i*(nx_/2)+(j+1)/2]+ur[i*(nx_/2)+(j-1)/2]) +
                                                    yst*(ur[(i-1)*((nx_)/2)+(j)/2]+ur[(1+i)*((nx_)/2)+(j)/2]) ;
                    r = r  + f_au*f_au;
                }
        }
        for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?1:2);j<nx__-1;j+=2){                   
                    
                    //f_au = f[j+i*nx_] - mst*u[j+i*nx_] + xst*(u[i*nx_+j+1]+u[i*nx_+j-1])+yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]);
                    int x = i*(nx_/2)+j/2;
                    f_au =f[i*nx_+j] - mst * ur[x] +
                                                   //xst*(ub[x]+ub[x+1]) +
                                                   //yst*(ub[x+nx_/2]+ub[x-nx_/2]) ;
                                                    xst*(ub[i*(nx_/2)+(j+1)/2]+ub[i*(nx_/2)+(j-1)/2]) +
                                                   yst*(ub[(i-1)*((nx_)/2)+(j)/2]+ub[(1+i)*((nx_)/2)+(j)/2]) ;
                    r = r  + f_au*f_au;
                }
        }
        return sqrt(r);
    }
    
    
    double solve_naive(){
       initialization_naive();
       // view(u);   
       //TODO OpenMP.    
       // std::cout<<nx_<<'\n';
       for(int it=0;it<c_;++it){//nr iterations
           
            for(int i=1;i<ny_-1;++i){//every black grid point
                for(int j=((1+i)&1)+1;j<nx_-1;j+=2){
                    /*TODO non-temporal writes*/
                    u[i*nx_+j] =(1.0/mst)*(f[i*nx_+j] + xst*(u[i*nx_+j+1]+u[i*nx_+j-1])+yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]));
                }
            } 
           
           
            for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=((i)&1)+1;j<nx_-1;j+=2){
                    /*TODO non-temporal writes*/
                    u[i*nx_+j] =(1.0/mst)*(f[i*nx_+j] + (xst*(u[i*nx_+j+1]+u[i*nx_+j-1])+yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_])));
                }
            }
          // if(it%100==0)std::cout<<residual_norm_naive()<<"\n";
           //view(u);
           //break;
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
        
       int last_row = (ny_-1)*nx_;  
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
    void view(double *a,double *b){
        for(int i=0;i<ny_;i++){
          //  std::cout<<(i*(nx_/2))<<'\n';
            for(int j=0;j<nx_;j++)
                if((i%2==0 && j%2==0)  || (i%2==1 && j%2==1))std::cout<<a[i*(nx_/2)+j/2]/*<<"b"<<(i*(nx_/2)+j/2)*/<<" ";
                else std::cout<<b[i*(nx_/2)+j/2]/*<<"r"<<(i*(nx_/2)+j/2)*/<<" ";
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
    
    public: void print_gnuplot(){
    //TODO print grid point in gnuplot format
    //do not forget about even nx_   
    //include boundaries!!!    
        std::ofstream out("solution.txt");
        int nx__ = nx_;
        if(!even)--nx__;
        for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?2:1);j<nx__-1;j+=2){                   
                    
                    int x = i*(nx_/2)+j/2;
                    out<<(j*hy_)<<" "<<(i*hx_)<<" "<<ub[x]<<"\n";
                }
        }
        for(int i=1;i<ny_-1;++i){//every red grid point
                for(int j=(i%2==0?1:2);j<nx__-1;j+=2){                   
                    int x = i*(nx_/2)+j/2;
                    out<<(j*hy_)<<" "<<(i*hx_)<<" "<<ur[x]<<"\n";
                }
        }
    }
    
    public: void print_gnuplot_naive(){
    //TODO print grid point in gnuplot format
    //do not forget about even nx_   
    //include boundaries!!!    
        std::ofstream out("solution.txt");
        
        for(int i=0;i<ny_;++i){//every red grid point
                for(int j=0;j<nx_;j+=1){                   
                    
                    out<<(j*hy_)<<" "<<(i*hx_)<<" "<<u[i*nx_+j]<<"\n";
                }
        }       
    }
   
    private:
    //we make nx_ even
    bool even;
    //optimize implemantation: red point, black points
    double *ur, * ub, * fr, * fb;
    //naive implementation:f points, u points
    double *f, *u;
    //delta x, delta y, left and right points stencil, top and bottom point stencil, middle point stencil
    double hx_,hy_, xst, yst, mst;
    //nr of iteration, nr of grid points x ax, nr of grid points y ax, #points grid, # threads
    int c_,nx_,ny_,pg,nrthreads_;
    //pi from maths
    const double pi;
};