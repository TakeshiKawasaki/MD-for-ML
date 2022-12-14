#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <omp.h>
#include "MT.h"

#define Np 4000 //# of the particles
#define rho 1.2 //# density
#define Nn 1000  //# of the neigbour lists
#define L pow(Np/rho,1./3.)
#define teq 10 //equilibration time
//#define tmax 1000 //production run time
#define dtbdhs 0.1
#define dtmd 0.001 //dt for molecular dynamics
#define dtbd 0.005 //dt for brownian dynamics
#define temp 1.0 // temperature
#define dim 3 //spatial dimension
#define cut 2.5 //potential cut off
#define skin 1.0// skin size for list update
#define d_gamma 0.001

const int nthread=10;

double unif_rand(double left, double right)
{
  return left + (right - left)*genrand_real1();
}

double gaussian_rand(void)
{
  static double iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if (iset == 0) {
    do {
      v1 = unif_rand(-1, 1);
      v2 = unif_rand(-1, 1);
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);

    gset = v1*fac;
    iset = 0.50;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}

void ini_coord_rand(double (*x)[dim]){
  for(int i = 0 ; i < Np; i++){
    x[i][0] = L*unif_rand(0.,1.);
    x[i][1] = L*unif_rand(0.,1.);
    x[i][2] = L*unif_rand(0.,1.);
  }
}


void set_diameter(int *a){
  for(int i=0;i<0.2*Np;i++){
    a[2*i]=2;
    a[2*i+1]=1;
  }
  for(int i=0.4*Np;i<Np;i++)
    a[i]=1;
}

void p_boundary(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]-=L*floor(x[i][j]/L);
}


void scale_sliding_blick(double (*x)[dim],double gamma)
{
  for(int i=0;i<Np;i++){
    if(x[i][0]<gamma*(x[i][1]))
      x[i][0]+=L;
    if(x[i][0]>L+gamma*(x[i][1]))
      x[i][0]-=L;

    if(x[i][1]<0.0){
      x[i][1]+=L;
      x[i][0]+=gamma*L;
    }
    if(x[i][1]>L){
      x[i][1]-=L;
      x[i][0]-=gamma*L;
    }
    if(x[i][2]<0.0)
      x[i][2]+=L;
    if(x[i][2]>L)
      x[i][2]-=L;
  }
}


void ini_matrix_dim(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
}

void ini_matrix_thread(double (*x)[nthread]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<nthread;j++)
      x[i][j]=0.0;
}

void ini_array(double *x){
  for(int i=0;i<Np;i++)
      x[i]=0.0;
}


int f(int i,int M)
{
  int k;
  k=i;
  if(k<0)
    k+=M;
  if(k>=M)
    k-=M;

  return k;
}

void cell_list(int (*list)[Nn],double (*x)[dim],int M,double gamma)
{
  int i,j,k;
  int nx,ny,nz;
  int l,m,n;
  double dx,dy,dy_temp,dz,r2;
  double thresh=cut+skin;
  
  int (*map)[Np+1]=new int[M*M*M][Np+1];


  for(i=0;i<M;i++)
    for(j=0;j<M;j++)
      for(k=0;k<M;k++)
	map[i+M*j+M*M*k][0]=0;

  //#pragma omp parallel for private(i,nx,ny,nz,n,m,l)
  for(i=0;i<Np;i++){
    nx=f((int)((x[i][0]-gamma*x[i][1])*M/L),M);
    ny=f((int)(x[i][1]*M/L),M);
    nz=f((int)(x[i][2]*M/L),M);
    for(n=nz-1;n<=nz+1;n++)
      for(m=ny-1;m<=ny+1;m++)
	for(l=nx-1;l<=nx+1;l++){
	  map[f(l,M)+M*f(m,M)+M*M*f(n,M)][map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0] +1]=i;
	  map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0]++;
	}  
  }

  #pragma omp parallel for private(i,j,k,nx,ny,nz,dx,dy,dy_temp,dz,r2)
  for(i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)((x[i][0]-gamma*x[i][1])*M/L),M);
    ny = f((int)(x[i][1]*M/L),M);
    nz = f((int)(x[i][2]*M/L),M);
    // std::cout<<nx<<std::endl;
    // printf("map=%d,%d\n",nx+M*ny+M*M*nz,map[nx+M*ny+M*M*nz][0]);

    for (k=1; k<=(map[nx+M*ny+M*M*nz][0]); k++){
      j = map[nx+M*ny+M*M*nz][k];
     
      if(j != i){
	dx =x[i][0] - x[j][0];
	dy =x[i][1] - x[j][1];
	dz =x[i][2] - x[j][2];
	dy_temp=dy;
	dy -= L*floor((dy+0.5*L)/L);
	dx -= gamma*L*floor((dy_temp+0.5*L)/L);
	dx -= L*floor((dx+0.5*L)/L);
	dz -= L*floor((dz+0.5*L)/L);

	r2 = dx*dx + dy*dy + dz*dz;

	if(r2<thresh*thresh){
	  list[i][0]++;
	  list[i][list[i][0]]=j;
	}
      }
    }
  }
  delete []map;
}

void calc_force_hs(double (*x)[dim],double (*f)[dim],int *a,double *U,int (*list)[Nn]){
  double dx,dy,dz,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,eij,dUrcut,Ucut,dr,t;
  ini_matrix_dim(f);
  *U=0;
  for(int i=0;i<Np;i++)
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dz=x[i][2]-x[list[i][j]][2];
     
      dx-=L*floor((dx+0.5*L)/L);
      dy-=L*floor((dy+0.5*L)/L);
      dz-=L*floor((dz+0.5*L)/L);
      dr2=dx*dx+dy*dy+dz*dz;
      if(a[i]+a[list[i][j]] == 2){
	aij=1.0;
	eij=1.0;
  //      1.0 0.8 0.88
  //	  1.0 1.5 0.5
      }
      if(a[i]+a[list[i][j]] == 3){
	aij=0.8;
	eij=1.5;
      }
      if(a[i]+a[list[i][j]] == 4){
	aij=0.88;
	eij=0.5;
      }

      if(dr2<aij*aij){
	//	printf("%f\n",t);
	t=sqrt(dr2/aij*aij);
	dr=sqrt(dr2);
	dUr=-(1.-t)/aij;
	f[i][0]-=dUr*dx/dr;
	//	f[list[i][j]][0]+=dUr*dx/dr;
	f[i][1]-=dUr*dy/dr;
	//	f[list[i][j]][1]+=dUr*dy/dr;
	f[i][2]-=dUr*dz/dr;
	//	f[list[i][j]][2]+=dUr*dz/dr;
      }
    }
}

void calc_force(double (*x)[dim],double (*f)[dim],int *a,double *U,double *rfxy,int (*list)[Nn],double gamma,double *stress){
  int i,j,n;
  double dx,dy,dy_temp,dz,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,eij,dUrcut,Ucut,dr;
  double V = L*L*L;
 
  ini_matrix_dim(f);
  ini_array(stress);
  // ini_matrix_thread(fx);
  // ini_matrix_thread(fy);
  // ini_matrix_thread(fz);
  // ini_matrix_thread(stress_thread);
  
  double U0=0,rfxy0=0.0;
  // omp_set_num_threads(nthread);

#pragma omp parallel for private(j,dx,dy,dy_temp,dz,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,eij,dUrcut,Ucut,dr) reduction(+:U0,rfxy0) 
  for(i=0;i<Np;i++)
    for(j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dz=x[i][2]-x[list[i][j]][2];
      dy_temp=dy;
      dy-=L*floor((dy+0.5*L)/L);
      dx-=gamma*L*floor((dy_temp+0.5*L)/L);
      dx-=L*floor((dx+0.5*L)/L);
      dz-=L*floor((dz+0.5*L)/L);
      
      dr2=dx*dx+dy*dy+dz*dz;
      
      if(a[i]+a[list[i][j]] == 2){
	aij=1.0;
	eij=1.0;
	//  1.0 0.8 0.88
	//  1.0 0.5 1.5
      }
      if(a[i]+a[list[i][j]] == 3){
	aij=0.8;
	eij=1.5;
      }
      if(a[i]+a[list[i][j]] == 4){
	aij=0.88;
	eij=0.5;
      }
      
      if(dr2<cut*cut*aij*aij){
	dr=sqrt(dr2);
	w2=aij*aij/dr2;
	w6=w2*w2*w2;
	w12=w6*w6;
	
	w2cut=1./cut/cut;
	w6cut=w2cut*w2cut*w2cut;
	w12cut=w6cut*w6cut;
	dUrcut=-48.*eij*w12cut/(cut*aij)+24.*eij*w6cut/(cut*aij);
	Ucut=4.*eij*w12cut-4.*eij*w6cut;
	
       	dUr=(-48.*eij*w12+24.*eij*w6)/dr2-dUrcut/dr;
	//	dUr=(-48.*eij*w12+24*eij*w6)/dr2;
	//      n=omp_get_thread_num();
	f[i][0] -=dUr*dx;
	//	fx[list[i][j]][n] +=dUr*dx;
	f[i][1] -=dUr*dy;
	//	fy[list[i][j]][n] +=dUr*dy;
	f[i][2] -=dUr*dz;
	//	fz[list[i][j]][n] +=dUr*dz;
      
 	U0 +=4.*eij*(w12-w6)-Ucut-dUrcut*(dr-cut*aij);
	//*U +=4.*eij*(w12-w6)-Ucut;
	rfxy0     += 0.5*dUr*dx*dy/V;
	stress[i] += 0.5*dUr*dx*dy/V;
	//stress_thread[list[i][j]][n] += 0.5*dUr*dx*dy/V;
	//std::cout<<"f_ij="<<dUr*dx<<std::endl;
      }
    }

  *U = 0.5*U0;
  *rfxy = rfxy0;
}

void eom_langevin_hs(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,double temp0,int (*list)[Nn],double *kine){
  int i,j;
  double zeta=1.0;
  double fluc=sqrt(2.*zeta*temp0*dt);
  *kine=0.0;
  calc_force_hs(x,f,a,&(*U),list);

#pragma omp parallel for private(j)
  for(i=0;i<Np;i++)
    for(j=0;j<dim;j++){
      v[i][j]+=-zeta*v[i][j]*dt+f[i][j]*dt+fluc*gaussian_rand();
      x[i][j]+=v[i][j]*dt;
      *kine +=v[i][j]*v[i][j]/3.0/Np;
    }
  p_boundary(x);
}


void eom_langevin(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,double temp0,int (*list)[Nn],double *kine,double *stress){
  int i,j;
  double zeta=1.0;
  double fluc=sqrt(2.*zeta*temp0*dt);
  double dummy;
  double kine0=0.0;
  calc_force(x,f,a,U,&dummy,list,0.0,stress);
#pragma omp parallel for private(j) reduction(+:kine0)
  for(i=0;i<Np;i++)
    for(j=0;j<dim;j++){
      v[i][j]+=-zeta*v[i][j]*dt+f[i][j]*dt+fluc*gaussian_rand();
      x[i][j]+=v[i][j]*dt;
      kine0 +=v[i][j]*v[i][j]/3.0/Np;
    }
  p_boundary(x);
  *kine=kine0;
}

void steepest_descent(double (*x)[dim],double (*f)[dim],double gamma,int (*list)[Nn],int *a,int M,double *U,double *stress){
  double dx,dy,dy_temp,dz,dt0=0.0001,zeta=0.0,sum_force =0.0,dr2max=0.0,dummy=0.0;
  double x0[Np][dim],v[Np][dim];
  int i,j;
  cell_list(list,x,M,gamma);
  int flag=0;
  
#pragma omp parallel for private(j)
  for(i=0;i<Np;i++)
    for(j=0;j<dim;j++){
      x0[i][j] = x[i][j];
      v[i][j] = 0.0;
    }

  for(;;){
    calc_force(x,f,a,U,&dummy,list,gamma,stress);
    sum_force=0.0;
    
#pragma omp parallel for private(i,j,dx,dy,dy_temp,dz) reduction(+:sum_force,flag)
    for(i=0;i<Np;i++){
      for(j=0;j<dim;j++){
	v[i][j] = f[i][j];
	x[i][j] += v[i][j]*dt0;
      }
      sum_force += sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]+f[i][2]*f[i][2])/Np;
      dy = x[i][1]-x0[i][1];
      dx = x[i][0]-x0[i][0];
      dz = x[i][2]-x0[i][2];
      dy_temp=dy;
      dy -= L*floor((dy+0.5*L)/L);
      dx -= gamma*L*floor((dy_temp+0.5*L)/L);
      dx -= L*floor((dx+0.5*L)/L);
      dz -= L*floor((dz+0.5*L)/L);
      
      if(dx*dx+dy*dy+dz*dz > skin*skin*0.25){   
	flag += 1;
	dr2max = dx*dx+dy*dy+dz*dz;
      }  
    }
    
    // printf("SD: sum_force=%.16f,x=%f\n",sum_force,x[0][0]); 
    
    if(flag != 0){
      printf("cell update: SD %f\n",dr2max);
      cell_list(list,x,M,gamma);
      flag=0;
      #pragma omp parallel for private(j)    
      for(int i=0;i<Np;i++)
        for(int j=0;j<dim;j++)
          x0[i][j]=x[i][j];
    }
    if(sum_force<1.e-1)
      break;
  }
}

void copy_array(double (*x)[dim],double (*x0)[dim]){
  int i,j;
#pragma omp parallel for private(j)
  for(i=0;i<Np;i++)
    for(j=0;j<dim;j++)
      x0[i][j]=x[i][j];
}

void norm_array(double *d,double (*x)[dim]){
  int i,j; 
  double d0=0.0;
  #pragma omp parallel for private(j) reduction(+:d0)
  for(i=0;i<Np;i++)
    for(j=0;j<dim;j++)
    d0 += x[i][j]*x[i][j];
  *d = sqrt(d0);
}

int FIRE(double (*x)[dim],double (*f)[dim],double gamma,int (*list)[Nn],int *a,int M,double *U,double *rfxy,double *stress)
{
  int i,j,imax;
  double alpha = 0.1,P,max;
  double dt0=0.0001;
  double v[Np][dim],x0[Np][dim],dx,dy,dy_temp,dz,v_nor,f_nor,dr2max;
  int count=0,count0=0;
  int conv=0,flag=0;
  double sum_force=0.0;
  double finc=1.1,falpha=0.99,fdec=0.5;
  
  cell_list(list,x,M,gamma); 
  
  #pragma omp parallel for private(j)
  for(i=0;i<Np;i++){
       for(j=0;j<dim;j++){
      v[i][j]  = 0.0;
      x0[i][j] = x[i][j];
    }
  }
  P=0.0;
  
  calc_force(x,f,a,U,rfxy,list,gamma,stress);
   
  for(;;){
    P=0.0;
    sum_force=0.0;
    count++;
    count0++;
    #pragma omp parallel for private(j)
    for(i=0;i<Np;i++)
      for(j=0;j<dim;j++){
	x[i][j]+=v[i][j]*dt0+0.5*f[i][j]*dt0*dt0;
	v[i][j]+=0.5*f[i][j]*dt0;
      }
    
    calc_force(x,f,a,U,rfxy,list,gamma,stress);
    
    #pragma omp parallel for private(j)
    for(i=0;i<Np;i++)
      for(j=0;j<dim;j++)
	v[i][j] += 0.5*f[i][j]*dt0;
    
    norm_array(&f_nor,f);
    norm_array(&v_nor,v);     
    
    
#pragma omp parallel for private(j,dx,dy,dy_temp,dz)  reduction(+:P,sum_force,flag)
    for(i=0;i<Np;i++){
      P +=v[i][0]*f[i][0]+v[i][1]*f[i][1]+v[i][2]*f[i][2];      
      for(j=0;j<dim;j++)
	v[i][j] = (1.0-alpha)*v[i][j]+alpha*f[i][j]/(f_nor+DBL_EPSILON)*v_nor;
    
      sum_force += sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]+f[i][2]*f[i][2])/Np;
      
      dy = x[i][1]-x0[i][1]; 
      dx = x[i][0]-x0[i][0]; 
      dz = x[i][2]-x0[i][2]; 
      dy_temp=dy;
      dy -= L*floor((dy+0.5*L)/L);
      dx -= gamma*L*floor((dy_temp+0.5*L)/L);
      dx -= L*floor((dx+0.5*L)/L);
      dz -= L*floor((dz+0.5*L)/L);

      if(dx*dx+dy*dy+dz*dz > skin*skin*0.25){
	flag += 1;
	imax = i;
	dr2max = dx*dx+dy*dy+dz*dz;
      }
    }
    //    #pragma omp barrier

    if(flag != 0){
      printf("cell update: FIRE %f,imax=%d,U=%f\n",dr2max,imax,*U/Np);
      cell_list(list,x,M,gamma);
      flag=0;
      #pragma omp parallel for private(j)
      for(i=0;i<Np;i++)
	for(j=0;j<dim;j++) 
	  x0[i][j]=x[i][j];      
    }
    
    //printf("FIRE: sum_force=%.16f,dt=%f,x=%f, alpha=%f,P=%f,gamma=%f \n",sum_force,dt0,x[0][0],alpha,P,gamma);

    if(P>=0){
      conv++;
      if(conv>5){
	dt0*=finc;
	if(dt0>0.002)
	  dt0=0.002;
	alpha*=falpha;
	conv=0;
      }
    }
    
    if(P<0){
      alpha=0.1;
      dt0*=fdec;
      conv=0.0;
      ini_matrix_dim(v);

    }
    if(sum_force <1e-12||count>1e9){
      scale_sliding_blick(x,gamma);
      std::cout<<"gamma="<<gamma<<"\t"<<"Iteration times ="<< count <<"\t"<<"energy="<<*U/(double)Np<<"\t"<<"stress ="<<*rfxy<<std::endl;
      break;  
    }
  }
  return 0;
}


void eom_md(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,int (*list)[Nn],double *stress){
  double dummy=0.0;
  int i,j;
 #pragma omp parallel for private(j) 
 for(i=0;i<Np;i++)
    for(j=0;j<dim;j++){
      x[i][j]+=v[i][j]*dt+0.5*f[i][j]*dt*dt;
      v[i][j]+=0.5*f[i][j]*dt;
    }
  calc_force(x,f,a,U,&dummy,list,0.0,stress);

  #pragma omp parallel for private(j)
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]+=0.5*f[i][j]*dt;
    }
  p_boundary(x);
}

void output_coord(double (*x)[dim],double(*x0)[dim],int *a,double gamma,double *stress){
  double dx,dy,dy_temp,dz; 
  char filename[128];
  std::ofstream file;
  sprintf(filename,"coord_disp_gamma%.3f.csv",gamma);
  file.open(filename);
  file<< std::setprecision(6)<<"# type[i],x[i],y[i],z[i],stress_xy[i]"<<std::endl;
  for(int i=0;i<Np;i++){
    /*  dy = x[i][1]-x0[i][1];
    dx = x[i][0]-x0[i][0];
    dz = x[i][2]-x0[i][2];
    dy_temp=dy;
    dy -= L*floor((dy+0.5*L)/L);
    dx -= gamma*L*floor((dy_temp+0.5*L)/L);
    dx -= L*floor((dx+0.5*L)/L);
    dz -= L*floor((dz+0.5*L)/L);*/
    if(a[i]==1)
      file<< std::setprecision(10)<<(int)0<<","<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<","<<stress[i]*Np<<std::endl;
    if(a[i]==2)
      file<< std::setprecision(10)<<(int)1<<","<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<","<<stress[i]*Np<<std::endl;
  }
  file.close();
}


void output(int k,double (*v)[dim],double U){
  char filename[128];
  double K=0.0;
  
  std::ofstream file;
  sprintf(filename,"energy.dat");
  file.open(filename,std::ios::app); //append
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      K+=0.5*v[i][j]*v[i][j];

  std::cout<< std::setprecision(6)<<k*dtmd<<"\t"<<K/Np*2./3.<<"\t"<<U/Np<<"\t"<<(K+U)/Np<<std::endl;
  file<< std::setprecision(6)<<k*dtmd<<"\t"<<K/Np*2./3.<<"\t"<<U/Np<<"\t"<<(K+U)/Np<<std::endl;
  file.close();
}

void update(double (*x_update)[dim],double (*x)[dim])
{
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim])
{
  int i;
  double dx,dy,dz;
  double disp,disp_max0;
  disp_max0 = *disp_max;
  //  #pragma omp parallel for private(i) reduction(max:disp_max0)
  for(i=0;i<Np;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dz=x[i][2]-x_update[i][2];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    dz-=L*floor((dz+0.5*L)/L);
    disp = dx*dx+dy*dy+dz*dz;
    if(disp > disp_max0)
      disp_max0 =disp;
  }
  *disp_max = disp_max0;
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update);
  if(*disp_max > skin*skin*0.25){
    cell_list(list,x,M,0.0);
    update(x_update,x);
    *disp_max=0.0;
    count=0;
  }
}

int main(){
  omp_set_num_threads(nthread);

#pragma omp parallel
  init_genrand((unsigned)time(NULL)*omp_get_thread_num());


  double x[Np][dim],x0[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],stress[Np];
  int list[Np][Nn],a[Np];
  double tout=0.0,U,rfxy,kine,disp_max=0.0,temp_anneal,gamma=0.0;
  int j=0;
  int M = (int)(L/(cut*sqrt(2.0)+skin));
  // std::cout<<L/(cut*sqrt(2.0)+skin)<<std::endl;
  set_diameter(a);
  ini_coord_rand(x);
 
  ini_matrix_dim(v);
  cell_list(list,x,M,0.0);
  std::cout<<"L="<<L<<" "<<"M="<<M<<std::endl;
  
  while(j*dtbdhs < 100.){
    j++;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin_hs(v,x,f,a,&U,dtbdhs,0.0,list,&kine);
    //  std::cout<<x[0][2]<<" "<<list[0][0]<<" "<<kine<<std::endl;
  }
  
  j=0;
  while(j*dtbd < teq){
    j++;
    temp_anneal = 4.0-j*dtbd*(4.0-temp)/teq;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp_anneal,list,&kine,stress);
    // std::cout<<kine<<" "<<U/Np<<std::endl;
     // std::cout<<x[0][2]<<" "<<f[0][0]<<" "<<kine<<std::endl;
  }
  
  j=0;
  while(j*dtbd < teq){
    j++;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp,list,&kine,stress);
    // std::cout<<kine<<" "<<U/Np<<std::endl;
    //  std::cout<<U/Np<<std::endl;
  }

  double  start =omp_get_wtime(); 
  steepest_descent(x,f,gamma,list,a,M,&U,stress);
  FIRE(x,f,gamma,list,a,M,&U,&rfxy,stress);

  double  end =omp_get_wtime();
 
  std::cout << "duration = " << (double)(end - start) << "sec.\n"; 
  copy_array(x,x0);
  int count=0;

  while(gamma < 0.5){
    count++;
    gamma += d_gamma;
    for(int i=0;i<Np;i++)
      x[i][0] += d_gamma*x[i][1];
    steepest_descent(x,f,gamma,list,a,M,&U,stress);
    FIRE(x,f,gamma,list,a,M,&U,&rfxy,stress);
    if(count == 1){
      output_coord(x,x0,a,gamma,stress);
      count = 0;
    }
    std::cout<<"gamma ="<< gamma <<std::endl;
  }
  return 0;
}
