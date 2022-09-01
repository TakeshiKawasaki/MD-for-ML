#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>
using namespace std;

#define Npm 30000
#define Pm 300

#define max 1000
#define Nf 1000
#define Ndm 1500
double pi= 3.14159265359;

double unif_rand(double left, double right) 
{
  return left + (right - left)*rand()/RAND_MAX;
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


int init_grv(char *iname,double *x,double *y,int Np,int Count,double L)
{
  int i;
  char filename[256];
  ifstream file;
  srand((unsigned) time(NULL));
  for(i=0;i<Np;i++){
    x[i]=L*rand()/(double)RAND_MAX;
    y[i]=L*rand()/(double)RAND_MAX;
  }

  return 0;
}

int output(char *Name,int l,double *x,double *y,
	   double N,int Np)
{
  int i;
  FILE *fp;
  char FileName[64];
  double x0[Npm];
  sprintf(FileName,"%s_%d.grv",Name,l);
  fp=fopen(FileName,"w");
  
  
  for(i=0;i<Np;i++){
   
    fprintf(fp,"%g\t%g\n",x[i],y[i]);
  }
  
  fclose(fp);
  printf("\n File %s_%d.bd OK\n",Name,l);
  return 0;
}

int scale_sliding_blick_x(double *x,double *y,double N,int Np,double gamma)
{
  int i;
  
  for(i=0;i<Np;i++)
    {
      if(x[i]<gamma*(y[i]))
	x[i]+=N;
      if(x[i]>N+gamma*(y[i]))
	x[i]-=N;
    }
  
  return 0;
}

int scale_sliding_blick_y(double *x,double *y,double N,int Np, double gamma)
{
  int i;
  for(i=0;i<Np;i++){
    if(y[i]<0.0){
      y[i]+=N;
      x[i]+=gamma*N;
      
    }
    if(y[i]>N){
      y[i]-=N;
      x[i]-=gamma*N;
    }
  }
  return 0;}


int scale_z(double *x,double N,int Np){
  int i;
  for(i=0;i<Np;i++){
    if(x[i]<0.0)
      x[i]+=N;
    if(x[i]>N)
      x[i]-=N;
  }
  return 0;
}


int output_grv(char *name,double *x,double *y,double t,
	       double N,int Np)
{
  int i;
  char filename[64];
  FILE *fp;
  for(i=0;i<1;i++)
    {
      sprintf(filename,"%s.%d",name,i);
      fp=fopen(filename,"a");
      fprintf(fp,"%g\t%g\t\t%g\n",t,x[i],y[i]);
      fclose(fp);
    }
  return 0;
}

int output_gr2(char *name,int l,double *x,double *y,int Np,int (*list)[Pm],double N){	
  int i,j;
  char filename[64];
  ofstream file;
  double x0[Npm];
  sprintf(filename,"%s_%d.gr2",name,l);
  file.open(filename);
  
 
  for(i=0;i<Np;i++)
    file << x[i] << "   " << y[i]  <<endl;
  
  file.close();
  return 0;
}

int ini_coord(double *x,double *y,double *x0,double *y0,int Np){
for(int j=0;j<Np;j++){
  x0[j]=x[j];
  y0[j]=y[j];
 }
 return 0;
}


int brokenbond(double *x,double *y,double *x0,double *y0,int Np,double L, int (*map)[Ndm],double (*BB)[1000], double gamma, double *r1,int m,int k){
  double dx,dy,r2,rab;
  for(int j=0;j<Np;j++){
    for(int i=1;i<=map[j][0];i++){
      dx=x[map[j][i]]-x[j];
      dy=y[map[j][i]]-y[j];
      if(dy>(L/2.)){
	dy-=L;
	dx-=gamma*L;
      }
      if(dy<-(L/2.)){
	dy+=L;
	dx+=gamma*L;
      }
      
      if(dx>(L/2.))
	dx-=L;
      if(dx<-(L/2))
	dx+=L;
      r2=dx*dx+dy*dy;
      rab=0.5*(r1[j]+r1[map[j][i]]);
      if(r2 > 1.0*rab*rab)
	BB[m][k]+=1.0/Np;
    }
  }
  return 0;
}


int calc_f(double *x,double *y,double *U
	       ,double *rfxy ,double *rf, double L,int Np, double *r1,int (*list)[Pm],double gamma) 
{                                        
  int i,j;
  double r;
  double t,fo,Uo;
  double dx,dy;
  double cut;
  double rab;
  
  cut = 1.0;

  *rfxy=0.0;
  *rf=0.0;
  *U=0.0;
 
  for(j=0;j<Np;j++){
    for(i=1;i<=list[j][0];i++){
      dx=x[list[j][i]]-x[j];
      dy=y[list[j][i]]-y[j];

      if(dy>(L/2.)){
	dy-=L;
	dx-=gamma*L;
      }
      if(dy<-(L/2.)){
	dy+=L;
	dx+=gamma*L;
      }

      if(dx>(L/2.))
	dx-=L;
      if(dx<-(L/2))
	dx+=L;
      
      r=sqrt(dx*dx+dy*dy);
      rab=(r1[list[j][i]]+r1[j])/2.0;
      t=r/rab; 
      
      if(t<=cut)
	fo = -(1-t)*sqrt(1-t)/rab;
      else
	fo = 0.0;
      
      *rfxy +=dx*dy*fo/r/(L*L);
      *rf -=r*fo/(2.0*L*L);
      
      if(t<=cut)       
	Uo= (1-t)*(1-t)*sqrt(1-t)/2.5;
      else
	Uo= 0.0;
     
      // printf("U=%f\n",Uo);
      *U+=Uo;
    }
  }
  return 0;
}

int	calc_force(double *x,double *y,double N,int Np,
		   double *r1,double *kx,double *ky,int (*list)[Pm],double gamma,double *rfxy,double *rf)
{
  int i,j;                                  
  int k;
  
  double r;
  double t,f;
  double dx,dy;
  
  double rab;  
  double cut;
  
  cut=1.0;
  *rfxy=0.0;
  *rf=0.0;
    
  for(k=0;k<Np;k++){
    kx[k]=0.0;
    ky[k]=0.0;
  }
  
  for(j=0;j<Np;j++)
    {
      for(i=1;i<=list[j][0];i++)
	{
	  dx=x[list[j][i]]-x[j];
	  dy=y[list[j][i]]-y[j];

	  /*	if(dx>(N/2)+gamma*dy)
	    dx-=N;
	    if(dx<-(N/2)-gamma*dy)
	    dx+=N;	*/
	  
	  if(dy>(N/2.)){
	    dy-=N;
	    dx-=gamma*N;
	  }
	  if(dy<-(N/2.)){
	    dy+=N;
	    dx+=gamma*N;
	  }
	  	  
	  if(dx>(N/2.))
	    dx-=N;
	  if(dx<-(N/2.))
	    dx+=N;
	  
	 	  
	  rab=(r1[list[j][i]]+r1[j])/2.0;
	  r=sqrt(dx*dx+dy*dy);
	  t=r/rab;  
	  
	  if(t<cut)
	    f =-(1-t)*sqrt(1-t)/rab;
	  else
	    f =0.0;
	  //  continue;
	  //	
	  *rfxy +=dx*dy*f/r/(N*N);
	  *rf -=r*f/(2.0*N*N);
	  
	  kx[j]+=f*dx/r;
	  kx[list[j][i]]-=f*dx/r;
	  ky[j]+=f*dy/r;
	  ky[list[j][i]]-=f*dy/r;
	  // printf("alpha=%f,r=%f,f=%f,kx=%f\n",alpha,t,f,kx[j]);
	}
    }
  return 0;
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
void update(double N,int Np,double *x,double *y,int M,double RCHK,int (*list)[Pm],double gamma)
{
  int i,j,k;
  int nx,ny;
  int l,m,n;
  double dx,dy,r;
  
  int (*map)[Npm]=new int[M*M*M][Npm];
  
  for(k=0;k<M;k++)
    for(j=0;j<M;j++)
      for(i=0;i<M;i++)
     	map[i+M*j+M*M*k][0]=0;
  
  for(i=0;i<Np;i++){
    nx=f((int)((x[i]-gamma*y[i])*M/N),M);/*new*/
    ny=f((int)(y[i]*M/N),M);
    
      for(m=ny-1;m<=ny+1;m++){
	for(l=nx-1;l<=nx+1;l++){
	  map[f(l,M)+M*f(m,M)][ map[f(l,M)+M*f(m,M)][0] +1]=i;
	  map[f(l,M)+M*f(m,M)][0]++;
	  //printf("%d\n", map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0]);
	}
      }
  }
  
  for (i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)((x[i]-gamma*y[i])*M/N),M);
    ny = f((int)(y[i]*M/N),M);
    
    for (k=1; k<=(map[nx+M*ny][0]); k++){
      j = map[nx+M*ny][k];
      if(j>i)	{
	dx =x[i] - x[j];
	dy =y[i] - y[j];
	
	if(dy<-N/2.0){
	  dy+=N;
	  dx+=gamma*N;
	}
	else if(dy> N/2.0){
	  dy-=N;
	  dx-=gamma*N;
	}
	
	if(dx<-N/2.0)
	  dx+=N;
	else if(dx> N/2.0)
	  dx-=N;
	r = dx*dx + dy*dy;
	
	if(r<RCHK*RCHK){
	  list[i][0]++;
	  list[i][list[i][0]]=j;
	}
      }
    }
  }	
  delete []map;
}


int langevin_eqSS(double *x,double *y,double *vx,double *vy,double dt,double N,int Np,double *kx,double *ky,
		double gamma_dot,int (*list)[Pm], double *r1,double RCHK,int M,double *gamma,double Pext)
{
  int i,imax,k;
  double dx,dy,dxmax=0.0,dymax=0.0,x0[Npm],y0[Npm],rfxy,rf,V_dot,N0; 
  for(i=0;i<Np;i++){
    x0[i]=x[i];
    y0[i]=y[i];
  }

  for(k=0;k<30000;k++){
    dxmax=0.0;
    dymax=0.0;
    calc_force(x,y,N,Np,r1,kx,ky,list,*gamma,&rfxy,&rf);  
    for(i=0;i<Np;i++){
      vx[i]=gamma_dot*y[i]+kx[i];
      vy[i]=ky[i]; 
      x[i]=x[i]+vx[i]*dt; 
      y[i]=y[i]+vy[i]*dt;
      // printf("vx=%f,vy=%f,vz=%f,kx=%f,ky=%f,kz=%f \n",vx[i],vy[i],vz[i],kx[i],ky[i],kz[i]);
      
      dy=y[i]-y0[i]; if(dy>0.5*N){dy-=N; dx-=*gamma*N;} if(dy<-0.5*N){dy+=N; dx+=*gamma*N;}
      dx=x[i]-x0[i]; if(dx>0.5*N) dx-=N; if(dx<-0.5*N) dx+=N;
      
      if(dx*dx+dy*dy>dxmax*dxmax+dymax*dymax){
        dxmax=dx;
        dymax=dy;
        imax=i;
      }
    }
    
    gamma_dot = -N*rfxy;
    //    V_dot= 0.1**N**N*(rf-Pext);
    
    *gamma+=gamma_dot*dt;
  
    //  *N+= pow(V_dot,1.0/3.0)*dt;
    
    if(dxmax*dxmax+dymax*dymax > 0.5*0.5){
      printf("langevin %f,%f,imax=%d,t=%d\n",dxmax,dymax,imax,k);
      update(N,Np,x,y,M,RCHK,list,*gamma);
      for(i=0;i<Np;i++){         
	x0[i]=x[i];
	y0[i]=y[i];
      }
      dxmax=0.0;
      dymax=0.0;
    }
  }
  
  printf("langevin gamma=%.35f\n",*gamma);
  scale_sliding_blick_x(x,y,N,Np,*gamma);
  scale_sliding_blick_y(x,y,N,Np,*gamma);
  return 0;
}


int displacement(double *x,double *y,double *x0,double *y0,double N,int Np,double gamma,double *disp)
{
  int i;
  double dx,dy;
  *disp=0.0;
  for(i=0;i<Np;i++){
    dy=y[i]-y0[i];
    dx=x[i]-x0[i];
    if(dy> (N/2.)){
      dy-=N;
      dx-=gamma*N;
    }
    if(dy<-(N/2.)){
      dy+=N;
      dx+=gamma*N;
    }
   
    if(dx>(N/2.))
      dx-=N;
    if(dx<-(N/2.))
      dx+=N;
    *disp+=sqrt(dx*dx+dy*dy)/(double)Np;
  }
  return 0;
}

int  set_map(double *x,double *y,double *r1,int (*map)[Ndm],double N,int Np,double gamma)
{
  int i,j;
  int n;
  double r;
  double dx,dy;
  double r0; 
 
  for(i=0;i<Np;i++)
    {
      n=0;
      for(j=0;j<Np;j++)
	if(i!=j)
	  {
	    dy=y[i]-y[j];
	    dx=x[i]-x[j];
	    if(dy> (N/2.)){
	      dy-=N;
	      dx-=gamma*N;
	    }
	    if(dy<-(N/2.)){
	      dy+=N;
	      dx+=gamma*N;
	    }

	
	    if(dx>(N/2.))
	      dx-=N;
	    if(dx<-(N/2.))
	      dx+=N;
	    
	   
	   	    
	    r=sqrt(dx*dx+dy*dy);
	    r0=(r1[i]+r1[j])*0.5;
	    if(r<r0)
	      {
		n++;
		map[i][n]=j;
	      }
	  }
      map[i][0]=n;
    }
  
  return 0;
}

int CG(double *x,double *y,double N,int Np,double *kx,double *ky,double gamma,int (*list)[Pm], double *r1,double RCHK,int M)
{
  int i,imax;
  double alpha = 0.1;
  double hx[Npm],hy[Npm],hxb[Npm],hyb[Npm],kxb[Npm],kyb[Npm],x0[Npm],y0[Npm],dx,dy,dxmax=0.0,dymax=0.0;
  int count=0;
  double U,Ub=0.0; 
  double grad2=0.0,grad2b=0.0;
  double sum_force;
  double rfxy,rf;
  
  // double RCHK=3.0;
  // int  M=(int)(N/RCHK);
  int conv=0;
  
  calc_force(x,y,N,Np,r1,kx,ky,list,gamma,&rfxy,&rf);

  for(i=0;i<Np;i++){
    grad2b+=kx[i]*kx[i]+ky[i]*ky[i];
    hxb[i]=0.0;
    hyb[i]=0.0;
    kxb[i]= kx[i];
    kyb[i]= ky[i];
    x0[i]=x[i];
    y0[i]=y[i];
  }
  // update(N,Np,x,y,M,RCHK,list,gamma); 
  for(;;){
    sum_force=0.0;

        for(i=0;i<Np;i++)
      grad2b+=kx[i]*kx[i]+ky[i]*ky[i];   
      
    count++;
    calc_force(x,y,N,Np,r1,kx,ky,list,gamma,&rfxy,&rf);
   
    // printf("k[1]=%f,x[1]=%f \n", kx[1],x[1]);
    // printf("U=%f,count=%d,k[1]=%f,x[1]=%f \n",U,count, kx[1],x[1]);
    // grad2=0.0;
   
    for(i=0;i<Np;i++){
   
      sum_force+=(kx[i]*kx[i]+ky[i]*ky[i])/Np;
      hx[i] = kx[i] + grad2/(grad2b+DBL_EPSILON)*hxb[i];
      hy[i] = ky[i] + grad2/(grad2b+DBL_EPSILON)*hyb[i];
      x[i]+= alpha*hx[i]; 
      y[i]+= alpha*hy[i];
     
      hxb[i]=hx[i];
      hyb[i]=hy[i];
      kxb[i]= kx[i];
      kyb[i]= ky[i];

      dy=y[i]-y0[i]; if(dy>0.5*N){dy-=N; dx-=gamma*N;} if(dy<-0.5*N){dy+=N; dx+=gamma*N;} 
      dx=x[i]-x0[i]; if(dx>0.5*N) dx-=N; if(dx<-0.5*N) dx+=N;
      
      if(dx*dx+dy*dy>dxmax*dxmax+dymax*dymax){
	dxmax=dx;
	dymax=dy;
	imax=i;
      }
    }
        
    if(dxmax*dxmax+dymax*dymax > 0.7*0.7){
      printf("OK %f,%f,imax=%d,count=%d\n",dxmax,dymax,imax,count);
      update(N,Np,x,y,M,RCHK,list,gamma);        
      for(i=0;i<Np;i++){         
	x0[i]=x[i];
	y0[i]=y[i];
      } 
      dxmax=0.0;
      dymax=0.0;
    }


    
    //   printf("U=%.12f,count=%d,x=%.12f \n",U,count,x[1]);
    //    scale_sliding_blick_x(x,y,N,Np,gamma);
    //  scale_sliding_blick_y(x,y,N,Np,gamma);
   
    // if((U-Ub)*(U-Ub) < 0.00000000001 && (U-Ub)>0.0)
    //  conv++;

    if(count==1000){
      printf("sum_force=%.31f,count=%d,x=%.25f,y=%.25f,gamma=%.10f \n",sqrt(sum_force),count,x[100],y[100],gamma);
      count=0;
      break;
    }

    //  Ub=U;
     grad2b=grad2;
   
  }

  scale_sliding_blick_x(x,y,N,Np,gamma);
  scale_sliding_blick_y(x,y,N,Np,gamma);
  
  return 0;
}




int Fire(double *x,double *y,double N,int Np,double *kx,double *ky,double gamma,int (*list)[Pm], double *r1,double gamma_dot,double RCHK,int M)
{
  int i,imax;
  double alpha = 0.1,P;
  double dt0=0.005;
  double hx[Npm],hy[Npm],hxb[Npm],hyb[Npm],kxb[Npm],kyb[Npm],vx[Npm],vy[Npm],x0[Npm],y0[Npm],dx,dy,dxmax,dymax,rfxy,rf;
  int count=0,count0=0;
  int conv=0;
  double sum_force=0.0;
  double finc=1.01,falpha=0.99,fdec=0.5;
  
  
  update(N,Np,x,y,M,RCHK,list,gamma); 
  for(i=0;i<Np;i++){
    x[i]+=gamma_dot*y[i];
    vx[i]= 0.0;
    vy[i]= 0.0;
    x0[i]=x[i];
    y0[i]=y[i];
  }
  P=0.0;
  dxmax=0.0;dymax=0.0; //call list update
  
  calc_force(x,y,N,Np,r1,kx,ky,list,gamma,&rfxy,&rf);
  
  for(;;){
    
    P=0.0;
    sum_force=0.0;
    count++;
    count0++;
    calc_force(x,y,N,Np,r1,kx,ky,list,gamma,&rfxy,&rf);
    
    for(i=0;i<Np;i++){
      vx[i]+=kx[i]*dt0;
      vy[i]+=ky[i]*dt0;
      
      x[i]+=vx[i]*dt0;
      y[i]+=vy[i]*dt0;
      
      P+=vx[i]*kx[i]+vy[i]*ky[i];
      vx[i]= (1.0-alpha)*vx[i]+alpha*kx[i]/sqrt(kx[i]*kx[i]+ky[i]*ky[i]+DBL_EPSILON)*sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
      vy[i]= (1.0-alpha)*vy[i]+alpha*ky[i]/sqrt(kx[i]*kx[i]+ky[i]*ky[i]+DBL_EPSILON)*sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
      sum_force+= sqrt(kx[i]*kx[i]+ky[i]*ky[i])/Np;
      
      dy=y[i]-y0[i]; if(dy>0.5*N){dy-=N; dx-=gamma*N;} if(dy<-0.5*N){dy+=N; dx+=gamma*N;} 
      dx=x[i]-x0[i]; if(dx>0.5*N) dx-=N; if(dx<-0.5*N) dx+=N;
      
      if(dx*dx+dy*dy>dxmax*dxmax+dymax*dymax){
	dxmax=dx;
	dymax=dy;
	imax=i;
      }
    }
    
    if(dxmax*dxmax+dymax*dymax > 0.7*0.7){
      printf("OK %f,%f,imax=%d\n",dxmax,dymax,imax);
      update(N,Np,x,y,M,RCHK,list,gamma);
      for(i=0;i<Np;i++){         
	x0[i]=x[i];
	y0[i]=y[i];
      }
      dxmax=0.0;
      dymax=0.0;
    }
    
    //  printf("sum_force=%.31f,alpha=%f,x=%.25f,y=%.25f,P=%f,gamma=%f \n",sum_force,alpha,x[100],y[100],P,gamma);
    if(P>=0){
      conv++;
      if(conv>5){
	dt0*=finc;
	if(dt0>0.05)
	  dt0=0.05;
	alpha*=falpha;
	conv=0;
      }
    }
    
    if(P<0){
      //  update(N,Np,x,y,M,RCHK,list,gamma);           
      alpha=0.1;
      dt0*=fdec;
      if(dt0<=0.00001)
	dt0=0.00001;
      conv=0.0;
      for(i=0;i<Np;i++){
	vx[i]= 0.0;
	vy[i]= 0.0;
      }
    }
    if(sum_force <1e-14||count>1e9){
      printf("sum_force=%.31f,alpha=%f,x=%.25f,y=%.25f,P=%.30f,gamma=%.10f,c=%d \n",sum_force,alpha,x[100],y[100],P,gamma,count);
      
      scale_sliding_blick_x(x,y,N,Np,gamma);
      scale_sliding_blick_y(x,y,N,Np,gamma);
      break;  
    }
  }
  return 0;
}


int FireP(double *x,double *y,double *N,int Np,double *kx,double *ky,double gamma,int (*list)[Pm], double *r1,double gamma_dot,double RCHK,int M,double rfex)
{
  int i,imax;
  double alpha = 0.1,P;
  double dt0=0.001;
  double hx[Npm],hy[Npm],hxb[Npm],hyb[Npm],kxb[Npm],kyb[Npm],vx[Npm],vy[Npm],x0[Npm],y0[Npm],dx,dy,dxmax,dymax,rfxy,rf;
  int count=0,count0=0;
  int conv=0;
  double sum_force=0.0;
  double finc=1.01,falpha=0.99,fdec=0.5;
  double Vdot;
  double Nc;
  double delrf;
  
  update(*N,Np,x,y,M,RCHK,list,gamma); 

  P=0.0;
  dxmax=0.0;dymax=0.0; //call list update
  Vdot=0.0;
  
  calc_force(x,y,*N,Np,r1,kx,ky,list,gamma,&rfxy,&rf);
  // rfex=rf;

  for(i=0;i<Np;i++){
    x[i]+=gamma_dot*y[i];
    vx[i]= 0.0;
    vy[i]= 0.0;
    x0[i]=x[i];
    y0[i]=y[i];
  }
  
  for(;;){
    P=0.0;
    sum_force=0.0;
    count++;
    count0++;
    calc_force(x,y,*N,Np,r1,kx,ky,list,gamma,&rfxy,&rf);
    
    Nc=*N;
    delrf=rf-rfex;
    Vdot=Vdot+1000.*delrf*dt0;
    *N+=0.5/(Nc)*Vdot*dt0;
    
    P+=delrf*Vdot;
    
    Vdot=(1.0-alpha)*Vdot+alpha*delrf/sqrt(delrf*delrf+DBL_EPSILON)*sqrt(Vdot*Vdot);
   // sum_force+=sqrt(delrf*delrf)/Np;
    
    for(i=0;i<Np;i++){
      vx[i]+=kx[i]*dt0;
      vy[i]+=ky[i]*dt0;
      
      x[i]+=vx[i]*dt0;
      y[i]+=vy[i]*dt0;
      // scaling the coordinates    
      x[i]=x[i]**N/Nc; 
      y[i]=y[i]**N/Nc; 
      
      P+=vx[i]*kx[i]+vy[i]*ky[i];
      vx[i]= (1.0-alpha)*vx[i]+alpha*kx[i]/sqrt(kx[i]*kx[i]+ky[i]*ky[i]+DBL_EPSILON)*sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
      vy[i]= (1.0-alpha)*vy[i]+alpha*ky[i]/sqrt(kx[i]*kx[i]+ky[i]*ky[i]+DBL_EPSILON)*sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
      sum_force+= sqrt(kx[i]*kx[i]+ky[i]*ky[i])/Np;
      
      dy=y[i]-y0[i]; if(dy>0.5**N){dy-=*N; dx-=gamma**N;} if(dy<-0.5**N){dy+=*N; dx+=gamma**N;} 
      dx=x[i]-x0[i]; if(dx>0.5**N) dx-=*N; if(dx<-0.5**N) dx+=*N;
      
      if(dx*dx+dy*dy>dxmax*dxmax+dymax*dymax){
	dxmax=dx;
	dymax=dy;
	imax=i;
      }
    }
    
    if(dxmax*dxmax+dymax*dymax > 0.5*0.5){
      printf("OK %f,%f,imax=%d\n",dxmax,dymax,imax);
      update(*N,Np,x,y,M,RCHK,list,gamma);
      for(i=0;i<Np;i++){         
	x0[i]=x[i];
	y0[i]=y[i];
      }
      dxmax=0.0;
      dymax=0.0;
    }
    
    if(P>=0){
      conv++;
      if(conv>5){
	dt0*=finc;
	if(dt0>0.1)
	  dt0=0.1;


	alpha*=falpha;
	conv=0;
      }
    }
    
    if(P<0){
      alpha=0.1;
      dt0*=fdec;
      conv=0.0;
      if(dt0<0.00001)
	dt0=0.00001;
      for(i=0;i<Np;i++){
	vx[i]= 0.0;
	vy[i]= 0.0;
      }
    }
    scale_sliding_blick_x(x,y,*N,Np,gamma);
    scale_sliding_blick_y(x,y,*N,Np,gamma);
    if(sum_force <1e-14||count>1e9){
      printf("sum_force=%.31f,alpha=%f,x=%.25f,y=%.25f,P=%.30f,gamma=%.10f,c=%d \n",sum_force,alpha,x[100],y[100],P,gamma,count);
      break;  
    }
  }
  return 0;
}


int FireSS(double *x,double *y,double N,int Np,double *kx,double *ky,double *gamma,int (*list)[Pm], double *r1,double gamma_dot,double RCHK,int M)
{
  int i,imax;
  double alpha = 0.1,P;
  double dt0=0.005;
  double hx[Npm],hy[Npm],hxb[Npm],hyb[Npm],kxb[Npm],kyb[Npm],vx[Npm],vy[Npm],x0[Npm],y0[Npm],dx,dy,dxmax,dymax,rfxy,rf;
  int count=0,count0=0;
  int conv=0;
  double sum_force=0.0;
  double finc=1.01,falpha=0.99,fdec=0.5;
  double gdot;
 
  
  
  update(N,Np,x,y,M,RCHK,list,*gamma); 
  for(i=0;i<Np;i++){
    x[i]+=gamma_dot*y[i];
    vx[i]= 0.0;
    vy[i]= 0.0;
    x0[i]=x[i];
    y0[i]=y[i];
  }
  P=0.0;
  gdot=0.0;
  rfxy=0.0;
  dxmax=0.0;dymax=0.0; //call list update
  
  calc_force(x,y,N,Np,r1,kx,ky,list,*gamma,&rfxy,&rf);
  
  for(;;){
    P=0.0;
    sum_force=0.0;
    count++;
    count0++;
    calc_force(x,y,N,Np,r1,kx,ky,list,*gamma,&rfxy,&rf);
    
    gdot-=100.*N*rfxy*dt0;
    *gamma+=gdot*dt0;
    P+=-rfxy*gdot*N*N;
    gdot=(1-alpha)*gdot-alpha*rfxy/sqrt(rfxy*rfxy+DBL_EPSILON)*sqrt(gdot*gdot);
    sum_force+=sqrt(rfxy*rfxy)*N/Np;     
    for(i=0;i<Np;i++){
      vx[i]+=kx[i]*dt0;
      vy[i]+=ky[i]*dt0;
      
      x[i]+=vx[i]*dt0;
      y[i]+=vy[i]*dt0;
      
      P+=vx[i]*kx[i]+vy[i]*ky[i];
      vx[i]= (1.0-alpha)*vx[i]+alpha*kx[i]/sqrt(kx[i]*kx[i]+ky[i]*ky[i]+DBL_EPSILON)*sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
      vy[i]= (1.0-alpha)*vy[i]+alpha*ky[i]/sqrt(kx[i]*kx[i]+ky[i]*ky[i]+DBL_EPSILON)*sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
      sum_force+= sqrt(kx[i]*kx[i]+ky[i]*ky[i])/Np;
      
      dy=y[i]-y0[i]; if(dy>0.5*N){dy-=N; dx-=*gamma*N;} if(dy<-0.5*N){dy+=N; dx+=*gamma*N;} 
      dx=x[i]-x0[i]; if(dx>0.5*N) dx-=N; if(dx<-0.5*N) dx+=N;
      
      if(dx*dx+dy*dy>dxmax*dxmax+dymax*dymax){
	dxmax=dx;
	dymax=dy;
	imax=i;
      }
    }
    
    if(dxmax*dxmax+dymax*dymax > 0.5*0.5){
      printf("OK %f,%f,imax=%d\n",dxmax,dymax,imax);
      update(N,Np,x,y,M,RCHK,list,*gamma);
      for(i=0;i<Np;i++){         
	x0[i]=x[i];
	y0[i]=y[i];
      }
      dxmax=0.0;
      dymax=0.0;
    }
    
    //  printf("sum_force=%.31f,alpha=%f,x=%.25f,y=%.25f,P=%f,gamma=%f \n",sum_force,alpha,x[100],y[100],P,gamma);
    if(P>=0){
      conv++;
      if(conv>5){
	dt0*=finc;
	if(dt0>0.05)
	  dt0=0.05;
	alpha*=falpha;
	conv=0;
      }
    }
    
    if(P<0){
      //  update(N,Np,x,y,M,RCHK,list,gamma);           
      alpha=0.1;
      dt0*=fdec;
      conv=0.0;
      gdot=0.0;
      if(dt0<=0.00001)
	dt0=0.00001;
      for(i=0;i<Np;i++){
	vx[i]= 0.0;
	vy[i]= 0.0;
      }
    }
    if(sum_force <1e-14||count>1e9){
      printf("sum_force=%.31f,alpha=%f,x=%.25f,y=%.25f,P=%.30f,gamma=%.10f,c=%d \n",sum_force,alpha,x[100],y[100],P,*gamma,count);
      
      scale_sliding_blick_x(x,y,N,Np,*gamma);
      scale_sliding_blick_y(x,y,N,Np,*gamma);
      break;  
    }
  }
  return 0;
}



int main(int argc, char *argv[])
{
  static double x[Npm],y[Npm];
  static double vx[Npm],vy[Npm];
  static double kx[Npm],ky[Npm];
  static double x0[Npm],y0[Npm],x1[Npm],y1[Npm],N1,r1[Npm],r0[Npm],phiJgamma[10000],BB[500][1000],gamma0;
  double E[10000],U,Uave,rfxy,rfxy0[50000][1000],rf0[50000][1000],E0[50000][1000],ZZ0[50000][1000],countz0[50000][1000],rfxyave,P[50000],countz[50000],sigmaxy[10000],rf,rfave,rfxyave2,Pext,phiJave=0.0;
  int (*list)[Pm]=new int[Npm][Pm];
  int (*map)[Ndm]=new int[Npm][Ndm];
  int (*map0)[Ndm]=new int[Npm][Ndm];
  double Gp[10000];
  double t,t1;
  //  double omega=2.*pi/10000.;
  double Amp;
  
  double Z[10000],Psi1[10000],Z0,phiJ;
  double P0;
  double gamma_dot,gamma_dot2;
  double gamma,gamma_ave[100000],rfxy_ave[100000];
  double alpha;
  double tq;
  double rmax,RCHK,Rave;
  int M;
  int l;
  int lp;
  int Q;
  int i,j,k;
  int ns,Nstep,Count=0;
  
  double N;
  double Psi,Psi0,Psi2;
  double disp=  0.0;
  
  int Np;
  char *iname;
  int imode;
  double tmax,dt,t_int;
  
  double Fint[Npm];
  double sigma;
  char *oname;
  int step;
  int Ngrv=0,Nq=0;
  int Nb;
  int n;
  char buf[64];
  char filename[64];
  char filename3[64];
  char filename4[64];
  char filename5[64];
  
  ifstream file;
  ofstream file2;
  ofstream file3;
  
  FILE *fp;
  
  double BN[Nf];
  double Bamp[Nf];
  double BPsi[Nf];
  int BNp[Nf];
  char Biname[64][Nf];
  int Bimode[Nf];
  double Btmax[Nf],Bdt[Nf],Bt_int[Nf];
  int ens;
  double Bf0[Nf],Br0[Nf];
  double BP0[Nf];
  double Beta[Nf];
  double Bgamma_dot[Nf];
  double Bsigma[Nf];
  double Balpha[Nf];
  
  int Bl[Nf];
  double Btq[Nf];
    
  char Boname[64][Nf];
  int Bstep[Nf];

  for(k=0;k<500;k++){
    Gp[k]=0.0;
    Psi1[k]=0.0;
    Z[k]=0.0;
    P[k]=0.0;
    E[k]=0.0;
    sigmaxy[k]=0.0;
    countz[k]=0.0;
    for(i=0;i<1000;i++){
      rfxy0[k][i]=0.0;	
      rf0[k][i]=0.0;
      E0[k][i]=0.0;
      BB[k][i]=0.0;
      ZZ0[k][i]=0.0;
      countz0[k][i]=0.0;
    }	
  }
  for(k=0;k<1000;k++)
    phiJgamma[k]=0.0;

  double delphi=1.e-5;
  
  for(ens=1;ens<=1000;ens++){
    
    sprintf(filename,"ini.txt");
    file.open(filename);
    
    file >> buf >> Nb;
    
    for(n=0;n<Nb;n++){
      file >> buf >> BPsi[n];
      file >> buf >> Bamp[n];
      file >> buf >> BN[n];
      file >> buf >> Biname[n] >> Bimode[n];
      file >> buf >> Btmax[n] >> Bdt[n] >> Bt_int[n];
      file >> buf >> Bl[n] >> Btq[n];
      file >> buf >> BP0[n];
      file >> buf >> Balpha[n];
      file >> buf >> Bgamma_dot[n];
      file >> buf >> Boname[n];
      file >> buf >> Bstep[n];
    }
    
    file.close();
    
    for(n=0;n<Nb;n++){
      //N     =BN[n];
      Psi =BPsi[n];
      N    =BN[n];
      iname =Biname[n];
      Amp= Bamp[n];
      imode =Bimode[n];
      tmax  =Btmax[n];
      dt    =Bdt[n];
      t_int =Bt_int[n];
      l     =Bl[n];
      tq    =Btq[n];
      P0  =atof(argv[2]);
      alpha=Balpha[n];
      gamma_dot   =Bgamma_dot[n];
      oname =Boname[n];
      step  =Bstep[n];
      
          
      sprintf(filename,"%s.para",oname);
      fp=fopen(filename,"w");
      
      fprintf(fp,"ini.exe\n");
      fprintf(fp,"Lattice size N : %f\n",N);
      fprintf(fp,"Particle number Np : %d\n",Np);
      fprintf(fp,"Initialize file : %s\t",iname);
      if(imode==0)
	fprintf(fp,"(new) \n");
      else
	fprintf(fp,"(continue) \n");
      fprintf(fp,"t_max,dt,t_int : %g, %g, %g\n",tmax,dt,t_int);
      fprintf(fp,"t_0,t_q: %d, %g\n",l,tq);
      fprintf(fp," P, gamma_dot: %g, %g\n",P0,gamma_dot);
      
      fprintf(fp,"Output step %d\n",step);
      fclose(fp);
      
      Np=0;
      
      for(j=100;j<100000;j+=2){
	Psi0=0.0;
	for(i=0;i<j;i+=2){
	  r1[i]=1.0;
	  r1[i+1]=1.4;
	  Psi0+=pi*(r1[i]*r1[i])/4.0/(N*N);
	  Psi0+=pi*(r1[i+1]*r1[i+1])/4.0/(N*N);
	}
	
	if(Psi0>Psi){
	  Np=j;
	  printf("np=%d\n",Np);
	  break;
	}
      }
      
      sprintf(filename3,"np%d_0.dat",Np);
      file2.open(filename3);
      for(j=0;j<Np;j++){
	file2 << r1[j] <<endl;
      }
      file2.close();
      /* 
	 sprintf(filename3,"np%d_0.grv",Np);
	 file2.open(filename3);
	 for(j=0;j<Np;j++){
	 file2 << x[j]<<" "<<y[j]<<endl;
	 }
	 file2.close();
      */
            
      printf("%f %f %f %f %f L=%f, phi=%f \n",r1[0],r1[1],r1[2],r1[Np-2],r1[Np-1],N,Psi0);
      printf("np%d_0.dat\n",Np);
      rmax=r1[0];
      Rave=0.0;
      gamma=0;
      
      for(j=0;j<Np;j++){
	Rave+=r1[j]/Np;
	
	if(r1[j]>rmax)
	  rmax=r1[j];
      }
      
      RCHK=2.5;
      // Rave=Rave*2.0;
      
      M=(int)(N/RCHK);
      
      sprintf(filename3,"%s/%s/sigma_%s_%d.dat",argv[1],argv[2],oname,ens);
      file2.open(filename3);
      file2  <<" # "<< "Psi-PsiJ" <<"   "<< "Z[k]" <<"   "<< "P[k]"  <<"   "<< "sigmaxy[k]" <<"   "<< "Gp[k]" <<"   "<< "E[k]"<< endl;
      
      if(imode==0){
	init_grv(iname,x,y,Np,Count,N);
	update(N,Np,x,y,M,RCHK,list,gamma);     
	// update(N,Np,x,y,z,M,RCHK,list,gamma);
	printf("rmax=%f\n",rmax );
	// calc_force(x,y,z,N,Np,r1,kx,ky,kz,list,gamma,alpha);
      }
      
      for(i=0;i<Np;i++){
	x0[i]=x[i];
	y0[i]=y[i];
	kx[i]=0.0;
	ky[i]=0.0;
	vx[i]=0.0;
	vy[i]=0.0;
      }
      
      printf("Now Simulating!\n");
      srand((unsigned)time(NULL));
      ns=0;
      Nstep=0;
      Q=1;
      Uave=0.0;
      rfxyave=0.0;
      rfxyave2=0.0;
      gamma_dot2=0.0;
      lp=0;
      double TINT=0.0;
      double TQ=0.0;
      int Nout=0;
      double gamma=0.0;     
      double amp=0;
      gamma=0.0;
      gamma_dot=0.0;
      Fire(x,y,N,Np,kx,ky,gamma,list,r1,gamma_dot,RCHK,M);
      k=0;
      int countJ=0;
      delphi=1.e-4;
      // double P0=1.e-4;
      
      for(;;){
	k++;
	Psi=0.0;
	for(j=0;j<Np;j++)
	  Psi+=pi*(r1[j]*r1[j])/4.0/(N*N);
	
	FireSS(x,y,N,Np,kx,ky,&gamma,list,r1,0.0,RCHK,M); 
	//	output(oname,k,x,y,N,Np);
	calc_f(x,y,&U,&rfxy,&rf,N,Np,r1,list,gamma);

	////////////    
	set_map(x,y,r1,map,N,Np,gamma);
	
	Z0=0.0;
	for(j=0;j<Np;j++)
	  Z0+=map[j][0];
	///////
       	printf("sigma=%.30f,P=%.30f,U=%.30f,phi=%f,Z=%f\n",rfxy,rf,U,Psi,Z0/Np);
	if(Psi>atof(argv[1]))
	  delphi=-1.e-4;

	
	if(U<1.e-5 && delphi<0)
	  delphi=-1.e-6;
	
	if(rf<=P0 && delphi<0){
	  phiJ=Psi;
	  break;
	}
	
	N/=pow((Psi+delphi)/Psi,0.5);
	
	for(j=0;j<Np;j++){
	  x[j]=x[j]/pow((Psi+delphi)/Psi,0.5); 
	  y[j]=y[j]/pow((Psi+delphi)/Psi,0.5);  
	}
      }
    
      printf("phiJ=%f,k=%d\n",phiJ,k);
      phiJave+=phiJ;
      file2  <<" # phiJ= "<<  phiJave/(double)ens  << endl;
     
      int m=0;
      for(j=0;j<Np;j++){
	x0[j]=x[j]; 
	y0[j]=y[j]; 
      }
      FireSS(x,y,N,Np,kx,ky,&gamma,list,r1,0.0,RCHK,M); 
      Psi=0.0;
      for(j=0;j<Np;j++)
	Psi+=pi*(r1[j]*r1[j])/4.0/(N*N);
      
      calc_f(x,y,&U,&rfxy,&rf,N,Np,r1,list,gamma);
      printf("sigma=%.30f,P=%.30f,U=%.30f,phi=%f,Z=%f\n,delta phi=%f",rfxy,rf,U,Psi,Z0/Np,Psi-phiJ);
      
      gamma0=gamma;
      calc_f(x,y,&U,&rfxy,&rf,N,Np,r1,list,gamma);
      /////////////////////////////
      sigmaxy[m]+=rfxy;
      P[m]+=rf;
      E[m]+=U;
      ////////////////////////
      set_map(x,y,r1,map,N,Np,gamma);
      for(j=0;j<Np;j++){
	Psi1[m]+=pi*(r1[j]*r1[j])/4.0/(N*N);
	if(map[j][0]>DBL_EPSILON){
	  Z[m]+=map[j][0];
	  countz[m]++;	
	}
      }
      
      Psi1[m]-=phiJ;
      
      printf("initial:::sigma=%.30f,P=%.30f,phi=%f\n",rfxy,rf,Psi);
      sprintf(filename4,"%s/%s/stress_strain_%s_%d.dat",argv[1],argv[2],oname,ens);
      file3.open(filename4);
      file3  <<" # "<< "gamma" <<"   "<< "stress[m]" <<"   "<< "Pressure[m]"<<"   "<< "energy[m]"<<"   "<< "ZZ[m]/contactz[m]" << endl;
    
      N1=N;
      
      rfxy0[m][0]+=rfxy;
      rf0[m][0]+=rf;
      E0[m][0]+=U;

      ini_coord(x,y,x0,y0,Np);
      set_map(x,y,r1,map0,N1,Np,gamma);
      
      for(j=1;j<31;j++){
	gamma_dot=1.e-9*pow(10,0.2*(j-1));//20200210
	gamma+=gamma_dot;
	FireP(x,y,&N1,Np,kx,ky,gamma,list,r1,gamma_dot,RCHK,M,P0); 
	//	Fire(x,y,N,Np,kx,ky,gamma,list,r1,gamma_dot,RCHK,M);
	calc_f(x,y,&U,&rfxy,&rf,N1,Np,r1,list,gamma);
	rfxy0[m][j]+=rfxy;
	rf0[m][j]+=rf;
	E0[m][j]+=U;
	set_map(x,y,r1,map,N,Np,gamma);
	brokenbond(x,y,x0,y0,Np,N1,map0,BB,gamma,r1,m,j);
	for(l=0;l<Np;l++){
	  if(map[l][0]>DBL_EPSILON){
	    ZZ0[m][j]+=map[l][0];
	    countz0[m][j]++;
	  }
	}
		
	Psi2=0.0;

	//     calc_f(x,y,&U,&rfxy,&rf,N1,Np,r1,list,gamma);
	for(i=0;i<Np;i++)
	  Psi2+=pi*(r1[i]*r1[i])/4.0/(N1*N1);
	phiJgamma[j]+=Psi2;

      }
      
      for(j=31;j<500;j++){
	gamma_dot=1.e-3;
	gamma+=gamma_dot;
	FireP(x,y,&N1,Np,kx,ky,gamma,list,r1,gamma_dot,RCHK,M,P0); 
	//	Fire(x,y,N,Np,kx,ky,gamma,list,r1,gamma_dot,RCHK,M); 
	calc_f(x,y,&U,&rfxy,&rf,N1,Np,r1,list,gamma);
	brokenbond(x,y,x0,y0,Np,N1,map0,BB,gamma,r1,m,j);
	rfxy0[m][j]+=rfxy;
	rf0[m][j]+=rf;
	E0[m][j]+=U;
	set_map(x,y,r1,map,N,Np,gamma);
	for(l=0;l<Np;l++){
	  if(map[l][0]>DBL_EPSILON){
	    ZZ0[m][j]+=map[l][0];
	    countz0[m][j]++;	
	  }
	}
	//CG(x,y,N,Np,kx,ky,gamma,list,r1,RCHK,M);
	
	Psi2=0.0;
	for(i=0;i<Np;i++) 
	  Psi2+=pi*(r1[i]*r1[i])/4.0/(N1*N1);
	//	calc_f(x,y,&U,&rfxy,&rf,N1,Np,r1,list,gamma);
	
	phiJgamma[j]+=Psi2;
      }
      
      printf("deformed:::sigma=%.30f,P=%.30f,phi=%f\n",rfxy,rf,Psi);
      Gp[m]=(rfxy0[m][5]-rfxy0[m][0])/(5.*2.e-4); 
      
      for(j=0;j<Np;j++){
	x[j]=x0[j]; 
	y[j]=y0[j]; 
      }
      gamma=gamma0;
      
      file2.precision(10);  
      file2 << setprecision(10) << (double)m <<"  "<< Psi-phiJ <<"   "<< Z[m]/(countz[m]+DBL_EPSILON)  <<"   "<< P[m]/((double)ens)  <<"   "<< sigmaxy[m]/((double)ens)  <<"   "<< Gp[m]/((double)ens)  <<"   "<< E[m]/((double)ens) << endl;
      
      j=0;
      file3.precision(10);  
      file3 << setprecision(10)<<" #  "<< 0.0 <<"   "<< rfxy0[m][j]/((double)ens)  <<"   "<< rf0[m][j]/((double)ens) <<"   "<< E0[m][j]/((double)ens) <<"   "<< ZZ0[m][j]/(countz0[m][j]+DBL_EPSILON)  << endl;
      
      for(j=1;j<31;j++){
	file3.precision(10);  
	file3 << setprecision(10)  <<   1.e-9*(1.-pow(10,0.2*(j)))/(1.-pow(10,0.2))  <<"   "<< rfxy0[m][j]/((double)ens)  <<"   "<< rf0[m][j]/((double)ens) <<"   "<< E0[m][j]/((double)ens) <<"   "<< ZZ0[m][j]/(countz0[m][j]+DBL_EPSILON)<<"   "<< phiJgamma[j]/(ens+DBL_EPSILON)  <<"   "<< BB[m][j]/(ens+DBL_EPSILON)  << endl;
      }
      
      for(j=31;j<500;j++){
	file3.precision(10);  
	file3 << setprecision(10)  << 1.e-9*(1-pow(10,0.2*30.))/(1-pow(10,0.2)) + 1.e-3*(double)(j-30)  <<"   "<< rfxy0[m][j]/((double)ens)  <<"   "<< rf0[m][j]/((double)ens) <<"   "<< E0[m][j]/((double)ens) <<"   "<< ZZ0[m][j]/(countz0[m][j]+DBL_EPSILON)<<"   "<< phiJgamma[j]/(ens+DBL_EPSILON)<<"   "<< BB[m][j]/(ens+DBL_EPSILON)   << endl;
      }
      file3.close();   
      Gp[m]=0.0;
      m++;
    }	  
    file2.close();  
  }
  delete []list;
  return 0;
}
