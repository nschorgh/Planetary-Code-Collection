/* Creates Gaussian surface

   created:  2000 Norbert Schorghofer
   modified: 2015
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"

/* functions from Numerical Recipes */
float ran1(long *);
void rlft3(float ***data, float **speq, unsigned long nn1,
	   unsigned long nn2, unsigned long nn3, int isign);
float selip(unsigned long k, unsigned long n, float arr[]);

/* other functions */
void derivatives(float **h, float **hx, float **hy, int N1, int N2);
void output(float **h, FILE *tmp, int N1, int N2);
float meansquare(float **h);

const int N=128;

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


void derivatives(float **h, float **hx, float **hy, int N1, int N2) {
  int i,j,ip1,im1,jp1,jm1;
  const float dx=1., dy=1.;
  
  if (N1!=N2) {
    printf("error N1!=N2\n"); 
    exit(1);
  }
  for(i=1;i<=N1;i++) {
    ip1=i+1; im1=i-1; 
    for(j=1;j<=N2;j++) {
      jp1=j+1; jm1=j-1; 
      if (im1<1) im1+=N1;
      if (ip1>N1) ip1-=N1;
      if (jm1<1) jm1+=N2;
      if (jp1>N2) jp1-=N2;
      hx[i][j]=(h[ip1][j]-h[im1][j])/(2*dx);
      hy[i][j]=(h[i][jp1]-h[i][jm1])/(2*dy); 
    }
  }
}

void output(float **h, FILE *tmp_f, int N1, int N2) {
  int i,j;
  for(i=1;i<=N2;i++) {
    for(j=1;j<=N1;j++) 
      fprintf(tmp_f,"%g ",h[j][i]);
    fprintf(tmp_f,"\n");
  }
}


float energy_r(float ***data, int N) {
  // energy in real space
  int i,j;
  double en=0.;
  for(i=1;i<=N;i++) for(j=1;j<=N;j++)
    en+=SQR(data[1][i][j]);
  return((float)(en/(N*N)));
}

float energy_f(float ***data, int N) {
  // energy in fourier space
  int i,j;
  double en=0.;
  for(i=1;i<=N;i++) for(j=1;j<=N;j++)
    en+=2.*SQR(data[1][i][j]);
  for(i=1;i<=N;i++) //j=1,2
    en-=SQR(data[1][i][1])+SQR(data[1][i][2]);
  return((float)(en/(N*N)));
}


void create_random_landscape(float **h, long *idum, float hurst) {
  // creates Gaussian surface
  float **speq,***data,freq,buf,am,ph;
  int i,j;
  const float pi=3.1415926535;
  const int MAXFREQ=(N/2)/2;
  const int MINFREQ=4.;
  
  printf("# ... creating landscape ...\n");
  data=f3tensor(1,1,1,N,1,N);
  speq=matrix(1,1,1,N<<1);
  for(i=1;i<=N;i++) for(j=1;j<=N;j++) data[1][i][j]=0.;

  for(i=1;i<=N;i++) for(j=1;j<=N/2;j++) {
    am=(float)N;
    ph=2.*pi*ran1(idum);
    data[1][i][2*j-1]=am*cos(ph);
    data[1][i][2*j]=am*sin(ph);
  } 
  /* fix special symmetries */
  for(i=2;i<=N/2;i++) { 
    data[1][N-i+2][1]=data[1][i][1]; 
    data[1][N-i+2][2]=-data[1][i][2];
  }  
  data[1][1][1]=0.; data[1][1][2]=0.;  /* constant part */
  for(j=1;j<=N;j++) data[1][N/2+1][j]=0.; /* Nyquist */
  for(j=2;j<=N/2;j++) {   // i=1
    freq=(float)((j-1)*(j-1));
    buf=pow(freq,1.+hurst);
    //buf/=1.-0.9*freq/(MAXFREQ*MAXFREQ);
    //buf=exp(sqrt(freq)/64.)*sqrt(freq);
    data[1][1][2*j-1]/=buf; 
    data[1][1][2*j]/=buf;
    if (freq>(float)(MAXFREQ*MAXFREQ)) {
      data[1][1][2*j-1]=0.; 
      data[1][1][2*j]  =0.;
    }
  }
  for(j=1;j<=N/2;j++) {
    for(i=2;i<=N/2;i++) {
      freq=(float)((i-1)*(i-1)+(j-1)*(j-1));
      buf=pow(freq,1.+hurst);
      //buf/=1.-0.9*freq/(MAXFREQ*MAXFREQ);
      //buf=exp(sqrt(freq)/64.)*sqrt(freq);
      data[1][i][2*j-1]/=buf; data[1][N-i+2][2*j-1]/=buf;
      data[1][i][2*j]/=buf;   data[1][N-i+2][2*j]/=buf;
      if (freq>(float)(MAXFREQ*MAXFREQ)) {
	data[1][i][2*j-1]=0.; data[1][N-i+2][2*j-1]=0.;
	data[1][i][2*j]  =0.; data[1][N-i+2][2*j]=  0.;
      } 
    }
  } 

  // erase long wavelengths
  for(j=2;j<=N/2;j++) {   // i=1
    freq=(float)((j-1)*(j-1));
    if (freq<MINFREQ*MINFREQ) {
      data[1][1][2*j-1]=0.; 
      data[1][1][2*j]  =0.;
    }
  }
  for(j=1;j<=N/2;j++) {
    for(i=2;i<=N/2;i++) {
      freq=(float)((i-1)*(i-1)+(j-1)*(j-1));
      if (freq<MINFREQ*MINFREQ) {
	data[1][i][2*j-1]=0.; data[1][N-i+2][2*j-1]=0.;
	data[1][i][2*j]  =0.; data[1][N-i+2][2*j]=  0.;
      } 
    }
  } 


  printf("# Energy: %g\n",energy_f(data,N)); 
  for (i=1;i<=N<<1;i++) speq[1][i]=0.;
  rlft3(data,speq,1,N,N,-1);
  for(i=1;i<=N;i++) for(j=1;j<=N;j++) data[1][i][j]/=N/2;
  printf("# Energy: %g\n",energy_r(data,N)); 

  for(i=1;i<=N;i++) for(j=1;j<=N;j++) h[i][j]=data[1][i][j];
  free_matrix(speq,1,1,1,N<<1);
  free_f3tensor(data,1,1,1,N,1,N);
}


float meansquare(float **h) {
  int i,j;
  float sum=0.;
  for(i=1;i<=N;i++)
    for(j=1;j<=N;j++)
      sum = sum + h[i][j]*h[i][j];
  sum /= N*N;
  return(sum);
}


int count_extrema(float **h) {
  int i,j,nr_ext=0,nr_max=0,nr_min=0;
  float pt;
  
  for(i=1;i<=N;i++) for(j=1;j<=N;j++) {
    pt=h[i][j];
    if (h[i-1][j]>pt && h[i+1][j]>pt && h[i][j-1]>pt && 
	h[i][j+1]>pt && h[i-1][j-1]>pt && h[i+1][j-1]>pt && 
	h[i-1][j+1]>pt && h[i+1][j+1]>pt) {
      nr_min++;
      //printf("minimum: %d %d\n",i,j);
    }
    if (h[i-1][j]<pt && h[i+1][j]<pt && h[i][j-1]<pt && 
	h[i][j+1]<pt && h[i-1][j-1]<pt && h[i+1][j-1]<pt && 
	h[i-1][j+1]<pt && h[i+1][j+1]<pt) {
      nr_max++;
      //printf("maximum: %d %d\n",i,j);
    }
  }
  printf("%d minima,  %d maxima\n",nr_min,nr_max);
  nr_ext=nr_min+nr_max;
  return(nr_ext);
}


void structure(float **h) {
  // calculate structure function
  int i,j,k,ipR,jpR,power,R=1;
  float s;
  
  power=(int)floor(log((float)N)/log(2.));
  for(k=1;k<=power;k++) {
    s=0.;
    for(i=1;i<=N;i++) for(j=1;j<=N;j++) {
      ipR=i+R;
      if (ipR>N) ipR-=N;
      s+=SQR(h[ipR][j]-h[i][j]);

      jpR=j+R;
      if (jpR>N) jpR-=N;
      s+=SQR(h[i][jpR]-h[i][j]);
    }
    s=s/(2*N*N);
    printf("%d %g\n",R,s);
    R*=2;
  }
}

/*
  structure function <|h(x)-h(x+R)|^2> ~ R^(2*rho)  0<rho<1
  spectral amplitude ~ 1/|k|^(1+rho) 
  rho = Hurst exponent
*/


int main() {
  const float hurst=0.9;   // H>=0
  int i,j;
  float **h,**hx,**hy,RMS=0.15,oldrms;
  long idum=-23460;
  FILE *info_f,*land_f;
  
  land_f=fopen("land.dat","wt"); 
  h=matrix(0,N+1,0,N+1);
  hx=matrix(1,N,1,N); hy=matrix(1,N,1,N);
  info_f=fopen("info.dat","wt");
  fprintf(info_f,"N=%d\nH=%g\nseed=%ld\nrms=%f\n\n",N,hurst,idum,RMS);
  fclose(info_f);

  create_random_landscape(h,&idum,hurst);

  derivatives(h,hx,hy,N,N);
  oldrms = sqrt(meansquare(hx)+meansquare(hy));
  //printf("%f %f\n",meansquare(h),oldrms);
  for(i=1;i<=N;i++)
    for(j=1;j<=N;j++)
      h[i][j]*=RMS/oldrms;
  derivatives(h,hx,hy,N,N);
  printf("%f %f\n",meansquare(h),sqrt(meansquare(hx)+meansquare(hy)));

  output(h,land_f,N,N);
  count_extrema(h);
  //structure(h);

  fclose(land_f);
  return(0);
}

