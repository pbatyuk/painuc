#include "_ROOT.h"
#include "_STD_UNIX.h"
//#include "circle.h"
#include <cfortran.h>
#include <minuit.h>


extern float Xcircle[50000], Ycircle[50000], Zcircle[50000], Wpoint[50000];
extern float XC, YC, RC, AC, BC, chi2, DRC, DAC, phif, t[50000], bufmax, bufmin;
extern short int Signe;
//extern float XC1, YC1, RC1;
extern float fcn_print;
extern int Npoints_fit;
extern int error_flag, is;

// for MINUIT usage
/*void fcn(int npar, double grad[3], double * fcnval,
                double xval[3],int *iflag, void (*Dummy)());
void fcn2(int npar, double grad[5], double * fcnval,
              double xval[5],int *iflag, void (*Dummy)());
*/
void fcn(int npar, double grad[3], double * fcnval,
           double xval[3],int *iflag, void (*Dummy)())
{
  int i;
  double x0, y0, rad;
  //  double a0, b0;        // unused
  double x, y, w;           // current point ...........
  //  double z;             // unused
  double curr,sum2;
  double phi;
  //  double dx, dy, dz;    // unused
  double drho;
  
  x0 = xval[0];
  y0 = xval[1];
  rad = xval[2];

  sum2 = 0.;
  bufmin = 1000.;
  bufmax =-1000.;
  for (i=0; i<Npoints_fit; i++)
    {
      x =(double) Xcircle[i];
      y =(double) Ycircle[i];
      t[i] = atan2((y-y0),(x-x0));
      if(t[i]<0.) t[i] += TMath::Pi()*2.;
      if(t[i] < bufmin) bufmin = t[i];
      if(t[i] > bufmax) bufmax = t[i];
    }
  //check for the -Pi - +Pi bounds
  if((bufmax - bufmin) > TMath::Pi()) 
    {
      bufmin = 1000.;
      bufmax =-1000.;
      for(i=0;i<Npoints_fit; i++)
	{
	  if (t[i] < TMath::Pi()) t[i] += 2.*TMath::Pi();
	  if(*(t+i) < bufmin) bufmin = *(t+i);
	  if(*(t+i) > bufmax) bufmax = *(t+i);
	}
    }
  is = 1;
  if(bufmax - t[0] > t[0] - bufmin) is = -1;
  //if(t[0] - t[Npoints_fit/2] < 0.) is = -1;
  phif = bufmax;
  if(is == -1) phif = bufmin;
  if(is==Signe || Signe == 10)
    {
      for (i=0; i<Npoints_fit; i++)
	{
	  x =(double) Xcircle[i];
	  y =(double) Ycircle[i];
	  phi = t[i] - phif;
	  w =(double) Wpoint[i];
	  if (w > 0.)
	    {
	      drho = pow((sqrt(pow(x-x0,2) + pow(y-y0,2)) - rad),2);
	      curr = w*drho;
	      sum2 = sum2 + curr;
	    }
	}
      *fcnval = sum2;
    }
  else
    *fcnval = 1000000.;
  switch(*iflag){
  case 1:
    //   Initialaise.
    //printf(" fcn_c called to init..\n");
    break;
  case 2:
    // derivatives.....
    break;
  case 3:
    // exit
    {
      //printf(" is %d %f\n", is, phif);
      //*fcnval = sum2;
      fcn_print = sum2;
      //for (i=0; i<Npoints_fit; i++)
      //  printf("i %3d phi %f\n", i, t[i]);
    }
    break;
  default:
    break;
  }
}
void fcn2(int npar, double grad[5], double * fcnval,
           double xval[5],int *iflag, void (*Dummy)())
{
  int i;
  double x0, y0, rad, a0, b0;
  double x, y, z, w; // current point ...........
  double curr,sum2;
  double phi;
  //  double dx, dy;      // unused 
  double dz, drho;
  
  x0 = xval[0];
  y0 = xval[1];
  rad = xval[2];
  a0 = xval[3];
  b0 = xval[4];


  sum2 = 0.;
  bufmin = 1000.;
  bufmax =-1000.;
  for (i=0; i<Npoints_fit; i++)
    {
      x =(double) Xcircle[i];//*(xcircle+i);
      y =(double) Ycircle[i];//*(ycircle+i);
      t[i] = atan2((y-y0),(x-x0));
      if(t[i]<0.) t[i] += TMath::Pi()*2.;
      if(t[i] < bufmin) bufmin = t[i];
      if(t[i] > bufmax) bufmax = t[i];
    }
  //check for the -Pi - +Pi bounds
  if((bufmax - bufmin) > TMath::Pi()) 
    {
      bufmin = 1000.;
      bufmax =-1000.;
      for(i=0;i<Npoints_fit; i++)
	{
	  if (t[i] < TMath::Pi()) t[i] += 2.*TMath::Pi();
	  if(*(t+i) < bufmin) bufmin = *(t+i);
	  if(*(t+i) > bufmax) bufmax = *(t+i);
	}
    }
  is = 1;
  if(bufmax - t[0] > t[0] - bufmin) is = -1;
  //if(t[0] - t[Npoints_fit/2] < 0.) is = -1;
  phif = bufmax;
  if(is == -1) phif = bufmin;
  if(is==Signe)
    {
      for (i=0; i<Npoints_fit; i++)
	{
	  x =(double) Xcircle[i];//*(xcircle+i);
	  y =(double) Ycircle[i];//*(ycircle+i);
	  phi = t[i] - phif;
	  z =(double) Zcircle[i];
	  w =(double) Wpoint[i];//*(Wpoint+i);
	  if (w > 0.)
	    {
	      dz = z - a0*phi - b0;
	      drho = pow((sqrt(pow(x-x0,2) + pow(y-y0,2)) - rad),2)+dz*dz;
	      curr = w*drho;
	      sum2 = sum2 + curr;
	    }
	}
      *fcnval = sum2;
    }
  else
    *fcnval = 1.e+16;
  switch(*iflag){
  case 1:
    //   Initialaise.
    //printf(" fcn_c called to init..\n");
    break;
  case 2:
    // derivatives.....
    break;
  case 3:
    // exit
    {
      //printf(" is %d %f\n", is, phif);
      //*fcnval = sum2;
      fcn_print = sum2;
      //for (i=0; i<Npoints_fit; i++)
      //  printf("i %3d phi %f\n", i, t[i]);
    }
    break;
  default:
    break;
  }
}
void fcn_simple(int npar, double grad[3], double * fcnval,
           double xval[3],int *iflag, void (*Dummy)())
{
  double x0, y0, rad;
  double x, y, w; // current point ...........
  double curr,sum2;
  
  x0 = xval[0];
  y0 = xval[1];
  rad = xval[2];


  switch(*iflag){
  case 1:
    //   Initialaise.
    //printf(" fcn_c called to init..\n");
    int i;
    sum2 = 0.;
    for (i=0; i<Npoints_fit; i++)
      {
	x =(double) Xcircle[i];//*(xcircle+i);
	y =(double) Ycircle[i];//*(ycircle+i);
	w =(double) Wpoint[i];//*(Wpoint+i);
	if (w > 0.)
	  {
	    curr = w*pow((sqrt(pow(x-x0,2) + pow(y-y0,2)) - rad), 2.);
	    sum2 = sum2 + curr;
	  }
      }
    *fcnval = sum2;
    break;
  case 2:
    // derivatives.....
    //grad[0] =-pow((pow(x-x0,2)+pow(y-y0,2))/pow(rad,2.), -0.5)*(x-x0)/pow(rad,4.)*x0;
    //grad[1] =-pow((pow(x-x0,2)+pow(y-y0,2))/pow(rad,2.), -0.5)*(y-y0)/pow(rad,4.)*y0;
    //grad[3] =-pow((pow(x-x0,2)+pow(y-y0,2))/pow(rad,2.), -0.5)/pow(rad,3.);
    break;
  default:
    {
      int i;
      sum2 = 0.;
      for (i=0; i<Npoints_fit; i++)
	{
	  x =(double) Xcircle[i];//*(xcircle+i);
	  y =(double) Ycircle[i];//*(ycircle+i);
	  w =(double) Wpoint[i];//*(Wpoint+i);
	  if (w > 0.)
	    {
	      curr = w*pow((sqrt(pow(x-x0,2) + pow(y-y0,2)) - rad), 2.);
	      sum2 = sum2 + curr;
	    }
	}
      *fcnval = sum2;
      fcn_print = sum2;
    }
    break;
  }
}
//FCALLSCSUB6(fcn,FCN,fcn,INT,DOUBLEV,PDOUBLE,DOUBLEV,INT,ROUTINE);
//FCALLSCSUB6(fcn2,FCN2,fcn2,INT,DOUBLEV,PDOUBLE,DOUBLEV,INT,ROUTINE);
/*void fcn_simple(int npar, double grad[3], double * fcnval,
                double xval[3],int *iflag, void (*Dummy)());
FCALLSCSUB6(fcn_simple,FCN_SIMPLE,fcn_simple,INT,DOUBLEV,PDOUBLE,DOUBLEV,INT,ROUTINE);
*/
void circle_minuit (int N)
  //	      float *A, float *B, float *R, float *pchi)
{
  int error_buf;
  double a,b,c;
  double A,B,D;
  double x1,x2,y1,y2,x3, y3, yc[4], xc[4];
  double cosA, sinA, test0, test1, R;
  int imin;
  if (N ==3)  // analytically!
    {
      a = sqrt(pow(Xcircle[2]-Xcircle[0],2) + pow(Ycircle[2]-Ycircle[0],2));
      b = sqrt(pow(Xcircle[1]-Xcircle[0],2) + pow(Ycircle[1]-Ycircle[0],2));
      c = sqrt(pow(Xcircle[2]-Xcircle[1],2) + pow(Ycircle[2]-Ycircle[1],2));
      if(a*b*c == 0.)
	{
	  if(a==0.) cosA = (pow(b,2)+pow(c,2)-pow(a,2))/(2.*b*c);
	  if(b==0.) cosA = (pow(c,2)+pow(a,2)-pow(b,2))/(2.*c*a);
	  if(c==0.) cosA = (pow(a,2)+pow(b,2)-pow(c,2))/(2.*a*b);
	}
      else
	cosA = (pow(b,2)+pow(c,2)-pow(a,2))/(2.*b*c);
      if(fabs(cosA) > 1.) cosA = 1.;
      sinA = sqrt(1.-pow(cosA,2));
      RC = 0.;
      if(sinA != 0.) {RC = (float) (a/2./sinA); R = a/2./sinA;}
      if (RC != 0.)
	{
	  imin = 10;
	  y1 = (double) Ycircle[0];
	  y2 = (double) Ycircle[2];
	  x1 = (double) Xcircle[0];
	  x2 = (double) Xcircle[2];
	  x3 = (double) Xcircle[1];
	  y3 = (double) Ycircle[1];
	  if (x1 == x2) 
	    {
	      y1 = (double) Ycircle[0];
	      y2 = (double) Ycircle[1];
	      x1 = (double) Xcircle[0];
	      x2 = (double) Xcircle[1];
	      x3 = (double) Xcircle[2];
	      y3 = (double) Ycircle[2];
	      if(x1 == x2) {
		RC = 250.;
		XC = 200.;
		YC = 200.;
		goto MINUIT;}
	    }
	  A = (y2*y2 - y1*y1 + x2*x2 - x1*x1)/2./(x2-x1);
	  B = (y2-y1)/(x1-x2);
	  a = 1. + B*B;
	  b = -2.*((x1-A)*B+y1);
	  c = (x1-A)*(x1-A) + y1*y1 - R*R;
	  D = b*b - 4.*a*c;
	  yc[0] = (-b+sqrt(D))/2./a;
	  yc[1] = (-b-sqrt(D))/2./a;
	  xc[0] = A+B*yc[0];
	  xc[1] = A+B*yc[1];
	  test0 = pow(x3-xc[0],2) + pow(y3-yc[0],2) - R*R;
	  test1 = pow(x3-xc[1],2) + pow(y3-yc[1],2) - R*R;
	  if(fabs(test0) < fabs(test1)) imin = 0;
	  if(fabs(test0) > fabs(test1)) imin = 1;
	  //printf("test0 %f test1 %f\n %f %f  %f %f imin %d\n", 
	  // test0, test1, xc[0], yc[0], xc[1], yc[1], imin);
	  
	  if(imin != 10)
	    {
	      YC = (float) yc[imin];
	      XC = (float) xc[imin];
	      //	      printf("Analytically parameters of the circle are\n");
	      //printf("XC %f YC %f RC %f\n", XC, YC, RC);
	      chi2 = 0.;
	      return;
	    }
	}
      RC = 250.;
      XC = 200.;
      YC = 200.;
    }
 MINUIT:
  char name1[10],name2[10],name3[10];
  double val1, val2, val3;
  double eval1, eval2, eval3;
  double bnd1, bnd2;
  
  error_flag=0;
  int ivarbl;
  double arglst[2];
  double f_null =0.;
  MNINIT(5,16,7);
  //MNINIT(5,6,7);
  MNSETI((char *)" just fit ");
  //   set of parameters
  // printf("beg.par: XC %f YC %f RC %f\n", XC, YC, RC);
  MNPARM(1, (char *)"X", (double) XC, (double) (10.), f_null, f_null, error_flag);
  MNPARM(2, (char *)"Y", (double) YC, (double) (10.), f_null, f_null, error_flag);
  if (RC != 0.)
    MNPARM(3, (char *)"R", (double) RC, (double) (RC/10.), f_null, f_null, error_flag);
  else
    MNPARM(3, (char *)"R", (double) 1000., (double) (10.), f_null, f_null, error_flag);
  if(Signe != 10)
    {
      MNEXCM(fcn, (char *)"SET NOW", arglst, 0, error_flag, 0);  // now warnings
      arglst[0] = -1.;
      MNEXCM(fcn, (char *)"SET PRI", arglst, 1, error_flag, 0);  // set minimal output
      arglst[0] = 1.;
      MNEXCM(fcn, (char *)"CAL", arglst, 1, error_flag, 0);
      arglst[0] = 2.;
      MNEXCM(fcn, (char *)"SET STR", arglst, 1, error_flag, 0);
      // fitting algorithm
      arglst[0] = 100000.;
      arglst[1] = 0.00001;
      MNEXCM(fcn, (char *)"MIG", arglst, 2, error_flag, 0);
      //MNEXCM(fcn,"MIG",arglst,0,error_flag,0);
      MNEXCM(fcn, (char *)"MINOS", arglst, 0, error_flag, 0);
      
      printf("error flag is %d\n", error_flag); 
      // readout the result
      arglst[0] = 3.;
      MNEXCM(fcn, (char *)"CAL", arglst, 1, error_buf, 0);
    }
  else
    {
      MNEXCM(fcn_simple, (char *)"SET NOW", arglst, 0, error_flag, 0);  // now warnings
      arglst[0] = -1.;
      MNEXCM(fcn_simple, (char *)"SET PRI", arglst, 1, error_flag, 0);  // set minimal output
      arglst[0] = 1.;
      MNEXCM(fcn_simple, (char *)"CAL", arglst, 1, error_flag, 0);
      arglst[0] = 2.;
      MNEXCM(fcn_simple, (char *)"SET STR", arglst, 1, error_flag, 0);
      // fitting algorithm
      arglst[0] = 100000.;
      arglst[1] = 0.00001;
      MNEXCM(fcn_simple, (char *)"MIG", arglst, 2, error_flag, 0);
      //MNEXCM(fcn,"MIG",arglst,0,error_flag,0);
      MNEXCM(fcn_simple, (char *)"MINOS", arglst, 0, error_flag, 0);
      
      printf("error flag is %d\n", error_flag); 
    }
  MNPOUT(1,name1,val1,eval1,bnd1,bnd2,ivarbl);
  MNPOUT(2,name2,val2,eval2,bnd1,bnd2,ivarbl);
  MNPOUT(3,name3,val3,eval3,bnd1,bnd2,ivarbl);
  // printout
  //printf(" MINUIT circle output: %d points, chi2 %f \n", N, fcn_print/((float) N - 2));
  //printf("var - %s = %f, err= %f\n",name1,val1,eval1);
  //printf("var - %s = %f, err= %f\n",name2,val2,eval2);
  //printf("var - %s = %f, err= %f\n",name3,val3,eval3);
  //getc(stdin);
  if(error_flag == 0 || RC == 250.)
    {
      XC = val1;
      YC = val2;
      RC = val3;
      DRC = eval3;
    }
  chi2 = fcn_print/((float) N - 2);
}
void helix_minuit (int N)
  //	      float *A, float *B, float *R, float *pchi)
{
  // int imin;    // unused
  int error_buf;
  char name1[10],name2[10],name3[10], name4[10], name5[10];
  double val1, val2, val3, val4, val5;
  double eval1, eval2, eval3, eval4, eval5;
  double bnd1, bnd2;
  
  //double *pval1=&val1, *pval2=&val2, *pval3=&val3,
  //  *peval1=&eval1,*peval2=&eval2,*peval3=&eval3,
  // *pbnd1=&bnd1,*pbnd2=&bnd2;

  error_flag=0;
  int ivarbl;
  //int *pivarbl=&ivarbl;
  double arglst[2];
  //  void *fcn();
  double f_null =0.;
  MNINIT(5,16,7);
  MNSETI( (char *)" just fit ");
  // set of parameters
  //printf("beg.par: XC %f YC %f RC %f\n", XC, YC, RC);
  MNPARM(1, (char *)"X", (double) XC, (double) (1.), f_null, f_null, error_flag);
  MNPARM(2, (char *)"Y", (double) YC, (double) (1.), f_null, f_null, error_flag);
  MNPARM(3, (char *)"R", (double) RC, (double) (TMath::Min(RC/100.,50.)), f_null, f_null, error_flag);
  MNPARM(4, (char *)"A", (double) AC, (double) (TMath::Min(abs(AC)/1000.,100.)), f_null, f_null, error_flag);
  MNPARM(5, (char *)"B", (double) BC, (double) (1.), f_null, f_null, error_flag);
  MNEXCM(fcn2, (char *)"SET NOW", arglst, 0, error_flag, 0);  // now warnings
  arglst[0] = -1.;
  MNEXCM(fcn2, (char *)"SET PRI", arglst, 1, error_flag, 0);  // set minimal output
  arglst[0] = 1.;
  MNEXCM(fcn2, (char *)"CAL", arglst, 1, error_flag, 0);
  arglst[0] = 2.;
  MNEXCM(fcn2, (char *)"SET STR", arglst, 1, error_flag, 0);
  // fitting algorithm
  arglst[0] = 100000.;
  arglst[1] = 0.00001;
  MNEXCM(fcn2, (char *)"MIG", arglst, 2, error_flag, 0);
  //MNEXCM(fcn,"MIG",arglst,0,error_flag,0);
  MNEXCM(fcn2, (char *)"MINOS", arglst, 0, error_flag, 0);

  //printf("error flag is %d\n", error_flag); 
  arglst[0] = 3.;
  MNEXCM(fcn2, (char *)"CAL", arglst, 1, error_buf, 0);
  // readout the result
  MNPOUT(1,name1,val1,eval1,bnd1,bnd2,ivarbl);
  MNPOUT(2,name2,val2,eval2,bnd1,bnd2,ivarbl);
  MNPOUT(3,name3,val3,eval3,bnd1,bnd2,ivarbl);
  MNPOUT(4,name4,val4,eval4,bnd1,bnd2,ivarbl);
  MNPOUT(5,name5,val5,eval5,bnd1,bnd2,ivarbl);
  // printout
  //printf(" MINUIT helix output: %d points, chi2 %f \n", N, fcn_print/((float) N - 4));
  //printf("var - %s = %f, err= %f\n",name1,val1,eval1);
  //printf("var - %s = %f, err= %f\n",name2,val2,eval2);
  //printf("var - %s = %f, err= %f\n",name3,val3,eval3);
  //printf("var - %s = %f, err= %f\n",name4,val4,eval4);
  //printf("var - %s = %f, err= %f\n",name5,val5,eval5);
  //getc(stdin);
  if(error_flag == 0 || RC == 250.)
    {
      XC = val1;
      YC = val2;
      RC = val3;
      DRC = eval3;
      AC = val4;
      BC = val5;
      DAC = eval4;
    }
  chi2 = fcn_print/((float) N - 4);
}
void helix_minuit_fixR (int N)
  //	      float *A, float *B, float *R, float *pchi)
{
  // int imin;    // unused 
  int error_buf;
  char name1[10],name2[10],name3[10], name4[10], name5[10];
  double val1, val2, val3, val4, val5;
  double eval1, eval2, eval3, eval4, eval5;
  double bnd1, bnd2;
  
  //double *pval1=&val1, *pval2=&val2, *pval3=&val3,
  //  *peval1=&eval1,*peval2=&eval2,*peval3=&eval3,
  // *pbnd1=&bnd1,*pbnd2=&bnd2;

  error_flag=0;
  int ivarbl;
  //int *pivarbl=&ivarbl;
  double arglst[2];
  //  void *fcn();
  double f_null =0.;
  MNINIT(5,16,7);
  MNSETI( (char *)" just fit ");
  // set of parameters
  //printf("beg.par: XC %f YC %f RC %f\n", XC, YC, RC);
  MNPARM(1, (char *)"X", (double) XC, (double) (10.), f_null, f_null, error_flag);
  MNPARM(2, (char *)"Y", (double) YC, (double) (10.), f_null, f_null, error_flag);
  MNPARM(3, (char *)"R", (double) RC, (double) (0.), f_null, f_null, error_flag);
  MNPARM(4, (char *)"A", (double) AC, (double) (0.01), f_null, f_null, error_flag);
  MNPARM(5, (char *)"B", (double) BC, (double) (1.), f_null, f_null, error_flag);
  MNEXCM(fcn2, (char *)"SET NOW", arglst, 0, error_flag, 0);  // now warnings
  arglst[0] = -1.;
  MNEXCM(fcn2, (char *)"SET PRI", arglst, 1, error_flag, 0);  // set minimal output
  arglst[0] = 1.;
  MNEXCM(fcn2, (char *)"CAL", arglst, 1, error_flag, 0);
  arglst[0] = 2.;
  MNEXCM(fcn2, (char *)"SET STR", arglst, 1, error_flag, 0);
  // fitting algorithm
  arglst[0] = 100000.;
  arglst[1] = 0.00001;
  MNEXCM(fcn2, (char *)"MIG", arglst, 2, error_flag, 0);
  //MNEXCM(fcn,"MIG",arglst,0,error_flag,0);
  MNEXCM(fcn2, (char *)"MINOS", arglst, 0, error_flag, 0);

  //printf("error flag is %d\n", error_flag); 
  arglst[0] = 3.;
  MNEXCM(fcn2, (char *)"CAL", arglst, 1, error_buf, 0);
  // readout the result
  MNPOUT(1,name1,val1,eval1,bnd1,bnd2,ivarbl);
  MNPOUT(2,name2,val2,eval2,bnd1,bnd2,ivarbl);
  MNPOUT(3,name3,val3,eval3,bnd1,bnd2,ivarbl);
  MNPOUT(4,name4,val4,eval4,bnd1,bnd2,ivarbl);
  MNPOUT(5,name5,val5,eval5,bnd1,bnd2,ivarbl);
  // printout
  //printf(" MINUIT helix output: %d points, chi2 %f \n", N, fcn_print/((float) N - 4));
  //printf("var - %s = %f, err= %f\n",name1,val1,eval1);
  //printf("var - %s = %f, err= %f\n",name2,val2,eval2);
  //printf("var - %s = %f, err= %f\n",name3,val3,eval3);
  //printf("var - %s = %f, err= %f\n",name4,val4,eval4);
  //printf("var - %s = %f, err= %f\n",name5,val5,eval5);
  //getc(stdin);
  if(error_flag == 0 || RC == 250.)
    {
      XC = val1;
      YC = val2;
      RC = val3;
      DRC = eval3;
      AC = val4;
      BC = val5;
      DAC = eval4;
    }
  chi2 = fcn_print/((float) N - 4);
}
void circle_minuit_fixR(int N)
  //	      float *A, float *B, float *R, float *pchi)
{
  // double a,b,c;                             // unused
  // double A,B,D;                             // unused
  // double x1,x2,y1,y2,x3, y3, yc[4], xc[4];  // unused
  // double cosA, sinA, test0, test1, R;       // unused
  // int imin;                                 // unused
  char name1[10],name2[10],name3[10];
  double val1, val2, val3;
  double eval1, eval2, eval3;
  double bnd1, bnd2;
  
  error_flag=0;
  int ivarbl;
  double arglst[2];
  double f_null =0.;
  MNINIT(5,16,7);
  MNSETI( (char *)" just fit ");
  // set of parameters
  //printf("beg.par: XC %f YC %f RC %f\n", XC, YC, RC);
  MNPARM(1, (char *)"X", (double) XC, (double) (10.), f_null, f_null, error_flag);
  MNPARM(2, (char *)"Y", (double) YC, (double) (10.), f_null, f_null, error_flag);
  MNPARM(3, (char *)"R", (double) RC, (double) 0., f_null, f_null, error_flag);
  MNEXCM(fcn, (char *)"SET NOW", arglst, 0, error_flag, 0);  // now warnings
  arglst[0] = -1.;
  MNEXCM(fcn, (char *)"SET PRI", arglst, 1, error_flag, 0);  // set minimal output
  arglst[0] = 1.;
  MNEXCM(fcn, (char *)"CAL", arglst, 1, error_flag, 0);
  arglst[0] = 2.;
  MNEXCM(fcn, (char *)"SET STR", arglst, 1, error_flag, 0);
  // fitting algorithm
  arglst[0] = 100000.;
  arglst[1] = 0.00001;
  MNEXCM(fcn, (char *)"MIG", arglst, 2, error_flag, 0);
  //MNEXCM(fcn,"MIG",arglst,0,error_flag,0);
  MNEXCM(fcn, (char *)"MINOS", arglst, 0, error_flag, 0);

  printf("error flag is %d\n", error_flag); 
  // readout the result
  MNPOUT(1,name1,val1,eval1,bnd1,bnd2,ivarbl);
  MNPOUT(2,name2,val2,eval2,bnd1,bnd2,ivarbl);
  MNPOUT(3,name3,val3,eval3,bnd1,bnd2,ivarbl);
  // printout
  printf(" MINUIT circle output: %d points, chi2 %f \n", N, fcn_print/((float) N - 2));
  printf("var - %s = %f, err= %f\n",name1,val1,eval1);
  printf("var - %s = %f, err= %f\n",name2,val2,eval2);
  printf("var - %s = %f, err= %f\n",name3,val3,eval3);
  //getc(stdin);
  if(error_flag == 0 || RC == 250.)
    {
      XC = val1;
      YC = val2;
      RC = val3;
    }
  chi2 = fcn_print/((float) N - 2);
}
void circle_analytical(int j)
{
  double a,b,c;
  double A,B,D;
  double x1,x2,y1,y2,x3, y3, yc[4], xc[4];
  double cosA, sinA, test0, test1, R;
  //  double  Min_rc;    // unused
  //  int i;             // unused
  int imin;
  double Xc[3], Yc[3];
  Xc[0] = Xcircle[0];
  Xc[1] = Xcircle[j/2];
  Xc[2] = Xcircle[j-1];

  Yc[0] = Ycircle[0];
  Yc[1] = Ycircle[j/2];
  Yc[2] = Ycircle[j-1];
  a = sqrt(pow(Xc[2]-Xc[0],2) + pow(Yc[2]-Yc[0],2));
  b = sqrt(pow(Xc[1]-Xc[0],2) + pow(Yc[1]-Yc[0],2));
  c = sqrt(pow(Xc[2]-Xc[1],2) + pow(Yc[2]-Yc[1],2));
  cosA = (pow(b,2)+pow(c,2)-pow(a,2))/(2.*b*c);
  if(fabs(cosA) > 1.) cosA = 1.;
  sinA = sqrt(1.-pow(cosA,2));
  RC = 0.;
  if(sinA != 0.) {RC = (float) (a/2./sinA); R = a/2./sinA;}
  if (RC != 0.)
    {
      imin = 10;
      y1 = Yc[0];
      y2 =  Yc[2];
      x1 =  Xc[0];
      x2 =  Xc[2];
      x3 =  Xc[1];
      y3 =  Yc[1];
      if (x1 == x2) 
	{
	  y1 =  Yc[0];
	  y2 =  Yc[1];
	  x1 =  Xc[0];
	  x2 =  Xc[1];
	  x3 =  Xc[2];
	  y3 =  Yc[2];
	  if(x1 == x2) {
	    RC = 250.;
	    XC = 200.;
	    YC = 200.;}
	}
      A = (y2*y2 - y1*y1 + x2*x2 - x1*x1)/2./(x2-x1);
      B = (y2-y1)/(x1-x2);
      a = 1. + B*B;
      b = -2.*((x1-A)*B+y1);
      c = (x1-A)*(x1-A) + y1*y1 - R*R;
      D = b*b - 4.*a*c;
      yc[0] = (-b+sqrt(D))/2./a;
      yc[1] = (-b-sqrt(D))/2./a;
      xc[0] = A+B*yc[0];
      xc[1] = A+B*yc[1];
      test0 = pow(x3-xc[0],2) + pow(y3-yc[0],2) - R*R;
      test1 = pow(x3-xc[1],2) + pow(y3-yc[1],2) - R*R;
      if(fabs(test0) < fabs(test1)) imin = 0;
      if(fabs(test0) > fabs(test1)) imin = 1;
      printf("test0 %f test1 %f\n %f %f  %f %f imin %d\n", 
	     test0, test1, xc[0], yc[0], xc[1], yc[1], imin);
      
      if(imin != 10)
	{
	  YC = (float) yc[imin];
	  XC = (float) xc[imin];
	  printf("Analytically parameters of the circle are\n");
	  printf("XC %f YC %f RC %f\n", XC, YC, RC);
	  chi2 = 0.;
	  return;
	}
    }
  RC = 250.;
  XC = 200.;
  YC = 200.;
}

