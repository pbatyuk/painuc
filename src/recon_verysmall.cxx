//  Reading of images from data file of image, filling data plot.

#include <dirent.h>
#include <unistd.h>
#include "_ROOT.h"
#include "_STD_UNIX.h"

/*
	Program GEOMETRY for establishing geometry in DSCS

		 (Dubna  Streamer Chamber Spectrometer):
				 S T R E A M E R

*/

#define ESC 0x1b        /* Define the escape key         */
#define TRUE  1         /* Define some handy constants   */
#define FALSE 0         /* Define some handy constants   */
#define PI  3.14159     /* Define a value for PI         */
#define ON  1           /* Define some handy constants   */
#define OFF 0           /* Define some handy constants   */
#define FOCUS 16.       /* focus of lens for videocamera */

#define OUT_filename "dubtoev.dat"

#define SCALE_BG 3.
#define R_MEAN 1120.;//1070.
#define P_MEAN 210.//214;

typedef unsigned* PUINT;
typedef PUINT* PPUINT;

//extern PPUINT Alloc_Buffer(int hsize, int vsize);
//extern void Free_Buffer(PPUINT buffer[2], int vsize);

extern PPUINT pCDBuffer[2];
extern TObjectTable *fObjTbl;                                                          
extern int np[2][6];
extern float rc_int, xc_int, yc_int;
extern float xtrack[2][6][10000], ytrack[2][6][10000];
extern float Xcircle[50000], Ycircle[50000], Zcircle[50000], Wpoint[50000];
extern float XC, YC, RC, AC, BC, chi2, DRC, DAC, phif, t[50000], bufmax, bufmin;
extern int Npoints_fit;
extern Float_t x_mark[50000], y_mark[50000];
extern int error_flag, is;
extern int   particle[10];
extern Float_t *px_mark, *py_mark;
extern TPolyMarker *markp, *markp2;
extern float Xvertax, Yvertax, Dvertax;
extern float brrr[2][6];
extern char eventname[20];
extern char gifname[100];
extern TCanvas *c2,*ch,*c3, *c4, *cp, *cw, *csl;
extern float stx, sty, endx, endy;
extern TGFileInfo fi;                                                                  
extern TPad *pleft, *pright, *ptmp, *pMsg;
extern TPad *hp0;
extern TPad *plc3, *prc3, *pMc3;
extern TPad *plc4, *prc4, *pMc4;
extern TPad *cp0, *cp1, *hp0, *cw0, *csl0;
extern TProfile *hprof, *hprof2, *hprof3, *hprof4;
extern TPaveLabel *pl, *pl1;
extern TGraph *gr;
extern TPaveLabel *pl11, *pl2, *ptext;
extern TMarker *mark, *mark2;
extern TPolyMarker *markp, *markp2;
extern TArc *myarc;
extern char *tokenPtr;
extern TH2F *left, *right;
extern short int Signe;
extern TH1F *hwidth, *hslice;

int xtrack1[2][6][10000], ytrack1[2][6][10000];
float angle[2][6];
float xcentr, ycentr;
float deltaV, disp;

//Minuit's variables

//for Gil's graphics
//TText *t;
//TH1D *hprof;
//

/*                                                */
/**************************************************/
/*                                                */
/*         Declare general variables for          */
/*                                                */
/*               S T R E A M E R                  */
/*                                                */
/*                   OPTICS                       */
/*                                                */
/**************************************************/
/*                                                */
float ztop, zbot, zmed, base, focus, zfilm[2], depth, dglass, S;
float d_fiduc_top[6][6], d_fiduc_bot[4][4];
float d_top_leftflm[6][6], d_bot_leftflm[4][4];
float d_top_rghtflm[6][6], d_bot_rghtflm[4][4];
float mag[2], mag_top[2], mag_bot[2], bm[2], cm[2];
int top, bot;
float mag_vert[2], mag_dir[2][10];
float xaxis[2], yaxis[2];
double xaxis_draw[2], yaxis_draw[2];
float  xfilm_top[2][6],  xfilm_bot[2][4];
float  yfilm_top[2][6],  yfilm_bot[2][4];
double  xfilmd_top[2][6],  xfilmd_bot[2][4];
double  yfilmd_top[2][6],  yfilmd_bot[2][4];
float  xfilm_top0[2][6],  xfilm_bot0[2][4];
float  yfilm_top0[2][6],  yfilm_bot0[2][4];
float  rfilm_top[2][6],  rfilm_bot[2][4];
float  xtrue_top[2][6],  xtrue_bot[2][4];
float  ytrue_top[2][6],  ytrue_bot[2][4];
float xchamb_top[2][6],  xchamb_bot[2][4];
float ychamb_top[2][6],  ychamb_bot[2][4];
float rchamb_top[2][6],  rchamb_bot[2][4];
int view, event_type, track, tracks;
float CCDBuffer[2][6][10001];
float xax_pxl[2], yax_pxl[2];
float xprong[2][6][10000], yprong[2][6][10000];
float xvertex[3], yvertex[3], zvertex, dip[10], Radius[3][10];
float xvertex_r[3], yvertex_r[3], zvertex_r;
float xstart[3][10], ystart[3][10], xend[3][10], yend[3][10];
float xcenter[3][10], ycenter[3][10], zcenter[10];
float xcenter_r[3][10], ycenter_r[3][10], zcenter_r[10];
float cosx[3][10], cosy[3][10], cosz[10], zstart[10], zend[10];
float cost0[2][10], sint0[2][10], tgt0[2][10], tgalpha[2][10];
float cosbeta[2][10], sinbeta[2][10], cosmu[2][10], sinmu[2][10];
float R0[2][10], Rvert[2][10], tgmu[2][10], tgtheta[10];
float stangle[3][10], endangle[3][10], dangle[3][10], range[3][10];
float vtxangle[3][10], vtxrange[3][10];
float mtrack[3][10], ntrack[3][10], qtrack[3][10];
float hvertex, hdir[10], xdir[3][10], ydir[3][10], zdir[10];
float cosx_rho[3][10], cosy_rho[3][10], tg_rho[3][10];
float rho[3][10], rho_meas[3][10], xcenter0[3][10], ycenter0[3][10];
float xcenter01[3][10], ycenter01[3][10];
float Radius0[3][10];
float mass[10], E[10], p[10], T[10], pincid, Etot, ptot;
float exp_lens, Xvertex[2], Yvertex[2], Dvertex[2];
float direct[3][10], rho_R0[2][10];
float tanhalfC[2][10], tanC[2][10], tanA, tanB;
float xcprime[3][10], ycprime[3][10], xvprime[3], yvprime[3];
int iangles;
float deltar[10], deltarho[4][4][10], rho_aux[2][4][10];
const float pixel2mm  = 0.0068;
float x3D[6][50000], y3D[6][50000], z3D[6][50000], phi[50000], 
  XCD[6], YCD[6], RCD[6];
float xcD[2][10], ycD[2][10];
float rcD[2][10], drcD[2][10];
float a[2][10], b[2][10], da[2][10];
int isig[2][10];
float xinter[10], yinter[10], zinter[10];
int   counter[10];
char *Part[6] = {"pi", "p", "d", "t", "He3", "He4"};
short int charge[6] = {1, 1, 1, 1, 2, 2};
float mass_part[6] = {139.5675, 938.2795, 1875.6278, 2808.9433, 2808.4138, 3727.409};

float Loss[10], losstag;
short int change_part;
int  isig_tr[10];
float phi_first[2][10], z_first[2][10];
float wes[10][50000];
int   nright[10];
float xcross3D[20], ycross3D[20], zcross3D[20];
float A_lineD[10], B_lineD[10];
float *x3d, *y3d, *z3d, *phid, *A_line, *B_line, *chi2_2, *DA_line;
int   npt3D[10], npt3D_left[10], npt3D_right[10];
float aver_l, aver_r;
double xtransfer, ytransfer, ytransfer2, xtransfer_err, ytransfer_err, ytransfer2_err;
float deltapx, deltapy, deltapz, deltaT, deltaE, deltaP;
int ierr;
float *chi2_1;
float chi2_line[2][10], chi2_circle[2][10], chi_tr[10];
float p2[10], T2[10], E2[10], m2[10], cosdir_x[2][10], cosdir_y[2][10], cosdir_z[2][10];
float cosxD[10], cosyD[10], coszD[10];
float lxy[10], L[10];
float dp[10];
//char ofilename[100];
double IntBr[2][6];
FILE *for_reconstr;
FILE *print_results;
/*                                                */
/*               Function prototypes              */
/*                                                */

void line(int N, float *X, float *Y, float *A, float *B, float *DA, float *VAR);
void XY_Vertex_Net(int NTracks, float *a, float *b, float *r, 
		   float *Xvertex, float *Yvertex, 
		   float *Dvertex);
void XY_Vertex(void);
void AreCross(float, float, float,
	      float, float, float);
void Compare(float, float);
void Apollonius(void);
void SameDistCirc(int i, int j);
void Save_Results(void);
//void Kinematics_Gil(void);
void Set_data(void);
Bool_t ChangePart(void);

void circle_minuit (int N);
void circle_minuit_fixR (int N);
void helix_minuit (int N);
void helix_minuit_fixR (int N);

/* void Initial_coordinates(  );   */
/* void Magnification(  );         */
float S_triangle(float a, float b, float c, float mmm, float nnn,
		float hhh, float *mm_tri);
float det4(float *a4);
float det5(float *a5);
void read_ascii(void);
void read_data(void);
void TrueCoordinates(void);
void ReconstructProper(void);
void Kinematics(void);
void Angles(void);
void Rho(void);
void Initialize(void);
int  gprintf(int *xloc, int *yloc, char *fmt, ... );

/*							*/
/*	Begin main function				*/
/*							*/

//extern unsigned _stklen = 16384;
float ratio;
extern int vsize;
extern int hsize;
float Zvertax;
//float XvertInTarg, YvertInTarg, DvertInTarg;
//unsigned short int Cross;
//float X1, Y1, X2, Y2;
//float a1[3][3], b1[3][3], r1[3][3], a2[3][3], b2[3][3], r2[3][3];
float Zv[10];

float xx,yy;
int fp, lp, mp;

float p0[1000],p1[1000];
float XFit[1000],YFit[1000],WFit[1000];
int NFit;
float tracklength[2][6];//track length  on the image in PXLs
float rompture[2][6];
float abright[2][6]; // tangent of fit of the brightness of track along its track 
//(background substracted)
float fx, fy, lx, ly;

float rel_loss[10], Mpr_trit, ang3_12, angprtr, excess_m;

void recon_main(char filename[], int ArrSize)
  //int main(int argc, char* argv[])
//int gil_main(char ofilename[], int ArrSize)
{
  //fObjTbl = new TObjectTable(10000);
  
  int i, ii, j, k, l, kk;
  char ifilename[100];
  char ofilename[100];
  
  float pi;
  
  strcpy(ifilename, filename);
  //deltaV = atof(argv[2]); 
  //TApplication theApp("App", &argc, argv);

  exp_lens = 1.0;
  
  printf("Now open input file %s\n", ifilename);
  if ((for_reconstr = fopen(ifilename,"rt")) == NULL)
    {
      fprintf(stderr, "Cannot open input file %s\n", ifilename);
      return ;
    }
  char dirout[10] = {"res/"};
  char dirin[10] = {"dig/"};
  strcpy(ofilename, dirout);
  printf ("string lenth %d\n",(int)strlen(dirout)); 
  k = (int) strlen(dirin);
  for (i = (int) strlen(dirout); i < (int) strlen(dirout) +
	 (int)strlen(ifilename) - (int) strlen(dirin); i++)
    {
      ofilename[i] = ifilename[k];
      k++;
    }
  for (i = (int) strlen(dirout) + (int) strlen(ifilename) - (int) strlen(dirin); 
       i < 100; i++)
    ofilename[i] = '\0';
  strcat(ofilename,".r0");

  printf("Now open output file %s\n", ofilename);
  if ((print_results = fopen(ofilename,"w")) == NULL)
    { 
      fprintf(stderr,"\n Cannot open results file %s \n",ofilename);
      exit(1);
    }
  read_data();
  fclose(for_reconstr);
  //double base;
  pi = PI;                    focus = FOCUS;
  //depth = 160. + 2./3.*8.;               //base = 253.031682; //253.7; //
  depth = 160. + 1./3.*8.;               //base = 253.031682; //253.7; //
  base = 253.;//251.;
  //mag_top[0] = 69.349327;         mag_top[1] = 69.639442;
  //mag_bot[0] = 79.141396;         mag_bot[1] = 79.486115;
  
  mag_top[0]  = 69.824161; mag_top[1] = 70.107374;
  mag_bot[0]  = 79.184760; mag_bot[1] = 79.505941;
  //phi[0]  = -0.000149; phi[1] =  -0.000173 ;
  //xtransfer = -0.000115;
  zfilm[0] = depth/(mag_bot[0] - mag_top[0]);
  zfilm[1] = depth/(mag_bot[1] - mag_top[1]);
  ztop = (mag_top[0]*zfilm[0] + mag_top[1]*zfilm[1])/2.;
  zbot = (mag_bot[0]*zfilm[0] + mag_bot[1]*zfilm[1])/2.;
  zmed = 0.5*(ztop + zbot);
  printf("zfilm 0 %7.5f 1 %7.5f mm   ztop = %7.2f mm   zbot = %7.2f mm\n",
	 zfilm[0], zfilm[1], ztop, zbot);
  printf("%f mag_top[0]*zfilm[0] %f mag_top[1]*zfilm[1]\n", mag_top[0]*zfilm[0], mag_top[1]*zfilm[1]);
  printf("%f mag_bot[0]*zfilm[0] %f mag_bot[1]*zfilm[1]\n", mag_bot[0]*zfilm[0], mag_bot[1]*zfilm[1]);


  top = 1; bot = 3;
  double ztop_diff, zbot_diff, ztop_diff0;
  ztop_diff = 0.;
  do
    {
      ztop_diff0 = ztop_diff;
      mag_top[0] = (mag_top[0] + mag_top[1]*mag_bot[0]/mag_bot[1])/2.;
      mag_top[1] = (mag_top[1] + mag_top[0]*mag_bot[1]/mag_bot[0])/2.;
      mag_bot[0] = (mag_bot[0] + mag_bot[1]*mag_top[0]/mag_top[1])/2.;
      mag_bot[1] = (mag_bot[1] + mag_bot[0]*mag_top[1]/mag_top[0])/2.;
      zfilm[0] = depth/(mag_bot[0] - mag_top[0]);
      zfilm[1] = depth/(mag_bot[1] - mag_top[1]);
      ztop_diff = fabs(mag_top[0]*zfilm[0] - mag_top[1]*zfilm[1]);
      zbot_diff = fabs(mag_bot[0]*zfilm[0] - mag_bot[1]*zfilm[1]);
      printf("%18.10f mag_top[0] %18.10f mag_top[1]\n", mag_top[0], mag_top[1]);
      printf("%18.10f mag_bot[0] %18.10f mag_bot[1]\n", mag_bot[0], mag_bot[1]);
      printf("ztop_diff %f zbot_diff %f \n\n", ztop_diff, zbot_diff);
    }
  while(ztop_diff > 0.00001 && ztop_diff > ztop_diff0);
  zfilm[0] = depth/(mag_bot[0] - mag_top[0]);
  zfilm[1] = depth/(mag_bot[1] - mag_top[1]);
  ztop = (mag_top[0]*zfilm[0] + mag_top[1]*zfilm[1])/2.;
  zbot = (mag_bot[0]*zfilm[0] + mag_bot[1]*zfilm[1])/2.;
  zmed = 0.5*(ztop + zbot);
  printf("zfilm 0 %7.5f 1 %7.5f mm   ztop = %7.2f mm   zbot = %7.2f mm\n",
	 zfilm[0], zfilm[1], ztop, zbot);
  
  /*                                                            */
  /*  Input (define or calculate or whatever...) the primarily  */
  /*  set data (distances between fiducials, coordinates of     */
  /*  optical axes, and so on.                                  */
  /*                                                            */
  for(i=0; i<=1; i++)
    for(j=0; j<6; j++)
      for(k=0; k<np[i][j]; k++)
	{
	  xprong[i][j][k] = xtrack[i][j][k];
	  yprong[i][j][k] = ytrack[i][j][k];
	}
  
  Set_data();
  
  TrueCoordinates();
  
  for(iangles=0; iangles<4; iangles++)
    {
      Angles();
      ReconstructProper();
    }
  Rho();
  
  Kinematics();
  
  ratio = (float) vsize/(float) hsize;
  printf("vsize = %d, hsize = %d, ratio = %f\n",vsize,hsize,ratio);
  
  float yleft, yright, ysize, xsize, xleft, xright;
  float ylefttmp, yrighttmp, xlefttmp, xrighttmp; 
  int ixWindow, iyWindow;
  float xWindow, yWindow, Wratio;
  float rhsize, rvsize;
  rhsize = (float) hsize;//pixel2mm;
  rvsize = (float) vsize;//pixel2mm;
  xWindow = 2.*rhsize + 4.*rhsize/100.;
  yWindow = 1.*rvsize + 2.*rvsize/100.;
  Wratio = yWindow/xWindow;
  if (Wratio <= 1.0)
    {
      xleft  = rhsize/100./xWindow;
      xright = 0.5 - xleft;
      yleft  = rvsize/100./yWindow;
      yright = 1.0 - rvsize/100./yWindow;
      ysize  = yright - yleft;
      xsize  = xright - xleft;

      xWindow = 1000;
      yWindow = xWindow*Wratio;
      ixWindow = (int) xWindow;
      iyWindow = (int) yWindow;
      printf("ratio = %f,%f,%f,%f,%f\n",Wratio,xleft,xright,yleft,yright);
      
    }
  else
    {
      xleft  = 0.05/ratio;
      xright = 0.45/ratio;
      yleft  = 0.55;
      yright = 0.95;
      ysize  = 0.40;
      xsize  = 0.40/ratio;
      yWindow = 600.;
      xWindow = yWindow/ratio;
      ixWindow = (int) xWindow;
      iyWindow = (int) yWindow;
      xlefttmp = 0.5 - 0.225;
      xrighttmp = 0.5 + 0.225;
      ylefttmp = 0.0;
      yrighttmp = 1. - yleft;   // 0.45;
    }
  
  //cout <<"x , y ="<<ixWindow<<" "<<iyWindow<<endl;
  
  delete gROOT->GetListOfCanvases()->FindObject("c3");
  c3 = new TCanvas("c3","Kinematics parameters",20,600,ixWindow/4.,iyWindow/2.);
  //pl11 = new TPaveLabel(-0.95,0.1,-0.65,0.9,ifilename,"br");
  c3->SetFillColor(42);
  //pl11->Draw();
  char PStr[60];
  TPaveText *text = new TPaveText(0.05, 0.05, 0.95, 0.95 );
  text->SetTextSize(0.2);
  //text->SetFillColor( 42 );
  for(track=0; track<tracks; track++)
    {
      sprintf( PStr, "p%d = %8.1f", 
  	       track, p[track]); 
     text->AddText(PStr);
    }
  text->Draw();
  c3->Modified();
  c3->Update();
  
  for (i = 0; i < 20; i++)
    {
      if(eventname[i]=='\n') eventname[i]='\0'; 
    }
  strcpy(gifname,"gif/");
  strcat(gifname, eventname);
  strcat(gifname, "c3.gif");
  c3->Print(gifname,"gif");
  strcpy(gifname,"gif/");
  strcat(gifname, eventname);
  strcat(gifname, "c2.gif");
  c2->Print(gifname,"gif");
  
  //theApp.Run(kTRUE);
  //sleep(10);
  return ;
}
/*                                                            */
/*  Input (define or calculate or whatever...) the primarily  */
/*  set data (distances between fiducials, coordinates of     */
/*  optical axes, and so on.                                  */
/*                                                            */
void Set_data()
{
  int i, k;
  float *dtop, *dbot;

  dtop = d_fiduc_top[0];          dbot =  d_fiduc_bot[0];

  /**************************************************************/
  /*                                                            */
  /*  Input coordinates of fiducials and optical axes on CCD    */
  /*  matrices ("film") in pixels for calculating various       */
  /*  distances. Then find correspondence between the fiducials */
  /*  on "film" and in chamber space...                         */
  /*                                                            */
  /**************************************************************/
  
  xfilm_top[0][0] = 198;          yfilm_top[0][0] = 1035. - 797;
  xfilm_top[1][0] = 248;          yfilm_top[1][0] = 1035. - 837;
  xfilm_top[0][1] = 187;          yfilm_top[0][1] = 1035. - 271;
  xfilm_top[1][1] = 244;          yfilm_top[1][1] = 1035. - 313;
  xfilm_top[0][2] = 1000;         yfilm_top[0][2] = 1035. - 263;
  xfilm_top[1][2] = 1043;         yfilm_top[1][2] = 1035. - 300;
  xfilm_top[0][3] = 1000;         yfilm_top[0][3] = 1035. - 788;
  xfilm_top[1][3] = 1055;         yfilm_top[1][3] = 1035. - 825;
  xfilm_top[0][4] = 600;          yfilm_top[0][4] = 1035. - 970;
  xfilm_top[1][4] = 655;          yfilm_top[1][4] = 1035. - 1013;
  xfilm_top[0][5] = 592;          yfilm_top[0][5] = 1035. - 84;
  xfilm_top[1][5] = 641;          yfilm_top[1][5] = 1035. - 130;
  
  xfilm_bot[0][0] = 595;          yfilm_bot[0][0] = 1035. - 267;
  xfilm_bot[1][0] = 645;          yfilm_bot[1][0] = 1035. - 371;
  xfilm_bot[0][1] = 950;          yfilm_bot[0][1] = 1035. - 495;
  xfilm_bot[1][1] = 999;          yfilm_bot[1][1] = 1035. - 595;
  xfilm_bot[0][2] = 599;          yfilm_bot[0][2] = 1035. - 727;
  xfilm_bot[1][2] = 651;          yfilm_bot[1][2] = 1035. - 830;
  xfilm_bot[0][3] = 244;          yfilm_bot[0][3] = 1035. - 497;
  xfilm_bot[1][3] = 297;          yfilm_bot[1][3] = 1035. - 602;
  /* The following coordinates of the optical axes are   */
  /* measured values, not calculated ...                 */
  xaxis[0] = 600.*pixel2mm;             yaxis[0] = (1035. - 261.)*pixel2mm;
  xaxis[1] = 648.*pixel2mm;             yaxis[1] = (1035. - 830.)*pixel2mm;
  //xaxis[0] = 3.92263951757274611; xaxis[1] = 4.23810946673662236;
  //yaxis[0] = 5.60838999898403312; yaxis[1] = 2.06678414682288913;

  /**************************************************************/
  for(i=0; i<=1; i++)
    {
      xax_pxl[i] = xaxis[i]/pixel2mm;
      yax_pxl[i] = yaxis[i]/pixel2mm;
      //xaxis[i] = xaxis[i]*pixel2mm;
      //yaxis[i] = yaxis[i]*pixel2mm;
      
      for(k=0; k<=5; k++)
	{
	  if(k<=3)
	    {xfilm_bot0[i][k] = xfilm_bot[i][k]*pixel2mm;
	    yfilm_bot0[i][k] = yfilm_bot[i][k]*pixel2mm;
	    xfilm_bot[i][k] = xfilm_bot[i][k]*pixel2mm;
	    yfilm_bot[i][k] = yfilm_bot[i][k]*pixel2mm;
	    rfilm_bot[i][k] = sqrt(pow(xfilm_bot[i][k] - xaxis[i], 2)
				   + pow(yfilm_bot[i][k] - yaxis[i], 2));
	    xfilm_bot[i][k] = (xfilm_bot[i][k] - xaxis[i])*
	      pow(rfilm_bot[i][k], exp_lens - 1.) + xaxis[i];
	    yfilm_bot[i][k] = (yfilm_bot[i][k] - yaxis[i])*
	      pow(rfilm_bot[i][k], exp_lens - 1.) + yaxis[i];  }
	  
	  xfilm_top0[i][k] = xfilm_top[i][k]*pixel2mm;
	  yfilm_top0[i][k] = yfilm_top[i][k]*pixel2mm;
	  xfilm_top[i][k] = xfilm_top[i][k]*pixel2mm;
	  yfilm_top[i][k] = yfilm_top[i][k]*pixel2mm;
	  rfilm_top[i][k] = sqrt(pow(xfilm_top[i][k] - xaxis[i], 2)
				 + pow(yfilm_top[i][k] - yaxis[i], 2));
	  xfilm_top[i][k] = (xfilm_top[i][k] - xaxis[i])*
	    pow(rfilm_top[i][k], exp_lens - 1.) + xaxis[i];
	  yfilm_top[i][k] = (yfilm_top[i][k] - yaxis[i])*
	    pow(rfilm_top[i][k], exp_lens - 1.) + yaxis[i];
	}
    }
  /**************************************************************/
  /*                                                            */
  /*      Input distances between fiducials for determining     */
  /*                       magnification                        */
  /*     (measurements made by E.M.Andreev and A.S.Moiseenko    */
  /*                    on January 15, 1999)                    */
  /*                                                            */
  /**************************************************************/
  /*           upper (top)   glass 111 - 6 fiducials            */
  /*                                                            */
  i = 6;
  *(dtop + 2*i + 3) = 207.9;          *(dtop + 2*i + 0) = 419.3;
  *(dtop + 2*i + 1) = 384.6;          *(dtop + 2*i + 5) = 208.0;
  *(dtop + 2*i + 4) = 384.1;
  *(dtop + 3*i + 0) = 384.2;          *(dtop + 3*i + 1) = 454.3;
  *(dtop + 3*i + 5) = 379.4;          *(dtop + 3*i + 4) = 248.8;
  *(dtop + 0*i + 1) = 208.5;          *(dtop + 0*i + 5) = 384.1;
  *(dtop + 0*i + 4) = 207.8;
  *(dtop + 1*i + 5) = 248.0;          *(dtop + 1*i + 4) = 380.0;
  *(dtop + 5*i + 4) = 453.6;
  i = 6;
  *(dtop + 3*i + 2) = 207.9;          *(dtop + 0*i + 2) = 419.3;
  *(dtop + 1*i + 2) = 384.6;          *(dtop + 5*i + 2) = 208.0;
  *(dtop + 4*i + 2) = 384.1;
  *(dtop + 0*i + 3) = 384.2;          *(dtop + 1*i + 3) = 454.3;
  *(dtop + 5*i + 3) = 379.4;          *(dtop + 4*i + 3) = 248.8;
  *(dtop + 1*i + 0) = 208.5;          *(dtop + 5*i + 0) = 384.1;
  *(dtop + 4*i + 0) = 207.8;
  *(dtop + 5*i + 1) = 248.0;          *(dtop + 4*i + 1) = 380.0;
  *(dtop + 4*i + 5) = 453.6;
  /**************************************************************/
  /*            lower (bot) glass 333 - 4 fiducials             */
  /*                                                            */
  i = 4;
  *(dbot + 0*i + 1) = 225.7;          *(dbot + 0*i + 2) = 248.2;
  *(dbot + 0*i + 3) = 226.8;
  *(dbot + 1*i + 2) = 227.3;          *(dbot + 1*i + 3) = 379.4;
  *(dbot + 2*i + 3) = 226.9;
  
  i = 4;
  *(dbot + 1*i + 0) = 225.7;          *(dbot + 2*i + 0) = 248.2;
  *(dbot + 3*i + 0) = 226.8;
  *(dbot + 2*i + 1) = 227.3;          *(dbot + 3*i + 1) = 379.4;
  *(dbot + 3*i + 2) = 226.9;
  return;
}


/*							      */
/*                    S_triangle:                             */
/*                                                            */
/*  In a triangle with known sides a, b, and c and opposite   */
/*                          angles A, B, and C we have:       */
/*       p = 0.5*(a + b + c) - semiperimeter;                 */
/*       r - radius of inserted ("vpisannoi") circle;         */
/*       r = sqrt((p-a)*(p-b)*(p-c)/p);                       */
/*       S - area:   S = r*p,              etc. (p.166)       */
/*							      */

float S_triangle(float a, float b, float c, float m, float n,
					 float h, float *mm)
{
	float r, p, S, m0, n0, h0, cosA, sinA;
	p = 0.5*(a + b + c);
	S = 0.;         r = 0.;
	if(p>a && p>b && p>c)
		{
				r = sqrt((p-a)*(p-b)*(p-c)/p);
				S = r*p;
		}
	cosA = (b*b + c*c - a*a)/(2.*b*c);
	if(fabs(cosA) > 1.) cosA = cosA/fabs(cosA);
	sinA = sqrt(1. - cosA*cosA);
/*                                                             */
/* The following are for a rectanglar triangle only:           */
/*  m0 = a*a/c;     n0 = b*b/c;      h0 = sqrt(m0*n0);         */
/*                                                             */
	h0 = b*sinA;        n0 = b*cosA;        m0 = c - n0;
	*mm = m0;           *(mm+1) = n0;       *(mm+2) = h0;
	//	printf("m = %8.3f mm    n = %8.3f mm     h = %8.3f mm\n",m0,n0,h0);
	return S;
}

/*                                                            */
/*                                                            */
/*   Calculation of determinant of (4x4) matrix ...           */
/*                                                            */
/*                  */
float det4(float *a4)
{
 float det3, det;
 float *ptr4;
 ptr4 = a4;
 det = 0;
 det3 = (*(ptr4+4*1+1)) * (*(ptr4+4*2+2)) * (*(ptr4+4*3+3)) +
	(*(ptr4+4*1+2)) * (*(ptr4+4*2+3)) * (*(ptr4+4*3+1)) +
	(*(ptr4+4*1+3)) * (*(ptr4+4*2+1)) * (*(ptr4+4*3+2)) -
	(*(ptr4+4*1+3)) * (*(ptr4+4*2+2)) * (*(ptr4+4*3+1)) -
	(*(ptr4+4*1+1)) * (*(ptr4+4*2+3)) * (*(ptr4+4*3+2)) -
	(*(ptr4+4*1+2)) * (*(ptr4+4*2+1)) * (*(ptr4+4*3+3));
 det = det + (*(ptr4+4*0+0)) * det3;
 //printf("det3(1) = %10.0f, det4 = %10.0f\n",det3,det);

 det3 = (*(ptr4+4*1+0)) * (*(ptr4+4*2+2)) * (*(ptr4+4*3+3)) +
	(*(ptr4+4*1+2)) * (*(ptr4+4*2+3)) * (*(ptr4+4*3+0)) +
	(*(ptr4+4*1+3)) * (*(ptr4+4*2+0)) * (*(ptr4+4*3+2)) -
	(*(ptr4+4*1+3)) * (*(ptr4+4*2+2)) * (*(ptr4+4*3+0)) -
	(*(ptr4+4*1+0)) * (*(ptr4+4*2+3)) * (*(ptr4+4*3+2)) -
	(*(ptr4+4*1+2)) * (*(ptr4+4*2+0)) * (*(ptr4+4*3+3));
 det = det - (*(ptr4+4*0+1)) * det3;
 //printf("det3(2) = %10.0f, det4 = %10.0f\n",det3,det);

 det3 = (*(ptr4+4*1+0)) * (*(ptr4+4*2+1)) * (*(ptr4+4*3+3)) +
	(*(ptr4+4*1+1)) * (*(ptr4+4*2+3)) * (*(ptr4+4*3+0)) +
	(*(ptr4+4*1+3)) * (*(ptr4+4*2+0)) * (*(ptr4+4*3+1)) -
	(*(ptr4+4*1+3)) * (*(ptr4+4*2+1)) * (*(ptr4+4*3+0)) -
	(*(ptr4+4*1+0)) * (*(ptr4+4*2+3)) * (*(ptr4+4*3+1)) -
	(*(ptr4+4*1+1)) * (*(ptr4+4*2+0)) * (*(ptr4+4*3+3));
 det = det + (*(ptr4+4*0+2)) * det3;
 //printf("det3(3) = %10.0f, det4 = %10.0f\n",det3,det);

 det3 = (*(ptr4+4*1+0)) * (*(ptr4+4*2+1)) * (*(ptr4+4*3+2)) +
	(*(ptr4+4*1+1)) * (*(ptr4+4*2+2)) * (*(ptr4+4*3+0)) +
	(*(ptr4+4*1+2)) * (*(ptr4+4*2+0)) * (*(ptr4+4*3+1)) -
	(*(ptr4+4*1+2)) * (*(ptr4+4*2+1)) * (*(ptr4+4*3+0)) -
	(*(ptr4+4*1+0)) * (*(ptr4+4*2+2)) * (*(ptr4+4*3+1)) -
	(*(ptr4+4*1+1)) * (*(ptr4+4*2+0)) * (*(ptr4+4*3+2));
 det = det - (*(ptr4+4*0+3)) * det3;
 //printf("det3(4) = %10.0f, det4 = %10.0f\n",det3,det);

 return det;
}

/*                  */
/*                                                            */
/*   Calculation of determinant of (5x5) matrix ...           */
/*                                                            */
/*                  */

float det5(float *a5)
{
	float determ4[5], a4[4][4];
	float det;
/*                                                            */
/*  The first elements of the (4x4)- and the (5x5)-matrices   */
/*  are a4[0][0] and a5[0][0], respectively.                  */
/*                                                            */
	int i, j, jj, k;
	float *ptr4, *ptr5;
	ptr5 = a5;     ptr4 = a4[0];
	for(j=0; j<=4; j++)
	 { k=-1;
		 for(jj=0; jj<=4; jj++)
			{if(jj!=j)
				{k++;
				 for(i=0; i<=3; i++) *(ptr4+4*i+k) = *(ptr5+5*(i+1)+jj);
				}
			}
		 *(determ4+j) = det4(ptr4);
		}
	det = 0.;
	for(j=0; j<=4; j++) det = det + pow((-1.), j)*(*(determ4+j));
	return det;
}
void read_data()
{
  double buf;
  int *xcircle, *ycircle;
  int itrack, track, i, k, ii, jj, bright;
  int imin, imax;
  int n;
  int j;
  int ima;
  float fpx[6], fpy[6], lpx[6], lpy[6];
  float tr_len[2][6];
  float wpoint[10000];
  short int circle;
  
  fgets(eventname, 19, for_reconstr);
  fputs(eventname, print_results);
  fscanf(for_reconstr,"%d  %d\n",&event_type, &tracks);
  printf("Event_type = %d, \n ",event_type);
  printf("there are  %d  tracks:\n\n",tracks);
  for(view=0; view<2; view++)
	 {
		//fscanf(for_reconstr,"%d\n",&ima);
		for(itrack=0; itrack<tracks; itrack++)
		  {
			 fscanf(for_reconstr,"%d",&track);
			 //  printf(" view:  %d    image:  %d \n", view, ima);
			 //printf(" track number:   %d \n", track);
			 xcircle = &xtrack1[view][track][0];
			 ycircle = &ytrack1[view][track][0];
			 fscanf(for_reconstr,"%f %f %f %f %f %f %f \n",
					  &rc_int, &xc_int, &yc_int, &stx, &sty, &endx, &endy);
			 fscanf(for_reconstr,"%d\n", &n);
			 printf("Number of points along track is %d\n", n);
			 //
			 //   Because of binning (pairing of 2 pixels into 1 bin) all
			 //   linear quantities measured or calculated in pixels must
			 //   be multiplied by a factor of 2 (...!!!...) for run
			 //
			 //                    November 1999
			 //
			 rc_int = rc_int*2.;
			 xc_int = xc_int*2.;
			 yc_int = yc_int*2.;
			 stx = stx*2.;          endx = endx*2.;
			 sty = sty*2.;          endy = endy*2.;
			 // printf(" rc_int  xc_int  yc_int     stx     sty    endx    endy\n");
			  printf("initial %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f \n\n",
			 	rc_int, xc_int, yc_int, stx, sty, endx, endy);
			 double bufmin = 1000.;
			 double bufmax =-1000.;
			 float WpMin = 0.;
			 float WpMax = 0.;
			 for(i=0; i<n; i++)
				{
				  fscanf(for_reconstr," %d %d %d ", &jj, &ii, &bright);
				  ii *= 2;     ii = ii;    /*  ii is the y-coordinate  */
				  jj *= 2;
				  *(Xcircle+i) = jj;             *(Ycircle+i) = ii;
				  wpoint[i] = (float) bright;
				  buf = atan2((Ycircle[i]-yc_int),(Xcircle[i]-xc_int));
				  if(buf<0.) buf += TMath::Pi()*2.;
				  if(buf < bufmin)
					 { imin = i; bufmin = buf;}
				  if(buf > bufmax)
					 { imax = i; bufmax = buf;}
				  if (wpoint[i] > WpMax) WpMax = wpoint[i];
				  if (wpoint[i] < WpMin) WpMin = wpoint[i];
				}
			 //if(itrack==1) continue;
			 //check for the -Pi - +Pi bounds
			 if((bufmax - bufmin) > TMath::Pi()) 
				{
				  bufmin = 1000.;
				  bufmax =-1000.;
				  for(i=0; i<n; i++)
					 {
						buf = atan2((Ycircle[i]-yc_int),(Xcircle[i]-xc_int));
						if(buf < bufmin)
						  { imin = i; bufmin = buf;}
						if(buf > bufmax)
						  { imax = i; bufmax = buf;}
					 }
				}
			 Radius[view][track] = rc_int;
			 xcenter[view][track] = xc_int;
			 ycenter[view][track] = yc_int;
			 
			 printf(" rc_int  xc_int  yc_int     stx     sty    endx    endy\n");
			 printf("%7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f \n",
					  rc_int, xc_int, yc_int, stx, sty, endx, endy);
			 printf("%7.1f %7.1f %7.1f  \n\n",
					  RC, XC, YC);
			 /*  printf("\n");   */
			 /* end of profile*/
		  }             /* end of cycle over itrack */

      fscanf(for_reconstr," Vertex position: X: %f, Y: %f, Error: %f\n",
	     &Xvertex[view],&Yvertex[view],&Dvertex[view]);
      Xvertex[view] = Xvertex[view]*2.;
      Yvertex[view] = Yvertex[view]*2.;
      Dvertex[view] = Dvertex[view]*2.;
      printf(" Vertex position: X: %f, Y: %f, Error: %f\n",
      	     Xvertex[view], Yvertex[view], Dvertex[view]);
      for (i=0; i< tracks;i++)
	{
	  float sx, sy, ex, ey;
	  sx = fpx[i]; sy = fpy[i]; ex = lpx[i]; ey = lpy[i];
	  float dist_vert0 = sqrt(pow(sy - Yvertex[view], 2)+
				  pow(sx - Xvertex[view], 2));
	  float dist_vert1 = sqrt(pow(ey - Yvertex[view], 2)+
				  pow(ex - Xvertex[view], 2));
	  if(dist_vert1 < dist_vert0) {float buf = sx; sx = ex; ex = buf;}
	  if(dist_vert1 < dist_vert0) {float buf = sy; sy = ey; ey = buf;}
	  angle[view][i] = TMath::ATan2((Double_t) (ey - sy),
					(Double_t) (ex - sx));
	  if(angle[view][i] < 0.) angle[view][i] += 2.*TMath::Pi();
	  printf("Image %d track %d angle %5.3f\n", 
		 view, i, angle[view][i]);
	}
      if(view == 1) 
		  {
			 fscanf(for_reconstr," %f %f \n",
					  &aver_l, &aver_r);
			 float buf1, buf2;
			 for(i=0; i<tracks; i++)
				{
				  fscanf(for_reconstr," %f %f \n",
							&buf1, &buf2);
				  //brrr[0][i] = buf1;
				  //brrr[1][i] = buf2;
				  //brrr[0][i] = IntBr[0][i]/tr_len[0][i];
				  //brrr[1][i] = IntBr[1][i]/tr_len[1][i];
				  //printf("Integral brightness  is %f track length %f   %f %f\n", 
				  //	IntBr[0][i], tr_len[0][i], IntBr[1][i], tr_len[1][i]);
				  printf(" track brightness %f %f\n", 
					  brrr[0][i], brrr[1][i]);
				}

			 printf(" Average bright among images left %f and right %f \n", 
					  aver_l, aver_r);
		  }
      
    }                 /* end of cycle over view   */
}

/*							*/
/*	This routine puts all coordinates in system 	*/
/*  with the y-axis along the base of the video-  */
/*	cameras, then exchanges x- and y-axes so as   */
/*  to work in a right-hand reference system.     */
/*						*/
void TrueCoordinates()
{
  
  int i1min, i1max;
  int i, j, k, nn, view;
  double ksin, kcross;
  double mdumb, ndumb, hdumb, mm[3], *ptr_mm;
  double sinphi[2], cosphi[2];
  double pi;
  double xtrans_track, ytrans_track;
  double sumz;

  double minsig, summag0, summag1, dbuf0, dbuf1, dbuf2, ztop_diff, ztop_diff0,
	 zbot_diff;
  int count, ind;
  int bdelta;
  double delta;
  
  double temp_buff, tg, cst0, znak;
  
  double Xvertax0[2], Yvertax0[2];
  /*                                                  */
  /*  Calculate the fiducial coordinates in chamber   */
  /*  relative to the position of the optical axes    */
  /*  for determining the angle through which all     */
  /*  are to be turned for putting them (on the       */
  /*  CCD image) into the required reference system.  */
  /*                                                  */
  
  pi = M_PI;
  ptr_mm = &mm[0];
  mdumb = hdumb = ndumb = 0.;
  double sumphi0 = 0.;
  double sumphi1 = 0.;

  
      /*                                                  */
      /*  Put xfilm and yfilm in system with xax and yax  */
      /*  as origins...                                   */
      /*                                                  */
  for(j = 0; j < 2; j++)
	 for(i = 0; i <6; i++)
		{ 
		  xfilm_top[j][i] -= xaxis[j];
		  yfilm_top[j][i] -= yaxis[j];
		}
  for(j = 0; j < 2; j++)
	 for(i = 0; i <4; i++)
		{ 
		  xfilm_bot[j][i] -= xaxis[j];
		  yfilm_bot[j][i] -= yaxis[j];
		}
  int iSecondLoop = 1;
  for(j = 0; j < 2; j++)
    for(i = 0; i <6; i++)
      { 
	xchamb_top[j][i] = -xfilm_top[j][i]*mag_top[j];
	ychamb_top[j][i] = -yfilm_top[j][i]*mag_top[j];
	
	xfilm_top0[j][i] = xfilm_top[j][i];
	yfilm_top0[j][i] = yfilm_top[j][i];
      }
  for(j = 0; j < 2; j++)
    for(i = 0; i <4; i++)
      { 
	xchamb_bot[j][i] = -xfilm_bot[j][i]*mag_bot[j];
	ychamb_bot[j][i] = -yfilm_bot[j][i]*mag_bot[j];
	
	xfilm_bot0[j][i] = xfilm_bot[j][i];
	yfilm_bot0[j][i] = yfilm_bot[j][i];
      }
  
  /*                                                       */
  /*  Determine the angle transforming fiducials "chamb"   */
  /*  to fiducials  "true" (in chamber space, like on CCD  */
  /*  images).  Use only "top" fiducials, since they are   */
  /*  all (the first 4) far from the optical axes.         */
  /*                                                       */
  
  double sinph0_1, sinph1_1;
  double cosph0_1, cosph1_1;
  double x_top[2][6],  x_bot[2][4];
  double y_top[2][6],  y_bot[2][4];
  
  double x0, x1, y0, y1;
  sumphi0 = 0.;
  sumphi1 = 0.;
  double xmag0, xmag1, xdiff;
  ksin = 0.; 
  double phi0_min = 10000.;
  double phi1_min = 10000.;
  for(j=0; j<4; j++)
  //for(j=0; j<=5; j++)
    { 
      //now work in right coord. system !!!! X-to right; Y- to the UP!!! 
      //Z- from the bottom to top
      i = 0;
      double a, b, c, cosA, sinA;
      a = xchamb_top[i][j]*xchamb_top[i][j] +
	ychamb_top[i][j]*ychamb_top[i][j];
      a = sqrt(a);
      b = xchamb_top[i+1][j]*xchamb_top[i+1][j] +
	ychamb_top[i+1][j]*ychamb_top[i+1][j];
      b = sqrt(b);
      c = (double) base;
      cosA = (b*b + c*c - a*a)/(2.*b*c);
      sinA = sqrt(1. - cosA*cosA);
      x_top[0][j] = x_top[1][j] = b*sinA;
      y_top[0][j] =-b*cosA + c;
      y_top[1][j] =-b*cosA;
      
      double x_t[2];
      double phi0_i[4], phi1_i[4];
      double phi_min = 10000.;
      for(k=0;k<4;k++)
	{
	  double sign0 = (TMath::Sin(TMath::Pi()/4.+k*TMath::Pi()/2.)/fabs(TMath::Sin(TMath::Pi()/4.+k*TMath::Pi()/2.)));
	  double sign1 = (TMath::Cos(TMath::Pi()/4.+k*TMath::Pi()/2.)/fabs(TMath::Cos(TMath::Pi()/4.+k*TMath::Pi()/2.)));
	  x_t[0] =  sign0*x_top[0][j];
	  x_t[1] =  sign1*x_top[1][j];
	  i = 0;
	  sinph0_1 = (-y_top[i][j]*xchamb_top[i][j] +
		      x_t[i]*ychamb_top[i][j]) /
	    (xchamb_top[i][j]*xchamb_top[i][j] +
	     ychamb_top[i][j]*ychamb_top[i][j]);
	  
	  cosph0_1 = (x_t[i]*xchamb_top[i][j]  +
		      y_top[i][j]*ychamb_top[i][j]) /
	    (xchamb_top[i][j]*xchamb_top[i][j] +
	     ychamb_top[i][j]*ychamb_top[i][j]);
	  x0 =  xfilm_top[i][j]*cosph0_1 + yfilm_top[i][j]*sinph0_1;
	  y0 = -xfilm_top[i][j]*sinph0_1 + yfilm_top[i][j]*cosph0_1;
	  xmag0 = fabs(x_t[i]/x0); 
	  i = 1;
	  sinph1_1 = (-y_top[i][j]*xchamb_top[i][j]  +
		      x_t[i]*ychamb_top[i][j]) /
	    (xchamb_top[i][j]*xchamb_top[i][j] +
	     ychamb_top[i][j]*ychamb_top[i][j]);
	  
	  cosph1_1 = (x_t[i]*xchamb_top[i][j]  +
		      y_top[i][j]*ychamb_top[i][j]) /
	    (xchamb_top[i][j]*xchamb_top[i][j] +
	     ychamb_top[i][j]*ychamb_top[i][j]);
	  x1 =  xfilm_top[i][j]*cosph1_1 + yfilm_top[i][j]*sinph1_1;
	  y1 = -xfilm_top[i][j]*sinph1_1 + yfilm_top[i][j]*cosph1_1;
	  phi0_i[k] = TMath::ATan2(sinph0_1, cosph0_1);
	  phi1_i[k] = TMath::ATan2(sinph1_1, cosph1_1);
	  xmag1 = fabs(x_t[i]/x1);
	  xdiff = fabs(x1*xmag1/xmag0 - x0);
	  //printf("xmag 0 %f and 1 %f, xdiff is %f\n", xmag0, xmag1, xdiff);
	  //printf("j %d k %d xtrue_top %f %f x0 %f x1 %f phi0 %f phi1 %f\n", 
	  //	  j, k, x_t[0], x_t[1], x0, x1, phi0_i[k], phi1_i[k]);
	  if(fabs(phi0_i[k]) + fabs(phi1_i[k]) < phi_min)
	    {
	      phi_min = fabs(phi0_i[k]) + fabs(phi1_i[k]);
	      phi0_min = phi0_i[k];
	      phi1_min = phi1_i[k];
	    }
	}
      printf("best variant is phi0 %f phi1 %f\n", phi0_min, phi1_min);
      sumphi0 += phi0_min;
      sumphi1 += phi1_min;
      ksin += 1.;
    }
  phi[0] = sumphi0/ksin;
  phi[1] = sumphi1/ksin;
  printf(" for TOP crosses: phi[0] = %8.5f phi[1] = %8.5f\n",
	 phi[0], phi[1]);
  

  double ksin_bot = 0.;
  for(j=0; j<=3; j++)
    { if(j == 0 || j == 2) continue;
    //now work in right coord. system !!!! X-to right; Y- to the UP!!! 
    //Z- from the bottom to top
    i = 0;
    double a, b, c, cosA, sinA;
    a = xchamb_bot[i][j]*xchamb_bot[i][j] +
      ychamb_bot[i][j]*ychamb_bot[i][j];
    a = sqrt(a);
    b = xchamb_bot[i+1][j]*xchamb_bot[i+1][j] +
      ychamb_bot[i+1][j]*ychamb_bot[i+1][j];
    b = sqrt(b);
    c = (double) base;
    cosA = (b*b + c*c - a*a)/(2.*b*c);
    sinA = sqrt(1. - cosA*cosA);
    x_bot[0][j] = x_bot[1][j] = b*sinA;
    y_bot[0][j] =-b*cosA + c;
    y_bot[1][j] =-b*cosA;
    
    double x_t[2];
    double phi0_i[4], phi1_i[4];
    double phi_min = 10000.;
    for(k=0;k<4;k++)
      {
	double sign0 = (TMath::Sin(TMath::Pi()/4.+k*TMath::Pi()/2.)/fabs(TMath::Sin(TMath::Pi()/4.+k*TMath::Pi()/2.)));
	double sign1 = (TMath::Cos(TMath::Pi()/4.+k*TMath::Pi()/2.)/fabs(TMath::Cos(TMath::Pi()/4.+k*TMath::Pi()/2.)));
	x_t[0] =  sign0*x_bot[0][j];
	x_t[1] =  sign1*x_bot[1][j];
	i = 0;
	sinph0_1 = (-y_bot[i][j]*xchamb_bot[i][j] +
		    x_t[i]*ychamb_bot[i][j]) /
	  (xchamb_bot[i][j]*xchamb_bot[i][j] +
	   ychamb_bot[i][j]*ychamb_bot[i][j]);
	
	cosph0_1 = (x_t[i]*xchamb_bot[i][j]  +
		    y_bot[i][j]*ychamb_bot[i][j]) /
	  (xchamb_bot[i][j]*xchamb_bot[i][j] +
	   ychamb_bot[i][j]*ychamb_bot[i][j]);
	x0 =  xfilm_bot[i][j]*cosph0_1 + yfilm_bot[i][j]*sinph0_1;
	y0 = -xfilm_bot[i][j]*sinph0_1 + yfilm_bot[i][j]*cosph0_1;
	xmag0 = fabs(x_t[i]/x0); 
	i = 1;
	sinph1_1 = (-y_bot[i][j]*xchamb_bot[i][j]  +
		    x_t[i]*ychamb_bot[i][j]) /
	  (xchamb_bot[i][j]*xchamb_bot[i][j] +
	   ychamb_bot[i][j]*ychamb_bot[i][j]);
	
	cosph1_1 = (x_t[i]*xchamb_bot[i][j]  +
		    y_bot[i][j]*ychamb_bot[i][j]) /
	  (xchamb_bot[i][j]*xchamb_bot[i][j] +
	   ychamb_bot[i][j]*ychamb_bot[i][j]);
	x1 =  xfilm_bot[i][j]*cosph1_1 + yfilm_bot[i][j]*sinph1_1;
	y1 = -xfilm_bot[i][j]*sinph1_1 + yfilm_bot[i][j]*cosph1_1;
	phi0_i[k] = TMath::ATan2(sinph0_1, cosph0_1);
	phi1_i[k] = TMath::ATan2(sinph1_1, cosph1_1);
	xmag1 = fabs(x_t[i]/x1);
	xdiff = fabs(x1*xmag1/xmag0 - x0);
	//printf("xmag 0 %f and 1 %f, xdiff is %f\n", xmag0, xmag1, xdiff);
	//printf("j %d k %d xtrue_bot %f %f x0 %f x1 %f phi0 %f phi1 %f\n", 
	//	  j, k, x_t[0], x_t[1], x0, x1, phi0_i[k], phi1_i[k]);
	if(fabs(phi0_i[k]) + fabs(phi1_i[k]) < phi_min)
	  {
	    phi_min = fabs(phi0_i[k]) + fabs(phi1_i[k]);
	    phi0_min = phi0_i[k];
	    phi1_min = phi1_i[k];
	  }
      }
    printf("best variant is phi0 %f phi1 %f\n", phi0_min, phi1_min);
    sumphi0 += phi0_min;
    sumphi1 += phi1_min;
    ksin += 1.;
    ksin_bot +=1.;
    }
  phi[0] = sumphi0/ksin;
  phi[1] = sumphi1/ksin;
  printf(" for all crosses: phi[0] = %8.5f phi[1] = %8.5f\n",
	 phi[0], phi[1]);
  
  //int buf;
  //cin >> buf;
  
  cosphi[0] = cos(phi[0]); sinphi[0] = sin(phi[0]);
  cosphi[1] = cos(phi[1]); sinphi[1] = sin(phi[1]);
  /*                                                     */
  /*  Now put all coordinates of points on track images  */
  /*  into "true" reference systems (as yet related to   */
  /*  the actual "view").  Do the same with the fiducial */
  /*  images.                                            */
  /*                                                     */
  for(i=0; i<=1; i++)      /*   i  is the    view  */
    {
      Xvertax0[i] = Xvertex[i];
      Yvertax0[i] = Yvertex[i];
      for(j=0; j<tracks; j++)
	{
	  Radius0[i][j] = Radius[i][j];
	  xcenter0[i][j] = xcenter[i][j];
	  ycenter0[i][j] = ycenter[i][j];
	  for(k=0; k<np[i][j]; k++)
	    {
	      xtrack[i][j][k] =
		(xprong[i][j][k]-xax_pxl[i])*cosphi[i] +
		(yprong[i][j][k]-yax_pxl[i])*sinphi[i];
	      
	      ytrack[i][j][k] =
		-(xprong[i][j][k]-xax_pxl[i])*sinphi[i] +
		(yprong[i][j][k]-yax_pxl[i])*cosphi[i];
	    }     /* end of cycle over k */
	  
	  xcenter[i][j] = (xcenter0[i][j]-xax_pxl[i])*cosphi[i]
	    + (ycenter0[i][j]-yax_pxl[i])*sinphi[i];
	  ycenter[i][j] =-(xcenter0[i][j]-xax_pxl[i])*sinphi[i]
	    + (ycenter0[i][j]-yax_pxl[i])*cosphi[i];
	  xcenter01[i][j] = xcenter[i][j];
	  ycenter01[i][j] = ycenter[i][j];
	}        /* end of cycle over j (tracks) */
      Xvertex[i] = (Xvertax0[i]-xax_pxl[i])*cosphi[i]
	+ (Yvertax0[i]-yax_pxl[i])*sinphi[i];
      Yvertex[i] =-(Xvertax0[i]-xax_pxl[i])*sinphi[i]
	+ (Yvertax0[i]-yax_pxl[i])*cosphi[i];
    }		
  for(i=0; i<=1; i++)
    {
      for(j=0; j<=5; j++)
	{
	  xfilm_top[i][j] = xfilm_top0[i][j]*cosphi[i]
	    + yfilm_top0[i][j]*sinphi[i];
	  yfilm_top[i][j] = xfilm_top0[i][j]*sinphi[i]
	    + yfilm_top0[i][j]*cosphi[i];
	}
      for(j=0; j<=3; j++)
	{
	  xfilm_bot[i][j] = xfilm_bot0[i][j]*cosphi[i]
	    + yfilm_bot0[i][j]*sinphi[i];
	  yfilm_bot[i][j] = xfilm_bot0[i][j]*sinphi[i]
	    + yfilm_bot0[i][j]*cosphi[i];
	}
    }
  
  /*                                                    */
  /*  We must now put all the coordinates into one and  */
  /*  the same reference frame, for instance, related   */
  /*  to the upper crosses of view numb. 0. In doing    */
  /*  such a transfer make use of the upper crosses,    */
  /*  like before.                                      */
  /*                                                    */
  
  kcross = 0.;
  xtransfer = 0.;       ytransfer = 0.;
  /*                                                    */
  /*   xtransfer and ytransfer are in mm., like xfilm   */
  /*   and yfilm for fiducials, so  they have  to  be   */
  /*   divided by      pixel2mm,   when we deal with    */
  /*   we deal with track coordinates. A factor "4"     */
  /*   will be added to correspond to run measured last */
  /*   autumn:                                          */
  /*                      May1999                       */
  /*                                                    */
  
  //for Gil's picture of crosses!(When he superposes top crosses and 
  // looks to the bottom ones  
  for(i=0; i<=3; i++)
    {
      kcross = kcross + 1.;
      xtransfer = xtransfer + (xfilm_top[0][i]*mag_top[0]/mag_top[1] - xfilm_top[1][i]);
      ytransfer = ytransfer + (yfilm_top[0][i]*mag_top[0]/mag_top[1] - yfilm_top[1][i]);
    }
  ytransfer = ytransfer/kcross;
  xtransfer = xtransfer/kcross;
  xtrans_track = xtransfer/pixel2mm;
  ytrans_track = ytransfer/pixel2mm;
  //printf("xtransfer = %8.4f     ytransfer = %8.4f in pixels\n",
  //xtrans_track, ytrans_track);
  for(i=0; i<6; i++)
    {
      xfilmd_top[0][i] = xfilm_top[0][i];
      yfilmd_top[0][i] = yfilm_top[0][i];
    }
  for(i=0; i<4; i++)
    {
      xfilmd_bot[0][i] = xfilm_bot[0][i];
      yfilmd_bot[0][i] = yfilm_bot[0][i];
    }
  for(i=0; i<6; i++)
    {
      xfilmd_top[1][i] = xfilm_top[1][i] + xtransfer;
      yfilmd_top[1][i] = yfilm_top[1][i] + ytransfer;
    }
  for(i=0; i<4; i++)
    {
      xfilmd_bot[1][i] = xfilm_bot[1][i] + xtransfer;
      yfilmd_bot[1][i] = yfilm_bot[1][i] + ytransfer;
    }
  
  xaxis_draw[0] = xaxis[0];
  yaxis_draw[0] = yaxis[0];
  xaxis_draw[1] = xaxis[1] + xtransfer;
  yaxis_draw[1] = yaxis[1] + ytransfer;
  // to balance crosses and points on 2 images (knowing that their X must be  equals!
  kcross = 0.;
  xtransfer = ytransfer = 0.;
  for(i=0; i<6; i++)
    {
      kcross = kcross + 1.;
      xtransfer = xtransfer + (xfilm_top[0][i]*mag_top[0]/mag_top[1] - xfilm_top[1][i]);
      ytransfer = ytransfer + (yfilm_top[0][i]*mag_top[0]/mag_top[1] - yfilm_top[1][i]);
    }
  ytransfer = ytransfer/kcross;
  sumz = 0.;
  for (i = 0; i <6; i++) {
    sumz +=  
      pow(((yfilm_top[0][i]*mag_top[0]/mag_top[1] - 
	    yfilm_top[1][i]) -  ytransfer),2);
  }
  ytransfer_err = TMath::Sqrt(sumz/kcross);
  printf("Ytransfer_top %f +/- error %f \n",
	 ytransfer, ytransfer_err);

  ytransfer2 = 0.;
  for(i=0; i<4; i++)
    {
      //if(i == 0 || i == 2) continue;      
      kcross = kcross + 1.;
      xtransfer = xtransfer + (xfilm_bot[0][i]*mag_bot[0]/mag_bot[1] - xfilm_bot[1][i]);
      ytransfer2 = ytransfer2 + (yfilm_bot[0][i]*mag_bot[0]/mag_bot[1] - yfilm_bot[1][i]);
    }
  xtransfer = xtransfer/kcross;
  ytransfer2 = ytransfer2/4.;
  sumz = 0.;
  for (i = 0; i <4; i++) sumz +=  
			   pow((yfilm_bot[0][i]*mag_bot[0]/mag_bot[1] - 
				yfilm_bot[1][i]) -  ytransfer2,2);
  ytransfer2_err = sqrt(sumz/4.);
  printf("Ytransfer_bot %f +/- error %f \n",
	 ytransfer2, ytransfer2_err);
  
  sumz = 0.;
  for (i = 0; i <6; i++) sumz +=  
			   pow((xfilm_top[0][i]*mag_top[0]/mag_top[1] - 
				xfilm_top[1][i]) -  xtransfer,2);
  for (i = 0; i <4; i++) sumz +=  
			   pow((xfilm_bot[0][i]*mag_bot[0]/mag_bot[1] - 
				xfilm_bot[1][i]) -  xtransfer,2);
  xtransfer_err = sqrt(sumz/10.);
  printf("Xtransfer %f +/- error %f \n",
	 xtransfer, xtransfer_err);
  
  for(i=0; i<6; i++)
    {
      xfilm_top[1][i] = xfilm_top[0][i]*mag_top[0]/mag_top[1] - xtransfer;
      yfilm_top[1][i] = yfilm_top[0][i]*mag_top[0]/mag_top[1] - ytransfer;
      printf("1: xfilm_top %f yfilm_top %f\n\n", xfilm_top[1][i], yfilm_top[1][i]);
    }
  for(i=0; i<4; i++)
    {
      //printf("0: xfilm_bot %f yfilm_bot %f\n", xfilm_bot[0][i], yfilm_bot[0][i]);
      //printf("1: xfilm_bot %f yfilm_bot %f\n", xfilm_bot[1][i], yfilm_bot[1][i]);
      xfilm_bot[1][i] = xfilm_bot[0][i]*mag_bot[0]/mag_bot[1] - xtransfer;
      yfilm_bot[1][i] = yfilm_bot[0][i]*mag_bot[0]/mag_bot[1] - ytransfer2;
      printf("1: xfilm_bot %f yfilm_bot %f\n\n", xfilm_bot[1][i], yfilm_bot[1][i]);
    }
  
  //if(iSecondLoop == 1)  return;
  if(iSecondLoop == 1)  goto NEXT;
  iSecondLoop = 1;
  double d_top_leftflm, d_top_rghtflm;
  double d_bot_leftflm, d_bot_rghtflm;
  double magt0[15], magb0[15], magt1[15], magb1[15];
  int ip, jp;
  minsig = 1000000.; 
  double dleft_max, drght_max;
  double dd;
  count = 0;     
  int ipmax0, ipmax1, jpmax0, jpmax1;
  int ifl0[6][6], ifl1[6][6];
  bdelta = 1000; //random numbers from 0 to 1000000
  delta = (double) 1./(double) bdelta;
  double dleft[15], drght[15], dcam[15];
  ind = 0;
  for(i = 0; i <6; i++)
    for(j = 0; j <6; j++)
      {
	ifl0[i][j] = 0; ifl1[i][j] = 0;
      }
  do
    {
      count = 0;
      dleft_max = drght_max = 0.;
      for(ip = 0; ip <5; ip++)
	{
	  for(jp = ip+1; jp <=5; jp++)
	    {
	      if(ifl0[ip][jp]==1) continue;
	      d_top_leftflm =
		pow((xfilm_top[0][ip] - xfilm_top[0][jp]), 2) +
		pow((yfilm_top[0][ip] - yfilm_top[0][jp]), 2);
	      d_top_leftflm = sqrt(d_top_leftflm)/pixel2mm;
	      if(d_top_leftflm > dleft_max) {dleft_max = d_top_leftflm; ipmax0 = ip; jpmax0=jp; dd = d_top_leftflm;}
	      count++;
	    }
	}
      printf(" %8.3f %8.3f Maximal d for %d and %d crosses on left image\n", 
	     dleft_max, dd, ipmax0, jpmax0);
      ifl0[ipmax0][jpmax0] = 1;
      if(count>0) 
	{
	  dleft[ind] = dd;
	  ind++;
	}
    }
  while(count != 0);
  count = 0;
  ind = 0;
  do
    {
      count = 0;
      dleft_max = drght_max = 0.;
      for(ip = 0; ip <5; ip++)
	{
	  for(jp = ip+1; jp <=5; jp++)
	    {
	      if(ifl1[ip][jp]==1) continue;
			/*
	      double dcur = 0.;
	      srand(time(NULL));
	      for(i=0; i<bdelta; i++)
		{
		  double xcur1, xcur2, ycur1, ycur2;
		  xcur1 = xfilm_top[1][ip]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  ycur1 = yfilm_top[1][ip]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  xcur2 = xfilm_top[1][jp]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  ycur2 = yfilm_top[1][jp]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  dcur  += sqrt(pow(xcur2-xcur1,2)+pow(ycur2-ycur1,2));
		}
	      dcur = dcur*delta;
			*/
	      d_top_rghtflm =
		pow((xfilm_top[1][ip] - xfilm_top[1][jp]), 2) +
		pow((yfilm_top[1][ip] - yfilm_top[1][jp]), 2);
	      d_top_rghtflm = sqrt(d_top_rghtflm)/pixel2mm;
	      if(d_top_rghtflm > drght_max) {drght_max = d_top_rghtflm; ipmax1 = ip; jpmax1=jp;dd = d_top_rghtflm;}
	      //if(dcur > drght_max) {drght_max = dcur; ipmax1 = ip; jpmax1=jp; dd = d_top_rghtflm;}
	      count++;
	    }
	}
      printf(" %8.3f %8.3f maximal %d and %d on the right one\n", 
	     drght_max, dd, ipmax1, jpmax1);
      ifl1[ipmax1][jpmax1] = 1;
      if(count>0) 
	{
	  drght[ind] = dd;
	  ind++;
	}
    }
  while(count != 0);
  count = 0;
  ind = 0;
  double dmax;
  for(ip = 0; ip <5; ip++)
    for(jp = ip+1; jp <=5; jp++)
      ifl1[ip][jp]=0;
  do
    {
      count = 0;
      dmax = 0.;
      for(ip = 0; ip <5; ip++)
	{
	  for(jp = ip+1; jp <=5; jp++)
	    {
	      if(ifl1[ip][jp]==1) continue;
	      if(d_fiduc_top[ip][jp] > dmax) {dmax = d_fiduc_top[ip][jp]; ipmax1 = ip; jpmax1=jp;}
	      count++;
	    }
	}
      printf(" %8.3f maximal %d and %d on the top of camera\n", 
	     dmax, ipmax1, jpmax1);
      ifl1[ipmax1][jpmax1] = 1;
      if(count>0) 
	{
	  dcam[ind] = dmax;
	  ind++;
	}
    }
  while(count != 0);
  summag0 = 0.;
  summag1 = 0.;
  for(i = 0; i< 3; i++)
    {
      magt0[i] = dcam[i]/dleft[i]/pixel2mm;
      magt1[i] = dcam[i]/drght[i]/pixel2mm;
      summag0 += magt0[i];
      summag1 += magt1[i];
    }
  dbuf0 = 0.;
  dbuf1 = 0.;
  dbuf2 = 0.;
  for(i = 3; i< 9; i++)
    {
      dbuf0 += dleft[i];
      dbuf1 += drght[i];
      dbuf2 += dcam[i];
    }
  for(i = 3; i< 9; i++)
    {
      magt0[i] = dbuf2/dbuf0/pixel2mm;
      magt1[i] = dbuf2/dbuf1/pixel2mm;
      if(i==3) 
	{
	  summag0 += magt0[i];
	  summag1 += magt1[i];
	}
    }
  dbuf0 = 0.;
  dbuf1 = 0.;
  dbuf2 = 0.;
  for( i = 9; i< 11; i++)
    {
      dbuf0 += dleft[i];
      dbuf1 += drght[i];
      dbuf2 += dcam[i];
    }
  for( i = 9; i< 11; i++)
    {
      magt0[i] = dbuf2/dbuf0/pixel2mm;
      magt1[i] = dbuf2/dbuf1/pixel2mm;
      if(i==9) 
	{
	  summag0 += magt0[i];
	  summag1 += magt1[i];
	}
    }
  dbuf0 = 0.;
  dbuf1 = 0.;
  dbuf2 = 0.;
  for( i = 11; i< 15; i++)
    {
      dbuf0 += dleft[i];
      dbuf1 += drght[i];
      dbuf2 += dcam[i];
    }
  for( i = 11; i< 15; i++)
    {
      magt0[i] = dbuf2/dbuf0/pixel2mm;
      magt1[i] = dbuf2/dbuf1/pixel2mm;
      if(i==11) 
	{
	  summag0 += magt0[i];
	  summag1 += magt1[i];
	}
    }
  for( i = 0; i< 15; i++)
    {
      printf("%d mag_top [0] %f [1] %f\n", i, magt0[i], magt1[i]);
    }
  mag_top[0] = magt0[top];//summag0/6.;
  mag_top[1] = magt1[top];//summag1/6.;
  printf("mean mag_top [0] %13.8f   [1] %13.8f\n", mag_top[0], mag_top[1]);

  ind = 0;
  for(ip = 0; ip <3; ip++)
    for(jp = ip+1; jp <=3; jp++)
      ifl0[ip][jp]=0;
  do
    {
      count = 0;
      dleft_max = drght_max = 0.;
      for(ip = 0; ip <3; ip++)
		  {
			 for(jp = ip+1; jp <=3; jp++)
				{
				  if(ifl0[ip][jp]==1) continue;
				  /*
				  double dcur = 0.;
			     srand(time(NULL));
			     for(i=0; i<bdelta; i++)
			       {
						double xcur1, xcur2, ycur1, ycur2;
						xcur1 = xfilm_bot[0][ip]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
						ycur1 = yfilm_bot[0][ip]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
						xcur2 = xfilm_bot[0][jp]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
						ycur2 = yfilm_bot[0][jp]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
						dcur  += sqrt(pow(xcur2-xcur1,2)+pow(ycur2-ycur1,2));
			       }
			     dcur = dcur*delta;
				  */
			     d_bot_leftflm =
			       pow((xfilm_bot[0][ip] - xfilm_bot[0][jp]), 2) +
			       pow((yfilm_bot[0][ip] - yfilm_bot[0][jp]), 2);
			     d_bot_leftflm = sqrt(d_bot_leftflm)/pixel2mm;
			     if(d_bot_leftflm > dleft_max) {dleft_max = d_bot_leftflm; ipmax0 = ip; jpmax0=jp;dd = d_bot_leftflm;}
			     //if(dcur > dleft_max) {dleft_max = dcur; ipmax0 = ip; jpmax0=jp; dd = d_bot_leftflm;}
			     count++;
				}
		  }
      printf(" BOT: %8.3f %8.3f Maximal d for %d and %d crosses on left image\n", 
				 dleft_max, dd, ipmax0, jpmax0);
      ifl0[ipmax0][jpmax0] = 1;
      if(count>0) 
		  {
			 dleft[ind] = dd;
			 ind++;
		  }
    }
  while(count != 0);
  count = 0;
  ind = 0;
  for(ip = 0; ip <3; ip++)
    for(jp = ip+1; jp <=3; jp++)
      ifl1[ip][jp]=0;
  do
    {
      count = 0;
      dleft_max = drght_max = 0.;
      for(ip = 0; ip <3; ip++)
	{
	  for(jp = ip+1; jp <=3; jp++)
	    {
	      if(ifl1[ip][jp]==1) continue;
			/*
	      double dcur = 0.;
	      srand(time(NULL));
	      for(i=0; i<bdelta; i++)
		{
		  double xcur1, xcur2, ycur1, ycur2;
		  xcur1 = xfilm_bot[1][ip]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  ycur1 = yfilm_bot[1][ip]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  xcur2 = xfilm_bot[1][jp]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  ycur2 = yfilm_bot[1][jp]/pixel2mm-0.5+delta*((double) (rand()%bdelta));
		  dcur  += sqrt(pow(xcur2-xcur1,2)+pow(ycur2-ycur1,2));
		}
	      dcur = dcur*delta;
			*/
	      d_bot_rghtflm =
		pow((xfilm_bot[1][ip] - xfilm_bot[1][jp]), 2) +
		pow((yfilm_bot[1][ip] - yfilm_bot[1][jp]), 2);
	      d_bot_rghtflm = sqrt(d_bot_rghtflm)/pixel2mm;
	      if(d_bot_rghtflm > drght_max) {drght_max = d_bot_rghtflm; ipmax1 = ip; jpmax1=jp; dd = d_bot_rghtflm;}
	      //if(dcur > drght_max) {drght_max = dcur; ipmax1 = ip; jpmax1=jp; dd = d_bot_rghtflm;}
	      count++;
	    }
	}
      printf("BOT: %8.3f %8.3f maximal %d and %d on the right one\n", 
	     drght_max, dd, ipmax1, jpmax1);
      ifl1[ipmax1][jpmax1] = 1;
      if(count>0) 
	{
	  drght[ind] = dd;
	  ind++;
	}
    }
  while(count != 0);
  count = 0;
  ind = 0;
  for(ip = 0; ip <3; ip++)
    for(jp = ip+1; jp <=3; jp++)
      ifl1[ip][jp]=0;
  do
    {
      count = 0;
      dmax = 0.;
      for(ip = 0; ip <3; ip++)
	{
	  for(jp = ip+1; jp <=3; jp++)
	    {
	      if(ifl1[ip][jp]==1) continue;
	      if(d_fiduc_bot[ip][jp] > dmax) {dmax = d_fiduc_bot[ip][jp]; ipmax1 = ip; jpmax1=jp;}
	      count++;
	    }
	}
      printf(" %8.3f maximal %d and %d on the bottom of the camera\n", 
	     dmax, ipmax1, jpmax1);
      ifl1[ipmax1][jpmax1] = 1;
      if(count>0) 
	{
	  dcam[ind] = dmax;
	  ind++;
	}
    }
  while(count != 0);

  summag0 = 0.;
  summag1 = 0.;
  for(i = 0; i< 3; i++)
    {
      magb0[i] = dcam[i]/dleft[i]/pixel2mm;
      magb1[i] = dcam[i]/drght[i]/pixel2mm;
      summag0 += magb0[i];
      summag1 += magb1[i];
    }
  dbuf0 = 0.;
  dbuf1 = 0.;
  dbuf2 = 0.;
  for(i = 3; i< 5; i++)
    {
      dbuf0 += dleft[i];
      dbuf1 += drght[i];
      dbuf2 += dcam[i];
    }
  for(i = 3; i< 5; i++)
    {
      magb0[i] = dbuf2/dbuf0/pixel2mm;
      magb1[i] = dbuf2/dbuf1/pixel2mm;
      if(i==3) 
	{
	  summag0 += magb0[i];
	  summag1 += magb1[i];
	}
    }
  dbuf0 = 0.;
  dbuf1 = 0.;
  dbuf2 = 0.;
  for( i = 5; i< 6; i++)
    {
      dbuf0 += dleft[i];
      dbuf1 += drght[i];
      dbuf2 += dcam[i];
    }
  for( i = 5; i< 6; i++)
    {
      magb0[i] = dbuf2/dbuf0/pixel2mm;
      magb1[i] = dbuf2/dbuf1/pixel2mm;
      if(i==5) 
	{
	  summag0 += magb0[i];
	  summag1 += magb1[i];
	}
    }
  for( i = 0; i< 6; i++)
    {
      printf("%d mag_bot [0] %f [1] %f\n", i, magb0[i], magb1[i]);
    }
  mag_bot[0] = magb0[bot];//summag0/5.;
  mag_bot[1] = magb1[bot];//summag1/5.;
  printf("mean mag_bot [0] %13.8f   [1] %13.8f\n", mag_bot[0], mag_bot[1]);

  depth = 160. + 8./3.;
  zfilm[0] = depth/(mag_bot[0] - mag_top[0]);
  zfilm[1] = depth/(mag_bot[1] - mag_top[1]);
  double z_top0, z_top1;
  z_top0 = mag_top[0]*zfilm[0];
  z_top1 = mag_top[1]*zfilm[1];
  double z_bot0, z_bot1;
  z_bot0 = mag_bot[0]*zfilm[0];
  z_bot1 = mag_bot[1]*zfilm[1];
  ztop = (z_top0 +z_top1)/2.; zbot = (z_bot0 + z_bot1)/2.;
  printf("ztop0 = %8.3f ztop1 = %8.3f zbot0 %f zbot1 %f zfilm 0 = %f 1 = %f \n",
	 z_top0, z_top1, z_bot0, z_bot1, zfilm[0], zfilm[1]);
  printf("ztop = %8.3f zbot %f zmed %f depth %f \n",
	 ztop, zbot, (zbot - ztop)/2., zbot - ztop);

  //double ztop_diff, zbot_diff, ztop_diff0;
  ztop_diff = 0.;
  do
	 {
		ztop_diff0 = ztop_diff;
		mag_top[0] = (mag_top[0] + mag_top[1]*mag_bot[0]/mag_bot[1])/2.;
		mag_top[1] = (mag_top[1] + mag_top[0]*mag_bot[1]/mag_bot[0])/2.;
		mag_bot[0] = (mag_bot[0] + mag_bot[1]*mag_top[0]/mag_top[1])/2.;
		mag_bot[1] = (mag_bot[1] + mag_bot[0]*mag_top[1]/mag_top[0])/2.;
		zfilm[0] = depth/(mag_bot[0] - mag_top[0]);
		zfilm[1] = depth/(mag_bot[1] - mag_top[1]);
		ztop_diff = fabs(mag_top[0]*zfilm[0] - mag_top[1]*zfilm[1]);
		zbot_diff = fabs(mag_bot[0]*zfilm[0] - mag_bot[1]*zfilm[1]);
		printf("%18.10f mag_top[0] %18.10f mag_top[1]\n", mag_top[0], mag_top[1]);
		printf("%18.10f mag_bot[0] %18.10f mag_bot[1]\n", mag_bot[0], mag_bot[1]);
		printf("ztop_diff %f zbot_diff %f \n\n", ztop_diff, zbot_diff);
	 }
  while(ztop_diff > 0.00001 && ztop_diff > ztop_diff0);
  zfilm[0] = depth/(mag_bot[0] - mag_top[0]);
  zfilm[1] = depth/(mag_bot[1] - mag_top[1]);
  ztop = (mag_top[0]*zfilm[0] + mag_top[1]*zfilm[1])/2.;
  zbot = (mag_bot[0]*zfilm[0] + mag_bot[1]*zfilm[1])/2.;
  zmed = 0.5*(ztop + zbot);
  printf("zfilm 0 %7.5f 1 %7.5f mm   ztop = %7.2f mm   zbot = %7.2f mm\n",
	 zfilm[0], zfilm[1], ztop, zbot);
  /*                                                         */
  /*   Now calculate the parameters of the circular tracks,  */
  /*   including their direction cosines on the plane.       */
  /*                                                         */
  
   NEXT:
  for(view=0; view<2; view++)
    {
      for(track=0; track<tracks; track++)
	{
	  /*                                                           */
	  /* Assign Vertex coordinates to first (...[view][track][0])  */
	  /* coordinates on track images.                              */
	  /*                                                           */
	  xtrack[view][track][0] = Xvertex[view];
	  ytrack[view][track][0] = Yvertex[view];      /* from Vera */
	  CCDBuffer[0][track][0] = 0.;
	  CCDBuffer[1][track][0] = 0.;
	}
    }
  //
  // Now start calculating the true characteristics of the vertex in
  // chamber space.
  //
  
  yvertex[0] = Yvertex[0];
  yvertex[1] = Yvertex[1];
  hvertex = (yvertex[1] - yvertex[0])/
    (yfilm_bot[1][2] - yfilm_bot[0][2]);
  hvertex = depth*hvertex*pixel2mm;
  zvertex = ztop + hvertex;
  printf("hvertex = %7.2f      zvertex = %7.2f\n",
   hvertex, zvertex);
  
  //
  // From mag_bot = (ztop+depth)/zfilm, mag_top = ztop/zfilm,
  //      mag_vert = (ztop+hvertex)/zfilm
  // calculate mag_vert[view] by excluding zfilm, ztop and zbot...
  //
  
  for(view=0; view<2; view++)
    {  mag_vert[view] = mag_top[view]
	 + hvertex*(mag_bot[view] - mag_top[view])/depth;
    printf("mag_vert[%d] = %7.2f\n", view,mag_vert[view]);
    //
    //  Calculate coordinates of vertex in chamber space. Do not
    //  forget to take into account pixel2mm, since till now all
    //  track characteristics were in pixels on the "film" (CCD
    //  matrix). Note also that the coordinates are already in
    //  reference systems with optical axes at the origins. In
    //  passing to chamber space exchange X-axis and the Y-axis,
    //  so as to work normally in a right-handed reference system.
    //  For the time being do not transform initial left-handed
    //  xvertex and yvertex: use other variables xvertex_r and
    //  yvertex_r.
    //
    
    xvertex_r[view] = Yvertex[view]*pixel2mm;
    yvertex_r[view] = Xvertex[view]*pixel2mm;
    if(view==1)  { xvertex_r[view] = xvertex_r[view] - ytransfer;
    yvertex_r[view] = yvertex_r[view] - xtransfer;
    }
    xvertex_r[view] = -xvertex_r[view]*mag_vert[view];
    yvertex_r[view] = -yvertex_r[view]*mag_vert[view];
    //printf("In space:  xvertex_r[%d] = %8.2f mm,    "
    //   "yvertex_r[%d] = %8.2f mm\n",
    //   view, xvertex_r[view], view, yvertex_r[view]);
    //printf("           hvertex = %7.2f        zvertex = %7.2f\n",
  //hvertex, zvertex);
    }
  for(track=0; track<tracks; track++)
    {
      for(view=0; view<2; view++)
	{
	  //
	  //  Determine direction cosines for track images on the CCD
	  //  plane and calculate the corresponding "ends" of direction
	  //  vectors in chamber space so as to find the dip angles and
	  //  direction cosines in chamber space. (We are still in
	  //  left-handed reference system.)
	  //
	  //  Start with calculating the tangent of track image direction
	  //  (Bronstein & Semendyaev, pp. 204, 182).
	  //
	  tg = (Yvertex[view] - ycenter[view][track])
	    /(Xvertex[view] - xcenter[view][track]);  /* in pixels */
	  tg = -1./tg;
	  tg = fabs(tg);             /* we determine the signs later */
	  cosx[view][track] = 1./sqrt(1. + tg*tg);
	  cosy[view][track] = tg*cosx[view][track];
	  delta = xtrack[view][track][1] - xtrack[view][track][0];
	  znak = delta/fabs(delta);
	  cosx[view][track] = cosx[view][track]*znak;
	  delta = ytrack[view][track][1] - ytrack[view][track][0];
	  znak = delta/fabs(delta);
	  cosy[view][track] = cosy[view][track]*znak;
	  //
	  //  Now, the equation of a straight line passing through point
	  //  (xvertex, yvertex) with inclination tg = cosy/cosx is
	  //  y = tg*x + yvertex - tg*xvertex. So we can, for instance,
	  //  set xdir[0] = xdir[1] = xtrack[0][track][10] and calculate
	  //  ydir[0] and ydir[1] to use the difference ydir[1]-ydir[0]
	  //  in finding hdir = zdir - ztop and the dip angle and direction
	  //  cosines in space. We are "sempre" in left-hand reference frame
	  //  where the image pairs of upper fiducials coincide.
	  //
	  xdir[view][track] = xtrack[0][track][10];
	  tg = cosy[view][track]/cosx[view][track];
	  ydir[view][track] = tg*xdir[view][track] + Yvertex[view]
	    -tg*Xvertex[view];
	  //printf("xdir[%d][%d] = %7.2f,  ydir[%d][%d] = %7.2f\n",
	  // view, track, xdir[view][track],
	  // view, track, ydir[view][track]);
	}
      //printf("cosx[0][%d] = %6.3f       cosy[0][%d] = %6.3f\n",
      //	     track, cosx[0][track],
      //     track, cosy[0][track]);
      //printf("cosx[1][%d] = %6.3f       cosy[1][%d] = %6.3f\n\n",
	//     track, cosx[1][track],
	  //   track, cosy[1][track]);
    }
  for(track=0; track<tracks; track++)
    {
      hdir[track] = (ydir[1][track] - ydir[0][track])
	/(yfilm_bot[1][2] - yfilm_bot[0][2]);
      hdir[track] = depth*hdir[track]*pixel2mm;
      zdir[track] = ztop + hdir[track];
      //printf("On CCD (in pixels): ydir[0][%d] = %7.2f,   "
      //     "ydir[1][%d] = %7.2f\n",
      //     track, ydir[0][track], track, ydir[1][track]);
      //printf("hvertex = %7.2f      zvertex = %7.2f\n",
  // hvertex, zvertex);
      //printf("hdir[%d] = %7.2f      zdir[%d] = %7.2f\n",
	  //   track, hdir[track], track, zdir[track]);
    }
  
  for(view=0; view<2; view++)
    {
      //
      //    Now start fully working in right-hand reference
      //    system in chamber space, meaning that the X- and
      //    Y-axes must be exchanged with each other. Start
      //    with xvertex and yvertex.
      //
      xvertex[view] = xvertex_r[view];
      yvertex[view] = yvertex_r[view];
    }
  for(track=0; track<tracks; track++)
    {
      //
      //  Calculate magnification corresponding to z-coordinate
      //  zdir[track]:  mag-dir[view][track].
      //
      for(view=0; view<2; view++)
	{  mag_dir[view][track] = mag_top[view]
	     + hdir[track]*(mag_bot[view] - mag_top[view])/depth;
	//printf("mag_dir[%d][%d] = %7.2f\n",
	//     view, track, mag_dir[view][track]);
	//
	//  Calculate coordinates of dir_end in chamber space. Do not
	//  forget to take into account pixel2mm, since till now all
	//  track characteristics were in pixels on the "film" (CCD
	//  matrix).
	//
	
	xdir[view][track] = xdir[view][track]*pixel2mm;
	ydir[view][track] = ydir[view][track]*pixel2mm;
	if(view==1)  { xdir[view][track] = xdir[view][track] - xtransfer;
	ydir[view][track] = ydir[view][track] - ytransfer;
	}
	//
	//      Exchange the X- and Y- axes with each other
	//      (for xdir and ydir), when calculating x- and
	//      y-coordinates in chamber space.
	//
	temp_buff = xdir[view][track];
	xdir[view][track] = -ydir[view][track]*mag_dir[view][track];
	ydir[view][track] = -temp_buff*mag_dir[view][track];
	//printf("In space:  xdir[%d][%d] = %8.2f mm,    "
	//     "xvertex[%d] = %8.2f mm\n",
	//     view, track, xdir[view][track], view, xvertex[view]);
	//printf("           ydir[%d][%d] = %8.2f mm,    "
	//     "yvertex[%d] = %8.2f mm\n",
	//     view, track, ydir[view][track], view, yvertex[view]);
	//printf("           hvertex = %7.2f        zvertex = %7.2f\n",
	//     hvertex, zvertex);
	}
    }
  
  for(track=0; track<tracks; track++)
    {
      tg = 0.;
      range[2][track] = 0.;
      for(view=0; view<2; view++)
	{
	  //
	  //     We have everything to calculate the dip
	  //     angles from both lenses as required  in
	  //     article  [PTE, numb.4 (1964) 66]   (for
	  //     both optical axes).
	  //
	  temp_buff = (xdir[view][track] - xvertex[view])
	    *(xdir[view][track] - xvertex[view])
	    +(ydir[view][track] - yvertex[view])
	    *(ydir[view][track] - yvertex[view]);
	  //printf("temp_buff = %8.1f, xdir[%d][%d] - xv[%d] = %8.4f\n",
	  //	 sqrt(temp_buff), view, track, view,
	  // xdir[view][track] - xvertex[view]);
	  //printf("temp_buff = %8.1f, ydir[%d][%d] - yv[%d] = %8.4f\n",
      // sqrt(temp_buff), view, track, view,
      //	 ydir[view][track] - yvertex[view]);
	  //
	  //   Calculate the true direction cosines in horizontal plane.
	  //
	  cosx[view][track] = (xdir[view][track] - xvertex[view])
	    /sqrt(temp_buff);
	  cosy[view][track] = (ydir[view][track] - yvertex[view])
	    /sqrt(temp_buff);
	  range[2][track] = range[2][track] + temp_buff
	    +(hdir[track]-hvertex)*(hdir[track]-hvertex);
	  //printf("range[2][%d]**2 = %10.1f mm**2 (in space)\n",
	  //	 track, range[2][track]);
	  temp_buff = (float) sqrt(temp_buff);
	  tg = (hdir[track]-hvertex)/temp_buff + tg;
	  dip[track] = (float) atan((hdir[track]-hvertex)/temp_buff)*(180./pi);
	  //
	  //    As required in ref. [PTE, numb.4 (1964) 66], positive
	  //    dip angles are directed upward, toward the lens.
	  //
	  dip[track] = -dip[track];
	  if(view==0)
	    {
	      //printf("track %d  view %d:   dip = %6.2f deg,  "
//	     "sqrt(dx**2+dy**2) = %6.2f mm\n",
//	     track, view, dip[track], temp_buff);
	    }
	  if(view==1)
	    {
	      //printf("         view %d:   dip = %6.2f deg,  "
//	     "sqrt(dx**2+dy**2) = %6.2f mm\n",
//	     view, dip[track], temp_buff);
	    }
	}
      range[2][track] = (float) sqrt(range[2][track]/2.);
      //printf("range[2][%d] = %9.1f mm (direct. vector in space)\n",
      //     track, range[2][track]);
      tgtheta[track] = tg/2.;
      tgtheta[track] = -tgtheta[track];
      dip[track] = (float) atan(tgtheta[track])*(180./pi);
      //printf("tgtheta[%d] = %7.2f,   dip[%d] = %6.2f deg\n",
      //     track, tgtheta[track], track, dip[track]);
      //fprintf(print_results,"tgtheta[%d] = %7.2f,  dip[%d] = %6.2f deg\n",
      //	      track, tgtheta[track], track, dip[track]);
      //	      if(track==0 || track==4) tgtheta[track] = -tgtheta[track];
    }
  
  for(track=0; track<tracks; track++)
    {
      cosx[2][track] = 0.;
      cosy[2][track] = 0.;
      cosz[track] = 0.;
      for(view=0; view<2; view++)
	{
	  //
	  //  Calculate the direction cosines of tracks in space.
	  //
	  temp_buff = cosx[2][track];
	  cosx[2][track] = (xdir[view][track] - xvertex[view])
	    / range[2][track];
	  if(view==1)
	    {
	      //printf("track %d view 0 cosx[2] = %8.4f,  "
//	     "view 1 cosx[2] = %8.4f\n",
//	     track, temp_buff, cosx[2][track]);
	      cosx[2][track] = (temp_buff + cosx[2][track])/2.;
	    }
	  temp_buff = cosy[2][track];
	  cosy[2][track] = (ydir[view][track] - yvertex[view])
	    / range[2][track];
	  if(view==1)
	    {
	      //printf("track %d view 0 cosy[2] = %8.4f,  "
//	     "view 1 cosy[2] = %8.4f\n",
//	     track, temp_buff, cosy[2][track]);
	      cosy[2][track] = (temp_buff + cosy[2][track])/2.;
	    }
	  temp_buff = cosz[track];
	  cosz[track] = (hdir[track]-hvertex)/range[2][track];
	  if(view==1)
	    {
	      //printf("track %d view 0 cosz[2] = %8.4f,  "
//     "view 1 cosz[2] = %8.4f\n",
//	     track, temp_buff, cosz[track]);
	      cosz[track] = (temp_buff + cosz[track])/2.;
	    }
	}
      //printf("cosx[0][%d] = %8.4f,  cosy[0][%d] = %8.4f\n"
	//     "cosx[1][%d] = %8.4f,  cosy[1][%d] = %8.4f\n",
	  //   track, cosx[0][track], track, cosy[0][track],
	    // track, cosx[1][track], track, cosy[1][track]);
      //fprintf(print_results,"cosx[0][%d] = %8.4f,  "
      //      "cosy[0][%d] = %8.4f\n"
      //      "cosx[1][%d] = %8.4f,  cosy[1][%d] = %8.4f\n",
      //      track, cosx[0][track], track, cosy[0][track],
      //      track, cosx[1][track], track, cosy[1][track]);
      cosx[0][track] = (cosx[0][track] + cosx[1][track])/2.;
      cosx[1][track] = cosx[0][track];
      cosy[0][track] = (cosy[0][track] + cosy[1][track])/2.;
      cosy[1][track] = cosy[0][track];
      //printf("<cosx[1][%d]> = %8.4f,  <cosy[1][%d]> = %8.4f\n"
      // "cosx**2 + cosy**2 (in horiz. plane) = %8.4f\n",
      //     track, cosx[1][track], track, cosy[1][track],
      //     cosx[1][track]*cosx[1][track]
      //     + cosy[1][track]*cosy[1][track]);
      //printf("cosx[2][%d] = %8.4f\n"
      //     "cosy[2][%d] = %8.4f\n"
      //     "cosz[%d]    = %8.4f\n"
      //     "cosx**2 + cosy**2 + cosz**2 = %8.4f\n",
      //     track, cosx[2][track], track, cosy[2][track],
      //     track, cosz[track],
      //     cosx[2][track]*cosx[2][track]
      //     +cosy[2][track]*cosy[2][track]
      //     +cosz[track]*cosz[track]);
      //fprintf(print_results,"<cosx[1][%d]> = %8.4f,  "
      //      "<cosy[1][%d]> = %8.4f\n"
      //      "cosx**2 + cosy**2 (in horiz. plane) = %8.4f\n",
      //      track, cosx[1][track], track, cosy[1][track],
      //      cosx[1][track]*cosx[1][track]
      //      + cosy[1][track]*cosy[1][track]);
      //fprintf(print_results,"cosx[2][%d] = %8.4f\n"
      //      "cosy[2][%d] = %8.4f\n"
      //      "cosz[%d]    = %8.4f\n"
      //      "cosx**2 + cosy**2 + cosz**2 = %8.4f\n",
      //      track, cosx[2][track], track, cosy[2][track],
      //      track, cosz[track],
      //      cosx[2][track]*cosx[2][track]
      //      +cosy[2][track]*cosy[2][track]
      //      +cosz[track]*cosz[track]);
      //printf("--------------\n");
    }
  
  for(track=0; track<tracks; track++)
    {
      for(view=0; view<2; view++)
	{
	  temp_buff = xcenter[view][track];
	  xcenter_r[view][track] = ycenter[view][track]*pixel2mm;
	  ycenter_r[view][track] = temp_buff*pixel2mm;
	  if(view==1)  {xcenter_r[view][track] = xcenter_r[view][track]
			  - ytransfer;
	  ycenter_r[view][track] = ycenter_r[view][track]
	    - xtransfer;
	  }
	  xcenter_r[view][track] = -xcenter_r[view][track]
	    *mag_vert[view];
	  ycenter_r[view][track] = -ycenter_r[view][track]
	    *mag_vert[view];
	  
	  xcenter[view][track] =  xcenter_r[view][track];
	  ycenter[view][track] =  ycenter_r[view][track];
	  Radius[view][track]  =  Radius[view][track]*pixel2mm
	    *mag_vert[view];
	  
	  //printf("In chamber space (relative to each axis):\n");
	  //printf("xcenter_r[%d][%d] = %8.2f mm,	  "
	  // "ycenter_r[%d][%d] = %8.2f mm\n",
	  // view,track,xcenter_r[view][track],
	  // view,track,ycenter_r[view][track]);
	  //printf("Radius[%d][%d] = %8.2f mm\n",
	  // view, track, Radius[view][track]);
	}
    }
  
  for(view=0; view<2; view++)
    {
      for(track=0; track<tracks; track++)
	{
	  //	 xcircle = &ytrack[view][track][0];
	  //	 ycircle = &xtrack[view][track][0];
	  //	nn = 3;    /* in case there are only 3 points along the track */
	  nn = np[view][track];
	  //
	  //	CIRCLE(xcircle, ycircle, &xcentr, &ycentr, &rcircle, &mcir,
	  //					&ncir, &qcir, &nn);
	  //
	  rc_int = Radius[view][track];     xc_int = xcenter[view][track];
	  yc_int = ycenter[view][track];
	  stx = xvertex[view];              sty = yvertex[view];
	  double bufmin = 1000.;
	  double bufmax =-1000.;
	  double phivert = atan2((yvertex[view]-ycenter[view][track]),
				 (xvertex[view]-xcenter[view][track]));
	  if(phivert<0.) phivert += TMath::Pi()*2.;

	  for(i=0; i<nn; i++)
	    {
	      double buf = atan2((ytrack[view][track][i]-ycenter[view][track]),
				 (xtrack[view][track][i]-xcenter[view][track]));
	      if(buf<0.) buf += TMath::Pi()*2.;
	      if(buf < bufmin)
		{ i1min = i; bufmin = buf;}
	      if(buf > bufmax)
		{ i1max = i; bufmax = buf;}
	    }
	  //check for the -Pi - +Pi bounds
	  if((bufmax - bufmin) > TMath::Pi()) 
	    {
	      int Buf = i1min;
	      i1min = i1max;
	      i1max = Buf;
	      double Buf2 = bufmin;
	      bufmin = bufmax;
	      bufmax = Buf2;
	    }
	  if(fabs(bufmax - phivert) > TMath::Pi())
	    {
	      if((bufmax - phivert) > 0.) 
		phivert += TMath::Pi()*2.;
	      if((bufmax - phivert) < 0.) 
		phivert -= TMath::Pi()*2.;
	    }
	  if(fabs(bufmin - phivert) > TMath::Pi())
	    {
	      if((bufmin - phivert) > 0.) 
		phivert += TMath::Pi()*2.;
	      if((bufmin - phivert) < 0.) 
		phivert -= TMath::Pi()*2.;
	    }
	  int iend = i1min;
	  if(fabs(bufmax - phivert) > fabs(bufmin - phivert)) iend = i1max;
	  endx = ytrack[view][track][iend]*pixel2mm;
	  endy = xtrack[view][track][iend]*pixel2mm;
	  //printf("%d, %d:  endx = %8.2f mm,    endy = %8.2f mm\n",
	  // view, track, endx, endy);
	  endx = -endx*mag_vert[view];
	  endy = -endy*mag_vert[view];
	  
	  mtrack[view][track] = -xcenter[view][track];
	  ntrack[view][track] = -ycenter[view][track];
	  qtrack[view][track] = mtrack[view][track]*mtrack[view][track]
	    + ntrack[view][track]*ntrack[view][track]
	    - Radius[view][track]*Radius[view][track];
	  
	  //printf("%7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f \n",
	  // rc_int, xc_int, yc_int, stx, sty, endx, endy);
	  
	  xstart[view][track] = stx;           ystart[view][track] = sty;
	  xend[view][track] = endx;            yend[view][track] = endy;
	  
	  range[view][track] =
	    sqrt((endx-stx)*(endx-stx) + (endy-sty)*(endy-sty));
	  range[view][track] = 2.*asin(range[view][track]/(2.*rc_int));
	  //
	  // Here the range is angular and expressed in radians.
	  //
	  stangle[view][track] = acos( (stx - xc_int)/rc_int);
	  if (sty - yc_int < 0 && stangle[view][track] <= pi/2.)
	    stangle[view][track] = - stangle[view][track];
	  else if (sty - yc_int < 0 && stangle[view][track] > pi/2.)
	    stangle[view][track] = 2.*pi - stangle[view][track];
	  endangle[view][track] = acos( (endx - xc_int)/rc_int);
	  if (endy - yc_int < 0 && endangle[view][track] <= pi/2.)
	    endangle[view][track] = - endangle[view][track];
	  else if (endy - yc_int < 0 && endangle[view][track] > pi/2.)
	    endangle[view][track] = 2.*pi - endangle[view][track];
	  
	  if(endangle[view][track] < stangle[view][track])
	    range[view][track] = -range[view][track];
	  dangle[view][track] = range[view][track]/100.;
	  //printf("stangle[%d][%d] = %8.2f       "
	  // "endangle[%d][%d] = %8.2f\n",
	  // view,track,stangle[view][track],
	  // view,track,endangle[view][track]);
	  //printf("range[%d][%d]   = %8.2f       "
	  // "dangle[%d][%d]   = %8.2f\n",
	  // view,track,range[view][track],
	  // view,track,dangle[view][track]);
	  
	}      /* end of loop over 'track' */
    }          /* end of loop over 'view'  */
  
  for(view=0; view<2; view++)
    {
      for(track=0; track<tracks; track++)
	{
	  //
	  //  We can now calculate tgt0[view][track], since we know
	  //  tgalpha = 1./(cosy/cosx), and alpha = t0 + mu.
	  //
	  tgalpha[view][track] = cosx[view][track]/cosy[view][track];
	  tgt0[view][track] = (tgalpha[view][track] - tgmu[view][track])
	    /(1. + tgalpha[view][track]*tgmu[view][track]);
	  tg = tgt0[view][track];
	  cst0 = 1./(1. + tg*tg);       /* check with previously */
	  /* calculated value      */
	  //printf("view %d track %d: cst0 = %8.5f, cost0 = %8.5f\n",
	  // view, track, cst0, cost0[view][track]);
	  
	}     /* end of cycle over track */
    }       /* end of cycle over view  */
  
}
/*                                                     */
/*  Angles: calculates the various angles such as      */
/*  direction from vertex to center, sint0, cost0,     */
/*  etc.                                               */
/*                                                     */

void Angles()
  
{
  for(track=0; track<tracks; track++)
    {
      for(view=0; view<2; view++)
	{
	  //
	  //   For determining the direction cosines from center
	  //   toward vertex make use of table in Bronstein, p.182.
	  //
	  if(iangles==0)
	    {
	      cosy_rho[view][track] =  cosx[view][track];
	      cosx_rho[view][track] = -cosy[view][track];
	    }
	  if(iangles==1)
	    {
	      cosy_rho[view][track] = -cosx[view][track];
	      cosx_rho[view][track] =  cosy[view][track];
	    }
	  if(iangles==2)
	    {
	      cosy_rho[view][track] =  cosx[view][track];
	      cosx_rho[view][track] =  cosy[view][track];
	    }
	  if(iangles==3)
	    {
	      cosy_rho[view][track] = -cosx[view][track];
	      cosx_rho[view][track] = -cosy[view][track];
	    }
	  //
	  //  Calculate the "real" xcenter0 and ycenter0 in chamber space
	  //  using the measured rho_meas[view][track] = Radius[view][track]
	  //  values as first approximations.
	  //
	  xcenter0[view][track] =
	    xvertex[view] - 0.5*(Radius[0][track]+Radius[1][track])
	    *cosx_rho[view][track];
	  ycenter0[view][track] =
	    yvertex[view] - 0.5*(Radius[0][track]+Radius[1][track])
	    *cosy_rho[view][track];
	  //
	  //  Here, the minus sign corresponds to the direction from the
	  //  vertex to the center...
	  //
	  //printf("cosx_rho[%d][%d] = %8.5f     cosy (track) = %8.5f,\n"
	  //	 "cosy_rho[%d][%d] = %8.5f     cosx (track) = %8.5f\n",
	  // view, track, cosx_rho[view][track],
	  // cosy[view][track],
	  // view, track, cosy_rho[view][track],
	  // cosx[view][track]);
	  //printf("---------------\n");
	}
      //printf("xcenter0[0][%d] = %8.2f mm,   "
  // "ycenter0[0][%d] = %8.2f mm\n",
  //     track, xcenter0[0][track],
  //     track, ycenter0[0][track]);
      //printf("xcenter0[1][%d] = %8.2f mm,   "
  //     "ycenter0[1][%d] = %8.2f mm\n",
  //     track, xcenter0[1][track],
  //     track, ycenter0[1][track]);
      //printf("xcenter[0][%d]  = %8.2f mm,    "
  //     "ycenter[0][%d]  = %8.2f mm\n",
  //     track, xcenter[0][track],
  //     track, ycenter[0][track]);
      //printf("xcenter[1][%d]  = %8.2f mm,    "
  //     "ycenter[1][%d]  = %8.2f mm\n",
  //     track, xcenter[1][track],
  //     track, ycenter[1][track]);
    }
  for(track=0; track<tracks; track++)
    {
      for(view=0; view<2; view++)
	{
	  //
	  //  Calculate sinto and costo in vertex plane (the same as in
	  //  any other plane in chamber space) of each track at the
	  //  vertex point, as required in article [PTE, numb.4 (1964) 66].
	  //
	  Rvert[view][track] = sqrt(
				    xvertex[view]*xvertex[view]
				    +yvertex[view]*yvertex[view]);
	  rho[view][track] =  sqrt(
				   (xvertex[view]-xcenter0[view][track])
				   * (xvertex[view]-xcenter0[view][track])
				   + (yvertex[view]-ycenter0[view][track])
				   * (yvertex[view]-ycenter0[view][track]));
	  //printf("Rvert[%d] = %8.1f: dist(vert-origin)\n"
	  // "rho[%d][%d] = %8.5f: dist(vert-centr),  "
	  // "Radius[%d][%d] = %8.2f \n",
	  // view, Rvert[view][track],
	  // view, track, rho[view][track],
	  // view, track, Radius[view][track]);
	  
	  R0[view][track] = sqrt(
				 xcenter0[view][track]*xcenter0[view][track]
				 + ycenter0[view][track]*ycenter0[view][track]);
	  //printf("R0[%d][%d] = %8.1f: dist(centr-origin)\n",
	  //		 view, track, R0[view][track]);
	  
	  tgmu[view][track] = ycenter0[view][track]/xcenter0[view][track];
	  
	  //printf("tgmu[%d][%d] = %8.4f\n", view, track, tgmu[view][track]);
	  //printf("  \n");
	}
    }
  
  //
  //     Find sint0[view][track] and cost0[view][track]
  //     by first turning the reference system about an
  //     angle mu, i.e. the angle of the direction toward
  //     the center of the circle, then calculating sint0
  //     and cost0 directly.
  //
  for(track=0; track<tracks; track++)
    {
      for(view=0; view<2; view++)
	{
	  cosmu[view][track] = xcenter0[view][track]/
	    R0[view][track];
	  sinmu[view][track] = ycenter0[view][track]/
	    R0[view][track];
	  
	  xcprime[view][track] =
	    xcenter0[view][track]*cosmu[view][track]
	    + ycenter0[view][track]*sinmu[view][track];
	  //printf("xcenter_prime = %8.2f\n", xcprime[view][track]);
	  ycprime[view][track] =
	    - xcenter0[view][track]*sinmu[view][track]
	    + ycenter0[view][track]*cosmu[view][track];
	  //printf("ycenter_prime = %8.2f\n", ycprime[view][track]);
	  
	  xvprime[view] =
	    xvertex[view]*cosmu[view][track]
	    + yvertex[view]*sinmu[view][track];
	  //printf("xvertex_prime = %8.2f\n", xvprime[view]);
	  yvprime[view] =
	    - xvertex[view]*sinmu[view][track]
	    + yvertex[view]*cosmu[view][track];
	  //printf("yvertex_prime = %8.2f\n", yvprime[view]);
	  sint0[view][track] = (yvprime[view] - ycprime[view][track])
	    / rho[view][track];
	  cost0[view][track] = (xvprime[view] - xcprime[view][track])
	    / rho[view][track];
	  //printf("sint0[%d][%d] = %8.4f,  cost0[%d][%d] = %8.4f  "
	  //	 "sin**2 + cos**2 = %8.2f\n",
	  // view, track, sint0[view][track],
	  // view, track, cost0[view][track],
	  // sint0[view][track]*sint0[view][track]
	  // +cost0[view][track]*cost0[view][track]);
	  //printf("   \n");
	}
    }
  //          End of routine     Angles(void)
}

/*                                                      */
/*  RECONSTRUCT PROPER:  Actually reconstructs tracks   */
/*  and calculates the vertex point, the direction      */
/*  cosines for each track at the vertex point, and,    */
/*  then, the helix of each track (thus, its dip angle  */
/*  and radius).                                        */
/*                                                      */

void ReconstructProper()
{
  int view, track, i;
  float sum_cos2;
  float pi;
  float rcheck0, rcheck1;
  float rmeas, brack1, brack2;
  float c, C;
  
  pi = M_PI;
  //
  //		 zvertex = ztop + hvertex
  //
  //printf("hvertex = %7.2f      zvertex = %7.2f\n",
// hvertex, zvertex);
  
  //
  // Now apply formula (7) of ref. [PTE, numb.4 (1964) 66] for calculation
  // of true radius of track helix. As first approximation consider
  //
  //        Radius[2][track] = 0.5*(Radius[0][track]+Radius[1][track]).
  //
  
  for(track=0; track<tracks; track++)
    {
      rho[0][track] = Radius[0][track];
      rho[1][track] = Radius[1][track];
      //printf("Initial:  rho[0][%d] = %8.2f,     rho[1][%d] = %8.2f\n",
//	     track, rho[0][track],  track, rho[1][track]);
      for(view=0; view<2; view++)
	{
	  c = R0[view][track];          /*   distance:  origin - center  */
	  C = rho_R0[view][track];      /*   angle opposite to   c       */
	  Radius[2][track] = 0.5*(Radius[0][track]+Radius[1][track]);
	  rcheck0 = rho[view][track];
	  rmeas = Radius[view][track];
	  rho_aux[view][iangles][track] = 0.;
	  
	  for(i=0; i<np[view][track]; i++)
	    {
	      brack1 = (1. - (tgtheta[track]/zvertex)
			*R0[view][track]*sint0[view][track] );
	      brack2 = brack1;
	      brack1 = brack1*brack1+(tgtheta[track]*tgtheta[track]
				      /(zvertex*zvertex))
		*(rcheck0+R0[view][track]*cost0[view][track])
		*(rcheck0+R0[view][track]*cost0[view][track]);
	      brack1 = pow(sqrt(brack1), 3);
	      brack2 = brack2 + 2.*rcheck0*(tgtheta[track]*tgtheta[track]
					    /(zvertex*zvertex))*(rcheck0+R0[view][track]
								 *cost0[view][track]);
	      rcheck1 = rmeas*brack2/brack1;
	      //
	      //   Calculate the new xcenter, ycenter and consequently the
	      //   new R0 and new sinto, cost0, knowing that the direction
	      //   toward the center (cosx_rho, cosy_rho) remains the same.
	      //
	      //		 printf("vw %d tr %d: xcenter0 = %8.1f  ycenter0 = %8.1f\n",
	      //						 view, track, xcenter0[view][track],
	      //													ycenter0[view][track]);
	      xcenter0[view][track] = xvertex[view]
		-rcheck1*cosx_rho[view][track];
	      ycenter0[view][track] = yvertex[view]
		-rcheck1*cosy_rho[view][track];
	      //		 printf("vw %d tr %d: xcenter0 = %8.1f  ycenter0 = %8.1f\n",
	      //						 view, track, xcenter0[view][track],
	      //													ycenter0[view][track]);
	      R0[view][track] = sqrt(
				     xcenter0[view][track]*xcenter0[view][track]
				     + ycenter0[view][track]*ycenter0[view][track]);
	      //		 printf("R0[%d][%d] = %8.1f: dist(centr-origin)\n",
	      //						 view, track, R0[view][track]);
	      
	      tgmu[view][track] = ycenter0[view][track]/xcenter0[view][track];
	      
	      //		 printf("tgmu[%d][%d] = %8.4f\n", view, track, tgmu[view][track]);
	      //		 printf("  \n");
	      //
	      //     Find sint0[view][track] and cost0[view][track]
	      //     by first turning the reference system about an
	      //     angle mu, i.e. the angle of the direction toward
	      //     the center of the circle, then calculating sint0
	      //     and cost0 directly.
	      //
	      cosmu[view][track] = xcenter0[view][track]/
		R0[view][track];
	      sinmu[view][track] = ycenter0[view][track]/
		R0[view][track];
	      
	      xcprime[view][track] =
		xcenter0[view][track]*cosmu[view][track]
		+ ycenter0[view][track]*sinmu[view][track];
	      //						printf("xcenter_prime = %8.2f\n", xcprime[view][track]);
	      ycprime[view][track] =
		- xcenter0[view][track]*sinmu[view][track]
		+ ycenter0[view][track]*cosmu[view][track];
	      //						printf("ycenter_prime = %8.2f\n", ycprime[view][track]);
	      
	      xvprime[view] =
		xvertex[view]*cosmu[view][track]
		+ yvertex[view]*sinmu[view][track];
	      //						printf("xvertex_prime = %8.2f\n", xvprime[view]);
	      yvprime[view] =
		- xvertex[view]*sinmu[view][track]
		+ yvertex[view]*cosmu[view][track];
	      //						printf("yvertex_prime = %8.2f\n", yvprime[view]);
	      sint0[view][track] = (yvprime[view] - ycprime[view][track])
		/ rcheck1;
	      cost0[view][track] = (xvprime[view] - xcprime[view][track])
		/ rcheck1;
//						printf("sint0[%d][%d] = %8.4f,  cost0[%d][%d] = %8.4f  "
//									 "sin**2 + cos**2 = %8.2f\n",
//										view, track, sint0[view][track],
//										view, track, cost0[view][track],
//										sint0[view][track]*sint0[view][track]
//									 +cost0[view][track]*cost0[view][track]);
//						printf("   \n");
	      
	      if(fabs(rcheck1 - rcheck0) < 0.01)
		{
		  //printf("rcheck1 - rcheck0 = %8.2f\n",rcheck1-rcheck0);
		  rho_aux[view][iangles][track] = rcheck1;
		  //							 rho[view][track] = rcheck1;
		  break;
		}
	      rcheck0 = rcheck1;
	    }
	  //printf("rho_aux[view %d][iangles %d][track %d] = %8.2f\n",
//	 view, iangles, track, rho_aux[view][iangles][track]);
	  //printf("  \n");
	}
    }
  for(track=0; track<tracks; track++)
    {
      //printf("dip angle[%d] = %5.1f\n", track,dip[track]);
      //		 printf("range[2][%d]  = %5.1f deg.\n",
      //					track, range[2][track]*(180./pi));
      view = 2;
      //printf("for track in space:\n");
      //printf("cosx[%d][%d] =  %6.3f,      cosy[%d][%d] =  %6.3f\n",
//     view, track, cosx[view][track],
//     view, track, cosy[view][track]);
      //printf("cosz[%d] =  %6.3f\n", track, cosz[track]);
      sum_cos2 =  cosx[view][track]*cosx[view][track]
	+ cosy[view][track]*cosy[view][track]
	+ cosz[track]*cosz[track];
      //printf("cosx**2 + cosy**2 + cosz**2 = %6.3f\n", sum_cos2);
      //printf("  \n");
    }
}


/*                                         */
/*  A routine for determining the actual   */
/*  radius in space and the deviations.    */
/*                                         */

void Rho()
{
  int i, iaux, j, jaux;
  float  deltarmin, pi;
  
  pi = M_PI;
  
  for(track=0; track<tracks; track++)
    {
      for(view=0; view<2; view++)
	{
	  for(i=0; i<4; i++)
	    {
	      //printf("rho_aux[%d][%d][%d] = %8.2f\n",
//	     view, i, track, rho_aux[view][i][track]);
	    }
	  //printf(" \n");
	}
    }
  
  for(track=0; track<tracks; track++)
    {
      for(i=0; i<4; i++)
	{
	  for(j=0; j<4; j++)
	    {
	      deltarho[i][j][track] = fabs(rho_aux[0][i][track]
					   -rho_aux[1][j][track]);
	      //printf("deltarho[%d][%d][%d] = %8.2f\n",
//	     i, j, track, deltarho[i][j][track]);
	    }
	  //printf(" \n");
	}
    }
  
  for(track=0; track<tracks; track++)
    {
      deltarmin = deltarho[0][0][track];
      rho[0][track] = rho_aux[0][0][track];
      rho[1][track] = rho_aux[1][0][track];
      iaux = 0;
      jaux = 0;
      for(i=0; i<4; i++)
	{
	  for(j=0; j<4; j++)
	    {
	      if(deltarho[i][j][track] <= deltarmin)
		{ deltarmin = deltarho[i][j][track];
		rho[0][track] = rho_aux[0][i][track];
		rho[1][track] = rho_aux[1][j][track];
		iaux = i;
		jaux = j;
		}
	      //printf("track %d:  iaux %d, jaux %d, deltarmin = %8.2f\n",
//	     track, iaux, jaux, deltarmin);
	    }
	  //printf(" ");
	}
      Radius[2][track] = (rho[0][track] + rho[1][track])/2.;
      //printf("Radius[0][%d] = %8.2f,     "
//     "Radius[1][%d] = %8.2f\n",
//     track, Radius[0][track], track, Radius[1][track]);
      rho[2][track] = (rho[0][track] + rho[1][track])/2.;
      deltar[track] = fabs(rho[0][track] - rho[1][track]);
      //printf("iaux %d, jaux %d:    deltar[%d] = %8.2f\n",
//     iaux, jaux,  track, deltar[track]);
      //printf("rho[0][%d] = %8.2f,    rho[1][%d] = %8.2f,    "
//    "rho[2][%d] = %8.2f mm\n",
//     track, rho[0][track], track, rho[1][track],
//     track, rho[2][track]);
    }
  //printf("  \n");
}

/*                                                     */
/* This routine calculates the kinematical quantities  */
/* for the reconstructed event.                        */
/*                                                     */

void	Kinematics()
{
  int track;
  float massH3, massHe4, masspi, massproton, massneutron;
  float deltapx, deltapy, deltapz, deltaT, deltaE, pi;
  float percent;
  
  pi = M_PI;
  masspi = 139.56;
  massproton = 938.26;
  massneutron = 939.55;
  massHe4 = 3727.23;
  massH3 = 2808.53;
  deltapx = 0.;    deltapy = 0.;     deltapz = 0.;
  deltaE = 0.;
  deltaT = 0.;
  for(track=0; track<tracks; track++)
    {
      percent = 100.*deltar[track]/rho[2][track];
//printf(" rho[0][%d] = %7.1f,  rho[1][%d] = %7.1f,\n"
//						 "deltarho[%d] = %6.1f (%4.1f% ),\n"
//  " rho[2][%d] = %7.1f +-%7.1f (%4.1f%)\n\n",
//track, rho[0][track], track, rho[1][track],
//						 track, deltarho[track], percent,
//     track, rho[2][track], deltar[track]/2.,percent/2.);
      
      //fprintf(print_results," rho[0][%d] =%7.1f,  rho[1][%d] =%7.1f,\n"
	      //						 "deltarho[%d] = %6.1f (%4.1f% ),\n"
      //      " rho[2][%d] = %7.1f +-%7.1f (%4.1f%)\n\n",
      //      track, rho[0][track], track, rho[1][track],
	      //						 track, deltarho[track], percent,
      //      track, rho[2][track], deltar[track]/2.,percent/2.);
    }
  for(track=0; track<tracks; track++)
    {
      mass[track] = massproton;
      if(track==0 || track==1) mass[track] = masspi;
      //
      //   Set the magnetic field to 7000 Gauss.
      //6500Gauss
      //p[track] = 0.21*Radius[2][track]
      p[track] = 0.195*Radius[2][track]
	/ cos(dip[track]*pi/180.);  /* in MeV/c */
      E[track] = sqrt(p[track]*p[track] + mass[track]*mass[track]);
      T[track] = E[track] - mass[track];
      deltapx = deltapx + p[track]*cosx[2][track];
      deltapy = deltapy + p[track]*cosy[2][track];
      deltapz = deltapz + p[track]*cosz[track];
      deltaT = deltaT + T[track];
      if(track == 0) deltaT = deltaT - 2.*T[track];
      //printf("p[%d] = %8.2f   T[%d] = %8.2f\n",
//     track, p[track], track, T[track]);
      
      //fprintf(print_results,"p[%d] = %8.2f   T[%d] = %8.2f\n",
      //      track, p[track], track, T[track]);
    }
  //printf("deltapx = %8.2f,  deltapy = %8.2f,  "
// "deltapz = %8.2f (in MeV/c)\n",
// deltapx, deltapy, deltapz);
  //printf("deltaT = %8.2f MeV\n", deltaT);
  //fprintf(print_results,"deltapx = %8.2f,  deltapy = %8.2f,  "
  //  "deltapz = %8.2f (in MeV/c)\n",
  //  deltapx, deltapy, deltapz);
  //fprintf(print_results,"deltaT = %8.2f MeV\n", deltaT);
  //fclose(print_results);
  
  //printf("   \n");
}

