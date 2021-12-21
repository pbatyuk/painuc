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

PPUINT Alloc_Buffer(int hsize, int vsize);
void Free_Buffer(PPUINT buffer[2], int vsize);
PPUINT pCDBuffer[2];

int xtrack1[2][6][10000], ytrack1[2][6][10000];
int np[2][6];
float angle[2][6];
float rc_int, xc_int, yc_int;
float stx, sty, endx, endy;
float xtrack[2][6][10000], ytrack[2][6][10000];
int   particle[10];
//float Wpoints[20000], *Wpoint;
float xcentr, ycentr;
float xcent_vx[10], ycent_vx[10], rcent_vx[10];
float Xvertax, Yvertax, Dvertax;
float brrr[2][6];
float deltaV, disp;

//Minuit's variables
float Xcircle[50000], Ycircle[50000], Zcircle[50000], Wpoint[50000];
short int Signe;
float XC, YC, RC, AC, BC, chi2, DRC, DAC, phif, t[50000], bufmax, bufmin;
float fcn_print;
int error_flag, is;
int Npoints_fit;

//for Gil's graphics
TCanvas *c3, *c4, *cp, *c2, *cw, *csl;
TPad *plc3, *prc3, *pMc3;
TPad *plc4, *prc4, *pMc4;
TPad *cp0, *cp1, *hp0, *cw0, *csl0;
TProfile *hprof, *hprof2, *hprof3, *hprof4;
//TText *t;
//TH1D *hprof;
TGraph *gr;
TPaveLabel *pl, *pl1;
TPaveLabel *pl11, *pl2, *ptext;
TMarker *mark, *mark2;
TPolyMarker *markp, *markp2;
Float_t *px_mark, *py_mark;
Float_t x_mark[50000], y_mark[50000];
const char *filetypes[] = { 
  "Left images",    "[L,l]*",
  "left images",    "l*",
  "All files",     "*",
  0,               0 
};
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
char gifname[100];
char eventname[20];
char *tokenPtr;
char iniDir[100];
// = {"/data/koval/DubTo/work_2003/restore/data/"};
char WorkDir[100];
int minl=0, maxl=4096, minr=0, maxr=4096;
int xlmin=108, xlmax=500, ylmin=68, ylmax=443;
int xrmin=128, xrmax=520, yrmin=52, yrmax=414;
//char ofilename[100];
double IntBr[2][6];
FILE *for_reconstr;
FILE *print_results;
TCanvas *ch;
TH1F *htest,*hbr, *hwidth, *hslice;
TH2F *left, *right;
gzFile im_right, im_left;                                                       
TGFileInfo* fi = NULL;                                                                  
TObjectTable *fObjTbl;                                                          
TPad *pleft, *pright, *ptmp, *pMsg;
/*                                                */
/*               Function prototypes              */
/*                                                */

void SetVeryNewPalette(void);
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
void Find_Tails(int itr, int image);

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
void find_points_close_to_circle(TProfile *hprof,
				 float fmin, float fmax, float xc, float yc, float rc, 
				 float rc0,
				 float x1, float y1, float x2, float y2, int ima, int iwidth);
void find_br_for_circlefit(TProfile *hprof,
			   float fmin, float fmax, float xc, float yc, float rc, 
			   float rc0,
			   float x1, float y1, float x2, float y2, int ima, int iwidth);
short int find_the_same_point(int i, int j, int np);
void find_widthofthecorridor(double *fmin, double *fmax, float xc, float yc, float rc, 
			      float x1, float y1, float x2, float y2, int ima, int itr, float *sigma, float *shift);

// void Pause(void);
int  gprintf(int *xloc, int *yloc, char *fmt, ... );

/*							*/
/*	Begin main function				*/
/*							*/

//extern unsigned _stklen = 16384;
float ratio;
extern const int vsize = 517;
extern const int hsize = 658;
extern const float RTarget = 600.;
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

TROOT simple("simple","Test of histogramming and I/O");

enum ETestCommandIdentifiers {
   M_FILE_OPEN,
   M_FILE_SAVE,
   M_FILE_SAVEAS,
   M_FILE_EXIT,

   M_TEST_DLG,
   M_TEST_MSGBOX,
   M_TEST_SLIDER,

   M_HELP_CONTENTS,
   M_HELP_SEARCH,
   M_HELP_ABOUT,

   M_CASCADE_1,
   M_CASCADE_2,
   M_CASCADE_3,

   VId1,
   HId1,
   VId2,
   HId2,

   VSId1,
   HSId1,
   VSId2,
   HSId2
};


Int_t mb_button_id[9] = { kMBYes, kMBNo, kMBOk, kMBApply,
                          kMBRetry, kMBIgnore, kMBCancel,
                          kMBClose, kMBDismiss };

EMsgBoxIcon mb_icon[4] = { kMBIconStop, kMBIconQuestion,
                           kMBIconExclamation, kMBIconAsterisk };

class TestDialog : public TGTransientFrame {

private:
  TGCompositeFrame    *fFrame1, *fF1, *fF2, *fF3, *fF4, *fF5, *fF6, *fF7;
  TGButton            *fOkButton, *fCancelButton, *fStartB, *fStopB;
  TGButton            *fRad1, *fRad2, *fRad3, *fRad4, *fRad5;
  TGPictureButton     *fPicBut1;
  //   TGRadioButton       *fRadio1, *fRadio2;
  TGRadioButton       *fR1[6], *fR2[6], *fR3[6], *fR4[6], *fR5[6], *fR6[6];
  TGCheckButton       *fCheck1;
  TGCheckButton       *fCheckMulti;
  TGListBox           *fListBox;
  TGComboBox          *fCombo;
  TGTab               *fTab;
  TGTextEntry         *fTxt1, *fTxt2;
  TGLayoutHints       *fL1, *fL2, *fL3, *fL4;
  TRootEmbeddedCanvas *fEc1, *fEc2;
  Int_t                fFirstEntry;
  Int_t                fLastEntry;
  TH1F                *fHpx;
  TH2F                *fHpxpy;
  TList               *fCleanup;

public:
   TestDialog(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
               UInt_t options = kMainFrame | kVerticalFrame);
   virtual ~TestDialog();

   virtual void CloseWindow();
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
};


TestDialog::TestDialog(const TGWindow *p, const TGWindow *main, UInt_t w,
                       UInt_t h, UInt_t options)
    : TGTransientFrame(p, main, w, h, options)
{
   // Create a dialog window. A dialog window pops up with respect to its
   // "main" window.

   // Used to store GUI elements that need to be deleted in the ctor.
   fCleanup = new TList;

   fFrame1 = new TGHorizontalFrame(this, 60, 20, kFixedWidth);

   fOkButton = new TGTextButton(fFrame1, "&Ok", 1);
   fOkButton->Associate(this);
   fCancelButton = new TGTextButton(fFrame1, "&Cancel", 2);
   fCancelButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                           2, 2, 2, 2);
   fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1);

   fFrame1->AddFrame(fOkButton, fL1);
   fFrame1->AddFrame(fCancelButton, fL1);

   fFrame1->Resize(150, fOkButton->GetDefaultHeight());
   AddFrame(fFrame1, fL2);

   int i;
   //--------- create Tab widget and some composite frames for Tab testing

   fTab = new TGTab(this, 300, 300);
   fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

   TGCompositeFrame *tf = fTab->AddTab("Track 0");
   fF1 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
   fF1->AddFrame(fR1[0] = new TGRadioButton(fF1, "pion",    11), fL1);
   fF1->AddFrame(fR1[1] = new TGRadioButton(fF1, "proton",  12), fL1);
   fF1->AddFrame(fR1[2] = new TGRadioButton(fF1, "deuteron",13), fL1);
   fF1->AddFrame(fR1[3] = new TGRadioButton(fF1, "tritium", 14), fL1);
   fF1->AddFrame(fR1[4] = new TGRadioButton(fF1, "He3",     15), fL1);
   fF1->AddFrame(fR1[5] = new TGRadioButton(fF1, "He4",     16), fL1);
   tf->AddFrame(fF1, fL3);

   for (i = 0; i < 6; ++i) 
     {
       fR1[i]->Associate(this);
     }

   fR1[0]->SetState(kButtonDown);

   if(tracks >= 2) {
     tf = fTab->AddTab("Track 1");
     fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			     200, 2, 2, 2);
     fF2 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
     fF2->AddFrame(fR2[0] = new TGRadioButton(fF2, "pion",    21), fL1);
     fF2->AddFrame(fR2[1] = new TGRadioButton(fF2, "proton",  22), fL1);
     fF2->AddFrame(fR2[2] = new TGRadioButton(fF2, "deuteron",23), fL1);
     fF2->AddFrame(fR2[3] = new TGRadioButton(fF2, "tritium", 24), fL1);
     fF2->AddFrame(fR2[4] = new TGRadioButton(fF2, "He3",     25), fL1);
     fF2->AddFrame(fR2[5] = new TGRadioButton(fF2, "He4",     26), fL1);
     tf->AddFrame(fF2, fL3);
     
     for (i = 0; i < 6; ++i) {
       fR2[i]->Associate(this);
     }
     
     fR2[0]->SetState(kButtonDown);
   }
   
   if(tracks >= 3) {
     tf = fTab->AddTab("Track 2");
     fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			     200, 3, 3, 3);
     fF3 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
     fF3->AddFrame(fR3[0] = new TGRadioButton(fF3, "pion",    31), fL1);
     fF3->AddFrame(fR3[1] = new TGRadioButton(fF3, "proton",  32), fL1);
     fF3->AddFrame(fR3[2] = new TGRadioButton(fF3, "deuteron",33), fL1);
     fF3->AddFrame(fR3[3] = new TGRadioButton(fF3, "tritium", 34), fL1);
     fF3->AddFrame(fR3[4] = new TGRadioButton(fF3, "He3",     35), fL1);
     fF3->AddFrame(fR3[5] = new TGRadioButton(fF3, "He4",     36), fL1);
     tf->AddFrame(fF3, fL3);
     
     for (i = 0; i < 6; ++i) {
       fR3[i]->Associate(this);
     }
     
     fR3[5]->SetState(kButtonDown);
   }
   if(tracks >= 4) {
     tf = fTab->AddTab("Track 3");
     fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			     200, 4, 4, 4);
     fF4 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
     fF4->AddFrame(fR4[0] = new TGRadioButton(fF4, "pion",    41), fL1);
     fF4->AddFrame(fR4[1] = new TGRadioButton(fF4, "proton",  42), fL1);
     fF4->AddFrame(fR4[2] = new TGRadioButton(fF4, "deuteron",43), fL1);
     fF4->AddFrame(fR4[3] = new TGRadioButton(fF4, "tritium", 44), fL1);
     fF4->AddFrame(fR4[4] = new TGRadioButton(fF4, "He3",     45), fL1);
     fF4->AddFrame(fR4[5] = new TGRadioButton(fF4, "He4",     46), fL1);
     tf->AddFrame(fF4, fL3);
     
     for (i = 0; i < 6; ++i) {
       fR4[i]->Associate(this);
     }
     
     fR4[0]->SetState(kButtonDown);
   }
   if(tracks >= 5) {
     tf = fTab->AddTab("Track 4");
     fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			     200, 4, 4, 4);
     fF5 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
     fF5->AddFrame(fR5[0] = new TGRadioButton(fF5, "pion",    51), fL1);
     fF5->AddFrame(fR5[1] = new TGRadioButton(fF5, "proton",  52), fL1);
     fF5->AddFrame(fR5[2] = new TGRadioButton(fF5, "deuteron",53), fL1);
     fF5->AddFrame(fR5[3] = new TGRadioButton(fF5, "tritium", 54), fL1);
     fF5->AddFrame(fR5[4] = new TGRadioButton(fF5, "He3",     55), fL1);
     fF5->AddFrame(fR5[5] = new TGRadioButton(fF5, "He4",     56), fL1);
     tf->AddFrame(fF5, fL3);
     
     for (i = 0; i < 6; ++i) {
       fR5[i]->Associate(this);
     }
   }
   if(tracks >= 6) {
     tf = fTab->AddTab("Track 5");
     fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
			     200, 4, 4, 4);
     fF6 = new TGCompositeFrame(tf, 60, 20, kVerticalFrame);
     fF6->AddFrame(fR6[0] = new TGRadioButton(fF6, "pion",    61), fL1);
     fF6->AddFrame(fR6[1] = new TGRadioButton(fF6, "proton",  62), fL1);
     fF6->AddFrame(fR6[2] = new TGRadioButton(fF6, "deuteron",63), fL1);
     fF6->AddFrame(fR6[3] = new TGRadioButton(fF6, "tritium", 64), fL1);
     fF6->AddFrame(fR6[4] = new TGRadioButton(fF6, "He3",     65), fL1);
     fF6->AddFrame(fR6[5] = new TGRadioButton(fF6, "He4",     66), fL1);
     tf->AddFrame(fF6, fL3);
     
     for (i = 0; i < 6; ++i) {
       fR6[i]->Associate(this);
     }
   }
   //--- end of last tab

   TGLayoutHints *fL5 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1);
   AddFrame(fTab, fL5);

   MapSubwindows();
   Resize(GetDefaultSize());

   // position relative to the parent's window
   Window_t wdum;
   int ax, ay;
   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
                          (((TGFrame *) main)->GetWidth() - fWidth) >> 1,
                          (((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                          ax, ay, wdum);
   Move(ax, ay);

   SetWindowName("Dialog");

   MapWindow();
   gClient->WaitFor(this);    // otherwise canvas contextmenu does not work
}

TestDialog::~TestDialog()
{
   // Delete test dialog widgets.

   fCleanup->Delete();
   delete fCleanup;
   delete fOkButton;
   delete fCancelButton;
   delete fStartB; delete fStopB;
   delete fFrame1;
   delete fRad1; delete fRad2; delete fRad3; delete fRad4; delete fRad5;
   delete fR1[0]; delete fR1[1]; delete fR1[2]; 
   delete fR1[3]; delete fR1[4]; delete fR1[5];
   delete fR2[0]; delete fR2[1]; delete fR2[2]; 
   delete fR2[3]; delete fR2[4]; delete fR2[5];
   delete fR3[0]; delete fR3[1]; delete fR3[2]; 
   delete fR3[3]; delete fR3[4]; delete fR3[5];
   delete fR4[0]; delete fR4[1]; delete fR4[2]; 
   delete fR4[3]; delete fR4[4]; delete fR4[5];
   delete fR5[0]; delete fR5[1]; delete fR5[2]; 
   delete fR5[3]; delete fR5[4]; delete fR5[5];
   delete fR6[0]; delete fR6[1]; delete fR6[2]; 
   delete fR6[3]; delete fR6[4]; delete fR6[5];
   delete fF1; delete fF2; delete fF3; delete fF4; delete fF5;
   delete fF6; delete fF7;
   delete fL3; delete fL4;
   delete fL1; delete fL2;
}

void TestDialog::CloseWindow()
{
   // Called when window is closed via the window manager.

   delete this;
}

Bool_t TestDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   // Process messages coming from widgets associated with the dialog.

  //char tmp[20];
  //static int newtab = 0;
   int i, itr;

   switch (GET_MSG(msg)) {
      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch(parm1) {
                  case 1:
                  case 2:
                     printf("\nTerminating dialog: %s pressed\n",
                            (parm1 == 1) ? "OK" : "Cancel");
                     CloseWindow();
                     break;
                  default:
                     break;
               }
               break;
	 case kCM_RADIOBUTTON:
	   if(parm1>10 && parm1 <20)
	     {
	       for (i=0; i<6; ++i)
		 {
		   fR1[i]->SetState(kButtonUp); 
		   if((parm1 % 10) - 1 == i ){
		     fR1[i]->SetState(kButtonDown);
		     itr = parm1/10 - 1;
		     particle[itr] = i;
		   }
		 }
	     }
	   if(parm1>20 && parm1 <30)
	     {
	       for (i=0; i<6; ++i)
		 {
		   fR2[i]->SetState(kButtonUp); 
		   if((parm1 % 10) - 1 == i ){
		     fR2[i]->SetState(kButtonDown);
		     itr = parm1/10 - 1;
		     particle[itr] = i;
		   }
		 }
	     }
	   if(parm1>30 && parm1 <40)
	     {
	       for (i=0; i<6; ++i)
		 {
		   fR3[i]->SetState(kButtonUp); 
		   if((parm1 % 10) - 1 == i ){
		     fR3[i]->SetState(kButtonDown);
		     itr = parm1/10 - 1;
		     particle[itr] = i;
		   }
		 }
	     }
	   if(parm1>40 && parm1 <50)
	     {
	       for (i=0; i<6; ++i)
		 {
		   fR4[i]->SetState(kButtonUp); 
		   if((parm1 % 10) - 1 == i ){
		     fR4[i]->SetState(kButtonDown);
		     itr = parm1/10 - 1;
		     particle[itr] = i;
		   }
		 }
	     }
	   if(parm1>50 && parm1 <60)
	     {
	       for (i=0; i<6; ++i)
		 {
		   fR5[i]->SetState(kButtonUp); 
		   if((parm1 % 10) - 1 == i ){
		     fR5[i]->SetState(kButtonDown);
		     itr = parm1/10 - 1;
		     particle[itr] = i;
		   }
		 }
	     }
	   if(parm1>60 && parm1 <70)
	     {
	       for (i=0; i<6; ++i)
		 {
		   fR6[i]->SetState(kButtonUp); 
		   if((parm1 % 10) - 1 == i ){
		     fR6[i]->SetState(kButtonDown);
		     itr = parm1/10 - 1;
		     particle[itr] = i;
		   }
		 }
	     }
	   break;
            default:
               break;
         }
         break;

      default:
         break;
   }
   return kTRUE;
}



int main(int argc, char** argv)
  //int main(int argc, char* argv[])
//int gil_main(char ofilename[], int ArrSize)
{
  fObjTbl = new TObjectTable(10000);
  fi = new TGFileInfo;
  int i, ii, j, k, l, kk;
  char ifilename[100];
  char ofilename[100];
  
  float pi;
  
  strcpy(ifilename, argv[1]);
  //deltaV = atof(argv[2]); 
  TApplication theApp("App", &argc, argv);

  exp_lens = 1.0;
  
  printf("Now open input file %s\n", ifilename);
  if ((for_reconstr = fopen(ifilename,"rt")) == NULL)
    {
      fprintf(stderr, "Cannot open input file %s\n", ifilename);
      return 1;
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
  SetVeryNewPalette();
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
  
  std::cout <<"x , y ="<<ixWindow<<" "<<iyWindow<<std::endl;
  
  c3 = new TCanvas("c3","Event reconstruction",20,20,ixWindow,iyWindow);
  plc3 = new TPad("plc3","Gil's picture", xleft, yleft, xright, yright, 21);
  plc3->Range(0, 0, hsize, vsize);
  printf("Range: xmax = %f, ymax = %f\n", rhsize, rvsize);
  prc3 = new TPad("prc3","In the space", xleft + 0.5, yleft, xright + 0.5, yright, 21);
  c3->ToggleEventStatus();
  
  plc3->Draw();
  plc3->UseCurrentStyle();
  
  prc3->Draw();
  prc3->UseCurrentStyle();
  
  plc3->cd();
  plc3->Draw("col");
  
  prc3->cd();
  prc3->Draw("col");
  
  c3->Update();
  
  plc3->cd();
  
  float xmax, ymax, xmax_tr;
  xmax = ymax = 0.;
  float xmin, ymin, xmin_tr;
  xmin = ymin = 1000000.;
  float xl[2][6], xr[2][6];
  xcentr = xaxis[0]/(2.*pixel2mm);
  ycentr = yaxis[0]/(2.*pixel2mm); 
  if(xcentr > xmax) xmax = xcentr;
  if(xcentr < xmin) xmin = xcentr;
  if(ycentr > ymax) ymax = ycentr;
  if(ycentr < ymin) ymin = ycentr;
  
  TArc *myarc = new TArc(xcentr,ycentr,9.,0.,360.);
  myarc->Draw();
  plc3->Modified();
  plc3->Update();
  ycentr = yaxis[0]/(2.*pixel2mm);
  for(kk=0; kk<2; kk++)
    {
      for(i=0; i<tracks; i++)
	{
	  xmax_tr = -1000000.;
	  xmin_tr =  1000000.;
	  
	  for(j=0; j<np[kk][i]; j++)
	    {
	      xx = xtrack[kk][i][j]/2. + xcentr;
	      if(kk==1) xx += xtransfer/pixel2mm/2.;
	      yy = ytrack[kk][i][j]/2. + ycentr; 
	      if(kk==1) yy += ytransfer/pixel2mm/2.;
	      if(xx > xmax) xmax = xx;
	      if(xx < xmin) xmin = xx;
	      if(yy > ymax) ymax = yy;
	      if(yy < ymin) ymin = yy;
	      if(kk==0) {
		if(xtrack[kk][i][j] > xmax_tr) xmax_tr = xtrack[kk][i][j];
		if(xtrack[kk][i][j] < xmin_tr) xmin_tr = xtrack[kk][i][j];}
	      if(kk==1) {
		if(xtrack[kk][i][j] > xmax_tr) xmax_tr = xtrack[kk][i][j] + xtransfer;
		if(xtrack[kk][i][j] < xmin_tr) xmin_tr = xtrack[kk][i][j] + xtransfer;}

	      xl[kk][i] = xmin_tr;
	      xr[kk][i] = xmax_tr;
	      if(xx != 0. && yy != 0.)
		{
		  mark = new TMarker((Float_t) xx, (Float_t) yy, 5);
		  mark->SetMarkerColor(51);
		  mark->SetMarkerStyle(8);
		  mark->SetMarkerSize(0.5);
		  fObjTbl->Add(mark);
		  mark->DrawMarker((Float_t) xx, (Float_t) yy);
		}
	    }
	}
      plc3->Modified();
      plc3->Update();
      if(kk==0) mark = new TMarker((Float_t) (Xvertex[kk]/2. + xcentr),
				    (Float_t) (Yvertex[kk]/2. + ycentr), 5);
      if(kk==1) mark = new TMarker((Float_t) (Xvertex[kk]/2. + xcentr + xtransfer/2./pixel2mm),
				    (Float_t) (Yvertex[kk]/2. + ycentr + ytransfer/2./pixel2mm), 5);
      mark->SetMarkerColor(50);
      mark->SetMarkerStyle(8);
      mark->SetMarkerSize(1.);
      fObjTbl->Add(mark);
      if(kk==0) mark->DrawMarker((Float_t) (Xvertex[kk]/2. + xcentr),
				 (Float_t) (Yvertex[kk]/2. + ycentr));
      if(kk==1) mark->DrawMarker((Float_t) (Xvertex[kk]/2. + xcentr + xtransfer/2./pixel2mm),
				 (Float_t) (Yvertex[kk]/2. + ycentr + ytransfer/2./pixel2mm));
    }
  plc3->Modified();
  plc3->Update();
  for(kk=0; kk<2; kk++)
    {
      for(i=0; i<=5; i++)
	{
	  xx = (xfilmd_top[kk][i] + xaxis_draw[0])/(2.*pixel2mm);
	  yy = (yfilmd_top[kk][i] + yaxis_draw[0])/(2.*pixel2mm); 
	  if(xx > xmax) xmax = xx;
	  if(xx < xmin) xmin = xx;
	  if(yy > ymax) ymax = yy;
	  if(yy < ymin) ymin = yy;
	  
	  if(yy > 480.) yy = 460.;
	  mark = new TMarker((Float_t) xx, (Float_t) yy, 5);
	  mark->SetMarkerColor(2);
	  mark->SetMarkerStyle(2);
	  mark->SetMarkerSize(1.);
	  fObjTbl->Add(mark);
	  mark->DrawMarker((Float_t) xx, (Float_t) yy);
	}
      for(i=0; i<=3; i++)
	{
	  xx =  xfilmd_bot[kk][i]/(2.*pixel2mm) + xcentr;
	  yy =  yfilmd_bot[kk][i]/(2.*pixel2mm) + ycentr; 
	  if(xx > xmax) xmax = xx;
	  if(xx < xmin) xmin = xx;
	  if(yy > ymax) ymax = yy;
	  if(yy < ymin) ymin = yy;
	  mark = new TMarker((Float_t) xx, (Float_t) yy, 5);
	  mark->SetMarkerColor(46);
	  mark->SetMarkerStyle(2);
	  mark->SetMarkerSize(1.);
	  fObjTbl->Add(mark);
	  mark->DrawMarker((Float_t) xx, (Float_t) yy);
	}
    }
  printf("Min x %8.3f, y %8.3f, Max x %8.3f, y %8.3f\n", xmin, ymin, xmax, ymax);
  plc3->Draw();
  plc3->Modified();
  plc3->Update();
  
  prc3->cd();
  
  TView *vview = new TView(1);
  const double ychamb = 500.;
  const double xchamb = ychamb/ratio/2.;
  vview->SetRange(-xchamb,-100.,-(zbot+100.),xchamb,-100.+ychamb,-(ztop-100.));
  TPolyMarker3D *pm3d = new TPolyMarker3D(50000);
  
  double x0, y0, y1;
  int track;
  l = 0;
  k = 0;
  double base = 253.;
  //249.105174619086522;//249.053515258138968;//249.13062783922717;//251.;
  double xx, yy, phirotate;
  phirotate = 0.;
  //  reconstruct the positions of the crosses
  //double deltaV; = -1.5*pixel2mm;
  deltaV = 0.;//-0.024;//(double) (xtrack[1][0][0]-xtransfer-xtrack[0][0][0])*pixel2mm;
  for(i=0; i<=5; i++)
    {
      double azy_line[2];
      double axy_line[2];
      double zp1, zp2, xp1, xp2, yp1, yp2;
      zp1 = zp2 = xp1 = xp2 = yp1 = yp2 = 0.;
      double xtrack_mm[2], ytrack_mm[2];
      double xt, yt, a1, b1, a2, b2, z0;
      double x01, y01, z01, x02, y02, z02;
      x01 = y01 = z01 = x02 = y02 = z02 = x0 = y0 = z0 = 0.;
      xtrack_mm[0] = (double) xfilm_top[0][i];
      ytrack_mm[0] = (double) yfilm_top[0][i];
      // rotate all the points on phirotate angle 
      // to be in the system perpendicular 
      // to the base axis
      xt = xtrack_mm[0];
      yt = ytrack_mm[0];
      azy_line[0] = yt/zfilm[0];
      axy_line[0] = yt/xt;
      xtrack_mm[1] = (double) xfilm_top[1][i];
      ytrack_mm[1] = (double) (yfilm_top[1][i] + base);//- ytransfer) +base;
      //check for the magnifications
      printf("Cross %d xtrack1/xtrack2 %f mag2/mag1 %f zfilm2/zfilm1 %f \n",
	     i, xtrack_mm[1]/xtrack_mm[0], mag_top[0]/mag_top[1],
	     zfilm[1]/zfilm[0]);
      xt = xtrack_mm[1];
      yt = ytrack_mm[1];
      a1 = azy_line[0]; b1 = 0.;
      a2 = (yt-base)/zfilm[1]; b2 = base;
      z0 = (b2 - b1)/(a1 - a2);
      y0 = a1*z0 + b1;
      a1 = axy_line[0]; b1 = 0.;
      a2 = (yt-base)/xt; b2 = base;
      x0 = (b2 - b1)/(a1 - a2);
      y0 = a1*x0 + b1;
      
      if((x0 !=0. && y0 != 0. && z0 != 0.) )
	//&& 
	// ((float) fabs(x0) <= xchamb && (float) fabs(y0-base) <= ychamb))
	{
	  xcross3D[l] = (float) x0;
	  ycross3D[l] = (float) y0;
	  zcross3D[l] = (float) z0;
	  printf("Cross %d coordinates: x %f, y %f z %f\n", l, x0, y0, z0);
	  l++;
	}
      else
	printf("%d Cross is outside of chamber space:x %6.3f y%6.3f z%8.3f\n", i, x0, y0, z0);
      
    }
  for(i=0; i<=3; i++)
    {

      double azy_line[2];
      double axy_line[2];
      double zp1, zp2, xp1, xp2, yp1, yp2;
      zp1 = zp2 = xp1 = xp2 = yp1 = yp2 = 0.;
      double xtrack_mm[2], ytrack_mm[2];
      double xt, yt, a1, b1, a2, b2, z0;
      double x01, y01, z01, x02, y02, z02;
      x01 = y01 = z01 = x02 = y02 = z02 = x0 = y0 = z0 = 0.;
      xtrack_mm[0] = (double) xfilm_bot[0][i];
      ytrack_mm[0] = (double) yfilm_bot[0][i];
      // rotate all the points on phirotate angle 
      // to be in the system perpendicular 
      // to the base axis
      xt = xtrack_mm[0];
      yt = ytrack_mm[0];
      azy_line[0] = yt/zfilm[0];
      axy_line[0] = yt/xt;
      xtrack_mm[1] = (double) xfilm_bot[1][i];
      ytrack_mm[1] = (double) (yfilm_bot[1][i] + base);//- ytransfer);
      printf("i %d xtrack1/xtrack2 %f mag2/mag1 %f zfilm2/zfilm1 %f \n",
	     i, xtrack_mm[1]/xtrack_mm[0], mag_bot[0]/mag_bot[1],
	     zfilm[1]/zfilm[0]);
      xt = xtrack_mm[1];
      yt = ytrack_mm[1];
      a1 = azy_line[0]; b1 = 0.;
      a2 = (yt-base)/zfilm[1]; b2 = base;
      z0 = (b2 - b1)/(a1 - a2);
      y0 = a1*z0 + b1;
      a1 = axy_line[0]; b1 = 0.;
      a2 = (yt-base)/xt; b2 = base;
      x0 = (b2 - b1)/(a1 - a2);
      y0 = a1*x0 + b1;
      if((x0 !=0. && y0 != 0. && z0 != 0.))
	//&& 
	// ((float) fabs(x0) <= xchamb && (float) fabs(y0-base) <= ychamb))
	{
	  xcross3D[l] = (float) x0;
	  ycross3D[l] = (float) y0;
	  zcross3D[l] = (float) z0;
	  printf("Cross %d coordinates: x %f, y %f z %f\n", l, x0, y0, z0);
	  l++;
	}
      else
	printf("%d Cross is outside of chamber space:x %6.3f y%6.3f z%8.3f\n", i, x0, y0, z0);
      
    }
  int ncrosses = l;
  int no_fit[6];
  for(track=0; track<tracks; track++)
    {
      no_fit[track] = 0;
      npt3D[track] = 0;
      npt3D_left[track] = 0;
      npt3D_right[track] = 0;
      int l = 0;
      // go to the coord.system  aligned to the "left" image then to "right" one
      //radius, coordinates of center for each image in one system
      double rtr0, rtr1, xc0, xc1, yc0, yc1;

      rtr0 = (double) (Radius0[0][track]*pixel2mm);
      xc0 = (double) (xcenter01[0][track]*pixel2mm);
      yc0 = (double) (ycenter01[0][track]*pixel2mm);

      rtr1 = (double) (Radius0[1][track]*pixel2mm);
      xc1 = (double) (xcenter01[1][track]*pixel2mm + xtransfer);
      yc1 = (double) (ycenter01[1][track]*pixel2mm + base);

      
      printf("rtr0 %f, xc0 %f, yc0 %f, rtr1 %f, xc1 %f, yc1 %f\n", rtr0,
           xc0, yc0, rtr1, xc1, yc1);
      //printf(" xl[1][track] %f xr %f\n", xl[1][track], xr[1][track]);
      //int buf;
      //cin >> buf;
      int nbad = 0;
      for(i=0; i<np[0][track]; i++)
	{
	  //if point lies between min-max bounds on X direction
	  //and if track goes not very along Y-axis
	  if((xtrack[0][track][i] < xl[1][track] ||
	     xtrack[0][track][i] > xr[1][track]) &&
	     (xr[1][track] - xl[1][track]) > 50. && i != 0) continue;
	  double azy_line[2];
	  double axy_line[2];
	  double zp1, zp2, xp1, xp2, yp1, yp2;
	  double xpim1, xpim0, ypim1, ypim0;
	  zp1 = zp2 = xp1 = xp2 = yp1 = yp2 = 0.;
	  double y2;
	  double xtrack_mm[2], ytrack_mm[2];
	  double xt0, xt1, a1, b1, a2, b2, z0;
	  double x01, y01, z01, x02, y02, z02;
	  x01 = y01 = z01 = x02 = y02 = z02 = x0 = y0 = z0 = 0.;
	  xtrack_mm[0] = (double) xtrack[0][track][i]*pixel2mm;// +deltaV/2.;
	  ytrack_mm[0] = (double) ytrack[0][track][i]*pixel2mm;
	  xpim0 = xtrack_mm[0];
	  ypim0 = ytrack_mm[0];
	  if(fabs(xpim0) < 0.1 && i != 0) continue;

	  azy_line[0] = ypim0/zfilm[0];
	  axy_line[0] = ypim0/xpim0;
	  
	  xpim1 = xpim0*zfilm[1]/zfilm[0];
	  double xpim1_pxl = xpim1/pixel2mm;
	  double ypim1buf = ypim0*zfilm[1]/zfilm[0];
	  int npts = 0;
	  for (j = 0; j < np[1][track]; j++)
	    {
	      if(fabs(xpim1_pxl - xtrack[1][track][j]) < 0.5) 
		{
		  if (ytransfer < ypim1buf - ytrack[1][track][j]*pixel2mm &&
		      ypim1buf - ytrack[1][track][j]*pixel2mm < ytransfer2) 
		    {
		      //printf("ytransfer %f ypim1buf - ytrack[1][track][j] %f ytransfer2 %f\n",
		      //     ytransfer, ypim1buf - ytrack[1][track][j]*pixel2mm, ytransfer2);
		      xpim1 = xtrack[1][track][j]*pixel2mm;
		      if(fabs(xpim1) < 0.1) continue;
		      ypim1 = ytrack[1][track][j]*pixel2mm + base;
		      a1 = azy_line[0]; b1 = 0.;
		      a2 = (ypim1-base)/zfilm[1]; b2 = base;
		      if(fabs(a1-a2)<0.001) continue;
		      z0 = (b2 - b1)/(a1 - a2);
		      y0 = a1*z0 + b1;
		      a1 = axy_line[0]; b1 = 0.;
		      a2 = (ypim1-base)/xpim1; b2 = base;
		      if(fabs(a1-a2)<0.001) continue;
		      x0 = (b2 - b1)/(a1 - a2);
		      y0 = a1*x0 + b1;
		      pm3d->SetPoint( k, x0, y0, z0 );
		      x3D[track][l] = (float) x0;
		      y3D[track][l] = (float) y0;
		      z3D[track][l] = (float) z0;
		      wes[track][l] = (float) (CCDBuffer[0][track][i] + CCDBuffer[1][track][j])/2.;
		      k++; l++;
		      npt3D[track]++;
		      npt3D_left[track]++;
		      npts++;
		    }
		}
	    }			

	  if(npts == 0 && ((i < np[0][track] - 1 && i > 1) || i == 0)) {
	  //add one point on the circle
	    double phibuf;
	    if(fabs(xpim0 - xc0) > 0.001) 
	      phibuf = atan2(ypim0 - yc0, xpim0 - xc0);
	    else
	      phibuf = TMath::Pi()/2.;
	    xt0 = rtr0*cos(phibuf)+xc0;
	    xpim1 = xpim0*zfilm[1]/zfilm[0];
	    xt1 = xpim1; 
			//xt0*zfilm[1]/zfilm[0];
	    if(fabs(xt1 - xc1) <= (rtr1-10.*pixel2mm) && 
	       fabs(xt1 - xc1) > 0.001 )
	      {//there is an intersection
		y1 = yc1 + sqrt(rtr1*rtr1 - (xt1 - xc1)*(xt1 - xc1)); 
		y2 = yc1 - sqrt(rtr1*rtr1 - (xt1 - xc1)*(xt1 - xc1)); 
		//printf(" Image 2 xt %f, y1 %f, y2 %f\n", xt, y1, y2);
		ypim1 = y1;
		//(y1 - yc1)/(xt1 - xc1)*xpim1 + 
		//(y1 - (y1 - yc1)/(xt1 - xc1)*xt1);
		a1 = azy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/zfilm[1]; b2 = base;
		if(fabs(a2-a1) < 0.001) 
		  {if(i>0)continue;
		  goto VERT0;}
		z01 = (b2 - b1)/(a1 - a2);
		y01 = a1*z01 + b1;
		a1 = axy_line[0]; b1 = 0.;
		if(fabs(xpim1) < 0.001) 
		  {if(i>0)continue;
		  goto VERT0;}
		a2 = (ypim1-base)/xpim1; b2 = base;
		if(fabs(a2-a1) < 0.001) 
		  {if(i>0)continue;
		  goto VERT0;}
		x01 = (b2 - b1)/(a1 - a2);
		y01 = a1*x01 + b1;
		ypim1 = y2;
		//(y2 - yc1)/(xt1 - xc1)*xpim1 + 
		//(y2 - (y2 - yc1)/(xt1 - xc1)*xt1);
		a1 = azy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/zfilm[1]; b2 =base;
		if(fabs(a2-a1) < 0.0001) 
		  {if(i>0)continue;
		  goto VERT0;}
		z02 = (b2 - b1)/(a1 - a2);
		y02 = a1*z02 + b1;
		a1 = axy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/xpim1; b2 = base;
		if(fabs(a2-a1) < 0.0001) 
		  {if(i>0)continue;
		  goto VERT0;}
		x02 = (b2 - b1)/(a1 - a2);
		y02 = a1*x02 + b1;
		if(z01 < -(ztop-100.) && z01 > -(zbot+100.)) {x0 = x01; y0 = y01; z0 = z01;}
		if(z02 < -(ztop-100.) && z02 > -(zbot+100.)) {x0 = x02; y0 = y02; z0 = z02;}
		if(x0 == 0. && y0 == 0. && z0 == 0.)
		  {
		    float distz_1, distz_2, distz_1top, distz_1bot, distz_2top, distz_2bot;
		    distz_1top = fabs(z01 + ztop); distz_1bot = fabs(z01 + zbot);
		    distz_2top = fabs(z02 + ztop); distz_2bot = fabs(z02 + zbot);
		    if(distz_1top < distz_1bot) distz_1 = distz_1top;
		    if(distz_1top > distz_1bot) distz_1 = distz_1bot;
		    if(distz_2top < distz_2bot) distz_2 = distz_2top;
		    if(distz_2top > distz_2bot) distz_2 = distz_2bot;
		    if(distz_1 < distz_2) {x0 = x01; y0 = y01; z0 = z01;}
		    if(distz_1 > distz_2) {x0 = x02; y0 = y02; z0 = z02;}
		  }
		//printf("x01 %6.3f y01 %6.3f z01 %8.3f\n", x01, y01, z01);
		//printf("x02 %6.3f y02 %6.3f z02 %8.3f\n", x02, y02, z02);
		//printf("x0 %6.3f y0 %6.3f z0 %8.3f\n", x0, y0, z0);
	      }
	    else if (fabs(xt1 - xc1) <= 0.001 )
	      {
		y1 = yc1 + rtr1; 
		y2 = yc1 - rtr1; 
		//printf(" Image 2 xt %f, y1 %f, y2 %f\n", xt, y1, y2);
		ypim1 = y1;
		a1 = azy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/zfilm[1]; b2 = base;
		z01 = (b2 - b1)/(a1 - a2);
		y01 = a1*z01 + b1;
		a1 = axy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/xpim1; b2 = base;
		x01 = (b2 - b1)/(a1 - a2);
		y01 = a1*x01 + b1;
		ypim1 = y2;
		a1 = azy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/zfilm[1]; b2 =base;
		z02 = (b2 - b1)/(a1 - a2);
		y02 = a1*z02 + b1;
		a1 = axy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/xpim1; b2 = base;
		x02 = (b2 - b1)/(a1 - a2);
		y02 = a1*x02 + b1;
		if(z01 < -(ztop-100.) && z01 > -(zbot+100.)) {x0 = x01; y0 = y01; z0 = z01;}
		if(z02 < -(ztop-100.) && z02 > -(zbot+100.)) {x0 = x02; y0 = y02; z0 = z02;}
	      }
	    else
	      {
		printf("There isn't direct corresponding point %d for %d track\n", i, track);
		nbad++;
	      }
	  }
	  VERT0:
	    if(i==0 && x0+y0+z0 == 0.) 
	      { //the point of vertex must be found anyway!!!!!!
		xpim0 = (double) xtrack[0][track][i]*pixel2mm;
		ypim0 = (double) ytrack[0][track][i]*pixel2mm;
		azy_line[0] = ypim0/zfilm[0];
		axy_line[0] = ypim0/xpim0;
		xpim1 = ((double) xtrack[1][track][i]*pixel2mm+xtransfer);
		ypim1 = ((double) ytrack[1][track][i]*pixel2mm+base);//-ytransfer+base;
		a1 = azy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/zfilm[1]; b2 = base;
		z0 = (b2 - b1)/(a1 - a2);
		y0 = a1*z0 + b1;
		a1 = axy_line[0]; b1 = 0.;
		a2 = (ypim1-base)/xpim1; b2 = base;
		x0 = (b2 - b1)/(a1 - a2);
		y0 = a1*x0 + b1;
	      }
	    if(i==0) 
	      {
		printf("track %d vertex point x0 %6.3f y0 %6.3f z0 %8.3f\n", track, x0, y0, z0);
		printf("track %d vertex point xpim0 %6.3f ypim0 %6.3f xpim1 %6.3f ypim1 %6.3f\n", 
		       track, xpim0, ypim0, xpim1, ypim1);
	      }
	    if(x0+y0+z0 !=0.){
	      pm3d->SetPoint( k, x0, y0, z0 );
	      x3D[track][l] = (float) x0;
	      y3D[track][l] = (float) y0;
	      z3D[track][l] = (float) z0;
	      wes[track][l] = (float) CCDBuffer[0][track][i]/2.;
	      k++; l++;
	      npt3D[track]++;
	      npt3D_left[track]++;
	    }
	}
      nright[track] = l;
      npt3D_right[track] = 0;
      for(i=0; i<np[1][track]; i++)
	{
	  //if point lies between min-max bounds on X direction
	  //and if track goes not very along Y-axis
	  if((xtrack[1][track][i] < xl[0][track] ||
	      xtrack[1][track][i] > xr[0][track]) &&
	     (xr[0][track] - xl[0][track]) > 50. && i!=0 ) continue;
	  double azy_line[2];
	  double axy_line[2];
	  double zp1, zp2, xp1, xp2, yp1, yp2;
	  double xpim1, xpim0, ypim1, ypim0;
	  zp1 = zp2 = xp1 = xp2 = yp1 = yp2 = 0.;
	  double y2;
	  double xtrack_mm[2], ytrack_mm[2];
	  double xt0, xt1, a1, b1, a2, b2, z0;
	  double x01, y01, z01, x02, y02, z02;
	  x01 = y01 = z01 = x02 = y02 = z02 = x0 = y0 = z0 = 0.;
	  xtrack_mm[1] = (double) xtrack[1][track][i]*pixel2mm + xtransfer;// - deltaV/2.;
	  ytrack_mm[1] = ((double)ytrack[1][track][i]*pixel2mm +base);//- ytransfer);
	  xpim1 = xtrack_mm[1];
	  ypim1 = ytrack_mm[1];
	  if(fabs(xpim1) < 0.1 && i != 0) continue;
	  azy_line[0] = (ypim1 - base)/zfilm[1];
	  axy_line[0] = (ypim1 - base)/xpim1;
	  
	  xpim0 = xpim1*zfilm[0]/zfilm[1];
	  double xpim1_pxl = xpim0/pixel2mm;
	  double ypim1buf = (ypim1-base)*zfilm[0]/zfilm[1];
	  int npts = 0;
	  for (j = 0; j < np[0][track]; j++)
	    {
	      if(fabs(xpim1_pxl - xtrack[0][track][j]) < 0.5) 
		{
		  if (ytransfer < ytrack[0][track][j]*pixel2mm - ypim1buf &&
		      ytrack[0][track][j]*pixel2mm - ypim1buf < ytransfer2) 
		    {
		      //printf("ytransfer %f ypim1buf - ytrack[0][track][j] %f ytransfer2 %f\n",
		      //     ytransfer, ytrack[0][track][j]*pixel2mm - ypim1buf, ytransfer2);
		      xpim0 = xtrack[0][track][j]*pixel2mm;
		      ypim0 = ytrack[0][track][j]*pixel2mm;
		      a1 = azy_line[0]; b1 = base;
		      a2 = ypim0/zfilm[0]; b2 = 0.;
		      if(fabs(xpim0) < 0.1) continue;
		      if(fabs(a1-a2) < 0.001) continue;
		      z0 = (b2 - b1)/(a1 - a2);
		      y0 = a1*z0 + b1;
		      a1 = axy_line[0]; b1 = base;
		      a2 = ypim0/xpim0; b2 = 0.;
		      if(fabs(a1-a2) < 0.001) continue;
		      x0 = (b2 - b1)/(a1 - a2);
		      y0 = a1*x0 + b1;
		      pm3d->SetPoint( k, x0, y0, z0 );
		      x3D[track][l] = (float) x0;
		      y3D[track][l] = (float) y0;
		      z3D[track][l] = (float) z0;
		      wes[track][l] = (float) (CCDBuffer[1][track][i]+CCDBuffer[0][track][j])/2.;
		      k++; l++;
		      npt3D[track]++;
		      npt3D_right[track]++;
		      npts++;
		    }
		}
	    }
	  if(npts ==0 && ((i < np[1][track] - 1 && i > 1)|| i==0)){
	    double phibuf;
	    if(fabs(xpim1 - xc1) > 0.001) 
	      phibuf = atan2(ypim1 - yc1, xpim1 - xc1);
	    else
	      phibuf = TMath::Pi()/2.;
	    xt1 = rtr1*cos(phibuf)+xc1;
	    xpim0 = xpim1*zfilm[0]/zfilm[1];
	    xt0 = xpim0;
		 //xt1*zfilm[0]/zfilm[1];
	    if(fabs(xt0 - xc0) <= (rtr0-10.*pixel2mm) && 
	       fabs(xt0 - xc0) > 0.001 )
	      {//there is an intersection
		y1 = yc0 + sqrt(rtr0*rtr0 - (xt0 - xc0)*(xt0 - xc0)); 
		y2 = yc0 - sqrt(rtr0*rtr0 - (xt0 - xc0)*(xt0 - xc0)); 
		//printf(" Image 2 xt %f, y1 %f, y2 %f\n", xt, y1, y2);
		ypim0 = y1;
		//(y1 - yc0)/(xt0 - xc0)*xpim0 + 
		//(y1 - (y1 - yc0)/(xt0 - xc0)*xt0);
		a1 = azy_line[0]; b1 = base;
		a2 = ypim0/zfilm[0]; b2 = 0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		z01 = (b2 - b1)/(a1 - a2);
		y01 = a1*z01 + b1;
		a1 = axy_line[0]; b1 = base;
		a2 = ypim0/xpim0; b2 = 0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		x01 = (b2 - b1)/(a1 - a2);
		y01 = a1*x01 + b1;
		ypim0 = y2;
		//(y2 - yc0)/(xt0 - xc0)*xpim0 + 
		//(y2 - (y2 - yc0)/(xt0 - xc0)*xt0);
		a1 = azy_line[0]; b1 = base;
		a2 = ypim0/zfilm[0]; b2 =0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		z02 = (b2 - b1)/(a1 - a2);
		y02 = a1*z02 + b1;
		a1 = axy_line[0]; b1 = base;
		a2 = ypim0/xpim0; b2 = 0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		x02 = (b2 - b1)/(a1 - a2);
		y02 = a1*x02 + b1;
		if(z01 < -(ztop-100.) && z01 > -(zbot+100.)) {x0 = x01; y0 = y01; z0 = z01;}
		if(z02 < -(ztop-100.) && z02 > -(zbot+100.)) {x0 = x02; y0 = y02; z0 = z02;}
		if(x0 == 0. && y0 == 0. && z0 == 0.)
		  {
		    float distz_1, distz_2, distz_1top, distz_1bot, distz_2top, distz_2bot;
		    distz_1top = fabs(z01 + ztop); distz_1bot = fabs(z01 + zbot);
		    distz_2top = fabs(z02 + ztop); distz_2bot = fabs(z02 + zbot);
		    if(distz_1top < distz_1bot) distz_1 = distz_1top;
		    if(distz_1top > distz_1bot) distz_1 = distz_1bot;
		    if(distz_2top < distz_2bot) distz_2 = distz_2top;
		    if(distz_2top > distz_2bot) distz_2 = distz_2bot;
		    if(distz_1 < distz_2) {x0 = x01; y0 = y01; z0 = z01;}
		    if(distz_1 > distz_2) {x0 = x02; y0 = y02; z0 = z02;}
		  }
		//printf("x01 %6.3f y01 %6.3f z01 %8.3f\n", x01, y01, z01);
		//printf("x02 %6.3f y02 %6.3f z02 %8.3f\n", x02, y02, z02);
		//printf("x0 %6.3f y0 %6.3f z0 %8.3f\n", x0, y0, z0);
	      }
	    else if (fabs(xt0 - xc0) <= 0.001 )
	      {//there is an intersection
		y1 = yc0 + rtr0; 
		y2 = yc0 - rtr0;
		//printf(" Image 2 xt %f, y1 %f, y2 %f\n", xt, y1, y2);
		ypim0 = y1;
		a1 = azy_line[0]; b1 = base;
		a2 = ypim0/zfilm[0]; b2 = 0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		z01 = (b2 - b1)/(a1 - a2);
		y01 = a1*z01 + b1;
		a1 = axy_line[0]; b1 = base;
		a2 = ypim0/xpim0; b2 = 0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		x01 = (b2 - b1)/(a1 - a2);
		y01 = a1*x01 + b1;
		ypim0 = y2;
		a1 = azy_line[0]; b1 = base;
		a2 = ypim0/zfilm[0]; b2 =0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		z02 = (b2 - b1)/(a1 - a2);
		y02 = a1*z02 + b1;
		a1 = axy_line[0]; b1 = base;
		a2 = ypim0/xpim0; b2 = 0.;
		if(fabs(a1-a2) < 0.001) 
		  {if(i>0)continue;
		  goto VERT1;}
		x02 = (b2 - b1)/(a1 - a2);
		y02 = a1*x02 + b1;
		if(z01 < -(ztop-100.) && z01 > -(zbot+100.)) {x0 = x01; y0 = y01; z0 = z01;}
		if(z02 < -(ztop-100.) && z02 > -(zbot+100.)) {x0 = x02; y0 = y02; z0 = z02;}
	      }
	    else
	      {
		printf("There isn't direct corresponding point %d for %d track, image 1\n", i, track);
		nbad++;
	      }
	  }
	VERT1:
	  if(i==0 && x0+y0+z0 == 0.)
	    { //the point of vertex must be found anyway!!!!!!
	      xpim0 = (double) xtrack[0][track][i]*pixel2mm;
	      ypim0 = (double) ytrack[0][track][i]*pixel2mm;
	      azy_line[0] = ypim0/zfilm[0];
	      axy_line[0] = ypim0/xpim0;
	      xpim1 = ((double) xtrack[1][track][i]*pixel2mm)+xtransfer;
	      ypim1 = ((double) ytrack[1][track][i]*pixel2mm+base);//-ytransfer+base;
	      a1 = azy_line[0]; b1 = 0.;
	      a2 = (ypim1-base)/zfilm[1]; b2 = base;
	      z0 = (b2 - b1)/(a1 - a2);
	      y0 = a1*z0 + b1;
	      a1 = axy_line[0]; b1 = 0.;
	      a2 = (ypim1-base)/xpim1; b2 = base;
	      x0 = (b2 - b1)/(a1 - a2);
	      y0 = a1*x0 + b1;
	      printf("vertex point x0 %6.3f y0 %6.3f z0 %8.3f\n", x0, y0, z0);
	    }
	  if(i==0) 
	    {
	      printf("track %d vertex point x0 %6.3f y0 %6.3f z0 %8.3f\n", track, x0, y0, z0);
	      printf("track %d vertex point xpim0 %6.3f ypim0 %6.3f xpim1 %6.3f ypim1 %6.3f\n", 
		     track, xpim0, ypim0, xpim1, ypim1);
	    }
	  if(x0+y0+z0 !=0.){
	    pm3d->SetPoint( k, x0, y0, z0 );
	    x3D[track][l] = (float) x0;
	    y3D[track][l] = (float) y0;
	    z3D[track][l] = (float) z0;
	    wes[track][l] = (float) CCDBuffer[1][track][i]/2.;
	    k++; l++;
	    npt3D[track]++;
	    npt3D_right[track]++;
	  }
	}
      //npt3D[track] = l;
      printf("npt3D[%d] = %d, npt_left = %d, npt_right = %d\n", 
	     track, npt3D[track], npt3D_left[track],npt3D_right[track]);
      if ((float) nbad/(float) (np[0][track] + np[1][track]) > 0.25) {
	if(Radius[2][track] != 0.){
	  //no_fit[track] = 1;
	  printf("Too many bad points -> Take fit from Gil's calculation!!!!\n");}
      }
    }
  printf("k %d\n",k);
  // set marker size, color & style
  pm3d->SetMarkerSize( 0.5 );
  pm3d->SetMarkerColor( 52);
  pm3d->SetMarkerStyle( 8 );
  pm3d->Draw();
  
  pl11 = new TPaveLabel(-0.99,0.87,-0.65,0.99,ifilename,"br");
  pl11->SetFillColor(42);
  pl11->SetTextSize(0.5);
  pl11->Draw();
  char PStr[60];
  TPaveText *text = new TPaveText(-0.1, 0.51, 0.99, 0.97 );
  text->SetFillColor( 42 );
  for(track=0; track<tracks; track++)
    {
      sprintf( PStr, "p%d = %8.1f, T%d = %8.1f", 
  	       track, p[track], track, T[track]); 
     text->AddText(PStr);
    }
  text->Draw();
  prc3->Update();
  
  c4 = new TCanvas("c4"," Plot in 3D",40,40,ixWindow,iyWindow);
  plc4 = new TPad("plc4","XY view", xleft, yleft, xright, yright, 21);
  plc4->Range(-xchamb,-100.,xchamb,-100.+ychamb);
  //prc4 = new TPad("prc4","Phi, Z view", xleft + 0.5, yleft, xright + 0.5, yright, 21);
  prc4 = new TPad("prc4","Phi*R, Z view", xleft + 0.5, yleft, xright + 0.5, yright, 21);
  //prc4->Range(-pi/6.,-(zbot+100.), 2.1*pi,-(ztop -100.));
  prc4->Range(-xchamb,-(zbot+50.), xchamb,-(ztop -50.));
  c4->ToggleEventStatus();
  
  plc4->Draw();
  plc4->UseCurrentStyle();
  
  prc4->Draw();
  prc4->UseCurrentStyle();
  
  plc4->cd();
  plc4->Draw("col");
  
  prc4->cd();
  prc4->Draw("col");
  
  c4->Update();
  plc4->cd();
  for(track=0; track<tracks; track++)
    {
      for(i=0; i<npt3D[track]; i++)
	{
	  if(wes[track][i] > 0.) 
	    {
	      mark = new TMarker((Float_t) x3D[track][i], (Float_t) y3D[track][i], 5);
	      int col = 0;
	      switch (track) {
	      case 0: 
		col = 1; break;
	      case 1:
		col = 51; break;
	      case 2:
		col = 50; break;
	      case 3:
		col = 38; break;
	      case 4:
		col = 42; break;
	      case 5:
		col = 37; break;
	      default:
		col = 1; break;
	      }
	      mark->SetMarkerColor(col);
	      mark->SetMarkerStyle(8);
	      mark->SetMarkerSize(0.2);
	      fObjTbl->Add(mark);
	      mark->DrawMarker((Float_t) x3D[track][i], (Float_t) y3D[track][i]);
	    }
	}
    }
  for(i=0; i<ncrosses; i++)
    {
      mark = new TMarker((Float_t) xcross3D[i], (Float_t) ycross3D[i], 5);
      mark->SetMarkerColor(1);
      mark->SetMarkerStyle(2);
      mark->SetMarkerSize(1.5);
      fObjTbl->Add(mark);
      mark->DrawMarker((Float_t) xcross3D[i], (Float_t) ycross3D[i]);
    }
  plc4->Modified();
  plc4->Update();
  
  prc4->cd();


  //second time Kinematics but in 3D
  //  for(mypass = 0; mypass<2; mypass++)
  int k1, k2, kt;
  deltapx = deltapy = deltapz = 0.;
  deltaE = 0.;
  deltaT = 0.;
  float XC_help, YC_help;
  for(track = 0; track<tracks; track++)
    {
      //for(view=0;view<2;view++)
		view = 0;
	{
	  xcD[view][track] = 0.;
	  ycD[view][track] = ycD[1][track] = 0.;
	  rcD[view][track] = rcD[1][track] = 0.;
	  a[view][track] = b[view][track] = a[1][track] = b[1][track] = 0.;
	  //if(view == 0) {k1 = 0; k2 = npt3D_left[track];}
	  //if(view == 1) {k1 = npt3D_left[track]; k2 = npt3D[track];}
	  k1 = 0;
	  k2 = npt3D[track];
	  if(no_fit[track]==1 || k2-k1 < 10) {
	    //take radius from Gil's calculations
	    isig[view][track] = 10;
	    no_fit[track] = 1;
	    continue;
	  }
	  Find_Tails(track, 0);
	  Xcircle[0] = x3D[track][fp];
	  Xcircle[1] = x3D[track][mp];
	  Xcircle[2] = x3D[track][lp];
	  Ycircle[0] = y3D[track][fp];
	  Ycircle[1] = y3D[track][mp];
	  Ycircle[2] = y3D[track][lp];
	  Wpoint[0] = Wpoint[1] = Wpoint[2] = 1.;

	  //printf("Xcircle[0] %f Xcircle[1] %f Xcircle[2] %f\n",
	  // Xcircle[0], Xcircle[1], Xcircle[2]);
	  //printf("Ycircle[0] %f Ycircle[1] %f Ycircle[2] %f\n",
	  // Ycircle[0], Ycircle[1], Ycircle[2]);
	  XC = YC = RC = 0.;
	  Npoints_fit = 3;
	  circle_minuit(Npoints_fit);
	  k = 0;
	  float WpMax = 0.;
	  for(i=k1; i<k2; i++)  
	    {
	      Xcircle[k] = x3D[track][i];
	      Ycircle[k] = y3D[track][i];
	      Zcircle[k] = z3D[track][i];
	      //Wpoint[k] = 1.;
			Wpoint[k] = wes[track][i];
	      if (Wpoint[i] > WpMax) WpMax = Wpoint[i];
	      k++;
	      //if(i==0) Wpoint[i] = 0.;//don't take the vertex point!
	    }
	  //define signe of track by 3 points-fit in order to
	  //change first approximation if signe of the track is wrong
	  Npoints_fit = k;
	  bufmin = 1000.;
	  bufmax =-1000.;
	  for(i=0;i<Npoints_fit; i++)
	    {
	      if(i==0) Wpoint[i] = WpMax;//take the vertex point!
			Wpoint[i] = Wpoint[i]/WpMax;
	      phi[i] = atan2((Ycircle[i]-YC),(Xcircle[i]-XC));
	      if(phi[i]<0.) phi[i] += TMath::Pi()*2.;
	      if(*(phi+i) < bufmin) bufmin = *(phi+i);
	      if(*(phi+i) > bufmax) bufmax = *(phi+i);
	    }
	  //check for the -Pi - +Pi bounds
	  if((bufmax - bufmin) > (float)TMath::Pi()) 
	    {
	      bufmin = 1000.;
	      bufmax =-1000.;
	      for(i=0;i<Npoints_fit; i++)
		{
		  if (*(phi+i) < (float)TMath::Pi()) *(phi+i) += 2.*(float)TMath::Pi();
		  if(*(phi+i) < bufmin) bufmin = *(phi+i);
		  if(*(phi+i) > bufmax) bufmax = *(phi+i);
		}
	    }
	  isig_tr[track] = 1;
	  isig[view][track] = 1;
	  if(bufmax - phi[0] > phi[0] - bufmin)
	    //if(phi[0] - phi[Npoints_fit/2] < 0.) 
	    isig_tr[track] = isig[view][track] = -1;
	  phi_first[view][track] = bufmax;
	  if(isig[view][track] == -1)
	    phi_first[view][track] = bufmin;
	  Signe=10; //! for Pi^- we don't know the signe of result tracks=>not define it
	  //Signe = isig[view][track];
	  if(track==0) Signe = 1;
	  //printf("bufmin %f phi[0] %f bufmax %f\n", bufmin, phi[0], bufmax);
	  //printf("isig %d signe %d\n", isig[view][track], Signe);
	  //printf(" before fit XC %f YC %f \n", XC, YC);

	  if(isig[view][track]!=Signe&&track==0) 
	    {
	      //look for the line along which then make a mirrow flip
	      //of the circle coordinate XC & YC
	      phid = &Xcircle[0];
	      z3d = &Ycircle[0];
	      A_line = &a[view][track]; 
	      B_line = &b[view][track];
	      DA_line = &da[view][track];
	      chi2_2 = &chi2_line[view][track];
	      line(Npoints_fit, phid, z3d, A_line, B_line, DA_line, chi2_2);
	      if(*A_line !=0.) 
		{
		  double  xpro, ypro;
		  double SinPhi = sin(atan(*A_line));
		  double CosPhi = cos(atan(*A_line));
		  //coordinates of the circle center in system of the center rotated
		  // along line of the points
		  xpro =  XC*CosPhi + (YC-*B_line)*SinPhi;
		  ypro = -XC*SinPhi + (YC-*B_line)*CosPhi;
		  //in mirrow flip y=-y then
		  ypro = -ypro;
		  XC =   xpro*CosPhi - ypro*SinPhi;
		  YC =   xpro*SinPhi + ypro*CosPhi +*B_line;
		}
	      else
		{
		  YC = -(YC-*B_line)+*B_line;
		}
	      //printf(" before fit a %f b %f XC %f YC %f \n", *A_line, *B_line, XC, YC);
	    }
	  
	  XC_help = XC; YC_help = YC;
	  circle_minuit(Npoints_fit);
	  printf("D3 circle fit : track %d rad %f cent %f %f DRC %f chi2 %f \n",
	   track, RC, XC, YC, DRC, chi2);	 

	  xcD[view][track] = XC;
	  ycD[view][track] = YC;
	  rcD[view][track] = RC;
	  drcD[view][track] = DRC;
	  chi2_circle[view][track] = chi2;
	  if(chi2 <= 0.)
	    chi2_circle[view][track] = 1000000.;
	  bufmin = 1000.;
	  bufmax =-1000.;
	  for(i=0;i<Npoints_fit; i++)
	    {
	      phi[i] = atan2((Ycircle[i]-YC),(Xcircle[i]-XC));
	      if(phi[i]<0.) phi[i] += TMath::Pi()*2.;
	      if(*(phi+i) < bufmin) bufmin = *(phi+i);
	      if(*(phi+i) > bufmax) bufmax = *(phi+i);
	    }
	  //check for the -Pi - +Pi bounds
	  if((bufmax - bufmin) > (float)TMath::Pi()) 
	    {
	      bufmin = 1000.;
	      bufmax =-1000.;
	      for(i=0;i<Npoints_fit; i++)
		{
		  if (*(phi+i) < (float)TMath::Pi()) *(phi+i) += 2.*(float)TMath::Pi();
		  if(*(phi+i) < bufmin) bufmin = *(phi+i);
		  if(*(phi+i) > bufmax) bufmax = *(phi+i);
		}
	    }
	  isig_tr[track] = 1;
	  isig[view][track] = 1;
	  if(bufmax - phi[0] > phi[0] - bufmin)
	    //if(phi[0] - phi[Npoints_fit/2] < 0.) 
	    isig_tr[track] = isig[view][track] = -1;
	  phi_first[view][track] = bufmax;
	  if(isig[view][track] == -1)
	    phi_first[view][track] = bufmin;
	  lxy[track] = fabs((bufmax-bufmin)*RC);
	  //suppose all tracks go from the vertex
	  for(i=0; i<Npoints_fit; i++) 
	    {
	      x_mark[i] = (phi[i]-phi_first[view][track]);
	      y_mark[i] = Zcircle[i];
	    }
	  phid = &x_mark[0];
	  z3d = &y_mark[0];
	  A_line = &a[view][track]; 
	  B_line = &b[view][track];
	  DA_line = &da[view][track];
	  chi2_2 = &chi2_line[view][track];
	  line(Npoints_fit, phid, z3d, A_line, B_line, DA_line, chi2_2);
	  chi2_line[view][track] = *chi2_2;
	  if(*chi2_2 <= 0.)
	    chi2_line[view][track] = 1000000.;
	  da[view][track] = *DA_line;
	  printf("D3 line fit : tr %d A_line %f B_line %f, chi2 %6.3f DA %f\n",
	      	 track, *A_line, *B_line, *chi2_2, *DA_line);
	  AC = *A_line;//a[view][track];
	  BC = *B_line;//b[view][track];
	  DAC = *DA_line;
	  a[view][track] = AC;//*A_line;//AC;
	  b[view][track] = BC;//*B_line;//BC;
	  da[view][track] = DAC;//*DA_line;//DAC;
	  Signe=10;
	  //Signe = isig[view][track];
	  if(track==0)Signe=1;
	  helix_minuit(Npoints_fit);
	  int renew = 0;
	  
	RENEWFIT0:
	  printf("D3 helix fit: tr %d R %f XC %f YC %f DRC %f \n",
		 track, RC, XC, YC, DRC);	 
	  printf("D3 helix fit: tr %d A %f B %f DA %f SIGN %1d PHIF %f\n",
	    	 track, AC, BC, DAC, is, phif);
 	  printf("D3 helix fit: ERROR_FLAG is    %d chi2 = %f \n", error_flag, chi2);
	  if((error_flag == 0 && chi2 > 0. && chi2 < 100. && fabs(RC-rcD[view][track])/rcD[view][track] < 0.1) || renew==1)
		 // in order to don't spoil the circle fit if helix one is too bad
	    {
 			printf("Take helix fit as final result of fit\n");
	      xcD[view][track] = XC;
	      ycD[view][track] = YC;
	      rcD[view][track] = RC;
	      a[view][track] = AC;//*A_line;//AC;
	      b[view][track] = BC;//*B_line;//BC;
	      drcD[view][track] = DRC;
	      da[view][track] = DAC;//*DA_line;//DAC;
	      chi2_circle[view][track] = chi2;
			//start to order the points along the circle 
			//beginning from the Vertex point
			// suppose that start and end angles are between 180degrees
			for(i=0;i<Npoints_fit; i++)
			  phi[i] = t[i];
			lxy[track] = fabs((bufmax-bufmin)*RC);
			if(view==1)
			  lxy[track] = (fabs(bufmax-bufmin)*RC+lxy[track])/2.;
			isig_tr[track] = is;
			isig[view][track] = is;
			phi_first[view][track] = phif;
			z_first[view][track] = b[view][track];
	    }
	  else
	    {
 			printf("Take circle and line fit as final results of fit\n");
	      chi2_circle[view][track] += chi2_line[view][track];
			printf("D3 circle fit : track %d rad %f cent %f %f chi2 %f \n",
					 track, rcD[view][track], xcD[view][track], ycD[view][track], 
					 chi2_circle[view][track]);	 
			printf("D3 line fit : tr %d A_line %f B_line %f\n",
					 track, a[view][track], b[view][track]);
			// all the parameters are taken from two previous fits (circle & line)
	    }
	  L[track] = sqrt(pow(fabs(bufmax-bufmin)*a[view][track],2) + pow(lxy[track],2));
	  if(track == 0 && renew == 0) 
	    {
	      double p0 = 0.03*6.5*rcD[view][track]*sqrt(1.+pow(a[view][track]/rcD[view][track],2))
			  *charge[particle[track]];
	      if (p0 < P_MEAN - 50. || p0 > P_MEAN + 50.)
			  {
				 XC = XC_help;//XC/RC*R_MEAN;
				 YC = YC_help;//YC/RC*R_MEAN;
				 RC = R_MEAN;
				 if(track==0)Signe=1;
				 helix_minuit_fixR(Npoints_fit);
				 //circle_minuit_fixR(Npoints_fit);
				 renew = 1;
				 goto RENEWFIT0;
			  }
	    }
	  for(i=0; i<Npoints_fit; i++) 
	    {
	      x_mark[i] = (phi[i]-phi_first[view][track])*RC;
	      //sqrt(pow(Ycircle[i]-YC,2)+pow(Xcircle[i]-XC,2));
	      y_mark[i] = Zcircle[i];
	      //printf("x_mark[%d] %f y_mark %f\n", i, x_mark[i], y_mark[i]);
	    }
	  px_mark = &x_mark[0];
	  py_mark = &y_mark[0];
	  markp = new TPolyMarker((Int_t) Npoints_fit, px_mark, py_mark, "same");
	  fObjTbl->Add(markp);
	  markp->SetMarkerStyle(8);
	  markp->SetMarkerSize(0.2);
	  int col = 0;
	  switch (track) {
	  case 0: 
	    col = 1; break;
	  case 1:
	    col = 51; break;
	  case 2:
	    col = 50; break;
	  case 3:
	    col = 38; break;
	  case 4:
	    col = 42; break;
	  case 5:
	    col = 37; break;
	  default:
	    col = 1; break;
	  }
	  markp->SetMarkerColor(col);
	  markp->Draw("same");	    
	  for(i=0; i<Npoints_fit; i++) 
	    {
	      x_mark[i] = (phi[i]-phi_first[view][track])*RC;
	      y_mark[i] = a[view][track]*(phi[i]-phi_first[view][track])+
		b[view][track];
	    }
	  px_mark = &x_mark[0];
	  py_mark = &y_mark[0];
	  markp = new TPolyMarker((Int_t) Npoints_fit, px_mark, py_mark, "same");
	  fObjTbl->Add(markp);
	  markp->SetMarkerStyle(8);
	  markp->SetMarkerSize(0.5);
	  markp->SetMarkerColor(1);
	  if(track==0) markp->SetMarkerColor(2);
	  markp->Draw("same");	    
	}
    }
  for(i=0; i<ncrosses; i++)
    {
      mark = new TMarker((Float_t) 0., (Float_t) zcross3D[i], 5);
      mark->SetMarkerColor(1);
      mark->SetMarkerStyle(2);
      mark->SetMarkerSize(1.5);
      fObjTbl->Add(mark);
      mark->DrawMarker((Float_t) 0., (Float_t) zcross3D[i]);
    }
  prc4->Modified(); 
  prc4->Update();
  //}
  //if(mypass == 0 && tracks > 1)
  //{
  //for(track = 0; track<tracks; track++)
  //printf(" Old Vertex coordinates: track %d x %f, y %f\n",
  //   track, x3D[track][0], y3D[track][0]);
  Xvertax = Yvertax = Dvertax = 0.;
  if(tracks != 1)
    {
      float xc[20], yc[20], rc[20];
      int itr = 0;
      for(track = 0; track<tracks; track++)
		  {
			 if(rcD[0][track] > 10. && 
				 //chi2_circle[0][track]<500. && 
				 npt3D_left[track] > 2)
				{
				  xc[itr] = xcD[0][track];
				  yc[itr] = ycD[0][track];
				  rc[itr] = rcD[0][track];
				  itr++;
				  if(itr==1)
					 {
						Xvertax = x3D[track][0];
						Yvertax = y3D[track][0];
					 }
				}
		  }
      if (itr > 1)
		  {
			 XY_Vertex_Net(itr, &xc[0], &yc[0], 
								&rc[0], &Xvertax, &Yvertax, 
								&Dvertax);
			 plc4->cd();
			 
			 mark = new TMarker(Xvertax,Yvertax,5);
			 mark->SetMarkerColor(38);
			 mark->SetMarkerStyle(8);
			 mark->SetMarkerSize(1.5);
			 mark->DrawMarker(Xvertax, Yvertax);
			 plc4->Modified();
			 plc4->Update();
			 printf(" Vertex coordinates: x %f, y %f; error %f\n",
					  Xvertax, Yvertax, Dvertax );
		  }
    }
  prc4->cd();
  //store points along track(1cm) close to the vertex
  float sumzv = 0.;
  ii = 0;
  for(track=0;track<tracks;track++)
	 {
		Zv[track] = 0.;
		kt = 0;
		view = 0;
		//for(view = 0;view<2; view++)
		//{
		//calculate the direction cosinuses
		float tmp_phi, dx, dy, dz, d;
		float phi_vert;
		float xc_short, yc_short, rc_short, ac_short, bc_short;
		int isig_short;
		float phi_firstshort;
		phi_vert = atan2((Yvertax - ycD[view][track]),
							  (Xvertax - xcD[view][track]));
		if(phi_vert < 0. ) phi_vert += TMath::Pi()*2.;
		if(fabs(phi_vert - phi_first[view][track]) > TMath::Pi())
		  {
			 if(isig[view][track] < 0) phi_vert -= TMath::Pi()*2.;
			 if(isig[view][track] > 0) phi_vert += TMath::Pi()*2.;
		  }
		//if(view == 0) {k1 = 0; k2 = npt3D_left[track];}
		//if(view == 1) {k1 = npt3D_left[track]; k2 = npt3D[track];}
		k1 = 0;
		k2 = npt3D[track];
		
		xc_short = xcD[view][track];
		yc_short = ycD[view][track];
		rc_short = rcD[view][track];
		ac_short = a[view][track];
		bc_short = b[view][track];
		isig_short = isig[view][track];
		phi_firstshort = phi_first[view][track];
		//printf("k %d\n",k);
		phi_vert = atan2((Yvertax - yc_short),
							  (Xvertax - xc_short));
		tmp_phi = fabs(fabs(phi_vert) - (float)TMath::Pi()/2.);
		if(phi_vert < 0. ) phi_vert += TMath::Pi()*2.;
		if(Yvertax - yc_short >= 0. && 
			Xvertax - xc_short >= 0.) 
		  { //I quadrant
			 dx =-10.*cos(tmp_phi);
			 dy = 10.*sin(tmp_phi);
		  }
		if(Yvertax - yc_short >= 0. && 
			Xvertax - xc_short < 0.)
		  { //II quadrant
			 dx =-10.*cos(tmp_phi);
			 dy =-10.*sin(tmp_phi);
		  }
		if(Yvertax - yc_short < 0. && 
			Xvertax - xc_short < 0.)
		  { //III quadrant
			 dx = 10.*cos(tmp_phi);
			 dy =-10.*sin(tmp_phi);
		  }
		if(Yvertax - yc_short < 0. && 
			Xvertax - xc_short >= 0.)
		  { //IV quadrant
			 dx = 10.*cos(tmp_phi);
			 dy = 10.*sin(tmp_phi);
		  }
		//tmp_phi = atan(10./(*RC));
		tmp_phi = atan(10./
							sqrt((Xvertax - xc_short)*
								  (Xvertax - xc_short) +
								  (Yvertax - yc_short)*
								  (Yvertax - yc_short)));
		dz = tmp_phi*ac_short;
		d  = sqrt(dx*dx + dy*dy + dz*dz);
		cosdir_x[view][track] = dx/d; 
		cosdir_y[view][track] = dy/d; 
		cosdir_z[view][track] = dz/d;
		if((isig_short == 1 && track != 0) || 
			(isig_short ==-1 && track == 0)) 
		  {
			 cosdir_x[view][track] = -dx/d; 
			 cosdir_y[view][track] = -dy/d; 
			 cosdir_z[view][track] = -dz/d;
		  }
		printf("%2d %2d %1d : cosx %f cosy %f cosz %f\n", 
				 track, view, isig_short,
				 cosdir_x[view][track], cosdir_y[view][track], cosdir_z[view][track]);
		//printf("%2d %2d: dipx %f dipy %f dipz %f\n", 
		//   track, view, acos(cosdir_x[view][track])*180./pi, 
		//   acos(cosdir_y[view][track])*180./pi, 
		//   90.-acos(cosdir_z[view][track])*180./pi);
		if(fabs(phi_vert - phi_firstshort) > TMath::Pi())
		  {
			 if(isig_short < 0) phi_vert -= TMath::Pi()*2.;
			 if(isig_short > 0) phi_vert += TMath::Pi()*2.;
		  }
		printf("z vertex %f\n",
				 (phi_vert - phi_firstshort)*ac_short+bc_short);
		x_mark[view] = (phi_vert - phi_firstshort)*rc_short;
		y_mark[view] = (phi_vert - phi_firstshort)*ac_short+bc_short;
		//printf("phi_vert %f phi_first %f x_mark %f\n", phi_vert, phi_firstshort, x_mark[view]);
		if(fabs(y_mark[view] - (-zmed)) < 200.) {
		  Zv[track] += y_mark[view];
		  kt++;}
		// }
		if(kt>0) Zv[track] /=(float) kt;
		if(kt<=0) Zv[track] =0.;
		printf("Z coordinate of vertex is %f\n", Zv[track]);
		if(kt>0) {sumzv += Zv[track]; ii++;}
		px_mark = &x_mark[0];
		py_mark = &y_mark[0];
		markp = new TPolyMarker((Int_t) 1, px_mark, py_mark, "same");
		fObjTbl->Add(markp);
		markp->SetMarkerStyle(8);
		markp->SetMarkerSize(1.);
		int col = 0;
		switch (track) {
		case 0: 
		  col = 1; break;
		case 1:
		  col = 51; break;
		case 2:
		  col = 50; break;
		case 3:
		  col = 38; break;
		case 4:
		  col = 42; break;
		case 5:
		  col = 37; break;
		default:
		  col = 1; break;
		}
		markp->SetMarkerColor(col);
		markp->Draw("same");	    
	 }
  
  //   vertex improving
  Zvertax = 0.;
  if(ii > 0) Zvertax = sumzv/((float) ii);
  ii = 0;
  sumzv = 0.;
  for(i=0;i<tracks;i++)
	 {
		if(Zv[i] != 0.) {sumzv += pow(Zv[i]-Zvertax,2); ii++;}
	 }
  if(ii > 0) disp = sqrt(sumzv/((float) ii));
  mark = new TMarker(0.,Zvertax,5);
  mark->SetMarkerColor(38);
  mark->SetMarkerStyle(8);
  mark->SetMarkerSize(1.5);
  mark->DrawMarker(0.,Zvertax);
  printf("Zvertax is %f disp %f\n", Zvertax, disp);
  prc4->Modified();
  prc4->Update();
  
  deltapx = deltapy = deltapz = 0.;
  float deltapx0, deltapx1, deltapy0, deltapy1, deltapz0, deltapz1;
  float dp0, dp1;
  float max_loss = 0.;

  dp0 = dp1 = 1000.;
  deltapx0 = deltapy0 = deltapz0 = 0.;
  deltapx1 = deltapy1 = deltapz1 = 0.;
  deltaT = 0.;
  if(tracks == 6) {
    particle[0] = particle[4] = 0; 
    particle[1] = particle[2] = particle[3] = particle[5] = 1;
  }
  
  for(track=0; track<tracks; track++)
	 {
		//suppose that all particles are the pions!!!!!!!!!!!
		particle[track] = 0;
		//Kinematics part
		double cosdx, cosdy, cosdz;
		p2[track] = T2[track] = E2[track] = 0.;
		dp[track] =-1.;
		if(no_fit[track] == 1) 
		  {
			 printf(" Take momentum from Gil calculations!\n");
			 p2[track] = 0.195*Radius[2][track]
				/ cos(dip[track]*pi/180.);  /* in MeV/c */
			 cosdx = cosx[2][track];
			 cosdy = cosy[2][track];
			 cosdz =-cosz[track];
			 printf("cosx %f cosy %f cosz %f\n", cosdx, cosdy, cosdz);
			 if(p2[track] != 0.)
				{
				  //isig_tr[track] = -1;
				  if(track == 0) 
					 isig_tr[track] = 1;
				  chi2_circle[0][track] = 1000.;
				}
		  }
		else 
		  {
			 p2[track] = 0.03*6.5*rcD[0][track]*sqrt(1.+pow(a[0][track]/rcD[0][track],2))
				*charge[particle[track]];
			 chi_tr[track] = chi2_circle[0][track];
			 //error of momentum: 
			 // dP = sqrt((dP/dR*delta R)^2 + (dP/da*delta a)^2)
			 // dP/dR = k/sqrt(1+(a/R)^2)
			 // dP/da = k*a/R/sqrt(1+(a/R)^2)
			 dp[track] = sqrt(pow(0.03*6.5/sqrt(1.+pow(a[0][track]/rcD[0][track],2))
										 *charge[particle[track]]*drcD[0][track],2)+
									pow(0.03*6.5*a[0][track]/sqrt(1.+pow(a[0][track]/rcD[0][track],2))
										 *charge[particle[track]]/rcD[0][track]*da[0][track],2));
			 float chi2_0 = chi2_circle[0][track];
			 printf("track %d momentum %f chi2 %f \n", track, p2[track], chi2_0);
			 if( Radius[2][track] > 0. )
				{
				  float bufp;
				  bufp = 0.195*Radius[2][track]
					 / cos(dip[track]*pi/180.);
				  if(p2[track] > 3000.)
					 {
						p2[track] = TMath::Min(p2[track], bufp);
						if(p2[track] == bufp) chi2_0 = 100.;
						printf("track %d momentum %f chi2 %f \n", track, p2[track], chi2_0);
					 }
				  else if(p2[track] < 40. )
					 {
						p2[track] = TMath::Max(p2[track], bufp);
						if(p2[track] == bufp) chi2_0 = 100.;
						printf("track %d momentum %f chi2 %f \n", track, p2[track], chi2_0);
					 }
				}
			 // now the information about chi2 written in chi2_circle!!!!!!!! after helix fit!!!
			 isig_tr[track] = isig[0][track];
			 chi_tr[track] = chi2_0;
			 cosdx = cosdir_x[0][track];
			 cosdy = cosdir_y[0][track];
			 cosdz = cosdir_z[0][track];
			 deltapx0 = deltapx0 + p2[track]*cosdir_x[0][track];
			 deltapy0 = deltapy0 + p2[track]*cosdir_y[0][track];
			 deltapz0 = deltapz0 + p2[track]*cosdir_z[0][track];
			 if(track == 0) {
				deltapx0 = deltapx0 - 2.*p2[track]*cosdir_x[0][track];
				deltapy0 = deltapy0 - 2.*p2[track]*cosdir_y[0][track];
				deltapz0 = deltapz0 - 2.*p2[track]*cosdir_z[0][track];
			 }
			 printf("p2[%d] %f \n", track, p2[track]);
		  }
		// find the kind of the particle
		cosxD[track] = cosdx;
		cosyD[track] = cosdy;
		coszD[track] = cosdz;
		deltapx = deltapx + p2[track]*cosdx;
		deltapy = deltapy + p2[track]*cosdy;
		deltapz = deltapz + p2[track]*cosdz;
		if(track == 0) {
		  deltapx = deltapx - 2.*p2[track]*cosdx;
		  deltapy = deltapy - 2.*p2[track]*cosdy;
		  deltapz = deltapz - 2.*p2[track]*cosdz;
		}
		// skip this part, take from the picture the kind of the particle
		double avbright;
		//for (j=0; j<2; j++)
		{
		  //avbright = TMath::Max(brrr[0][track],brrr[1][track]);
		  avbright = (brrr[0][track]+brrr[1][track])/2.;
		  avbright = TMath::Max(avbright,1.);
		  //if (j == 0) avbright -= aver_l;
		  //if (j == 1) avbright -= aver_r;
		  Loss[track] = avbright*p2[track]*p2[track];
		}
		
		if(track==0) losstag = Loss[track];
		Loss[track] = Loss[track]/losstag;
		//printf("Loss, particle [%d]: %f %s\n",track,Loss[track],Part[particle[track]]); 
		//
		if(track != 0)
		  {
			 if(Loss[track]*sqrt(1. - coszD[track]*coszD[track])
				 > max_loss)
				max_loss = Loss[track]*
				  sqrt(1. - coszD[track]*coszD[track]);
		  }
		//if(Loss[track] < 10. ) {particle[track] = 0; goto LabPart;}//pion
		//if(Loss[track] < 80. ) {particle[track] = 1; goto LabPart;}//proton
		//if(Loss[track] < 250. ) {particle[track] = 2; goto LabPart;}//deuteron
		//if(Loss[track] < 900. ) {particle[track] = 3; goto LabPart;}//tritium
		//particle[track] = 4; //He3
		
		m2[track] = mass_part[particle[track]];
		//   Set the magnetic field to 7000 Gauss.
		//6500Gauss
		E2[track] = sqrt(p2[track]*p2[track] + m2[track]*m2[track]);
		T2[track] = E2[track] - m2[track];
		deltaT = deltaT + T2[track];
		if(track == 0) {
		  deltaT = deltaT - 2.*T2[track];
		}
		printf("p[%d] = %8.2f+/-%6.3f  chi_tr[%d] = %12.6f particle %3s\n",
				 track, p2[track], dp[track], track, chi_tr[track], Part[particle[track]]);
	 }

  //  printf("0:deltapx = %8.2f,  deltapy = %8.2f,  "
  //			"deltapz = %8.2f (in MeV/c)\n",
  //		deltapx0, deltapy0, deltapz0);
  // printf("1:deltapx = %8.2f,  deltapy = %8.2f,  "
  //		"deltapz = %8.2f (in MeV/c)\n",
  //		deltapx1, deltapy1, deltapz1);
  //printf("deltapx = %8.2f,  deltapy = %8.2f,  "
  //		"deltapz = %8.2f (in MeV/c)\n",
  //		deltapx, deltapy, deltapz);
  //printf("deltaT = %8.2f MeV\n", deltaT);
  //   vertex improving

  // try to determine what are the particles
  double maxbr = 0;
  double minbr = 1000000;
  int itrackmax = -1;
  int itrackmin = 10;
  double p2br[6];
  //first suppose that Z=2 correspond to the brightest track(if track number equal 3!)
  for(track=0; track<tracks; track++)
    {
      double avbright = (brrr[0][track]+brrr[1][track])/2.;
      avbright = TMath::Max(avbright,1.);
      if(avbright>maxbr) {maxbr = avbright; itrackmax = track;}
      if(avbright<minbr) {minbr = avbright; itrackmin = track;}
      if(tracks==3) 
	{
	}
    }
  double losspar = 0;
  for(track=0; track<tracks; track++)
    {
      //if(tracks==3) 
      //if(track != itrackmax && track != itrackmin) itracklight = track;
      double avbright = (brrr[0][track]+brrr[1][track])/2.;
      avbright = TMath::Max(avbright,1.);
      p2br[track] = p2[track]*p2[track]*avbright;
      if(tracks==3&&track==itrackmax)
	p2br[track] = p2[track]*p2[track]*4*avbright*sqrt(1. - coszD[track]*coszD[track]);
    }
  if(tracks==3)
    {
      
    }

  c4->cd();
  plc4->cd();
  text = new TPaveText(-xchamb, ychamb/2., xchamb, -100.+ychamb);
  text->SetFillColor( 0 );
  text->SetTextSize(0.05);
  for(track=0; track<tracks; track++)
    {
      double avbright = (brrr[0][track]+brrr[1][track])/2.;
      avbright = TMath::Max(avbright,1.);
      printf("p[%d] %6.0f br %5.0f P^2*br %f loss %f %f\n", track, p2[track], avbright, p2br[track], 
	     p2br[track]/p2br[itrackmin], (p2br[track]-p2br[itrackmin])/p2br[itrackmax]);
      sprintf( PStr, "p%d=%8.0f br=%8.0f LossRatio=%8.1f  %7.4f", 
  	       track, p2[track], avbright, p2br[track]/p2br[itrackmin], (p2br[track]-p2br[itrackmin])/(p2br[itrackmax]-p2br[itrackmin])); 
     text->AddText(PStr);
    }
  text->Draw();
  plc4->Draw();
  c4->Update();

  for(track=0; track<tracks; track++)
	 {
		rel_loss[track] = Loss[track]
		  *sqrt(1. - coszD[track]*coszD[track])/max_loss;
		printf(" rel_loss[%d] = %f - relative to maximum losses "
				 "(H3 or He3)\n",
				 track,rel_loss[track]);
	 }
  //29.12.2002 VK: the call to Kinematics_Gil CANCELLED!!!
  //Kinematics_Gil();
  
  plc4->cd();
  pl2 = new TPaveLabel(-0.99*xchamb,0.9*(base+ychamb),-0.05*xchamb,0.99*(base+ychamb),
							  ifilename,"br");
  pl2->SetFillColor(42);
  pl2->SetTextSize(0.5);
  pl2->Draw();
  char PStr2[80];
  TPaveText *text2 = new TPaveText(0.*xchamb, 0.65*(base+ychamb), 0.99*xchamb, 0.99*(base+ychamb));
  text2->SetFillColor( 42 );
  for(track=0; track<tracks; track++)
	 {
		for(i=0; i<npt3D[track]; i++) 
		  {
			 phi[i] = atan2((y3D[track][i]-ycD[0][track]),
								 (x3D[track][i]-xcD[0][track]));
			 xx = xcD[0][track] + rcD[0][track]*cos(phi[i]);
			 yy = ycD[0][track] + rcD[0][track]*sin(phi[i]);
			 mark = new TMarker((Float_t) xx, (Float_t) yy, 5);
			 mark->SetMarkerColor(1);
			 mark->SetMarkerStyle(8);
			 mark->SetMarkerSize(0.1);
			 fObjTbl->Add(mark);
			 mark->DrawMarker((Float_t) xx, (Float_t) yy);
		  }
		sprintf( PStr2, "p%d %6.2f+/-%6.4f  is %3s", 
					track, p2[track], dp[track], Part[particle[track]]); 
		text2->AddText(PStr2);
	 }
  sprintf( PStr2, "--------------------------------------------");
  text2->AddText(PStr2);
  sprintf( PStr2, "deltaT = %8.2f MeV", deltaT);
  text2->AddText(PStr2);
  //text2->Draw();
  plc4->Modified();
  plc4->Update();
  
  c4->cd();
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
  strcpy(gifname,"gif/");
  strcat(gifname, eventname);
  strcat(gifname, "cp.gif");
  cp->Print(gifname,"gif");
  strcpy(gifname,"gif/");
  strcat(gifname, eventname);
  strcat(gifname, "cw.gif");
  cw->Print(gifname,"gif");
  strcpy(gifname,"gif/");
  strcat(gifname, eventname);
  strcat(gifname, "csl.gif");
  csl->Print(gifname,"gif");
  strcpy(gifname,"gif/");
  strcat(gifname, eventname);
  strcat(gifname, "c4.gif");
  c4->Print(gifname,"gif");
  change_part = 0;
  Save_Results();

  //theApp.Run(kTRUE);
  //sleep(10);
  return(0);
}
void Save_Results()
{
  int i, j;
  fprintf(print_results," Event type: %d\n",event_type);
  fprintf(print_results, " Number of tracks: %d\n", tracks);
  fprintf(print_results, " Direction cosines for [left-to-right, right-to-left way]: X Y Z directions\n");
  for(i = 0; i<tracks; i++)
    for(j = 0; j<2; j++) 
      fprintf(print_results, " Track %d way %d: %f %f %f\n", 
	      i, j, cosxD[i], cosyD[i], coszD[i]);
  fprintf(print_results, " Vertex position: %f %f %f  error in XY plane %f\n", 
	  Xvertax, Yvertax, Zvertax, Dvertax);
  for(j = 0; j< 2; j++)
    for(i = 0; i<tracks; i++)
      fprintf(print_results, 
	      " Measured radiuses image %d track %d: %f center %f %f\n", 
	      j, i, Radius0[j][i], xcenter01[j][i], 
	      ycenter01[j][i]);
	      //write errors of Radius, a_line and of momentum instead!!!
  //j, i, dp[i], drcD[j][i], da[j][i]);
  
  fprintf(print_results, 
	  " Reconstructed radii, X,Y center positions, chi2 [left-to-right right-to-left ways]\n:");
  for(i = 0; i<tracks; i++)
    for(j = 0; j<2; j++) 
      fprintf(print_results, "Track %d way %d: r %f x %f y %f chi %f\n", 
	      i, j, rcD[j][i],  xcD[j][i],  ycD[j][i], chi2_circle[j][i]);
  fprintf(print_results, 
	  " Reconstructed line in (PHI;Z) view A_line B_line chi2 [left-to-right right-to-left ways]\n");
  for(i=0; i<tracks; i++) 
    for(j = 0; j<2; j++) 
      fprintf(print_results, 
	      " Track %d way %d:  a %f b %f chi %f\n",
	      //	      i, j, a[j][i], chi_tr[i], chi2_line[j][i]);
	      i, j, abright[j][i], tracklength[j][i], rompture[j][i]);

  //now z coordinate of the vertex is saved
  for(i = 0; i<tracks; i++)
    fprintf(print_results, 
	    " Track length: track %d in XY plane %f in space %f \n", 
				i, lxy[i], Zv[i]);//L[i]
  fprintf(print_results, 
	  " Average image brightness: on image 0 %f on image 1 %f\n", deltaV, disp);

  for(i=0; i<tracks; i++) 
    fprintf(print_results, 
	    " Average brightness along track %d on image 0 %f on image 1 %f\n",
	    i, brrr[0][i], brrr[1][i]);

  for(i = 0; i<tracks; i++)
    fprintf(print_results, 
	    " track %d particle %3s P %f, E %f, T %f, Loss %f\n", 
	    i, Part[particle[i]], p2[i], E2[i], (float) isig_tr[i], Loss[i]);
  fprintf(print_results, 
	  "deltapx = %8.2f,  deltapy = %8.2f,  "
	  "deltapz = %8.2f (in MeV/c)\n",
	  deltapx, deltapy, deltapz);
  fprintf(print_results, "deltaT = %8.2f MeV\n", deltaT);
  fprintf(print_results, 
	  " After Kinematics routine :\n");
  for(i = 0; i<tracks; i++)
    fprintf(print_results, 
	    " track %d Relative Loss %f\n", 
	    i, rel_loss[i]);
  fprintf(print_results, 
	  " Mass proton_H3 %f\n", Mpr_trit-mass_part[5]);
  fprintf(print_results, 
	  " Complanar angle 3_12 %f deg.\n", ang3_12);
  fprintf(print_results," Angle proton_tritium %f\n", angprtr);
  fprintf(print_results," Excess of mass %f\n", excess_m);




  for(i = 0; i<tracks; i++)
    fprintf(print_results,
	    "Track %2d hits counter %1d in x %7.1f y %7.1f z %7.1f\n", 
	    i, counter[i], xinter[i], yinter[i], 
	    zinter[i]);
  
  fclose(print_results);

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
void read_ascii()
{
  int i, j, image;
  char ctmp[10000];
  char lfile[1000], rfile[1000];
  char filename[100], tfilename[100];
  char *cptmp1, *cptmp2;

  fObjTbl = new TObjectTable(10000);
  //fi.fFileTypes = (char **)filetypes;
  fi->fFileTypes = filetypes;
  
  strcpy(WorkDir,gSystem->WorkingDirectory());
  
  strcpy(iniDir,WorkDir);
  strcat(iniDir,"/data");
  printf(" iniDir %s\n", iniDir);
  
  gSystem->ChangeDirectory(iniDir);
  /* print files in current directory in reverse order */
  struct dirent **namelist;
  int n;
  int find_file = 0;
  
  n = scandir(".", &namelist, 0, alphasort);
  for(int iname = 0; iname <n; iname++)
    {
      fi->fFilename = namelist[iname]->d_name;
      printf("%s \n", fi->fFilename);
      strcpy(lfile,gSystem->BaseName(fi->fFilename));
      if(lfile[0] != 'L' && lfile[0] != 'l') continue;
      /*    forming whole filenames islower  */
      strcpy(lfile,gSystem->BaseName(fi->fFilename));
		// take only file which fits to eventname
      gSystem->ChangeDirectory(WorkDir);
      find_file = 0;
      // take the first part of event name up to the first point
      tokenPtr = strtok(lfile,".");
      if(tokenPtr == NULL) continue;
      strcpy(filename, lfile);
      tokenPtr = strtok(filename,".");
      if(tokenPtr == NULL) continue;
      strcpy(tfilename, eventname);
      tokenPtr = strtok(tfilename,".");
      int imax = strlen(filename);
      int imax_gif = strlen(tfilename);
      //if(imax-imax_gif != 1) continue;
      if(imax-imax_gif != 0) continue;
		printf("filename %s eventname %s\n", filename, eventname);
      for(i = 1; i < imax; i++)
		  {
			 if(i-1 < imax_gif)
				{
				  if(eventname[i-1]!=filename[i] ) 
					 {find_file = 0; break;}
				  else
					 find_file++;
				}
		  }
      if(find_file > 0) {
		  strcpy(lfile,gSystem->BaseName(fi->fFilename));
		  strcpy(rfile,lfile);
		  printf("%s event file \n", eventname); 
		  break;}
	 }
  gSystem->ChangeDirectory(iniDir);
  if (isupper(rfile[0]))
	 rfile[0] = 'R';
  else
	 rfile[0] = 'r';
  if(find_file==0) 
	 {
		printf("program stops because of bad fit with ascii data, no corresponding file found!\n");
		//stop program because of the absence  of the ascii file in data directory
		exit(1);
	 }
  std::cout<< "lfile = "<<lfile<<std::endl;
  std::cout<< "rfile = "<<rfile<<std::endl;

  cptmp1 = gSystem->ConcatFileName(iniDir,(char *)lfile);
  cptmp2 = gSystem->ConcatFileName(iniDir,(char *)rfile);
  /****************   open data files   ****************/
  if ((im_left = gzopen(cptmp1,"r")) == NULL)
	 {
		fprintf(stderr,"\n Cannot open input file %s\n",cptmp1);
		exit(1);
	 }
  if ((im_right = gzopen(cptmp2,"r")) == NULL)
	 {
		fprintf(stderr,"\n Cannot open input file %s\n",cptmp2);
		exit(1);
	 }
  gSystem->ChangeDirectory(WorkDir);
  //  printf("current folder  - %s \n",gSystem->WorkingDirectory());
  
  gzgets(im_left, ctmp,10000);
  
  int ntmp;
  std::istringstream Inbuffer(ctmp);
  
  int length = 0;
  while ( Inbuffer >> ntmp )
	 length++;
  
  printf("Numbers of INTs = %d\n",length);
  gzrewind(im_left);
  
  int length1 = 0;
  while( gzgets(im_left, ctmp,10000) != NULL)
	 length1++;
  
  gzrewind(im_left);
  //vsize = length1; hsize = length;
  printf("vsize = %d, hsize = %d\n",length1,length);
  
  if ((pCDBuffer[0] = Alloc_Buffer(vsize,hsize)) == NULL)
	 {
		printf("memory allocation error - aborting\n");
		exit(1);
	 }
  if ((pCDBuffer[1] = Alloc_Buffer(vsize,hsize)) == NULL)
	 {
		printf("memory allocation error - aborting\n");
		exit(1);
	 }
  
  unsigned Bright;
  for ( image=0; image<2; image++ )
	 {
		for ( i=0; i<vsize; i++ )
		  {
			 if (image == 0)
				gzgets(im_left,ctmp,10000);
			 else
				gzgets(im_right,ctmp,10000);
			 
			 std::istringstream Inbuffer(ctmp);
			 /*    j - x coordinate ;  i - y coordinate   */ 
			 for ( j=0; j<hsize; j++ )
				{
				  Inbuffer >> Bright;
				  //	      printf(" i,j = %i \n",Bright);
				  pCDBuffer[image][i][j] = Bright;
				  //if ( Bright > maxBright ) 
				  //maxBright = Bright;
				  ////	      getc(stdin);
				}
		  }
		if (image == 0) gzclose(im_left);
		else gzclose(im_right);
	 }
  left = new TH2F("left","left image", 
						(Int_t) length,0,(Axis_t) length,(Int_t) length1,0,(Axis_t) length1);
  
  right = new TH2F("right","right image",
						 (Int_t) length,0,(Axis_t) length,(Int_t) length1,0,(Axis_t) length1);

  Axis_t xx,yy;
  Stat_t ww;
  for (image=0; image<2; image++)
    {
      for (i=0; i<vsize; i++)
		  {
			 for (j=0; j<hsize; j++)
				{
				  ww = (Stat_t) pCDBuffer[image][i][j]; ///maxBright;
				  xx = (Axis_t) (j + 0.5);
				  yy = (Axis_t) (vsize - i - 0.5);  /* histogram filling in oppozitive*/ 
				  //   cout << "brightness ="<<pCDBuffer[image][i][j]<<" weight = "
				  //		     << ww << endl;
				  //	  getc(stdin);		
				  if (image == 0) 
					 {
						left->Fill(xx,yy,ww);
						if (xx>xlmin&&xx<xlmax&&yy>ylmin&&yy<ylmax)
						  {
							 if ((int)ww>maxl) maxl=(int)ww;
							 if ((int)ww<minl) minl=(int)ww;
						  }
					 }
				  else
					 {
						right->Fill(xx,yy,ww);
						if (xx>xrmin&&xx<xrmax&&yy>yrmin&&yy<yrmax)
						  {
							 if ((int)ww>maxr) maxr=(int)ww;
							 if ((int)ww<minr) minr=(int)ww;
						  }
					 }
				}
		  }
    }
  float ratio;
  ratio = (float) vsize/(float) hsize;
  
  float yleft, yright, ysize, xsize, xleft, xright;
  float ylefttmp, yrighttmp, xlefttmp, xrighttmp; 
  int ixWindow, iyWindow;
  float xWindow, yWindow, Wratio;
  float rhsize, rvsize;
  rhsize = (float) hsize;
  rvsize = (float) vsize;
  xWindow = 2.*rhsize + 4.*rhsize/100.;
  yWindow = 1.45*rvsize + 4.*rvsize/100.;
  Wratio = yWindow/xWindow;
  if (Wratio <= 1.0)
	 {
		xleft  = rhsize/100./xWindow;
		xright = 0.5 - xleft;
		yleft  = (0.45*rvsize + rvsize/100.)/yWindow;
		yright = 1.0 - rvsize/100./yWindow;
		ysize  = yright - yleft;
		xsize  = xright - xleft;
		ylefttmp = rvsize/100./yWindow;
		yrighttmp = (0.45*rvsize - rvsize/100.)/yWindow;
		xlefttmp = 0.5 - (yrighttmp - ylefttmp)/2.*Wratio;
		xrighttmp = 0.5 + (yrighttmp - ylefttmp)/2.*Wratio;
		
		xWindow = 1000;
		yWindow = xWindow*Wratio;
		ixWindow = (int) xWindow;
		iyWindow = (int) yWindow;
		//      printf("ratio = %f,%f,%f,%f,%f\n",Wratio,xleft,xright,yleft,yright);
		
	 }
  else
	 {
	 }
  c2 = new TCanvas("c2","Event digitizing",10,10,ixWindow,iyWindow);
  pleft = new TPad("pleft","Left image",    xleft, yleft, xright, yright, 21);
  pright = new TPad("pright","Right image", xleft + 0.5, yleft, xright + 0.5, yright, 21);
  //  ptmp = new TPad("ptmp"," ",xlefttmp, ylefttmp, xrighttmp, yrighttmp, 21);
  
  pMsg = new TPad("pMsg"," ",0.01,0.01,xlefttmp-0.01,yrighttmp,28);
  
  pleft->Draw();
  pright->Draw();
  
  pleft->UseCurrentStyle();
  pright->UseCurrentStyle();
  
  //  ptmp->Draw();
  //ptmp->UseCurrentStyle();
  
  pMsg->Draw();

  pMsg->cd();
  pl1 = new TPaveLabel(0.07,0.15,0.93,0.45,eventname,"br");
  pl1->SetFillColor(45);
  pl1->SetTextSize(0.5);
  pl1->Draw();
  left->SetContour(200);
  right->SetContour(200);
  pleft->cd();
  left->Draw("col");
  pright->cd();
  right->Draw("col");
  
  c2->Modified(); 
  c2->Update();
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
  
  cp = new TCanvas("cp","Profile br along lenght",30,30,1200,800);
  cp->ToggleEventStatus();
  cp->Divide(2,1);
  
  fgets(eventname, 19, for_reconstr);
  read_ascii();
  cw = new TCanvas("hwidth","br along radius",25,25,1000,800);
  cw->Divide(2,1);

  gStyle->cd();
  cw->SetFillColor(0);

  csl = new TCanvas("hslice","br along radius-background",20,20,1000,800);
  csl->Divide(2,1);

  gStyle->cd();
  csl->SetFillColor(0);

  fputs(eventname, print_results);
  fscanf(for_reconstr,"%d  %d\n",&event_type, &tracks);
  printf("Event_type = %d, \n ",event_type);
  printf("there are  %d  tracks:\n\n",tracks);
  //cin >> buf;
  for(view=0; view<2; view++)
	 {
		cp->cd(view+1);
		cp0 = (TPad*)gPad;
		gPad->UseCurrentStyle();
		cp0->Divide(1,tracks); 
		cp0->UseCurrentStyle();
		
		cw->cd();
		cw->cd(view+1);
 		cw0 = (TPad*)gPad;
		gPad->UseCurrentStyle();
		cw0->Divide(1,tracks); 
		cw0->UseCurrentStyle();
		gStyle->SetPadTopMargin(0.12);
		gStyle->SetPadLeftMargin(0.12);
		gStyle->SetPadRightMargin(0.12);
		gStyle->SetPadBottomMargin(0.12);
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		
		csl->cd();
		csl->cd(view+1);
 		csl0 = (TPad*)gPad;
		gPad->UseCurrentStyle();
		csl0->Divide(1,tracks); 
		csl0->UseCurrentStyle();
		
		//fscanf(for_reconstr,"%d\n",&ima);
		for(itrack=0; itrack<tracks; itrack++)
		  {
			 fscanf(for_reconstr,"%d",&track);
			 cp0->cd(track+1);
			 gPad->SetTopMargin(0.005);
			 gPad->SetLeftMargin(0.005);
			 gPad->SetRightMargin(0.005);
			 gPad->SetBottomMargin(0.005);
			 //			 gPad->SetFillColor();
				
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
 			 fpx[track] = Xcircle[imin];      
			 fpy[track] = Ycircle[imin];
			 lpx[track] = Xcircle[imax];     
			 lpy[track] = Ycircle[imax];
			 fx = Xcircle[imin];
			 fy = Ycircle[imin];
			 lx = Xcircle[imax];
			 ly = Ycircle[imax];
			 double delta_angle = (bufmax-bufmin)/TMath::Pi()*180./100./
			   sqrt(pow(lx-fx,2) + pow(ly-fy,2));
			 //make a variable definding if track is a circle or the line
			 circle = 1;
			 if(delta_angle < 0.00001) circle=0;
			 if( circle) tracklength[view][track] = (bufmax-bufmin)*rc_int/2.;
			 else tracklength[view][track] = sqrt(pow(lx-fx,2) + pow(ly-fy,2));
			 /*
			 c2->cd();
			 if(view==0)  pleft->cd();
			 if(view==1)  pright->cd();
			 if(circle)
			   {
			     TEllipse *myarc = new TEllipse(xc_int/2., yc_int/2., rc_int/2., 0, 
							    bufmin*180./TMath::Pi(), bufmax*180./TMath::Pi(), 0);
			     
			     fObjTbl->Add(myarc);
			     myarc->SetLineWidth(0.02);
			     myarc->SetLineColor(2);
			     myarc->SetNoEdges(1);
			     myarc->Draw("same");	    
			   }
			 else
			   {
			     TLine *myline = new TLine(fx/2.,fy/2.,lx/2.,ly/2.);
			     fObjTbl->Add(myline);
			     //myarc->SetMarkerStyle(20);
			     myline->SetLineWidth(0.02);
			     myline->SetLineColor(2);
			     myline->Draw("same");	    
			   }
			 c2->Modified(); 
			 c2->Update();
			 */
			 char ctmp1[3],ctmp2[3], name[10], name2[10];
			 int nchannels = (int) (tracklength[view][track]);
			 sprintf(ctmp1,"%u",view);
			 sprintf(ctmp2,"%u",track);
			 strcpy(name,"h");
			 strcat(name,ctmp1);
			 strcat(name,"tr");
			 strcat(name,ctmp2);
			 float YMIN = 0.;
			 float YMAX = 4096.;
			 //float b1, b2, b3;
			 cw->cd();
			 cw->cd(view+1);
			 cw0->cd(track+1);
			 cw0->UseCurrentStyle();
			 gStyle->SetPadTopMargin(0.12);
			 gStyle->SetPadLeftMargin(0.12);
			 gStyle->SetPadRightMargin(0.12);
			 gStyle->SetPadBottomMargin(0.12);
			 gStyle->SetHistFillColor(4);
			 //gStyle->SetOptFit(1);
			 //hwidth = new TH1F("hwidth","hw", 101,-50.,51.);
			 hwidth = new TH1F("hwidth","hw", 41,-20.,21.);
			 float iw;
			 float shift=0.;
			 find_widthofthecorridor(&bufmin, &bufmax,
						 xc_int/2., yc_int/2.,rc_int/2.,
						 fx/2., fy/2., lx/2., ly/2., view, itrack, &iw, &shift);
			 printf("SHIFT on radius correction was found as %7.2f pixels\n", shift);
			 /*
			 c2->cd();
			 if(view==0)  pleft->cd();
			 if(view==1)  pright->cd();
			 if(circle)
			   {
			     TEllipse *myarc = new TEllipse(xc_int/2., yc_int/2., rc_int/2., 0, 
							    bufmin*180./TMath::Pi(), bufmax*180./TMath::Pi(), 0);
			     
			     fObjTbl->Add(myarc);
			     myarc->SetLineWidth(0.02);
			     myarc->SetLineColor(3);
			     myarc->SetNoEdges(1);
			     myarc->Draw("same");	    
			   }
			 else
			   {
			     TLine *myline = new TLine(fx/2.,fy/2.,lx/2.,ly/2.);
			     fObjTbl->Add(myline);
			     //myarc->SetMarkerStyle(20);
			     myline->SetLineWidth(0.02);
			     myline->SetLineColor(3);
			     myline->Draw("same");	    
			   }
			 c2->Modified(); 
			 c2->Update();
			 */
			 cp->cd();
			 cp->cd(view+1);
			 cp0->cd(track+1);
			 cp1 = (TPad*)gPad;
			 gPad->UseCurrentStyle();
			 cp1->Divide(4,1); 
			 cp1->UseCurrentStyle();
			 gStyle->SetOptStat(0);
			 gStyle->SetOptFit(0);
			 delta_angle = (bufmax-bufmin)/TMath::Pi()*180./100./
			   sqrt(pow(lx-fx,2) + pow(ly-fy,2));
			 //make a variable definding if track is a circle or the line
			 circle = 1;
			 if(delta_angle < 0.00001) circle=0;
			 if( circle) tracklength[view][track] = (bufmax-bufmin)*rc_int/2.;
			 else tracklength[view][track] = (sqrt(pow(lx-fx,2) + pow(ly-fy,2)))/2.;
			 nchannels = (int) (tracklength[view][track]);
			 float step = 2.*iw; 
			 //step to shift in pxls equal TWO sigma value found 
			 //from previous find_widthofthecorridor
			 cp->cd();
			 cp->cd(view+1);
			 cp0->cd(track+1);
			 cp1->cd(2);
			 strcpy(name2,name);
			 strcat(name2,"2");
			 hprof2  = new TProfile(name2,name2,nchannels,0,
						tracklength[view][track]+1.,YMIN,YMAX);
 			 float fx1 = fx+step*cos(bufmin);      
			 float fy1 = fy+step*sin(bufmin);
			 float lx1 = lx+step*cos(bufmax);     
			 float ly1 = ly+step*sin(bufmax);
			 int iwfixed = 4;
			 find_points_close_to_circle(hprof2,bufmin,bufmax,
						     xc_int/2.,yc_int/2.,rc_int/2.+step/2., rc_int/2.,
						     fx1/2., fy1/2., lx1/2., ly1/2., view, iwfixed);
			 hprof2->Draw();
			 cp->cd();
			 cp->cd(view+1);
			 cp0->cd(track+1);
			 cp1->cd(3);
 			 fx1 = fx-step*cos(bufmin);      
			 fy1 = fy-step*sin(bufmin);
			 lx1 = lx-step*cos(bufmax);     
			 ly1 = ly-step*sin(bufmax);
			 strcpy(name2,name);
			 strcat(name2,"3");
			 hprof3  = new TProfile(name2,name2,nchannels,0,
						tracklength[view][track]+1.,YMIN,YMAX);
			 find_points_close_to_circle(hprof3,bufmin,bufmax,
						     xc_int/2.,yc_int/2.,rc_int/2.-step/2., rc_int/2.,
						     fx1/2., fy1/2., lx1/2., ly1/2., view, iwfixed);
			 hprof3->Draw();
			 cp->cd();
			 cp->cd(view+1);
			 cp0->cd(track+1);
			 cp1->cd(1);
			 strcpy(name2,name);
			 strcat(name2,"1");
			 hprof  = new TProfile(name2,name2,nchannels,0,
					       tracklength[view][track]+1.,YMIN,YMAX);
			 find_points_close_to_circle(hprof, bufmin, bufmax,
						     xc_int/2., yc_int/2.,rc_int/2.,rc_int/2.,
						     fx/2., fy/2., lx/2., ly/2., view, (int) iw);
			 hprof->Draw();

			 //brrr[view][track] = TMath::Max(1., b1-(b2+b3)/2.);

			 float bincont, binentr, binerr;
			 float bincont2, binerr2;
			 float bincont3, binerr3;
			 strcpy(name2,name);
			 strcat(name2,"4");
			 cp->cd();
			 cp->cd(view+1);
			 cp0->cd(track+1);
			 cp1->cd(4);
			 TH1F *th4  = new TH1F(name2,name2,nchannels,0,
					       tracklength[view][track]+1.);

			 int n0 = 0;
			 int Nentr = 0;
			 for(i=1;i<nchannels+1;i++)
				{
				  if(hprof->GetBinContent(i)==0.)
				    {
				      bincont = hprof->GetBinContent(i-1);
				      binentr = hprof->GetBinEntries(i-1);
				      binerr  = hprof->GetBinError(i-1);
				    }
				  else
				    {
				      bincont = hprof->GetBinContent(i);
				      binentr = hprof->GetBinEntries(i);
				      binerr  = hprof->GetBinError(i);
				    }

				  if(hprof2->GetBinContent(i)==0.)
				    {
				      bincont2 = hprof2->GetBinContent(i-1);
				      bincont -= bincont2/2.;
				      binerr2 = hprof2->GetBinError(i-1);
				    }
				  else
				    {
				      bincont2 = hprof2->GetBinContent(i);
				      bincont -= bincont2/2.;
				      binerr2 = hprof2->GetBinError(i);
				    }

				  if(hprof3->GetBinContent(i)==0.)
				    {
				      bincont3 = hprof3->GetBinContent(i-1);
				      bincont -= bincont3/2.;
				      binerr3 = hprof3->GetBinError(i-1);
				    }
				  else
				    {
				      bincont3 = hprof3->GetBinContent(i);
				      bincont -= bincont3/2.;
				      binerr3 = hprof3->GetBinError(i);
				    }
				  //th4->SetBinContent(i, 0.);
				  //th4->SetBinError(i, 100000.);
				  binerr  =  sqrt(pow(binerr,2) + pow(binerr2,2) + pow(binerr3,2));
				  //printf(" bincont %f binentr %f binerr %f\n", bincont, binentr, binerr);
				  if(bincont > 0) 
				    {n0++;
					 Nentr += (int) binentr;
					 //xxx = hprof->GetBinCenter(i);
					 //printf("i %d x %f\n", i, xxx);
					 th4->SetBinContent(i, bincont);
				         th4->SetBinError(i, binerr);
					 }
				  //define here slope of background for given slice
				  if(circle) 
				    {
				      p1[i] = (bincont2 - bincont3)/step;
				      p0[i] = bincont2 - p1[i]*(rc_int/2.+step/2.);
				    }
				  else
				    {
				      p1[i] = (bincont2 - bincont3)/step;
				      p0[i] = bincont2 - p1[i]*step/2.;
				    }
				}
			 //hprof4->Draw();
			 th4->Fit("pol1");
			 TF1 *fit = th4->GetFunction("pol1");
			 Double_t ch = fit->GetChisquare();
			 abright[view][track] = fit->GetParameter(1);
			 Double_t bbright = fit->GetParameter(0);
			 Double_t ea = fit->GetParError(1);
			 Double_t eb = fit->GetParError(0);
			 printf("chi2 %f a %f+/-%f b %f+/-%f\n", ch, abright[view][track], ea, bbright, eb);
			 float brightness = 0.5*(bbright+abright[view][track]*tracklength[view][track]+bbright)*
			   ((float) n0)/tracklength[view][track];
			 rompture[view][track] = ((float) n0)/tracklength[view][track];
			 TText *t = new TText();
			 t->SetTextFont(62);
			 t->SetTextColor(2);
			 t->SetTextSize(0.10);
			 t->SetTextAlign(12);
			 char text[100];
			 sprintf(text,"%6.0f %4.2f\n", brightness,rompture[view][track]);
			 float profmax = th4->GetMaximum();
			 t->DrawText((double) nchannels/4., (double) 0.95*profmax,text);
			 brrr[view][track] = TMath::Max(1., (double) brightness);
			 printf("brightness for track %d is %f\n", track, brrr[view][track]);
			 //here balance the brigness around circle by substracting background 
			 //along each slice along R
			 csl->cd();
			 csl->cd(view+1);
			 csl0->cd(track+1);
			 //csl0->UseCurrentStyle();
			 gStyle->SetPadTopMargin(0.1);
			 gStyle->SetPadLeftMargin(0.1);
			 gStyle->SetPadRightMargin(0.1);
			 gStyle->SetPadBottomMargin(0.1);
			 gStyle->SetHistFillColor(3);
			 int iwsl = 4;
			 hslice = new TH1F("hslice","slice", iwsl*2+1,(Float_t) -iwsl,(Float_t) (iwsl+1));
			 find_br_for_circlefit(hprof, bufmin, bufmax,
					       xc_int/2., yc_int/2.,rc_int/2.,rc_int/2.,
					       fx/2., fy/2., lx/2., ly/2., view, iwsl);
			 cp->Modified();
			 cp->Update();
			 c2->Modified();
			 c2->Update();
			 csl->Modified();
			 csl->Update();
			 //

			 printf(" Npoints_fit %d\n", Npoints_fit);
			 //redefinition of the points to make a circle fit (XFit-> To Xcircle and etc)
			 /*
			 for(i = 0; i<Npoints_fit; i++) 
			   {
			     x_mark[i] = Xcircle[i];
			     y_mark[i] = Ycircle[i];
			     wes[0][i] = Wpoint[i];
			   }
			 int Nbuffer = Npoints_fit;
			 for(i = 0; i<NFit; i++) 
			   {
			     Xcircle[i] = XFit[i];
			     Ycircle[i] = YFit[i];
			     Wpoint[i] = WFit[i];
			   }
			 Npoints_fit = NFit;
			 */

			 Signe=10;
			 XC = xc_int;
			 YC = yc_int;
			 RC = rc_int;
			 circle_minuit(Npoints_fit);
			 printf("ASCII image circle results: %7.1f %7.1f %7.1f  \n\n",
			 	  RC, XC, YC);
			 xc_int = XC;
			 yc_int = YC;
			 rc_int = RC;
			 //int buf1;
			 //cin>>buf1;
			 bufmin = 1000.;
			 bufmax =-1000.;
			 for(i=0;i<Npoints_fit;i++)
			   {
			     buf = atan2((Ycircle[i]-yc_int),(Xcircle[i]-xc_int));
			     if(buf<0.) buf += TMath::Pi()*2.;
			     if(buf < bufmin)
			       { imin = i; bufmin = buf;}
			     if(buf > bufmax)
			       { imax = i; bufmax = buf;}
			   }
			 //check for the -Pi - +Pi bounds
			 if((bufmax - bufmin) > TMath::Pi()) 
			   {
			     bufmin = 1000.;
			     bufmax =-1000.;
			     for(i=0; i<Npoints_fit; i++)
			       {
				 buf = atan2((Ycircle[i]-yc_int),(Xcircle[i]-xc_int));
				 if(buf < bufmin)
				   { imin = i; bufmin = buf;}
				 if(buf > bufmax)
				   { imax = i; bufmax = buf;}
			       }
			   }
			 xc_int = XC;
			 yc_int = YC;
			 rc_int = RC;
			 /*
			 c2->cd();
			 if(view==0)  pleft->cd();
			 if(view==1)  pright->cd();
			 if(circle)
			   {
			     TEllipse *myarc = new TEllipse(xc_int/2., yc_int/2., rc_int/2., 0, 
							    bufmin*180./TMath::Pi(), bufmax*180./TMath::Pi(), 0);
			     
			     fObjTbl->Add(myarc);
			     myarc->SetLineWidth(0.02);
			     myarc->SetLineColor(5);
			     myarc->SetNoEdges(1);
			     myarc->Draw("same");	    
			   }
			 else
			   {
			     TLine *myline = new TLine(fx/2.,fy/2.,lx/2.,ly/2.);
			     fObjTbl->Add(myline);
			     //myarc->SetMarkerStyle(20);
			     myline->SetLineWidth(0.02);
			     myline->SetLineColor(5);
			     myline->Draw("same");	    
			   }
			 c2->Modified(); 
			 c2->Update();
			 */
			 tr_len[view][track] = tracklength[view][track]*2.;//(bufmax-bufmin)*RC;
			 /*
			 WpMin = 1000.;
			 WpMax =-1000.;
			 for(i = 0; i<Npoints_fit; i++) 
			   {
			     if(Wpoint[i]>WpMax) WpMax = Wpoint[i];
			     if(Wpoint[i]<WpMin) WpMin = Wpoint[i];
			   }
			 */
			 //WpMin -= 1;
			 //WpMax += 1;
			 k = 0;
			 IntBr[view][track] = 0.;
			 //			 Npoints_fit = n;
			 for(i = 0; i<Npoints_fit; i++) 
			   {
			     //if(WpMin != WpMax) Wpoint[i]  = (Wpoint[i]-WpMin)/(WpMax - WpMin);
			     //if(WpMin == WpMax) Wpoint[i]  = Wpoint[i]/WpMax;
			     double dist = sqrt(pow((Ycircle[i]-YC),2) + 
						pow((Xcircle[i]-XC),2)) - RC;
			     IntBr[view][track] += Wpoint[i];
			     //if(fabs(dist) < 1.55) 
			     if(fabs(dist) < 2.55 && k<30000-1) 
			       {
				 CCDBuffer[view][track][k] = Wpoint[i];
				 //  printf(" %3d %3d %4d ", jj, ii, bright);
				 xtrack[view][track][k] = (float) Xcircle[i]; 
				 ytrack[view][track][k] = (float) Ycircle[i];
				 k++;
			       }
			   }
			 np[view][track] = k;
			 printf("points number along track CIRCLE is %d\n", 
				k);
			 Radius[view][track] = RC;
			 xcenter[view][track] = XC;
			 ycenter[view][track] = YC;
			 
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
  
  cp->Update();
  c2->cd();
  TPaveLabel *pbr = new TPaveLabel(0.51,0.01,0.99,0.3,eventname,"br");
  pbr->SetFillColor(42);
  pbr->SetTextSize(0.3);
  pbr->Draw();
  char PStr[60];
  TPaveText *textbr = new TPaveText(0.51, 0.01, 0.99, 0.3 );
  textbr->SetFillColor( 42 );
  for(track=0; track<tracks; track++)
    {
		sprintf( PStr, "track %1d brightness %6.0f %4.2f %6.0f %4.2f",
					track, brrr[0][track],rompture[0][track],brrr[1][track],rompture[1][track]); 
		textbr->AddText(PStr);
    }

  textbr->Draw();
  c2->Update();


  //int buf1;
  //cin>>buf1;
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

void line(int N, float *X, float *Y, float *A, float *B, float *DA, float *VAR)
{
  double S, Sx, Sy, Sxx, Sxy, delta, xt, yt, wt, buf;
  //call norm_weights(n,w,w_norm)
  float RES;
  //     clear all accumulators
  float SUMRES = 0.;
  int i;
  //     accumulate sums
  S   = Sx  = Sy  = Sxx = Sxy = (double) 0.0;
  *A  = 0.0;
  int NPOS = 0;

  for(i = 0; i<N; i++)
    {
      wt = (double) Wpoint[i];
      xt = (double) *(X+i);
      yt = (double) *(Y+i);
      S  += wt;
      Sx += xt*wt;
      Sy += yt*wt;
      Sxx += (xt*xt)*wt;
      Sxy += xt*yt*wt;
    }
  delta = S*Sxx-Sx*Sx;
  if(delta == 0.) delta = 0.000001;
  buf = (S*Sxy-Sx*Sy)/delta;
  *A = (float) buf;
  *DA = (float) sqrt(S/delta);
  buf  = (Sxx*Sy-Sx*Sxy)/delta;
  *B = (float) buf;
  //     compute the residuals
  
  for(i=0; i<N; i++)
    {
      RES = *(Y+i)-*A*(*(X+i))-*B;
      if (Wpoint[i] > 0.0) 
	{
	  SUMRES+=pow(RES,2)*(Wpoint[i]);
	  NPOS += 1;
	}
    }
  //c     compute the estimate of variance of statistical errors
  //c     (mean square deviation)
  //c     vk it was an error here previously
  *VAR=SUMRES/(float)(NPOS-2);
}

void Find_Tails(int itr, int image)
{
  float x, y;
  float x_min, x_max, y_min, y_max;
  float x_mid, y_mid;
  int jc1, jc2, jr1, jr2;
  int LIMIT0, LIMIT1;
  x_min = y_min = 10000.;
  x_max = y_max = -10000.;
  x_mid = y_mid = 0.;
  jc1 = jc2 = jr1 = jr2 = 0;
  if (image == 0) {LIMIT0 = 0; LIMIT1 = npt3D_left[itr];}
  if (image == 1) {LIMIT0 = npt3D_left[itr]; LIMIT1 = npt3D[itr];}
  for(int i = LIMIT0; i < LIMIT1; i++)
    {
      x = x3D[itr][i];
      y = y3D[itr][i];
      if(x < x_min) {x_min = x; jc1 = i;}
      if(y < y_min) {y_min = y; jr1 = i;}
      if(x > x_max) {x_max = x; jc2 = i;}
      if(y > y_max) {y_max = y; jr2 = i;}
    }
  if(x_max-x_min > y_max - y_min) 
    {
      fp = jc1;
      lp = jc2;
      x_mid = (x_max+x_min)/2.;
    }
  else
    {
      fp = jr1;
      lp = jr2;
      y_mid = (y_max+y_min)/2.;
    }
  float dist_min = 10000.;
  for (int i = LIMIT0; i < LIMIT1; i++)
    {
      x = x3D[itr][i];
      y = y3D[itr][i];
      if(x_mid != 0. && fabs(x - x_mid) < dist_min) 
	{dist_min = fabs(x - x_mid); jc1 = i;}
      if(y_mid != 0. && fabs(y - y_mid) < dist_min) 
	{dist_min = fabs(y - y_mid); jr1 = i;}
    }
  if(x_mid != 0.) 
    mp = jc1;
  else
    mp = jr1;
  }
  
    
PPUINT Alloc_Buffer(int vsize, int hsize)
{
  PPUINT ret_val;
  int i;

  ret_val = (PPUINT) malloc(vsize*sizeof(PUINT));
  if (ret_val == NULL)
    return NULL;                                                    
  for(i=0; i<vsize; i++)
    { 
      ret_val[i] = (PUINT) malloc(hsize*sizeof(unsigned));
      if(ret_val[i] == NULL)
	return NULL;
    }
  return ret_val;
}

void Free_Buffer(PPUINT buffer[2], int vsize)
{
  int i;
  for(i = 0; i < vsize; i++)
    {
      free(buffer[0][i]);
      free(buffer[1][i]);
    }
  free(buffer[0]);
  free(buffer[1]);
}
/* new routine of SetNewPalette - root version 3.02 */
void SetVeryNewPalette()
{
  //  TColor *color;
  const int ncolors = 202;
  int itmp;
  Float_t light;
  Int_t Apalet[200];

  /*
  TStyle *palett = new TStyle; //->SetPalette(500);
  */
  /*
  for ( itmp=0; itmp<50; itmp++) {
    Apalet[itmp] = gStyle->GetColorPalette(itmp+1);
    //    cout<<"palett... "<<itmp<<" = "<<Apalet[itmp]<<endl;
    }
  */
  //  TColor *newcolor = new TColor();
 
  for ( itmp=0; itmp<ncolors; itmp++)
    {
      //      TColor *color = gROOT->GetColor(itmp+1);
      light = (Float_t)itmp/(ncolors-1);
      TColor *color = new TColor(230+itmp,light,light,light,"");

      //      color->SetRGB(light,light,light);
      //      cout<<itmp<<"    =   "<<light<<endl;
      Apalet[itmp]=230 + itmp + 1;
      //      Apalet[itmp]=0;
    }

  //  palett->SetPalette(ncolors+1,Apalet);
  gStyle->SetPalette(ncolors+1,Apalet);
  //  gStyle->SetPalette(200,Apalet);
  //  palett->SetPalette(1,0);
   //  cc->cd();  
   //  cc->Close();
  //  newcolor->DisplayColorTable();

  // theApp.Run(kTRUE);
  /*  c3 = new TCanvas("c3","right image",580,10,570,500); */


  //  TStyle *palett1 = new TStyle; //->SetPalette(500);
  //  palett->SetPalette(10);
  //  TStyle *palett_old = new TStyle;


  /*
   for ( itmp=0; itmp<50; itmp++)
    cout << "color number "<<itmp<< "palett =" 
  	 <<palett->GetColorPalette(itmp)<<endl;
  */

  //  cout << "Stat flag.."<< palett->GetOptStat()<<endl;
  //  TPaveStats *stat = new TPaveStats();
  //  left->SetStats(0);
  ////  palett->SetOptStat(0);
  gStyle->SetOptStat(0);
  //  stat->SaveStyle();


  //  left->SetFillColor(20);
  //  TStyle *mystyle=new Tstyle();
  gStyle->SetPadTopMargin(0.005);
  gStyle->SetPadBottomMargin(0.005);
  gStyle->SetPadLeftMargin(0.005);
  gStyle->SetPadRightMargin(0.005);

  gStyle->SetHistFillColor(38);
  gStyle->SetHistFillStyle(1001);
  gStyle->SetOptLogz(1);
  gStyle->cd();
  //  palett->cd();
}
void find_points_close_to_circle(TProfile *pro,
				  float fmin, float fmax, float xc, float yc, float rc, 
				  float rc0,
				  float x1, float y1, float x2, float y2, int ima, 
				  int iwidth)
{
  int i, j;
  float bright;
  double delta_angle;
  double range = fmax-fmin;
  int curx = 0; int cury = 0;
  double dr;
  int ip; int ii;
  int np = 0;
  double pibytwo = TMath::Pi()/2.;
  
  //prevent cicling in stright line tracks
  dr = 0.1*sqrt(2.); //step to move along r
  double cur_angle;
  delta_angle = range/TMath::Pi()*180./100./sqrt(pow(x2-x1,2) + pow(y2-y1,2));
  //printf("delta_angle %f\n",delta_angle);
  if(delta_angle > 0.00001) 
	 {
      cur_angle = fmin;
      int binxprev_circle = 0; int binyprev_circle = 0;
      int epx = (int) (xc + rc*cos(cur_angle+range));
      int epy = (int) (yc + rc*sin(cur_angle+range));
      //printf("epx %d epy %d, endangle %f cur %f\n", epx, epy, fmax, cur_angle);
      for(ii = 0; ii<1000000; ii++)
	{
	  ip = 0; 
	  int binxprev = 0;
	  int binyprev = 0;
	  
	  curx = (int) (xc + rc*cos(cur_angle));
	  cury = (int) (yc + rc*sin(cur_angle));
	  if (binxprev_circle != curx || binyprev_circle != cury)
	    {
	      //store this pont 
	      x_mark[np] = curx;
	      y_mark[np] = cury;
	      i = (int) (vsize - cury - 1);
	      j = (int) (curx);
	      if(i>=0&&j>=0&&i<vsize&&j<hsize)
		{
		  bright = pCDBuffer[ima][i][j];
		  pro->Fill((cur_angle-fmin)*rc0,bright,1);
		  np++;
		}
	      binxprev_circle = curx; binyprev_circle = cury;
	      binxprev = curx; binyprev = cury;
	      ip++;
	      float rcur = rc + dr;
	      do {
		curx = (int) (xc + rcur*cos(cur_angle));
		cury = (int) (yc + rcur*sin(cur_angle));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      pro->Fill((cur_angle-fmin)*rc0,bright,1);
		      np++;
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur += dr;
	      } while (ip < iwidth+1);
	      ip = 1;
	      rcur = rc - dr;
	      binxprev = binxprev_circle; 
	      binyprev = binyprev_circle;
	      do {
		curx = (int) (xc + rcur*cos(cur_angle));
		cury = (int) (yc + rcur*sin(cur_angle));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      pro->Fill((cur_angle-fmin)*rc0,bright,1);
		      np++;
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur -= dr;
	      } while(ip < iwidth+1);
	    }
	  cur_angle += delta_angle;
	  if(curx == epx && cury == epy) break;
	  if(cur_angle > (double) fmax) break;
	}
    }
  else
    {
      printf("There is a line, no circle!\n");
      double delta_l = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
      //printf("delta_l %f\n",delta_l);
      cur_angle = atan2(y2-y1,x2-x1);
      
      int binxprev_line = 0; int binyprev_line = 0;
      i = 0;
      double dl = 0.;
      float epx = x2;
      float epy = y2;
      //printf("epx %f epy %f, cur_angle %f\n", epx, epy, cur_angle);
      for(ii = 0; ii<1000000; ii++)
	{
	  ip = 0; 
	  int binxprev = 0;
	  int binyprev = 0;
	  
	  curx = (int) (x1 + dl*cos(cur_angle));
	  cury = (int) (y1 + dl*sin(cur_angle));
	  if (binxprev_line != curx || binyprev_line != cury)
	    {
	      //store this pont 
	      x_mark[np] = curx;
	      y_mark[np] = cury;
	      
	      i = (int) (vsize - cury - 1);
	      j = (int) (curx);
	      if(i>=0&&j>=0&&i<vsize&&j<hsize)
		{
		  bright = pCDBuffer[ima][i][j];
		  pro->Fill(dl,bright,1);
		  np++;
		}
	      binxprev_line = curx; binyprev_line = cury;
	      binxprev = curx; binyprev = cury;
	      ip++;
	      float rcur = dr;
	      do {
		curx = (int) (x1 + dl*cos(cur_angle) + rcur*cos(cur_angle+pibytwo));
		cury = (int) (y1 + dl*sin(cur_angle) + rcur*sin(cur_angle+pibytwo));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      pro->Fill(dl,bright,1);
		      np++;
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur += dr;
	      } while (ip < iwidth+1);
	      ip = 1;
	      rcur = -dr;
	      binxprev = binxprev_line; 
	      binyprev = binyprev_line;
	      do {
		curx = (int) (x1 + dl*cos(cur_angle) + rcur*cos(cur_angle+pibytwo));
		cury = (int) (y1 + dl*sin(cur_angle) + rcur*sin(cur_angle+pibytwo));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      pro->Fill(dl,bright,1);
		      np++;
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur -= dr;
	      } while(ip < iwidth+1);
	    }
	  dl += 0.1*sqrt(2.); //step to move along r
	  if(curx == epx && cury == epy) break;
	  if(dl >= delta_l) break;
	}
    }
  //pro->Draw();
  /*
  float bincont, binentr, binerr;
  int n0 = 0;
  int Nentr = 0;
  for(i=0;i<np;i++)
    {
      bincont = pro->GetBinContent(i);
      binentr = pro->GetBinEntries(i);
      binerr  = pro->GetBinError(i);
      //printf(" bincont %f binentr %f binerr %f\n", bincont, binentr, binerr);
      if(binentr==0) n0++;
      Nentr += (int) binentr;
    }
  pro->Fit("pol1");
  TF1 *fit = pro->GetFunction("pol1");
  Double_t ch = fit->GetChisquare();
  Double_t a = fit->GetParameter(1);
  Double_t b = fit->GetParameter(0);
  Double_t ea = fit->GetParError(1);
  Double_t eb = fit->GetParError(0);
  printf("chi2 %f a %f+/-%f b %f+/-%f\n", ch, a, ea, b, eb);
  float tracklength[view][track] = (fmax-fmin)*rc0;
  float brightness = 0.5*(b+a*tracklength[view][track]+b)*Nentr/tracklength[view][track];
  TText *t = new TText();
  t->SetTextFont(62);
  t->SetTextColor(36);
  t->SetTextSize(0.08);
  t->SetTextAlign(12);
  char text[100];
  sprintf(text," %6.0f \n", brightness);
  float profmax = pro->GetMaximum();
  int nchannels = (int) (tracklength[view][track]);
  t->DrawText((double) nchannels/4., (double) 0.95*profmax,text);
  
  printf("br %f \n", brightness);
  */
  //cp->Update();

//   c2->cd();
//   printf("npoints is %d\n", np);
//   if(ima==0)  pleft->cd();
//   if(ima==1)  pright->cd();
//   px_mark = &x_mark[0];
//   py_mark = &y_mark[0];
//   markp = new TPolyMarker((Int_t) np, px_mark, py_mark, "same");
//   fObjTbl->Add(markp);
//   markp->SetMarkerStyle(20);
//   markp->SetMarkerSize(0.2);
//   markp->SetMarkerColor(50);
//   markp->Draw("same");	    
  //c2->Modified(); 
  //c2->Update();
}
short int find_the_same_point(int i, int j, int np)
{
  for(int k = 0; k < np; k++)
    {
      if((int) x_mark[k]==i && (int) y_mark[k]==j) 
	{
	  return 1;
	}
    }
  return 0;
}
void find_widthofthecorridor(double *fmin, double *fmax, float xc, float yc, float rc, 
			      float x1, float y1, float x2, float y2, int ima, int itr, float *sigma, float *shift)
{
  int i, j;
  float bright;
  double delta_angle;
  double range = *fmax-*fmin;
  int curx = 0; int cury = 0;
  int iwidth = 20;//50.?????????/
  double dr;
  int ip; int ii;
  int np = 0;
  double pibytwo = TMath::Pi()/2.;
  
  dr = 0.1*sqrt(2.); //step to along r
  double cur_angle;
  //prevent cicling in stright line tracks
  iwidth = (int) (TMath::Min((float) iwidth, (float)(rc-10.)));
  //step along track -> delta phi
  delta_angle = range/TMath::Pi()*180./100./sqrt(pow(x2-x1,2) + pow(y2-y1,2));
  //printf("delta_angle %f\n",delta_angle);
  if(delta_angle > 0.00001) 
    {
      cur_angle = *fmin;
      int binxprev_circle = 0; int binyprev_circle = 0;
      int epx = (int) (xc + rc*cos(cur_angle+range));
      int epy = (int) (yc + rc*sin(cur_angle+range));
      //printf("epx %d epy %d, endangle %f cur %f\n", epx, epy, *fmax, cur_angle);
      for(ii = 0; ii<1000000; ii++)
	{
	  ip = 0; 
	  int binxprev = 0;
	  int binyprev = 0;
	  
	  curx = (int) (xc + rc*cos(cur_angle));
	  cury = (int) (yc + rc*sin(cur_angle));
	  if (binxprev_circle != curx || binyprev_circle != cury)
	    {
	      //store this point 
	      
	      i = (int) (vsize - cury - 1);
	      j = (int) (curx);
	      if(i>=0&&j>=0&&i<vsize&&j<hsize)
		{
		  bright = pCDBuffer[ima][i][j];
		  hwidth->Fill((Axis_t) ip, (Stat_t) (bright/(*fmax-*fmin)/(rc+ip)));
		}
	      
	      binxprev_circle = curx; binyprev_circle = cury;
	      binxprev = curx; binyprev = cury;
	      ip++;
	      float rcur = rc + dr;
	      do {
		curx = (int) (xc + rcur*cos(cur_angle));
		cury = (int) (yc + rcur*sin(cur_angle));
		if(curx != binxprev || cury != binyprev) {
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      hwidth->Fill((Axis_t) ip, (Stat_t) (bright/(*fmax-*fmin)/(rc+ip)));
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur += dr;
	      } while (ip < iwidth+1);
	      ip = 1;
	      rcur = rc - dr;
	      binxprev = binxprev_circle; 
	      binyprev = binyprev_circle;
	      do {
		curx = (int) (xc + rcur*cos(cur_angle));
		cury = (int) (yc + rcur*sin(cur_angle));
		if(curx != binxprev || cury != binyprev) {
		  
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      hwidth->Fill((Axis_t) -ip, (Stat_t) (bright/(*fmax-*fmin)/(rc-ip)));
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur -= dr;
	      } while(ip < iwidth+1);
	    }
	  cur_angle += delta_angle;
	  if(curx == epx && cury == epy) break;
	  if(cur_angle > (double) *fmax) break;
	}
    }
  else
    {
      printf("There is a line, no circle!\n");
      double delta_l = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
      //printf("delta_l %f\n",delta_l);
      cur_angle = atan2(y2-y1,x2-x1);
      
      int binxprev_line = 0; int binyprev_line = 0;
      i = 0;
      double dl = 0.;
      float epx = x2;
      float epy = y2;
      //printf("epx %f epy %f, cur_angle %f\n", epx, epy, cur_angle);
      for(ii = 0; ii<1000000; ii++)
	{
	  ip = 0; 
	  int binxprev = 0;
	  int binyprev = 0;
	  
	  curx = (int) (x1 + dl*cos(cur_angle));
	  cury = (int) (y1 + dl*sin(cur_angle));
	  if (binxprev_line != curx || binyprev_line != cury)
	    {
	      //store this pont 
	      
	      i = (int) (vsize - cury - 1);
	      j = (int) (curx);
	      if(i>=0&&j>=0&&i<vsize&&j<hsize)
		{
		  bright = pCDBuffer[ima][i][j];
		  hwidth->Fill((Axis_t) ip, (Stat_t) (bright/delta_l));
		}
	      
	      binxprev_line = curx; binyprev_line = cury;
	      binxprev = curx; binyprev = cury;
	      ip++;
	      float rcur = dr;
	      do {
		curx = (int) (x1 + dl*cos(cur_angle) + rcur*cos(cur_angle+pibytwo));
		cury = (int) (y1 + dl*sin(cur_angle) + rcur*sin(cur_angle+pibytwo));
		if(curx != binxprev || cury != binyprev) {
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      hwidth->Fill((Axis_t) ip, (Stat_t) (bright/delta_l));
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur += dr;
	      } while (ip < iwidth+1);
	      ip = 1;
	      rcur = -dr;
	      binxprev = binxprev_line; 
	      binyprev = binyprev_line;
	      do {
		curx = (int) (x1 + dl*cos(cur_angle) + rcur*cos(cur_angle+pibytwo));
		cury = (int) (y1 + dl*sin(cur_angle) + rcur*sin(cur_angle+pibytwo));
		if(curx != binxprev || cury != binyprev) {
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      hwidth->Fill((Axis_t) -ip, (Stat_t) (bright/delta_l));
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur -= dr;
	      } while(ip < iwidth+1);
	    }
	  dl += 0.1*sqrt(2.); //step to move along r
	  if(curx == epx && cury == epy) break;
	  if(dl >= delta_l) break;
	}
    }
  hwidth->Draw();
  //TF1 *myfit = new TF1("myfit","[0]*exp(-0.5*pow((x-[4])/[1],2))+[2]+[3]*x", (float) -(iwidth), (float) iwidth);
TF1 *myfit = new TF1("myfit","[0]*exp(-0.5*pow((x-[4])/[1],2))+[2]+[3]*x+[5]*exp(-0.5*pow((x-[6])/[7],2))", (float) -(iwidth), (float) iwidth);
  myfit->SetParName(0,"C");
  myfit->SetParName(1,"Sigma");
  myfit->SetParName(2,"B");
  myfit->SetParName(3,"A");
  myfit->SetParName(4,"X0");
  myfit->SetParName(5,"A2");
  myfit->SetParName(6,"x02");
  myfit->SetParName(7,"Sigma2");
  myfit->SetParameter(0, 100);
  myfit->SetParLimits(0, 10, 2000.);
  myfit->SetParameter(1, 10);
  myfit->SetParLimits(1, 5., 25.);
  myfit->SetParameter(2, 200);
  myfit->SetParameter(3, 0);
  myfit->SetParameter(4, 0);
  myfit->SetParLimits(4,-10.,10.);
  myfit->SetParameter(5, 100);
  myfit->SetParLimits(5, 10, 4000.);
  myfit->SetParLimits(6,-2., 2.);
  myfit->SetParameter(7, 1.);
  myfit->SetParLimits(7,1.,3.);
  hwidth->Fit("myfit");
  //hwidth->Draw("without fit");
  TF1 *fit = hwidth->GetFunction("myfit");
  Double_t ch = fit->GetChisquare();
  *sigma = (float) fit->GetParameter(1);
  float sigma2 = (float) fit->GetParameter(7);
  float x0_2 = (float) fit->GetParameter(6);
  float x0 = (float) fit->GetParameter(4);
  Double_t Br = fit->GetParameter(0);
  Double_t Br2 = fit->GetParameter(5);
  Double_t B = fit->GetParameter(2);
  Double_t A = fit->GetParameter(3);
  //Double_t sigma = 5.; Double_t ch =0; Double_t Br=0;
  //printf("chi2 %f sigma %f br %f\n", ch, *sigma, Br);
  //if(itr==0&&sigma>4.) sigma=4;
  TText *t = new TText();
  t->SetTextFont(62);
  t->SetTextColor(36);
  t->SetTextSize(0.08);
  t->SetTextAlign(12);
  char text[100],text2[100];
  sprintf(text,"A %5.0f sigma %5.2f X0 %5.2f",Br, *sigma, x0);
  sprintf(text2,"A2%4.0f sigma2 %5.2f X02 %5.2f shift %5.3f pixels", 
	  Br2, sigma2, x0_2, *shift);
  float profmax = hwidth->GetMaximum();
  t->DrawText((double) -15., (double) 0.25*profmax,text);
  t->DrawText((double) -15., (double) 0.15*profmax,text2);
  cw->Modified();
  cw->Update();
  rc = (rc + x0_2);
  cur_angle = atan2(y2-y1,x2-x1);
  
  x1 = x1 + x0_2*cos(cur_angle+pibytwo);
  y1 = y1 + x0_2*sin(cur_angle+pibytwo);
  
  x2 = x2 + x0_2*cos(cur_angle+pibytwo);
  y2 = y2 + x0_2*sin(cur_angle+pibytwo);
  
  if(sigma2 > 2.*TMath::Abs(x0_2))
    {
      *sigma= sigma2;
      if (TMath::Abs(*shift) < 2.)
	{
	  rc_int = rc*2.;
	  fx = x1*2.;
	  lx = x2*2.;
	  fy = y1*2.;
	  ly = y2*2.;
	}
      return;
    }
  else if(TMath::Abs(*shift) != 0.)
    {
      *sigma= 5.;
      return;
    }
  //*shift>2.5) return;
  //else do correction on radius or line position)(if not circle but line)!!
  hwidth->Delete();
  //hwidth = new TH1F("hwidth","hw", 101,-50.,51.);
  hwidth = new TH1F("hwidth","hw", 41,-20.,21.);
  //*fmin = atan2((y1-yc_int),(x1-xc_int));
  //if(*fmin<0.) *fmin += TMath::Pi()*2.;
  //*fmax = atan2((y2-yc_int),(x2-xc_int));
  //if(*fmax<0.) *fmax += TMath::Pi()*2.;
  //if(TMath::Abs(*fmax - *fmin) > TMath::Pi()) 
  //  {
  //   *fmin += TMath::Pi()*2;
  //    float buf = *fmin;
  //    *fmin = *fmax;
  //    *fmax = buf;
  //  }
  //printf("fmax %7.2f fmin %7.2f\n", *fmax, *fmin);
  *shift += x0_2;
  find_widthofthecorridor(fmin, fmax, xc, yc, rc, 
			  x1, y1, x2, y2, ima, itr, sigma, shift);
  //sigma = TMath::Min(sigma, sigma2);
  //return sigma;
}
void find_br_for_circlefit(TProfile *pro,
			   float fmin, float fmax, float xc, float yc, float rc, 
			   float rc0,
			   float x1, float y1, float x2, float y2, int ima, 
			   int iwidth)
{
  int i, j;
  float bright;
  double delta_angle;
  double range = fmax-fmin;
  int curx = 0; int cury = 0;
  double dr;
  int ip; int ii;
  int np = 0;
  double pibytwo = TMath::Pi()/2.;
  float True_bright;
  float XM, YM;
  float Mb;
  int ifit;
  
  //prevent cicling in stright line tracks
  dr = 0.1*sqrt(2.); //step to move along r
  double cur_angle;
  ifit = 0;
  delta_angle = range/TMath::Pi()*180./100./sqrt(pow(x2-x1,2) + pow(y2-y1,2));
  //printf("delta_angle %f\n",delta_angle);
  if(delta_angle > 0.00001) 
	 {
      cur_angle = fmin;
      int binxprev_circle = 0; int binyprev_circle = 0;
      int epx = (int) (xc + rc*cos(cur_angle+range));
      int epy = (int) (yc + rc*sin(cur_angle+range));
      //printf("epx %d epy %d, endangle %f cur %f\n", epx, epy, fmax, cur_angle);
      for(ii = 0; ii<1000000; ii++)
	{
	  ip = 0; 
	  int binxprev = 0;
	  int binyprev = 0;
	  
	  curx = (int) (xc + rc*cos(cur_angle));
	  cury = (int) (yc + rc*sin(cur_angle));
	  Mb=0;
	  if (binxprev_circle != curx || binyprev_circle != cury)
	    {
	      //store this pont 
	      x_mark[np] = curx;
	      y_mark[np] = cury;
	      i = (int) (vsize - cury - 1);
	      j = (int) (curx);
	      int bin = hprof->FindBin((Axis_t) (cur_angle-fmin)*rc0);
	      if(i>=0&&j>=0&&i<vsize&&j<hsize)
		{
		  bright = pCDBuffer[ima][i][j];
		  True_bright = bright - (p0[bin] + p1[bin]*rc0);
		  hslice->Fill((Axis_t) ip, (Stat_t)(True_bright/(fmax-fmin)/(rc+ip)));
		  if(True_bright > Mb) 
		    {
		      Mb = True_bright; XM = 2.*x_mark[np]; YM = 2.*y_mark[np];
		    }
		  if(find_the_same_point(curx,cury,np)==0)
		    {
		      Xcircle[np] = 2.*curx;
		      Ycircle[np] = 2.*cury;
		      Wpoint[np] = True_bright;
		      np++;
		    }
		}
	      binxprev_circle = curx; binyprev_circle = cury;
	      binxprev = curx; binyprev = cury;
	      ip++;
	      float rcur = rc + dr;
	      do {
		curx = (int) (xc + rcur*cos(cur_angle));
		cury = (int) (yc + rcur*sin(cur_angle));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      True_bright = bright - (p0[bin] +p1[bin]*rcur);
		      hslice->Fill((Axis_t) ip, (Stat_t)(True_bright/(fmax-fmin)/(rc+ip)));
		      if(True_bright >= Mb) 
			{
			  if(True_bright == Mb) 
			    {
			      XM = (XM + 2*curx)/2.; YM = (YM + 2*cury)/2.;
			    }
			  else
			    {
			      Mb = True_bright; XM = 2.*curx; YM = 2.*cury;
			    }
			}
		      if(find_the_same_point(curx,cury,np)==0)
			{
			  Xcircle[np] = 2.*x_mark[np];
			  Ycircle[np] = 2.*y_mark[np];
			  Wpoint[np] = True_bright;
			  np++;
			}
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur += dr;
	      } while (ip < iwidth+1);
	      ip = 1;
	      rcur = rc - dr;
	      binxprev = binxprev_circle; 
	      binyprev = binyprev_circle;
	      do {
		curx = (int) (xc + rcur*cos(cur_angle));
		cury = (int) (yc + rcur*sin(cur_angle));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      True_bright = bright - (p0[bin] +p1[bin]*rcur);
		      hslice->Fill((Axis_t) -ip, (Stat_t)(True_bright/(fmax-fmin)/(rc-ip)));
		      if(True_bright >= Mb) 
			{
			  if(True_bright == Mb) 
			    {
			      XM = (XM + 2*curx)/2.; YM = (YM + 2*cury)/2.;
			    }
			  else
			    {
			      Mb = True_bright; XM = 2.*curx; YM = 2.*cury;
			    }
			}
		      if(find_the_same_point(curx,cury,np)==0)
			{
			  Xcircle[np] = 2.*x_mark[np];
			  Ycircle[np] = 2.*y_mark[np];
			  Wpoint[np] = True_bright;
			  np++;
			}
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur -= dr;
	      } while(ip < iwidth+1);
	      XFit[ifit] = XM;
	      YFit[ifit] = YM;
	      WFit[ifit] = Mb;
	      ifit++;
	    }
	  cur_angle += delta_angle;
	  if(curx == epx && cury == epy) break;
	  if(cur_angle > (double) fmax) break;
	}
    }
  else
    {
      printf("There is a line, no circle!\n");
      double delta_l = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
      //printf("delta_l %f\n",delta_l);
      cur_angle = atan2(y2-y1,x2-x1);
      
      int binxprev_line = 0; int binyprev_line = 0;
      i = 0;
      double dl = 0.;
      float epx = x2;
      float epy = y2;
      //printf("epx %f epy %f, cur_angle %f\n", epx, epy, cur_angle);
      for(ii = 0; ii<1000000; ii++)
	{
	  ip = 0; 
	  int binxprev = 0;
	  int binyprev = 0;
	  
	  curx = (int) (x1 + dl*cos(cur_angle));
	  cury = (int) (y1 + dl*sin(cur_angle));
	  int bin = hprof->FindBin((Axis_t) dl);
	  if (binxprev_line != curx || binyprev_line != cury)
	    {
	      //store this pont 
	      x_mark[np] = curx;
	      y_mark[np] = cury;
	      
	      i = (int) (vsize - cury - 1);
	      j = (int) (curx);
	      if(i>=0&&j>=0&&i<vsize&&j<hsize)
		{
		  bright = pCDBuffer[ima][i][j];
		  True_bright = bright - (p0[bin] +p1[bin]*ip);
		  hslice->Fill((Axis_t) ip, (Stat_t) (True_bright/delta_l));
		  if(True_bright > Mb) {Mb = True_bright; XM = 2.*x_mark[np]; YM = 2.*y_mark[np];}
		  if(find_the_same_point(curx,cury,np)==0)
		    {
		      Xcircle[np] = 2.*x_mark[np];
		      Ycircle[np] = 2.*y_mark[np];
		      Wpoint[np] = True_bright;
		      np++;
		    }
		}
	      binxprev_line = curx; binyprev_line = cury;
	      binxprev = curx; binyprev = cury;
	      ip++;
	      float rcur = dr;
	      do {
		curx = (int) (x1 + dl*cos(cur_angle) + rcur*cos(cur_angle+pibytwo));
		cury = (int) (y1 + dl*sin(cur_angle) + rcur*sin(cur_angle+pibytwo));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      True_bright = bright - (p0[bin] +p1[bin]*ip);
		      hslice->Fill((Axis_t) ip, (Stat_t) (True_bright/delta_l));
		      if(True_bright > Mb) {Mb = True_bright; XM = 2.*x_mark[np]; YM = 2.*y_mark[np];}
		      if(find_the_same_point(curx,cury,np)==0)
			{
			  Xcircle[np] = 2.*x_mark[np];
			  Ycircle[np] = 2.*y_mark[np];
			  Wpoint[np] = True_bright;
			  np++;
			}
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur += dr;
	      } while (ip < iwidth+1);
	      ip = 1;
	      rcur = -dr;
	      binxprev = binxprev_line; 
	      binyprev = binyprev_line;
	      do {
		curx = (int) (x1 + dl*cos(cur_angle) + rcur*cos(cur_angle+pibytwo));
		cury = (int) (y1 + dl*sin(cur_angle) + rcur*sin(cur_angle+pibytwo));
		if(curx != binxprev || cury != binyprev) {
		  x_mark[np] = curx;
		  y_mark[np] = cury;
		  i = (int) (vsize - cury - 1);
		  j = (int) (curx);
		  if(i>=0&&j>=0&&i<vsize&&j<hsize)
		    {
		      bright = pCDBuffer[ima][i][j];
		      True_bright = bright - (p0[bin] +p1[bin]*ip);
		      hslice->Fill((Axis_t) -ip, (Stat_t) (True_bright/delta_l));
		      if(True_bright > Mb) {Mb = True_bright; XM = 2.*x_mark[np]; YM = 2.*y_mark[np];}
		      if(find_the_same_point(curx,cury,np)==0)
			{
			  Xcircle[np] = 2.*x_mark[np];
			  Ycircle[np] = 2.*y_mark[np];
			  Wpoint[np] = True_bright;
			  np++;
			}
		    }
		  binxprev = curx; binyprev = cury;
		  ip++;
		}
		rcur -= dr;
	      } while(ip < iwidth+1);
	    }
	  dl += 0.1*sqrt(2.); //step to move along r
	  if(curx == epx && cury == epy) break;
	  if(dl >= delta_l) break;
	}
      XFit[ifit] = XM;
      YFit[ifit] = YM;
      WFit[ifit] = Mb;
      ifit++;
    }
  hslice->Draw();
  Npoints_fit = np;
  NFit = ifit;
}
