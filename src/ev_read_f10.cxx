//  Reading of images from data file of image, filling data plot.

#include "TPostScript.h"
#include <dirent.h>
#include "_ROOT.h"
#include "_STD_UNIX.h"
#include "TPolyMarker.h"

#include "TColor.h"
//#include "TQObject.h"
#include "RQ_OBJECT.h"
#include "vertex.h"

#include <cfortran.h>
#include <minuit.h>

#include "EnergySelection.h"

typedef unsigned* PUINT;
typedef PUINT* PPUINT;

PPUINT pCDBuffer[2];

//void exec_event(Int_t event, Int_t x, Int_t y, TObject *selected);
void exec_event(int,int ,int , TObject *);
void usage(char* argv[]);
int end_input();

PPUINT Alloc_Buffer(int hsize, int vsize);
void Free_Buffer(PPUINT buffer[2], int vsize);
   
void GetNewEvent();
void RemoveCreated();
Bool_t CheckOutFile();
void SetPart(void);
void Draw_region(int xpix, int ypix, int w);
    
void Smooth(int image, int ipass);
void LifeFirst(int image, int iter, int ipass);
void LifeSecond(int image, int iter, int ipass);
void SortPointsByBrightness(int image);
void indexx(unsigned long n, int arr[], unsigned long indx[]);

void MakeMainSegments(int image, int ipass);
void SelectTracks(int version, int image, int ipass);
float D_to_Line(int xp, int yp, float Slope, float Y_Origin);
void EraseNearestToTrackPoints(int image, int Start, int End, int ipass, int what);

void FindVertex(int image, int ipass);
void LookForWeakTracks(int image, int ipass);
void XY_Vertex_Net(int NTracks, float *a, float *b, float *r, 
		   float *Xvertex, float *Yvertex, 
		   float *Dvertex);
void FindPairsOnImages(void);
void RebuildTracks(int);
int  FindFirstPion(void);
void FindNewVertex(int image);
void DigitizeEvent(void);
int AddPointsToTrack(int, int, int, int);
void ProlongtoVertex(int image);
void JoinSegTracks(int image);
int SelectPointsOnTracks(int image, int N, int *xcbuf, int *ycbuf);
int DeleteWeakPoints(int image, int ibr);
double FindRcrit(int id);

/*
void find_neib(int image, int ireg, int i, int j, int ip);
void find_neib2(int image, int ireg, int i, int j, int i, int j);
void find_neib3(int image, int ireg, int i, int j, int i, int j, int columnmax, int rowmax,int n, int ncrit, int br_temp);
*/

void circle_minuit (int N);
void circle_minuit_fixR (int N);
void two_circles_minuit (int N);
void N_circles_fit(int N);
void N_lines_fit(int N);
void circle_analytical(int j);

extern const char* ReconMainSign();
extern void recon_main(char [], int);

  //double DistMax(float *x, float *y);
#define OUT_filename "dubtoev.dat"
FILE *out_file;   /* output file handler */
FILE *in_file;   /* output file handler */

TCanvas *c2,*ch,*c3, *c4, *cp, *cw, *csl;
TPad *pleft, *pright, *ptmp, *pMsg;
TPad *pbleft, *pbright;
TPad *plc3, *prc3, *pMc3;
TPad *plc4, *prc4, *pMc4;
TPad *cp0, *cp1, *hp0, *cw0, *csl0;
TProfile *hprof, *hprof2, *hprof3, *hprof4;

TH1F *htest,*hbr;
TH1F *hwidth, *hslice;
TH1F *bleft, *bright;
TH2F *left, *right;
//^egcs  gzFile *im_right, *im_left;
gzFile im_right, im_left;
TGFileInfo fi;
TPaveLabel *pl, *pl1;
TGraph *gr;
TPaveLabel *pl11, *pl2, *ptext;
TMarker *mark, *mark2;
TObjectTable *fObjTbl;
TPolyMarker *markp, *markp2;
TArc *myarc;
TMarker *marke;
Float_t *px_mark, *py_mark;

TStopwatch timer;

Double_t xmax_br_l=4096,xmin_br_l=0,xmax_br_r=4096,xmin_br_r=0;
int hsize, vsize;
int cur_image;
int start_image = 0;
int event_kind=1, num_tracks=1, ipoint[2]={0,0}, itrack[2];

const char *filetypes[] = { 
  "Left images",    "[L,l]*",
  "left images",    "l*",
  "All files",     "*",
  0,               0 
};

//char iniDir[100] = {"/data/koval/Dubto/automaton/test_data/"};
//char iniDir[100] = {"/data/koval/DubTo/automaton/test_data/"};
//char iniDir[100] = {"/data2/koval/DubTo/automaton/data/"};
char iniDir[100] = {"/data/koval/DubTo/automaton_2003/data/"};
//char iniDir[100] = {"/data/dubto/data99_1/cdrom/nov8/"};
//char iniDir[100] = {"/mnt/cdrom/nov12"};
//char iniDir[100] = {"/data/dubto/data99_1/cdrom/nov9"};
//char iniDir[100] = {"/home/gilpont/DubTo/Vera/031100/nov14"};
int minl=0, maxl=4096, minr=0, maxr=4096;
int xlmin=108, xlmax=500, ylmin=68, ylmax=443;
int xrmin=128, xrmax=520, yrmin=52, yrmax=414;
char WorkDir[100];
char eventname[200];
char ofilename[100], gifname[100];
char test_outfile[100];
char *tokenPtr;

//block for segments & tracks
#define NAT_LOG 2.718281828
#define  MAXNUMBEROFPOINTS 10000
#define POINTSLIMIT 350.
#define RADIUSLIMIT 30.
#define DIST_CRIT  15.
#define BIG 10000000000.
#define NEIB 10
int tailx[2][2], taily[2][2], iw;
float rad_tail[2][3][4];
short int sostn[8];
short int sosts[658][517];
int rsost[658][517];
short int flag[517][658];
float aver[2], main_aver[2];

float avermax;
int xrest[2][MAXNUMBEROFPOINTS], yrest[2][MAXNUMBEROFPOINTS];
int xmain[2][MAXNUMBEROFPOINTS], ymain[2][MAXNUMBEROFPOINTS];
int Npoint[2], Nreg[2], Npoint_pass, Npoint_main, Not_Used;
float eff;
#define MAXNUMBEROFSEGMENTS 50000
#define MAXNP 2000
#define MAXNUMBEROFTRACKS 40
int NumberOfSegments;
int StartOfSegment[MAXNUMBEROFSEGMENTS], EndOfSegment[MAXNUMBEROFSEGMENTS];
double Dist[MAXNUMBEROFSEGMENTS],  CosTheta[1000];
double cos_a[MAXNUMBEROFSEGMENTS],  cos_b[MAXNUMBEROFSEGMENTS];
//float Br1[MAXNUMBEROFSEGMENTS], Br2[MAXNUMBEROFSEGMENTS];
double MaxCosPhi;
int NumberOfTracks, Level, NumberOfTrackPoints;
int LengthOfBestTrack[MAXNUMBEROFTRACKS];
int CurrentTrack[1000];
int BestTrack[MAXNUMBEROFTRACKS][1000];
int LevelOfSegment[MAXNUMBEROFSEGMENTS];
double CosPhiMax;// = 0.9948;//0.9948; //10 degrees
//const double SinPhiMax = 0.0;/1745;//1 degrees
short int TrackPoint[MAXNUMBEROFPOINTS];
int X[2][MAXNUMBEROFTRACKS][MAXNUMBEROFPOINTS], Y[2][MAXNUMBEROFTRACKS][MAXNUMBEROFPOINTS];
int ntr[2], np[2][MAXNUMBEROFTRACKS], Br_tr[2][MAXNUMBEROFTRACKS];
float  Av_Br[2][MAXNUMBEROFTRACKS];
int longtrack[2][MAXNUMBEROFTRACKS];
int fp_longtrack[2][MAXNUMBEROFTRACKS], lp_longtrack[2][MAXNUMBEROFTRACKS];
int mp_longtrack[2][MAXNUMBEROFTRACKS];
float Ltr[2][MAXNUMBEROFTRACKS], Dens[2][MAXNUMBEROFTRACKS];
float elx[10000], ely[10000];

//float xcD[2][MAXNUMBEROFTRACKS], ycD[2][MAXNUMBEROFTRACKS], rcD[2][MAXNUMBEROFTRACKS];
//float chi2_circle[2][MAXNUMBEROFTRACKS];
float X_Vert[2], Y_Vert[2], Min_Vertex_Error[2];
//int mark[2][MAXNUMBEROFTRACKS];

  int xmask[2][10][10], ymask[2][10][10];
  float xc_int, yc_int, rc_int;
  float stx, sty, endx, endy;
  float brrr[2][6];
  float xcent_vx[10], ycent_vx[10], rcent_vx[10];
  float Xvertax, Yvertax, Dvertax;
int nomer_tr[2][MAXNUMBEROFTRACKS];
//float angle[2][MAXNUMBEROFTRACKS];
//float Angle[2][MAXNUMBEROFTRACKS];
int GoodEvent, NGoodTracks[2];

float xtrack[2][6][10000], ytrack[2][6][10000];
float xt[2][6][10000], yt[2][6][10000];
float wt[2][6][10000], xp[10000], yp[10000];;
int   particle[10]={0,0,0,0,0,0,0,0,0,0};

float Xcircle[50000], Ycircle[50000], Zcircle[50000], Wpoint[50000];
float XC, YC, RC, AC, BC, chi2, DRC, DAC, phif, t[50000], bufmax, bufmin;
short int Signe;
//float XC, YC, RC, chi2, XC1, YC1, RC1;
float XCN[6], YCN[6], RCN[6];

int NTracks[2], pair[MAXNUMBEROFTRACKS], Track[2][MAXNUMBEROFTRACKS];
int PionTrack[2];
float radiuslimit, dist_crit1, dist_crit2;

int nneib[10000], brneib[10000];
int dlinamax, xmax, ymax, square;
float meanneib;
int maxneib;
float sumbr;
float meangrad, maxgrad;
int NumberJoinedTracks;
// globals for MINUIT
int Npoints_fit, ncircles;
float fcn_print;
int error_flag, is;
Float_t x_mark[50000], y_mark[50000];
short int not_join[MAXNUMBEROFTRACKS][MAXNUMBEROFTRACKS];

#define POINTS 3
//float xe[2][5000], ye[2][5000], we[2][5000];
//float elx[5000], ely[5000];
//int ne[2];
//int niter;
#include "TPolyLine.h"
//#include "tsp_my.cxx"
//TSP *eleft, *eright;
TPolyLine *poly1, *poly2;




TROOT simple("simple","Test of histogramming and I/O");

class MyChoiseFrame : public TGMainFrame {

  TGGroupFrame    *fG1, *fG2;
  TGLabel         *fLabel;
  TGTextButton    *fButton1, *fButton2;
  TGRadioButton   *fRBut1, *fRBut2, *fRBut3, *fRBut4, *fRBut5, *fRBut6;
  TGRadioButton   *fRBut21, *fRBut22, *fRBut23, *fRBut24, *fRBut25, *fRBut26;
  TGLayoutHints   *fLayout, *fLayout1;

public:
  MyChoiseFrame(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~MyChoiseFrame()  ;
  virtual void DisableWindow();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
  virtual Bool_t HandleKey(Event_t*);
};

MyChoiseFrame::~MyChoiseFrame()
{
  delete fG1; delete fG2;
  delete fLabel;
  delete fLayout; delete fLayout1;
  delete fButton1; delete fButton2;
  delete fRBut1; delete fRBut2; delete fRBut3; 
  delete fRBut4; delete fRBut5; delete fRBut6;
  delete fRBut21; delete fRBut22; delete fRBut23; delete fRBut24;
  delete fRBut25; delete fRBut26; 
}

void MyChoiseFrame::DisableWindow()
{
  // Called when window is closed via the window manager.
  //   delete fButton1; delete fButton2;
  //   delete fRBut1; delete fRBut2; delete fRBut3; delete fRBut4; delete fRBut5;
  //   delete fLayout;
     delete this;
   //   TGMainFrame::CloseWindow();
  //  gApplication->Terminate(0);
}

MyChoiseFrame::MyChoiseFrame(const TGWindow *p, UInt_t w, UInt_t h)
  : TGMainFrame(p, w, h)
{
  //------------------------------
  // Create new graphics context for drawing of TGlabel and TGButtons
  // with text in red.
  FontStruct_t labelfont;
  labelfont = fClient->GetFontByName(gEnv->GetValue("Gui.BoldFont","-adobe-helvetica-bold-r-*-*-14-*-*-*-*-*-iso8859-1"));
  // Define new graphics context. Only the fields specified in
  // fMask will be used (for default TGLabel context see
  // http://root.cern.ch/root/html/src/TGClient.cxx.html).
  GCValues_t   gval;
  GContext_t            fRedTextGC;
  gval.fMask = kGCForeground | kGCFont;
  gval.fFont = gVirtualX->GetFontHandle(labelfont);
  fClient->GetColorByName("red", gval.fForeground);
  
  fRedTextGC = gVirtualX->CreateGC(fClient->GetRoot()->GetId(), &gval);
  //________________________________________________________________
  fG1 = new TGGroupFrame(this, new TGString("Tracks"),kVerticalFrame);
  fG2 = new TGGroupFrame(this, new TGString("Kind of event"),kVerticalFrame);
  
  // Create a main frame with a number of different buttons.
  fButton1 = new TGTextButton(this, "Info", 7);
  // Change background of fButton to yellow
  ULong_t yellow;
  fClient->GetColorByName("yellow", yellow);
  fButton1->ChangeBackground(yellow);

  fButton2 = new TGTextButton(this, "Open new event", 8);
   // Change background of fButton to red
  ULong_t red;
  fClient->GetColorByName("green", red);
  fButton2->ChangeBackground(red);


  fLabel = new TGLabel(this,"Set number of tracks",fRedTextGC,labelfont);
  fRBut1 = new TGRadioButton(fG1, "One track   ", 1);
  fRBut2 = new TGRadioButton(fG1, "Two tracks  ", 2);
  fRBut3 = new TGRadioButton(fG1, "Three tracks", 3);
  fRBut4 = new TGRadioButton(fG1, "Four tracks ", 4);
  fRBut5 = new TGRadioButton(fG1, "Five tracks ", 5);
  fRBut6 = new TGRadioButton(fG1, "Six tracks  ", 6);
  fLayout1 = new TGLayoutHints(kLHintsTop | kLHintsExpandX,
			       2, 2, 3, 0);
  fLayout = new TGLayoutHints(kLHintsLeft|kLHintsTop,
			      6, 2, 5, 0);

  fRBut21 = new TGRadioButton(fG2, "1 --- elastic", 21);
  fRBut22 = new TGRadioButton(fG2, "2 --- Pi0+p+He3", 22);
  fRBut23 = new TGRadioButton(fG2, "3 --- p+p+p+n", 23);
  fRBut24 = new TGRadioButton(fG2, "4 --- without origin Pi+", 24);
  fRBut25 = new TGRadioButton(fG2, "5 --- ", 25);
  fRBut26 = new TGRadioButton(fG2, "6 --- ", 26);

  AddFrame(fLabel,fLayout1);
  AddFrame(fButton1, fLayout1);
  AddFrame(fButton2, fLayout1);

  AddFrame(fG1,fLayout1);
  AddFrame(fG2,fLayout1);

  fG1->AddFrame(fRBut1, fLayout);
  fG1->AddFrame(fRBut2, fLayout);
  fG1->AddFrame(fRBut3, fLayout);
  fG1->AddFrame(fRBut4, fLayout);
  fG1->AddFrame(fRBut5, fLayout);
  fG1->AddFrame(fRBut6, fLayout);

  fG2->AddFrame(fRBut21, fLayout1);
  fG2->AddFrame(fRBut22, fLayout1);
  fG2->AddFrame(fRBut23, fLayout1);
  fG2->AddFrame(fRBut24, fLayout1);
  fG2->AddFrame(fRBut25, fLayout1);
  fG2->AddFrame(fRBut26, fLayout1);

  fRBut1->Associate(this);
  fRBut2->Associate(this);
  fRBut3->Associate(this);
  fRBut4->Associate(this);
  fRBut5->Associate(this);
  fRBut6->Associate(this);

  fRBut21->Associate(this);
  fRBut22->Associate(this);
  fRBut23->Associate(this);
  fRBut24->Associate(this);
  fRBut25->Associate(this);
  fRBut26->Associate(this);

  MapSubwindows();

  Layout();

  SetWindowName("How many tracks?");
  SetIconName("How many tracks?");
  SetWMPosition(805,350);
  fRBut1->SetState(kButtonDown);
  fRBut21->SetState(kButtonDown);
  MapWindow();
}

MyChoiseFrame *maiWin;


Bool_t MyChoiseFrame::HandleKey(Event_t *event)
{
  printf("Choise's Key pressed\n");
  return kTRUE;
}

Bool_t MyChoiseFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{

  int retval;
  //*icons:kMBIconStop,kMBIconQuestion,kMBIconExclamation,kMBIconAsterisk*/
  //*{ kMBYes, kMBNo, kMBOk, kMBApply,
  //kMBRetry, kMBIgnore, kMBCancel,
  //  kMBClose, kMBDismiss }; */
  EMsgBoxIcon icontype = kMBIconAsterisk;  
  // Process events generated by the buttons in the frame.
  //  cout << "mess = "<<GET_MSG(msg)<<endl;
  //  cout << "submessage = "<<GET_SUBMSG(msg)<<endl;

  switch (GET_MSG(msg)) {
  case kC_COMMAND:
    switch (GET_SUBMSG(msg)) {
    case kCM_BUTTON:
      //      printf("text button id %ld pressed\n", parm1);
      if (parm1 == 7)
	new TGMsgBox(fClient->GetRoot(), this,
		     "Info",  "You have to put 3 points per each track",
		     icontype, kMBOk, &retval,kVerticalFrame);
      if (parm1 == 8)
	{
	  maiWin->DisableWindow();
	  RemoveCreated();
	  GetNewEvent();
	}
      break;
    case kCM_RADIOBUTTON:
      if ( parm1 >= 1 && parm1 <= 6)
	{
	  if (parm1 == 1)
	    {
	      fRBut2->SetState(kButtonUp);
	      fRBut3->SetState(kButtonUp);
	      fRBut4->SetState(kButtonUp);
	      fRBut5->SetState(kButtonUp);
	      fRBut6->SetState(kButtonUp);
	    }
	  if (parm1 == 2)
	    {
	      fRBut1->SetState(kButtonUp);
	      fRBut3->SetState(kButtonUp);
	      fRBut4->SetState(kButtonUp);
	      fRBut5->SetState(kButtonUp);
	      fRBut6->SetState(kButtonUp);
	    }
	  if (parm1 == 3)
	    {
	      fRBut1->SetState(kButtonUp);
	      fRBut2->SetState(kButtonUp);
	      fRBut4->SetState(kButtonUp);
	      fRBut5->SetState(kButtonUp);
	      fRBut6->SetState(kButtonUp);
	    }
	  if (parm1 == 4)
	    {
	      fRBut1->SetState(kButtonUp);
	      fRBut2->SetState(kButtonUp);
	      fRBut3->SetState(kButtonUp);
	      fRBut5->SetState(kButtonUp);
	      fRBut6->SetState(kButtonUp);
	    }
	  if (parm1 == 5)
	    {
	      fRBut1->SetState(kButtonUp);
	      fRBut2->SetState(kButtonUp);
	      fRBut3->SetState(kButtonUp);
	      fRBut4->SetState(kButtonUp);
	      fRBut6->SetState(kButtonUp);
	    }
	  if (parm1 == 6)
	    {
	      fRBut1->SetState(kButtonUp);
	      fRBut2->SetState(kButtonUp);
	      fRBut3->SetState(kButtonUp);
	      fRBut4->SetState(kButtonUp);
	      fRBut5->SetState(kButtonUp);
	    }
	  num_tracks = parm1;
	  //	  printf("number of tracks in event %ld \n",parm1);
	}
      if (parm1 >=21 && parm1<=24)
	{
	  if (parm1 == 21)
	    {
	      fRBut22->SetState(kButtonUp);
	      fRBut23->SetState(kButtonUp);
	      fRBut24->SetState(kButtonUp);
	      fRBut25->SetState(kButtonUp);
	      fRBut26->SetState(kButtonUp);
	    }	  
	  if (parm1 == 22)
	    {
	      fRBut21->SetState(kButtonUp);
	      fRBut23->SetState(kButtonUp);
	      fRBut24->SetState(kButtonUp);
	      fRBut25->SetState(kButtonUp);
	      fRBut26->SetState(kButtonUp);
	    }	  
	  if (parm1 == 23)
	    {
	      fRBut21->SetState(kButtonUp);
	      fRBut22->SetState(kButtonUp);
	      fRBut24->SetState(kButtonUp);
	      fRBut25->SetState(kButtonUp);
	      fRBut26->SetState(kButtonUp);
	    }	  
	  if (parm1 == 24)
	    {
	      fRBut21->SetState(kButtonUp);
	      fRBut22->SetState(kButtonUp);
	      fRBut23->SetState(kButtonUp);
	      fRBut25->SetState(kButtonUp);
	      fRBut26->SetState(kButtonUp);
	    }	  
	  if (parm1 == 25)
	    {
	      fRBut21->SetState(kButtonUp);
	      fRBut22->SetState(kButtonUp);
	      fRBut23->SetState(kButtonUp);
	      fRBut24->SetState(kButtonUp);
	      fRBut26->SetState(kButtonUp);
	    }	  
	  if (parm1 == 26)
	    {
	      fRBut21->SetState(kButtonUp);
	      fRBut22->SetState(kButtonUp);
	      fRBut23->SetState(kButtonUp);
	      fRBut24->SetState(kButtonUp);
	      fRBut25->SetState(kButtonUp);
	    }	  
	  event_kind = parm1 - 20;
	}
      //      printf("radio button id %ld pressed\n", parm1);

      break;
    default:
      break;
    }
  default:
    break;
  }
  return kTRUE;
}

class WaitFrame : public TGMainFrame {
  TGLabel         *fLabel;
  TGLayoutHints   *fLayout1;

public:
  WaitFrame(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~WaitFrame();
  virtual void CloseWaitFrame();
};

WaitFrame::~WaitFrame()
{
  delete fLabel;
  delete fLayout1;
}

void WaitFrame::CloseWaitFrame()
{
  delete this;
}

WaitFrame::WaitFrame(const TGWindow *p, UInt_t w, UInt_t h)
  :  TGMainFrame(p, w, h)
{
  //------------------------------
  // Create new graphics context for drawing of TGlabel and TGButtons
  // with text in red.
  FontStruct_t labelfont;
  labelfont = fClient->GetFontByName(gEnv->GetValue("Gui.ProportionalFont","-adobe-courier-bold-r-*-*-16-*-*-*-*-*-iso8859-1"));
  // Define new graphics context. Only the fields specified in
  // fMask will be used (for default TGLabel context see
  // http://root.cern.ch/root/html/src/TGClient.cxx.html).
  GCValues_t   gval;
  GContext_t            fRedTextGC;
  gval.fMask = kGCForeground | kGCFont;
  gval.fFont = gVirtualX->GetFontHandle(labelfont);
  fClient->GetColorByName("red", gval.fForeground);
  
  fRedTextGC = gVirtualX->CreateGC(fClient->GetRoot()->GetId(), &gval);
  //________________________________________________________________

  fLabel = new TGLabel(this,"Wait a while",fRedTextGC,labelfont);
  fLayout1 = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY);

  AddFrame(fLabel,fLayout1);

  MapSubwindows();

  Layout();

  SetWindowName("Wait!!");
  //  SetIconName("Wait!!");
  SetWMPosition(300,300);
  MapWindow();
}
// set particle type
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
class MyTestDialog : public TGTransientFrame {

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
   MyTestDialog(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
               UInt_t options = kMainFrame | kVerticalFrame);
   virtual ~MyTestDialog();

  //virtual void CloseWindow();
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
};


MyTestDialog::MyTestDialog(const TGWindow *p, const TGWindow *main, UInt_t w,
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

   if(num_tracks >= 2) {
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
     
     //fR2[0]->SetState(kButtonDown);
   }
   
   if(num_tracks >= 3) {
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
     
     //fR3[5]->SetState(kButtonDown);
   }
   if(num_tracks >= 4) {
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
     
     //fR4[0]->SetState(kButtonDown);
   }
   if(num_tracks >= 5) {
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
   if(num_tracks >= 6) {
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

   SetWindowName("Set Particles for each track");
   SetWMPosition(400,600);

   MapWindow();
   gClient->WaitFor(this);    // otherwise canvas contextmenu does not work
}

MyTestDialog::~MyTestDialog()
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

/*
void MyTestDialog::CloseWindow()
{
   // Called when window is closed via the window manager.

   delete this;
}
*/
Bool_t MyTestDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   // Process messages coming from widgets associated with the dialog.

   //   char tmp[20];               // unused
   //   static int newtab = 0;      // unused
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

/* new routine of SetNewPalette - root version 3.02 */
void SetVeryNewPalette()
{
  //  TColor *color;
  const int ncolors = 551;
  int itmp;
  Float_t light;
  Int_t Apalet[551];

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
 
  for ( itmp=0; itmp<ncolors+1; itmp++)
    {
      //      TColor *color = gROOT->GetColor(itmp+1);
      light = (Float_t)itmp/(ncolors-1);
      TColor *color = new TColor(1230+itmp,light,light,light,"");  // ignore "unused" compilation warning !!!

      //      color->SetRGB(light,light,light);
      //      cout<<itmp<<"    =   "<<light<<endl;
      Apalet[itmp]=1230 + itmp + 1;
      //      Apalet[itmp]=0;
    }

  //  palett->SetPalette(ncolors+1,Apalet);
  gStyle->SetPalette(ncolors,Apalet);
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

Bool_t CheckOutFile()
{
  char ctmp[3];
  int itemp;

  if (!gSystem->AccessPathName((char *)ofilename))
    {
      int retval;
      new TGMsgBox(gClient->GetRoot(),gClient->GetRoot(),
		   "Output File Exists!!",
		   "Would you like to create another version for the same event?",
		   kMBIconQuestion,kMBYes|kMBNo,&retval);
      if (retval == kMBNo) 
	{
	  RemoveCreated();
	  GetNewEvent();
	  return 0;
	}
    }

  /*   searching for empty filename   */
  itemp = 0;
  strcpy(test_outfile,ofilename);
  while (!gSystem->AccessPathName((char *)test_outfile))
    {
      itemp +=1;
      sprintf(ctmp,"%u",itemp);
      strcpy(test_outfile,ofilename);
      strcat(test_outfile,".");
      strcat(test_outfile,ctmp);
    }
  /*  new filename */
  //  strcpy(ofilename,test_outfile);  
  printf("ofilename %s\n test_outfile %s\n", ofilename, test_outfile);
  if ((out_file = fopen(test_outfile,"w")) == NULL)
    { 
      fprintf(stderr,"\n Cannot open output file %s \n",test_outfile);
      exit(1);
    }
  fprintf(out_file,"%s\n",eventname);
  return 1;
}

void RemoveCreated()
{

  left->~TH2F();  right->~TH2F();

  if (pl != NULL)
    {
      pl->~TPaveLabel();
    }
  if (pl1 != NULL)
    {
      pl1->~TPaveLabel();
    }

  fObjTbl->Print();
  fObjTbl->Delete();
  /*
  c4->Delete();
  c3->Delete();
  cp->Delete();
  cw->Delete();
  csl->Delete();
  hwidth->Delete();
  hslice->Delete();
  hprof->Delete();
  hprof2->Delete();
  hprof3->Delete();
  hprof4->Delete();
  */
  event_kind = 1;
  num_tracks=1;
  start_image = 0;
  ipoint[1] = 0; ipoint[2] = 0;
 
  // fObjTbl->Print();
  gObjectTable->Print();
}

WaitFrame *WaitMsg;


void GetNewEvent()
{
  int length = 0;
  char lfile[1000], rfile[1000];
  char *cptmp1, *cptmp2;
  char ctmp[10000];
  char fout_file[1000];
  // int maxl=0, minl=4096, maxr=0, minr=4096;   // unused

  //^egcs  fi.fFileTypes = (char **)filetypes;
  fi.fFileTypes = filetypes;

  strcpy(WorkDir,gSystem->WorkingDirectory());

  strcpy(iniDir,WorkDir);
  strcat(iniDir,"/data");
  printf(" iniDir %s\n", iniDir);
  gSystem->ChangeDirectory(iniDir);
       /* print files in current directory in reverse order */
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen,&fi);

  timer.Reset();
  timer.Start(kTRUE);

  if (fi.fFilename == NULL) exit(0);
  WaitMsg = new WaitFrame(gClient->GetRoot(),300,70);


  strcpy(iniDir,gSystem->DirName(fi.fFilename));
  printf("selected file  - %s \n",fi.fFilename);

  strcpy(lfile,gSystem->BaseName(fi.fFilename));
  
  /*    forming whole filenames islower  */
  strcpy(rfile,lfile);
  if (isupper(rfile[0]))
	 rfile[0] = 'R';
  else
	 rfile[0] = 'r';
  //cout<< "lfile = "<<lfile<<endl;
  //cout<< "rfile = "<<rfile<<endl;
  
  cptmp1 = gSystem->ConcatFileName(iniDir,(char *)lfile);
  cptmp2 = gSystem->ConcatFileName(iniDir,(char *)rfile);
  
  //cout<< "clfile = "<<cptmp1<<endl;
  //cout<< "crfile = "<<cptmp2<<endl;
  
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
  
  int iitmp = 0; 
  cptmp1 = strcpy(fout_file,cptmp1);
  
  if ((cptmp2 = strrchr(fout_file,'.')) != NULL) // looking for "."
	 *cptmp2 = '\0';
  
  if ((cptmp2 = strrchr(fout_file,'/')) != NULL) 
	 iitmp = ((intptr_t)cptmp2 - (intptr_t)cptmp1)/sizeof(char) + 1;
  
  iitmp += 1;     /*  removing first char from name  */
  iitmp *=sizeof(char);
  
  strcpy(eventname,fout_file + iitmp);
  //fprintf(stderr,"\n event : %s \n",eventname);
  char dirout[10] = {"dig/"};
  strcpy(ofilename, dirout);
  int i;
  int k = 0;
  for (i = (int) strlen(dirout); 
		 i<(int) strlen(dirout) + (int)strlen(eventname); i++)
	 {ofilename[i] = eventname[k]; k++;}
  for (i= (int) strlen(dirout) + (int)strlen(eventname);i < 100;i++)
	 ofilename[i] = '\0';
  strcat(ofilename,".dig");
  
  //  fprintf(stderr,"\n event filename : %s \n",ofilename);
  
  gzgets(im_left, ctmp,10000);
  
  int ntmp;
  //istrstream Inbuffer(ctmp,10000);
  std::istringstream Inbuffer(ctmp);
  
  while ( Inbuffer >> ntmp )
	 length++;
  
  printf("Numbers of INTs = %d\n",length);
  gzrewind(im_left);
  
  int length1 = 0;
  while( gzgets(im_left, ctmp,10000) != NULL)
	 length1++;
  
  gzrewind(im_left);
  vsize = length1; hsize = length;
  printf("vsize = %d, hsize = %d\n",vsize,hsize);
  
  
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
  int image, j;
  //unsigned maxBright = 0;
  for ( image=0; image<2; image++ )
	 {
		itrack[image] = 1;
		ipoint[image] = 0;
		for ( i=0; i<vsize; i++ )
		  {
			 if (image == 0)
				gzgets(im_left,ctmp,10000);
			 else
				gzgets(im_right,ctmp,10000);
			 
			 //istrstream Inbuffer(ctmp,10000);
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
  
  TPaveLabel *a = 0;
  a = (TPaveLabel*)gROOT->FindObject("pbr"); if (a) a->Delete(); a = 0;  
  TH2F *c = 0;                                                               
  c = (TH2F*)gROOT->FindObject("left"); if (c) c->Delete(); c = 0;  
  left = new TH2F("left","bitmap for left image", 
						(Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);
  
  c = (TH2F*)gROOT->FindObject("right"); if (c) c->Delete(); c = 0;  
  right = new TH2F("right","bitmap for right image",
						 (Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);
  
  TH1F *bc = 0;
  bc = (TH1F*)gROOT->FindObject("bleft"); if (bc) bc->Delete(); bc = 0;
  bleft = new TH1F("bleft","left brightness",4096,-0.5,4095.5);
  bc = (TH1F*)gROOT->FindObject("bright"); if (bc) bc->Delete(); bc = 0;
  bright = new TH1F("bright","right brightness",4096,-0.5,4095.5);


  Axis_t xx,yy;
  Stat_t ww;
  //  Stat_t ww1;    // unused
  //  float sumw;    // unused
  //  int inum;    // unused
  for (image=0; image<2; image++)
    {
      for (i=0; i<vsize; i++)
	{
	  for (j=0; j<hsize; j++)
	    {
	      ww = (Stat_t) pCDBuffer[image][i][j];
	      xx = (Axis_t) j;//(j + 0.5);
	      yy = (Axis_t) (vsize - i - 1);//0.5);  
	      if (image == 0) 
		{
		  left->Fill(xx,yy,ww);
		  bleft->Fill((Float_t)ww);
		}
	      else
		{
		  right->Fill(xx,yy,ww);
		  bright->Fill((Float_t)ww);
		}
	    }
	}
    }
    WaitMsg->CloseWaitFrame();
  //  getc(stdin);
  /*
  free(pCDBuffer[0]);
  free(pCDBuffer[1]);
  */
	 //Free_Buffer(pCDBuffer,vsize);


  //  Int_t logz=1;
  //  pleft->SetLogz(logz);
  //  pright->SetLogz(logz);
  pleft->UseCurrentStyle();
  pright->UseCurrentStyle();
  //SetVeryNewPalette();

  left->SetContour(550);
  right->SetContour(550);
  pleft->cd();
  //pleft->SetLogz(1);
  left->Draw("colz");
  pright->cd();
  //  pright->SetLogz(1);
  right->Draw("colz");

  //left->SetMaximum(maxl);
  //left->SetMinimum(minl);
  //right->SetMaximum(maxr);
  //right->SetMinimum(minr);
  
  pMsg->cd();
  pl1 = new TPaveLabel(0.07,0.15,0.93,0.45,eventname,"br");
  pl1->SetFillColor(45);
  pl1->SetTextSize(0.5);
  pl1->Draw();

  pbleft->cd();
  bleft->Draw();
  pbright->cd();
  bright->Draw();

  c2->Update();

  //  cout << "before close"<<endl;

  maiWin = new MyChoiseFrame(gClient->GetRoot(), 160, 320);
  
}

int end_input()
{
  int sample = 2 * num_tracks * POINTS;
  int itmp = 0;
  for (int i=0; i<2; i++)   // loop over 
    itmp += ipoint[i];
  if ( itmp == sample )
    return -1;
  else
    return itmp;
}

void TH1::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  TH2F *htmp;
  

  Float_t xx,yy;
  float sumbright_all;
  int binxmax, binymax;
  float rmax;
  int bright;
  
  //  t0 = time(t0);
  gPad->SetCursor(kCross);
  
  switch (event)
    {
    case kButton1Double:
      {
	//	printf("Botton 1 double pressed\n");
	xx = gPad->AbsPixeltoX(px);  
	yy = gPad->AbsPixeltoY(py); 

	TAxis *xaxis = GetXaxis();
	TAxis *yaxis = GetYaxis();
	Int_t binx, biny;
	binx = xaxis->FindBin(xx);
	biny = yaxis->FindBin(yy);
	
	if (strcmp(this->GetName(), "bleft") == 0 )
	  {
	    Double_t xmin_br_l = xaxis->GetBinCenter(binx);
	    left->SetMinimum(xmin_br_l);
	    left->SetMaximum(xmax_br_l);
	    pleft->Modified();
	    pleft->Update();
	    return;
	  }
	if (strcmp(this->GetName(), "bright") == 0 )
	  {
	    Double_t xmin_br_r = xaxis->GetBinCenter(binx);
	    right->SetMinimum(xmin_br_r);
	    right->SetMaximum(xmax_br_r);
	    pright->Modified();
	    pright->Update();
	    return;
	  }
      }
      break;
    case kButton1Down:
      //      printf("Botton 1 pressed\n");      



      //      t0 = second();
      //      cout << "time 0 ="<<t0<<endl;
      
      xx = gPad->AbsPixeltoX(px);  
      yy = gPad->AbsPixeltoY(py); 
      //      cout << "coordinates : "<<"px="<<px<<"  py="<<py<<endl;
      //      cout << "XY coordinates ="<<xx<< "  "<<yy<<endl;
      //      cout << "XY coordinates ="<<gPad->PadtoX(xx)<< "  " 
      //	   <<gPad->PadtoY(yy)<<endl;
      //      cout << "XY coordinates ="<<gPad->GetEventX()<< "  " 
      //	   <<gPad->GetEventY()<<endl;
      {
	TAxis *xaxis = GetXaxis();
	TAxis *yaxis = GetYaxis();
	Int_t binx, biny;
	binx = xaxis->FindBin(xx);
	biny = yaxis->FindBin(yy);
	//	cout << "bins ...."<<binx<< "   "<<biny<<endl;
	Int_t bin_arr = this->GetBin(binx,biny);
	//	cout << "content ..."<< GetBinContent(bin_arr)<<endl;
	Int_t nxbins = this->GetNbinsX();
	Int_t nybins = this->GetNbinsY();
	Float_t widthx = (Float_t)hsize/nxbins;
	Float_t widthy = (Float_t)vsize/nybins;
	Int_t nside = 10;
	Int_t nnewbins = 2*nside + 1;

	/*	
	std::cout<<"object :"<<gPad->GetName()<<"  left ="<<
	  strcmp(this->GetName(), "bleft")<<"  right ="<<
	  strcmp(this->GetName(), "bright")<<std::endl;
	*/
	if (strcmp(this->GetName(), "bleft") == 0 )
	  {
	    Double_t xmax_br_l = xaxis->GetBinCenter(binx);
	    Double_t xmin_br_l = left->GetMinimum();
	    left->SetMaximum(xmax_br_l);
	    left->SetMinimum(xmin_br_l);
	    pleft->Modified();
	    pleft->Update();
	    return;
	  }
	if (strcmp(this->GetName(), "bright") == 0 )
	  {
	    Double_t xmax_br_r = xaxis->GetBinCenter(binx);
	    Double_t xmin_br_r = right->GetMinimum();
	    right->SetMaximum(xmax_br_r);
	    right->SetMinimum(xmin_br_r);
	    pright->Modified();
	    pright->Update();
	    return;
	  }
	 
	if (strcmp(this->GetName(), "htmp") != 0)
	  {
	    if (start_image == 0)
	      {
		//		TH2F *tmp;
		//		int retval;
		char msg[100];
		sprintf(msg,"%u tracks were choosen for the event\n",num_tracks);

		maiWin->DisableWindow();
		pMsg->cd();
		//		TAttText::SetTextFont(63);
		pl = new TPaveLabel(0.07,0.55,0.93,0.8,msg,"tl");
		pl->SetFillColor(44);
		pl->SetTextSize(0.5);
		pl->Draw();
		
		if (!CheckOutFile()) return;
		//SetPart();
	      }
	    start_image +=1;
	  }
	
	if (strcmp(this->GetName(),"left") == 0)
	  { cur_image = 0; pleft->cd(); }
	else if (strcmp(this->GetName(), "right") == 0)
	  { cur_image = 1; pright->cd(); }
	else
	  {
	    if (strcmp(this->GetName(), "htmp") == 0)
	      {
		//  Float_t xx1,xx2,yy1,yy2;    // unused
		timer.Reset();
		timer.Start(kTRUE);
		delete this;
		
		if (cur_image == 0)
		  pleft->cd();
		else if (cur_image == 1)
		  pright->cd();
		
		int iipoint;
		++ipoint[cur_image];
		iipoint = ipoint[cur_image] % POINTS;
		if ( iipoint == 0 )  iipoint = 3;
		itrack[cur_image] = (ipoint[cur_image]-1)/POINTS + 1;
		
		int itr;
		itr = itrack[cur_image] - 1 ;
		
		TMarker *mark = new TMarker(xx,yy,5);
		fObjTbl->Add(mark);
		mark->SetMarkerColor(46);
		//		cout << gPad->GetX1()<<endl;
		//		cout << gPad->GetX2()<<endl;
		//		cout << gPad->GetY1()<<endl;
		//		cout << gPad->GetY2()<<endl;
		mark->SetMarkerStyle(8);
		mark->SetMarkerSize(1.);
		mark->DrawMarker(xx,yy);
		//cout << "real time0 ="<<timer.RealTime()<<endl;
		//cout << "CPU time0 ="<<timer.CpuTime()<<endl;
		timer.Continue();
		//		mark->Draw("same");
		//		gPad->GetPadPar(xx1,yy1,xx2,yy2);
		//		gPad->SetPad(-0.12,-0.12,1.12,1.12);
		//		cout <<xx1<<yy1<<xx2<<yy2<<endl;
						
		if (cur_image == 0)
		  { 
		    pleft->Modified();
		    pleft->Update();
		  }
		else if (cur_image == 1)
		  {
		    pright->Modified();
		    pright->Update();
		  }
		timer.Stop();
		//cout << "real time ="<<timer.RealTime()<<endl;
		//cout << "CPU time ="<<timer.CpuTime()<<endl;
		
		{
		  /*
		    xaxis = GetXaxis();
		    yaxis = GetYaxis();
		    binx = xaxis->FindBin(xx);
		    biny = yaxis->FindBin(yy);*/
		  //		    cout <<"Bins..."<< xx<<" "<<yy<<endl;
		  //		    cout << ->GetName();
		  
		  //		    cout << "IMAGE = "<<cur_image<<endl;
		  //		    cout << "point "<<ipoint[cur_image]<<" "<<iipoint<<endl;
		  //		    cout << "track = "<<itrack[cur_image]<<endl;
		  
		  if (itrack[cur_image] <= num_tracks)
		    {
		      xmask[cur_image][itr][iipoint-1]=(int) xx;
		      ymask[cur_image][itr][iipoint-1]=(int) yy;
		      printf(" %d %d %d x %d y %d\n",cur_image,itr,iipoint, 
			     xmask[cur_image][itr][iipoint-1], 
			     ymask[cur_image][itr][iipoint-1]);
		      /*  disable of number of tracks panel */
		      
		      if (end_input() == -1)
			{
			  //			    cout << " F I N I S H " << endl;
			  //simple.ProcessLine(".q");
			  fprintf(out_file," %d      %d \n",event_kind,num_tracks);
			  int ima;
			  int track;
			  float aver_l, aver_r;
			  for ( ima = 0; ima<2; ima++)
			    {
			      int sum = 0;
			      
				//  printout to the output file 
			      //fprintf(out_file," %d \n",ima);
			      if ( ima == 0 )
				{ 
				  pleft->cd(); 
				  xaxis = left->GetXaxis();
				  yaxis = left->GetYaxis();
				  for (biny=100; biny<=400; biny++)
				    {
				      for (binx=100; binx<=500; binx++)
					{
					  bin_arr = left->GetBin(binx,biny);
					  sum += (int) left->GetBinContent(bin_arr);
					}
				    }
				  aver_l = (float) sum/120000.;
				}
			      if ( ima == 1 )
				{ 
				  pright->cd(); 
				  xaxis = right->GetXaxis();
				  yaxis = right->GetYaxis();
				  for (biny=100; biny<400; biny++)
				    {
				      for (binx=100; binx<500; binx++)
					{
					  bin_arr = right->GetBin(binx,biny);
					  sum += (int) right->GetBinContent(bin_arr);
					}
				    }
				  aver_r = (float) sum/120000.;
				}
			      //for (int mypass=0; mypass<2; mypass++)
			      for (int mypass=0; mypass<1; mypass++) {
				  // cycle for passing the vertex point in the second tour
				  for (track=0; track<num_tracks; track++) {
				      //if(mypass == 1)
				      {
				//  printout to the output file 
					fprintf(out_file," %d \n",track);
				      }
				      int i;
				      float buf;
				      float bufmin, bufmax;
				      int imin, imax;
				      float x_cur, y_cur;
				      
				      for(i=0; i<POINTS; i++) {
					  Xcircle[i] = (float) xmask[ima][track][i];
					  Ycircle[i] = (float) ymask[ima][track][i];
					  Wpoint[i]=1.;
				      }
				      XC = YC = RC = 0.;
				      Npoints_fit = POINTS;
				      Signe=10;
				      circle_minuit(POINTS);
				      // suppose that start and end angles are between 180degrees
				      bufmin = 1000.;
				      bufmax =-1000.;
				      rc_int = RC; xc_int = XC; yc_int = YC;
				      for(i=0;i<POINTS; i++) {
					  buf = atan2((Ycircle[i]-YC),(Xcircle[i]-XC));
					  if(buf < bufmin) {
					      imin = i; bufmin = buf;
					  }
					  if(buf > bufmax) {
					      imax = i; bufmax = buf;
					  }
				      }
				      //check for the -Pi - +Pi bounds
				      if((bufmax - bufmin) > (float)TMath::Pi()) {
					  bufmin = 1000.;
					  bufmax =-1000.;
					  for(i=0;i<POINTS; i++) {
					      buf = atan2((Ycircle[i]-YC),(Xcircle[i]-XC));
					      if (buf < 0.) buf += 2.*(float)TMath::Pi();
					      if(buf < bufmin) {
						  imin = i; bufmin = buf;
					      }
					      if(buf > bufmax) {
						  imax = i; bufmax = buf;
					      }
					  }
				      }
				      stx = Xcircle[imin];      sty = Ycircle[imin];
				      endx = Xcircle[imax];     endy = Ycircle[imax];
				      printf("imin %d imax %d\n", imin, imax);
				      //printf("Circle 2 Image=%d track=%d radius=%f points=%f %f, chi2 = %f\n",
				      //     ima,track,rc_int,xc_int,yc_int, qchi2);

				      double phi1, phi2, range, stangle, endangle;
				      /*                                                     */
				      /*  Find the angle supporting the chord between start  */
				      /*  and end of the track, i.e. between the first and   */
				      /*  third points measured on the mask... (see p. 169   */
				      /*  of Br.& Sm): chord = 2.*rad*sin(alpha/2).          */
				      /*                                                     */
				      /*             WORK WITH RADIANS !!!                   */
				      /*                                                     */
				      /*           but for drawing of arcs use degrees       */
				      
				      
				      range = sqrt((endx-stx)*(endx-stx) + (endy-sty)*(endy-sty));
				      range = 2.*asin(range/(2.*rc_int));
				      
				      //  simplest way to get angle in 2Pi
				      stangle = atan2((sty-yc_int),(stx-xc_int));
				      endangle = atan2((endy-yc_int),(endx-xc_int));
				      if (stangle > endangle) endangle += 2.*TMath::Pi();

				      phi1 = stangle*180./TMath::Pi(); 
				      phi2 = endangle*180./TMath::Pi();
				      int epx = (int) endx;
				      int epy = (int) endy;

				      double delta_angle;
				      //prevent cicling in stright line tracks
				      if(endangle-stangle > 0.) {
					  delta_angle = range/2./sqrt(pow(endx-stx,2) + pow(endy-sty,2));
					  printf("delta_angle %f cells %f \n",delta_angle, 
						 2.*sqrt(pow(endx-stx,2) + pow(endy-sty,2)));
					  if(stangle+delta_angle == stangle) printf("delta range is too small!!!\n");
				      } else {
					printf("There is a line, no circle!\n");
				      }

				      double cur_angle = stangle;
				      double nat_log=2.718281828459045;
				      /*
				      char ctmp1[3],ctmp2[3], name[10];
				      sprintf(ctmp1,"%u",ima);
				      sprintf(ctmp2,"%u",track);
				      strcpy(name,"h");
				      strcat(name,ctmp1);
				      strcat(name,"tr");
				      strcat(name,ctmp2);
				      hbr = new TH1F(name,name, 41,-20.,20.);
						*/
				      sumbright_all = 0.;
				      int iwidth = 0;//may be 5? was 3
				      int xpix; int ypix;
				      int wp[100];
				      int curx = 0; int cury = 0;
				      int binxprev_circle = 0; int binyprev_circle = 0;
				      i = 0;
				      double dr;
				      dr = 0.1*sqrt(2.); //step to move along r
				      int npoints = 0;
				      int npp = 0;
				      int npoints_bright = 0;
				      int ip; int ii;
				      epx = (int) (xc_int + rc_int*cos(cur_angle+range));
				      epy = (int) (yc_int + rc_int*sin(cur_angle+range));
				      printf("epx %d epy %d, endangle %f cur %f\n", epx, epy, endangle, cur_angle);
				      for(ii = 0; ii<1000000; ii++) {
					  ip = 0; 
					  int wpmax = 0; int xpmax = 0; int ypmax = 0;
					  int binxprev = 0;
					  int binyprev = 0;
					  
					  x_cur = xc_int + rc_int*cos(cur_angle);
					  y_cur = yc_int + rc_int*sin(cur_angle);
					  xpix = (int) x_cur; ypix = (int) y_cur;
					  curx = xpix; cury = ypix;
					  //if (track == 1) printf("curx %d, epx %d,  cury %d, epy %d\n",
					  // curx, epx, cury, epy);
					  Int_t binx, biny;
					  binx = xaxis->FindBin(xpix);
					  biny = yaxis->FindBin(ypix);
					  if (binxprev_circle != binx || binyprev_circle != biny) {
					      // now go along the radius
					      binxprev_circle = binx; binyprev_circle = biny;
					      
					      if ( ima == 0) {
						  bin_arr = left->GetBin(binx,biny);
						  bright = (int) left->GetBinContent(bin_arr);
					      } else {
						  bin_arr = right->GetBin(binx,biny);
						  bright = (int) right->GetBinContent(bin_arr);
					      }
					      int ifl = 0;
					      if(ifl == 0 && bright > wpmax) {
						wpmax = bright; xpmax = xpix; 
						ypmax = ypix; rmax = rc_int;
						binxmax = binx; binymax = biny; 
					      }
					      binxprev = binx; binyprev = biny;
					      ip++;
					      float rcur = rc_int + dr;
					      do {
						bright = 0;
						xpix = (int) (xc_int + rcur*cos(cur_angle));
						ypix = (int) (yc_int + rcur*sin(cur_angle));
						binx = xaxis->FindBin(xpix);
						biny = yaxis->FindBin(ypix);
						if(binx != binxprev || biny != binyprev) {
						  if ( ima == 0) {
						      bin_arr = left->GetBin(binx,biny);
						      bright = (int) left->GetBinContent(bin_arr);
						  } else {
						      bin_arr = right->GetBin(binx,biny);
						      bright = (int) right->GetBinContent(bin_arr);
						  }
						  int ifl = 0;
						  if(ifl == 0 && bright > wpmax) {
						    wpmax = bright; xpmax = xpix; 
						    ypmax = ypix; rmax = rcur;
						    binxmax = binx; binymax = biny; 
						  }
						  binxprev = binx; binyprev = biny;
						  ip++;
						}
						rcur += dr;
					      } while (ip < iwidth);
					      ip = 1;
					      rcur = rc_int - dr;
					      binxprev = binxprev_circle; 
					      binyprev = binyprev_circle;
					      do {
						bright = 0;
						xpix = (int) (xc_int + rcur*cos(cur_angle));
						ypix = (int) (yc_int + rcur*sin(cur_angle));
						binx = xaxis->FindBin(xpix);
						biny = yaxis->FindBin(ypix);
						if(binx != binxprev || biny != binyprev) {
						  if ( ima == 0) {
						      bin_arr = left->GetBin(binx,biny);
						      bright = (int) left->GetBinContent(bin_arr);
						  } else {
						      bin_arr = right->GetBin(binx,biny);
						      bright = (int) right->GetBinContent(bin_arr);
						  }
						  int ifl = 0;
						  if(ifl == 0 && bright > wpmax) {
						    wpmax = bright; xpmax = xpix; 
						    ypmax = ypix; rmax = rcur;
						    binxmax = binx; binymax = biny; 
						  }
						  binxprev = binx; binyprev = biny;
						  ip++;
						}
						rcur -= dr;
					      } while(ip < iwidth);
					      //now go from point with maximal brightness 
					      //along radius in both sides before brightness
					      //becomes less in e times
					      float background_cadr;
					      float gradient_prev = 0.;
					      if(ima == 0) background_cadr = 0.75*aver_l;
					      if(ima == 1) background_cadr = 0.75*aver_r;
					      if( npoints > 10000) break;
					      ip = 0;
					      iwidth = 2;
					      //hbr->Fill((Axis_t) ip, (Stat_t) wpmax);
					      wt[ima][track][npoints] = 
						(((float) wpmax) -  background_cadr);
					      wpmax = wp[ip] = (int) wt[ima][track][npoints];
					      sumbright_all += (float) wpmax;
					      //xp[npoints_bright] = 
					      xt[ima][track][npoints] = (float) xpmax;
					      //yp[npoints_bright] = 
					      yt[ima][track][npoints] = (float) ypmax;
					      if(wpmax < 0 ) {wt[ima][track][npoints]  = 0.; wp[ip]  = wpmax = 0;}
					      x_mark[npoints] = (Float_t) xpmax;
					      y_mark[npoints] = (Float_t) ypmax;
					      //Draw_region(image, xpmax, ypmax, (int) wt[image][track][npoints]);
					      npoints++;
					      npoints_bright++; ip++;
					      npp++;
					      binxprev = binxmax; binyprev = binymax;
					      rcur = rmax + dr;
					      if (wpmax > 0) {
						do {
						  bright = 0;
						  xpix = (int) (xc_int + rcur*cos(cur_angle));
						  ypix = (int) (yc_int + rcur*sin(cur_angle));
						  binx = xaxis->FindBin(xpix);
						  biny = yaxis->FindBin(ypix);
						  if(binx != binxprev || biny != binyprev) {
						    if ( ima == 0) {
							bin_arr = left->GetBin(binx,biny);
							bright = (int) left->GetBinContent(bin_arr);
						    } else {
							bin_arr = right->GetBin(binx,biny);
							bright = (int) right->GetBinContent(bin_arr);
						    }
						    //hbr->Fill((Axis_t) ip, (Stat_t) bright);
						    binxprev = binx; binyprev = biny;
						    int ifl = 0;
						    if(ifl == 0) {
						      wp[ip] = (int) (((float) bright) - background_cadr);
						      sumbright_all += wp[ip];;
						      if(wp[ip] < 0) wp[ip] = 0;
						      if(ip < iwidth && npoints < 10000 - 1) {
							xt[ima][track][npoints] = (float) xpix; 
							yt[ima][track][npoints] = (float) ypix;
							wt[ima][track][npoints]  = (float) wp[ip];
							x_mark[npoints] = (Float_t) xpix;
							y_mark[npoints] = (Float_t) ypix;
							//Draw_region(image, xpix, ypix, (int) wt[image][track][npoints]);
							npoints++;
						      }
						      //else
						      //Draw_region(image, xpix, ypix, (int) wp[ip]);
						      npoints_bright++; ip++;
						      if(wp[ip - 1] == 0) break;
						      if(((float) wpmax)/((float) wp[ip - 1]) > 2.*nat_log) break;
						      float gradient = (float) (wp[ip - 1] - wp[ip - 2]);
						      gradient = gradient / (float) wp[ip - 2];
						      // if grad was negative but it becames to be positive- break
						      if(ip > 2 && gradient*gradient_prev < 1. && gradient > 0.) break;
						      // if gradient changes very slowly - break
						      if(fabs(gradient) < 0.1 && ip > 10) break;
						      gradient_prev = gradient;
						    }
						  }
						  rcur += dr;
						} while (ip < 50);
						rcur = rmax - dr;
						binxprev = binxmax; 
						binyprev = binymax;
						ip = 1;
						do {
						  bright = 0;
						  xpix = (int) (xc_int + rcur*cos(cur_angle));
						  ypix = (int) (yc_int + rcur*sin(cur_angle));
						  binx = xaxis->FindBin(xpix);
						  biny = yaxis->FindBin(ypix);
						  if(binx != binxprev || biny != binyprev) {
						    if ( ima == 0) {
							bin_arr = left->GetBin(binx,biny);
							bright = (int) left->GetBinContent(bin_arr);
						    } else {
							bin_arr = right->GetBin(binx,biny);
							bright = (int) right->GetBinContent(bin_arr);
						    }
						    binxprev = binx; binyprev = biny;
						    //hbr->Fill((Axis_t) -ip, (Stat_t) bright);
						    int ifl = 0;
						    if(ifl == 0) {
						      wp[ip] = (int) (((float) bright) - background_cadr);
						      sumbright_all += wp[ip];;
						      if(wp[ip] < 0) wp[ip] = 0;
						      if(ip < iwidth && npoints < 10000 - 1) {
							xt[ima][track][npoints] = (float) xpix; 
							yt[ima][track][npoints] = (float) ypix;
							wt[ima][track][npoints]  = (float) wp[ip];
							x_mark[npoints] = (Float_t) xpix;
							y_mark[npoints] = (Float_t) ypix;
							//Draw_region(image, xpix, ypix, (int) wt[image][track][npoints]);
							npoints++;
						      }
						      //else
						      //Draw_region(image, xpix, ypix, (int) wp[ip]);
						      npoints_bright++; ip++;
						      if(wp[ip - 1] == 0) break;
						      if(((float) wpmax)/((float) wp[ip - 1]) > 2.*nat_log) break;
						      float gradient = (float) (wp[ip - 1] - wp[ip - 2]);
						      gradient = gradient / (float) wp[ip - 2];
						      // if grad was negative but it becames to be positive- break
						      if(ip > 2 && gradient*gradient_prev < 1. && gradient > 0.) break;
						      // if gradient changes very slowly - break
						      if(fabs(gradient) < 0.1 && ip > 10) break;
						      gradient_prev = gradient;
						    }
						  }
						  rcur -= dr;
						} while (ip < 50);
					      }
					  }
					  cur_angle += delta_angle;
					  if(curx == epx && cury == epy) break;
					  if(cur_angle > endangle) break;
				      }
				      printf("npoints %d is %d, in a corridor %d of average width %6.3f\n", 
					     track, npoints, npoints_bright, (float) npoints_bright/(float) npp);
				      //hbr->Scale(1./(float) npp);
				      float bright_av_all;
				      //bright_av_all = sumbright_all/(float) npoints_bright;
				      bright_av_all = sumbright_all/(float) npp;
				      brrr[ima][track] = bright_av_all*(float) npoints_bright/(float) npp;
				      printf("average brightness on track %d is %f\n", 
					     track, bright_av_all);
				      np[ima][track] = npoints;// 100;
				      float WpMin = 0.;
				      float WpMax = 4096.;
				      for(i=0; i<np[ima][track]; i++) {
					  Xcircle[i] = (float) xt[ima][track][i];
					  Ycircle[i] = (float) yt[ima][track][i];
					  Wpoint[i] = ((float)  wt[ima][track][i]);
					  if (Wpoint[i] > WpMax) WpMax = Wpoint[i];
					  if (Wpoint[i] < WpMin) WpMin = Wpoint[i];
				      }
				      for(i = 0; i< np[ima][track]; i++) {
					  if(WpMin != WpMax) Wpoint[i]  = (Wpoint[i]-WpMin)/(WpMax - WpMin);
					  if(WpMin == WpMax) Wpoint[i]  = Wpoint[i]/WpMax;
				      }
				      px_mark = &x_mark[0];
				      py_mark = &y_mark[0];
				      markp = new TPolyMarker((Int_t) npoints, px_mark, py_mark, "same");
				      fObjTbl->Add(markp);
				      markp->SetMarkerStyle(8);
				      markp->SetMarkerSize(0.2);
				      markp->SetMarkerColor(41);
				      if(ima == 0) pleft->cd();
				      if(ima == 1) pright ->cd();
				      markp->Draw("same");	     
				      
				      Npoints_fit = npoints;
				      Signe=10;
				      circle_minuit(npoints);
				      
				      // suppose that start and end angles are between 180degrees
				      rc_int = RC; xc_int = XC; yc_int = YC;
				      // suppose that start and end angles are between 180degrees
				      bufmin = 1000.;
				      bufmax =-1000.;
				      for(i=0;i<np[ima][track]; i++) {
					  buf = atan2((Ycircle[i]-YC),(Xcircle[i]-XC));
					  //x_cur = XC + RC*cos(buf);
					  //y_cur = YC + RC*sin(buf);
					  //TMarker *mark = new TMarker(x_cur,y_cur,7);
					  //fObjTbl->Add(mark);
					  //mark->SetMarkerStyle(8);
					  //mark->SetMarkerSize(0.1);
					  //mark->SetMarkerColor(51+mypass);
					  //mark->Draw("same");	     
					  //draw circle for the corridor
					  x_cur = XC + RC*cos(buf);
					  y_cur = YC + RC*sin(buf);
					  
					  if(buf < bufmin) {
					      imin = i; bufmin = buf;
					  }
					  if(buf > bufmax) {
					      imax = i; bufmax = buf;
					  }
				      }
				      //					  gPad->Modified();
				      //					  gPad->Update();
				      
				      //check for the -Pi - +Pi bounds
				      if((bufmax - bufmin) > (float)TMath::Pi()) {
					  bufmin = 1000.;
					  bufmax =-1000.;
					  for(i=0;i<np[ima][track]; i++) {
					      buf = atan2((Ycircle[i]-YC),(Xcircle[i]-XC));
					      if (buf < 0.) buf += 2.*(float)TMath::Pi();
					      if(buf < bufmin) {
						  imin = i; bufmin = buf;
					      }
					      if(buf > bufmax) {
						  imax = i; bufmax = buf;
					      }
					  }
				      }
				      TEllipse *myarc = new TEllipse(XC, YC, RC, 0,
								     bufmin*180./TMath::Pi(), bufmax*180./TMath::Pi(), 0);
				      fObjTbl->Add(myarc);
				      myarc->SetLineWidth(0.02);
				      myarc->SetLineColor(51);
				      myarc->SetNoEdges(1);
				      myarc->Draw("same");	    
				      
				      stx = Xcircle[imin];   sty = Ycircle[imin];
				      endx = Xcircle[imax];     endy = Ycircle[imax];
				      // vertex store
				      
				      xcent_vx[track] = xc_int;
				      ycent_vx[track] = yc_int;
				      rcent_vx[track] = rc_int;
				      //if (mypass == 1)
				      //{
				      //  printout to the output file 	 
				      fprintf(out_file,"%7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f \n",
					      rc_int, xc_int, yc_int, stx, sty, endx, endy);
				      fprintf(out_file,"%7d\n", np[ima][track]);
				      int k = 0;
				      for ( int ic=0; ic<np[ima][track]; ic++) {
					  k++;
					  bright = 0;
					  binx = (Int_t) Xcircle[ic]; biny = (Int_t) Ycircle[ic];
					  if ( ima == 0) {
					      bin_arr = left->GetBin(binx,biny);
					      bright = (int) left->GetBinContent(bin_arr);
					      //bright -= aver_l;
					  } else {
					      bin_arr = right->GetBin(binx,biny);
					      bright = (int) right->GetBinContent(bin_arr);
					      //bright -= aver_r;
					  }
					  //if (bright < 0) bright = 0;
					  //write out sumbright instead of bright!!  
					  //fprintf(out_file," %lu %lu %d ", binx, biny, (int)bright);
					  fprintf(out_file," %u %u %d ", binx, biny, (int) wt[ima][track][ic]);
					  if(k == 5) {
					    k = 0;
					    fprintf(out_file,"\n"); 
					  }
				      }
				      //}
				      //end over tracks loop
				  }
				  //   vertex finding
				  if(mypass == 0) {
				      Xvertax = Yvertax = Dvertax = 0.;
				      // ==============  ADDED MISSING INITIALIZATION TO 0.0 FOR THE
				      // ==============  FOLLOWING THREE VARIABLES. flv 25 March 2010
				      float XVertex_3=0.0, YVertex_3=0.0, DVertex_3=0.0;
				      if(num_tracks != 1) {
					XY_Vertex_Net(num_tracks, &xcent_vx[0], &ycent_vx[0], 
						      &rcent_vx[0], &Xvertax, &Yvertax, 
						      &Dvertax);
				        TMarker *mark = new TMarker(Xvertax,Yvertax,5);
					mark->SetMarkerColor(38);
					mark->SetMarkerStyle(5);
					mark->SetMarkerSize(1.5);
					mark->DrawMarker(Xvertax, Yvertax);
					//gPad->Modified();
					//					  gPad->Update();
					  
					printf(" Vertex coordinates: x %f, y %f; error %f\n",
						 Xvertax, Yvertax, Dvertax );
				      }
				      if(num_tracks == 3 && DVertex_3 < Dvertax &&
					 XVertex_3 < hsize && YVertex_3 < vsize) {
					  Dvertax = DVertex_3;
					  Xvertax = XVertex_3;
					  Yvertax = YVertex_3;
				      }
				  }
				  //if (mypass == 1)
				  fprintf(out_file,
					  " Vertex position: X: %f, Y: %f, Error: %f \n", 
					  Xvertax, Yvertax, Dvertax);
			      }
			      gPad->Update();
			    }
			    fprintf(out_file, " %f %f \n", aver_l, aver_r);
			    int i;
			    for(i = 0; i< num_tracks; i++) {
			      fprintf(out_file, " %f %f \n", brrr[0][i], brrr[1][i]);
			    }
			  
			  
			  // end loop over images
			  //			    gApplication->Terminate(0);
			  fclose(out_file);
			  // run the Gil Pontecorvo's programme
			  recon_main(test_outfile, 100);
			  int retval;
			  new TGMsgBox(gClient->GetRoot(),gClient->GetRoot(),
				       "Open new event?",
				       "Would you like to open new event or close application?",
				       kMBIconStop , kMBYes|kMBDismiss|kMBClose, &retval);
			  //			  printf(" return code = %i \n",retval);
			  
			  if (retval == kMBClose) gApplication->Terminate(0);
			  if (retval == kMBYes)
			    {
			      RemoveCreated();
			      GetNewEvent();
			    }
			}
		    }
		  else
		    {
		      --ipoint[cur_image];
		      //cout << " All tracks for this image are enumerated\n";
		      int retval;
		      new TGMsgBox(gClient->GetRoot(),gClient->GetRoot(),
				   "msg",
				   "All tracks for this image are enumerated",
				   kMBIconStop , kMBDismiss, &retval);
		    }
		}
		return;
	      }
	    else
	      {
		//cout << "histogram does not exist\n";
		exit(1);
	      }
	  }
	/* at this point only "real" images can occur */
	
	if (ipoint[cur_image] == num_tracks * POINTS) 
	  {
	    /* Are all points for the image enumerated? */
	    
	    int retval;
	    new TGMsgBox(gClient->GetRoot(),gClient->GetRoot(),
			 "msg",
			 "All tracks for this image are enumerated",
			 kMBIconStop , kMBDismiss, &retval);
	    return;
	  }
	/*
	htmp = new TH2F("htmp","Image ZOOM", 
			nnewbins,(Axis_t)(widthx*(binx-nside-1)),
			(Axis_t)(widthx*(binx+nside)),
			nnewbins,(Axis_t)(widthy*(biny-nside-1)),
			(Axis_t)(widthy*(biny+nside)));
	*/
	TH2F *c = 0;                                                               
	c = (TH2F*)gROOT->FindObject("htmp"); if (c) c->Delete(); c = 0;  
	htmp = new TH2F("htmp","ZOOM", 
			nnewbins,(Axis_t)(widthx*(binx-nside)),
			(Axis_t)(widthx*(binx+nside+1)),
			nnewbins,(Axis_t)(widthy*(biny-nside)),
			(Axis_t)(widthy*(biny+nside+1)));
	{
	  
	  //	  cout << this->GetName()<<"   image ="<<cur_image<<endl;
	  int itmpx = 0; 
	  //int lowx = binx-nside-1;
	  int lowx = binx-nside;
	  //	int maxx = binx+nside;
	  //int lowy = biny-nside-1;
	  int lowy = biny-nside;
	  //	int maxy = biny+nside;
	  int iix,iiy;
	  Stat_t weight;
	  Int_t bin_arr1;
	  //	  itmpx +=1; itmpy +=1; /* from next bin */
	  //	  cout << " ZOOM ----- "<<nnewbins<<endl;
	  while (++itmpx<=nnewbins)
	    {
	      int itmpy = 0;
	      iix = itmpx + lowx;
	      while ( ++itmpy <= nnewbins )
		{
		  iiy = itmpy + lowy;
		  /*		  cout << itmpx<<"="<<iix 
				  << " " << itmpy<<"="<<iiy<<endl;
				  //		  printf(" %d  %d  %f  \n",iix,iiy,weight);
		  */
		  bin_arr = htmp->GetBin(itmpx,itmpy);
		  if ( (iix > nxbins || iix < 1) || 
		       (iiy > nybins || iiy < 1))
		    weight = 0;
		  else
		    //		    weight = (Stat_t)pCDBuffer[cur_image][iiy][iix];
		    {
		      bin_arr1 = this->GetBin(iix,iiy);
		      weight = this->GetBinContent(bin_arr1);
		      htmp->SetBinContent(bin_arr,weight);
		    }
		}
	    }
	  //	  cout << widthx*(binx-nside-1) << " " <<widthx*(binx+nside)<<endl;
	  //	  cout << widthy*(biny-nside-1) << " " <<widthy*(biny+nside)<<endl;
	  //	  ct = new TCanvas("ct","tmp canv",150,150,200,200);
	  /*	  ct->SetFillColor(42);
		  ct->GetFrame()->SetFillColor(21);
		  ct->GetFrame()->SetBorderSize(6);
		  ct->GetFrame()->SetBorderMode(-1);*/
	  //	  ct->SetLogz();
	  //	  ct->ToggleEventStatus();
	  //	  ct->RangeAxis(0.25,0.25,0.75,0.75);
	  //	  ct->cd();
	  //	  gPad->SetPad(-0.12,-0.12,1.12,1.12);
	  ptmp->cd();
	  ptmp->SetLogz();
	  
	  //	  cout<<"DIMENSION = "<<ct->GetUxmax()<<endl;
	  htmp->SetMaximum(htmp->GetMaximum());
	  htmp->SetMinimum(htmp->GetMinimum());
	  htmp->SetContour(550);
	  //	  htmp->UseCurrentStyle();
	  htmp->Draw("col");
	  //	  gPad->Modified();
	  gPad->Update();
	  //	  t1 =  second();
	  //	  cout<< "Time = "<<(t1-t0)<<endl;
	  //	  ptmp->Update();
	  
	}
      }
      break;

    case kButton3Down:
      printf("Botton 3 pressed\n");
      break;
    }
}
//^MAIN
int main(int argc, char** argv)
{
  // flv - added lines have "// flv " comment on the right
  // these are to save the original command line parameters  // flv
  int  flv_argc;                                             // flv
  char *flv_argv[10];                                        // flv
  int flv_i;                                                 // flv
  int flv_i_up;                                              // flv
  char flv_file[1000];                                       // flv
                                                             // flv
  flv_argc = argc;                                           // flv
  if (flv_argc <10 ) {                                       // flv
    flv_i_up = flv_argc;                                     // flv
  } else {                                                   // flv
    flv_i_up = 10;                                           // flv
  }                                                          // flv
  for (flv_i = 0; flv_i < flv_i_up; flv_i++) {               // flv
    flv_argv[flv_i] = argv[flv_i];                           // flv
  }                                                          // flv
  // flv - end of added lines
  
  //cout<<"before paleee"<<endl;
  TApplication theApp("App", &argc, argv);

  int length = 0;
  char ctmp[10000];
  char lfile[1000], rfile[1000];
  char fout_file[1000];
  char *cptmp1, *cptmp2;

  //cout<<"before palette"<<endl;
  SetVeryNewPalette();
  //cout<<"after palette"<<endl;

  fObjTbl = new TObjectTable(10000);
  //^egcs  fi.fFileTypes = (char **)filetypes;
  fi.fFileTypes = filetypes;

  strcpy(WorkDir,gSystem->WorkingDirectory());

  strcpy(iniDir,WorkDir);
  strcat(iniDir,"/data");
  printf(" iniDir %s\n", iniDir);
  gSystem->ChangeDirectory(iniDir);
       /* print files in current directory in reverse order */
  //  WaitMsg = new WaitFrame(fClient->GetRoot(),300,70);

  // flv - added lines have "// flv " comment on the right
  // these are to use the possible filename given as command line parameter                         // flv
  if (flv_argc == 1) {                                                                              // flv
    new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);    // this line was in the original code
  } else if (flv_argc == 2) {                                                                       // flv  
    // check that the name of the executable matches with the charge of incident pion in the        // flv
    // event name given as command line parameter                                                   // flv 
    if ( (strstr(flv_argv[0],"positive") != NULL && strstr(flv_argv[1],"positive") != NULL) ||
         (strstr(flv_argv[0],"negative") != NULL && strstr(flv_argv[1],"negative") != NULL)   ) {   // flv
      // get iniDir from command line parameter                                                     // flv
      strcpy(iniDir, gSystem->DirName(flv_argv[1]));                                                // flv
      gSystem->ChangeDirectory(iniDir);                                                             // flv
      printf(" *** - iniDir changed to %s (from command line parameter) \n", iniDir);               // flv
      // set the filename from command line parameter                                               // flv
      strcpy(flv_file,gSystem->BaseName(flv_argv[1]));                                              // flv
      fi.fFilename = flv_file;                                                                      // flv
      // print few help lines                                                                       // flv
      printf(" *** - filename preset to %s (from command line parameter) \n", fi.fFilename);        // flv
      printf(" *** -  In the 'Open' dialog window you simply have to click the 'Open' button \n");  // flv
      printf(" *** -  WARNING. Check the event number for safety. If the preset event number\n");   // flv
      printf(" *** -           is not the one you want to measure click the 'Cancel' button).\n");  // flv
      // new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen, &fi);                       // flv
    } else {                                                                                        // flv
      printf("WARNING: the sign (positive or negative) of input filename (%s) \n",  flv_argv[1]);   // flv
      printf("         does not match with the program name (%s) \n", flv_argv[0]);                 // flv
      printf("         Please select a different input filename or click the 'Cancel' button.\n");  // flv
      new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen,&fi);                        // flv
    }                                                                                               // flv
  } else {                                                                                          // flv
    printf("ERROR: unexpected value of argc (= %d). Terminating %s \n", flv_argc, flv_argv[0] );    // flv
    exit(0);                                                                                        // flv
  }                                                                                                 // flv
  printf("================================================================\n");                     // flv
  printf("Running the program for %s pions of mean momentum %f MeV/c\n", ReconMainSign(), P_MEAN);  // flv
  printf("Magnetic field set at %f kGauus\n", MFIELD);                                              // flv
  printf("================================================================\n");                     // flv

  // flv - end of added lines

  timer.Reset();
  timer.Start(kTRUE);

  if (fi.fFilename == NULL) exit(0);

  strcpy(iniDir,gSystem->DirName(fi.fFilename));
  printf("selected file  - %s \n",fi.fFilename);

  strcpy(lfile,gSystem->BaseName(fi.fFilename));
  
  /*    forming whole filenames islower  */
  strcpy(rfile,lfile);
  if (isupper(rfile[0]))
	 rfile[0] = 'R';
  else
	 rfile[0] = 'r';
  //cout<< "lfile = "<<lfile<<endl;
  //cout
  //cout<< "rfile = "<<rfile<<endl;

  cptmp1 = gSystem->ConcatFileName(iniDir,(char *)lfile);
  cptmp2 = gSystem->ConcatFileName(iniDir,(char *)rfile);
  
  //cout<< "clfile = "<<cptmp1<<endl;
  //cout<< "crfile = "<<cptmp2<<endl;
  
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
  
  int iitmp = 0; 
  cptmp1 = strcpy(fout_file,cptmp1);
  
  if ((cptmp2 = strrchr(fout_file,'.')) != NULL) // looking for "."
	 *cptmp2 = '\0';
  
  if ((cptmp2 = strrchr(fout_file,'/')) != NULL) 
	 iitmp = ((intptr_t)cptmp2 - (intptr_t)cptmp1)/sizeof(char) + 1;
  
  iitmp += 1;     /*  removing first char from name  */
  iitmp *=sizeof(char);
  
  strcpy(eventname,fout_file + iitmp);
  //fprintf(stderr,"\n event : %s \n",eventname);
  char dirout[10] = {"dig/"};
  strcpy(ofilename, dirout);
  int i;
  int k = 0;
  for (i = (int) strlen(dirout); 
		 i<(int) strlen(dirout) + (int)strlen(eventname); i++)
	 {ofilename[i] = eventname[k]; k++;}
  for (i= (int) strlen(dirout) + (int)strlen(eventname);i < 100;i++)
	 ofilename[i] = '\0';
  strcat(ofilename,".dig");
  
  //  fprintf(stderr,"\n event filename : %s \n",ofilename);
  
  gzgets(im_left, ctmp,10000);
  
  int ntmp;
  //istrstream Inbuffer(ctmp,10000);
  std::istringstream Inbuffer(ctmp);
  
  while ( Inbuffer >> ntmp )
    length++;
  
  printf("Numbers of INTs = %d\n",length);
  gzrewind(im_left);
  
  int length1 = 0;
  while( gzgets(im_left, ctmp,10000) != NULL)
    length1++;
  
  gzrewind(im_left);
  vsize = length1; hsize = length;
  printf("vsize = %d, hsize = %d\n",vsize,hsize);
  
  
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
  int image, j;
  //unsigned maxBright = 0;
  for ( image=0; image<2; image++ )
    {
      itrack[image] = 1;
      ipoint[image] = 0;
      for ( i=0; i<vsize; i++ )
	{
	  if (image == 0)
	    gzgets(im_left,ctmp,10000);
	  else
	    gzgets(im_right,ctmp,10000);
	  
	  std::istringstream Inbuffer(ctmp);
	  //istrstream Inbuffer(ctmp,10000);
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
  
  TH2F *c = 0;                                                               
  c = (TH2F*)gROOT->FindObject("left"); if (c) c->Delete(); c = 0;  
  left = new TH2F("left","bitmap for left image", 
		  (Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);
  
  c = (TH2F*)gROOT->FindObject("right"); if (c) c->Delete(); c = 0;  
  right = new TH2F("right","bitmap for right image",
		   (Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);
  
  TH1F *bc = 0;
  bc = (TH1F*)gROOT->FindObject("bleft"); if (bc) bc->Delete(); bc = 0;
  bleft = new TH1F("bleft","left brightness",4096,-0.5,4095.5);
  bc = (TH1F*)gROOT->FindObject("bright"); if (bc) bc->Delete(); bc = 0;
  bright = new TH1F("bright","right brightness",4096,-0.5,4095.5);

  Axis_t xx,yy;
  Stat_t ww;
  //  Stat ww1;    // unused
  //  float sumw;    // unused
 //   int inum;    // unused
  for (image=0; image<2; image++)
    {
      for (i=0; i<vsize; i++)
	{
	  for (j=0; j<hsize; j++)
	    {
	      ww = (Stat_t) pCDBuffer[image][i][j];
	      xx = (Axis_t) j;//(j + 0.5);
	      yy = (Axis_t) (vsize - i - 1);//0.5);  
	      if (image == 0) 
		{
		  left->Fill(xx,yy,ww);
		  bleft->Fill((Float_t)ww);
		}
	      else
		{
		  right->Fill(xx,yy,ww);
		  bright->Fill((Float_t)ww);
		}
	    }
	}
    }
  float ratio;
  ratio = (float) vsize/(float) hsize;
  //  printf("vsize = %d, hsize = %d, ratio = %f\n",vsize,hsize,ratio);
  //  ptmp = new TPad("ptmp"," ",0.37,0.01,0.63,0.34,21);
  
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
  
  //cout <<"x , y ="<<ixWindow<<" "<<iyWindow<<endl;
  
  delete gROOT->GetListOfCanvases()->FindObject("c2");
  c2 = new TCanvas("c2","Event digitizing",10,10,ixWindow,iyWindow);
  pleft = new TPad("pleft","Left image",    xleft, yleft, xright, yright, 21);
  pright = new TPad("pright","Right image", xleft + 0.5, yleft, xright + 0.5, yright, 21);
  ptmp = new TPad("ptmp"," ",xlefttmp, ylefttmp, xrighttmp, yrighttmp, 21);
  
  pMsg = new TPad("pMsg"," ",0.01,0.01,xlefttmp-0.01,yrighttmp,28);

  float xrighttmp2 = (0.99 - xrighttmp-0.01)/2.+ xrighttmp+0.01;
  pbleft = new TPad("pbleft","",xrighttmp+0.01,0.01,xrighttmp2,yrighttmp,24);
  pbright = new TPad("pbright","",xrighttmp2,0.01,0.99,yrighttmp,24);
  //  std::cout <<"xrighttmp = "<<xrighttmp<<std::endl;
  //  std::cout <<"xrighttmp2 = "<<xrighttmp2<<std::endl;
  
  pleft->Draw();
  pright->Draw();
  
  pleft->UseCurrentStyle();
  pright->UseCurrentStyle();
  
  ptmp->Draw();
  ptmp->UseCurrentStyle();
  
  pMsg->Draw();
  
  pMsg->cd();
  pl1 = new TPaveLabel(0.07,0.15,0.93,0.45,eventname,"br");
  pl1->SetFillColor(45);
  pl1->SetTextSize(0.5);
  pl1->Draw();

  left->SetContour(550);
  right->SetContour(550);
  pleft->cd();
  left->Draw("colz");
  pright->cd();
  right->Draw("colz");
  
  c2->cd();
  pbleft->Draw();
  c2->cd();
  pbright->Draw();

  pbleft->SetLogy();
  pbright->SetLogy();
  
  pbleft->cd();
  bleft->Draw();
  pbright->cd();
  bright->Draw();
  

  c2->Modified(); 
  c2->Update();

  /*
  c2->LoadRQ_OBJECT();
  c2->cd();

  std::cout<<"connect ="<<c2->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "", 0,
			 "execevent(int,int,int,TObject*)")<<std::endl;
  */
  /*  pbright->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
      "exec_event(Int_t,Int_t,Int_t,TObject*)");
  */

  
  //new canvas for GIF file
  //strcpy(gifname, eventname);
  //strcat(gifname, ".gif");
  //c2->Print(gifname,"gif");
  maiWin = new MyChoiseFrame(gClient->GetRoot(), 200, 380);
  
  //cout << "Real time after drawing ="<<timer.RealTime()<< endl;
  //cout << "CPU time after drawing ="<<timer.CpuTime()<< endl;
  printf("There was running event number %s\n\a", eventname);

  //gGXW->Bell(100);
  //gClient->Bell(100);
  theApp.Run(kTRUE);
 
}
#include "TClass.h"
/*^execevent  */
void execevent(int event, int x, int y, TObject *selected)
{
  TCanvas *c = (TCanvas*) gTQSender;
  printf("Pad %s: event=%d, x=%d, y=%d, selected=%s\n", c->GetName(),
	 event, x, y, selected->IsA()->GetName());
}
/***************  Usage ******************/
void usage(char* argv[])
{

  fprintf(stderr,"\n\n Usage: %s left_file right_file\n\n",argv[0]);
}
/***************** Argument line parsing ***********************************/
//
//int parse_options(int argc, char **argv)
//{
//  int c;                                                      /* option char */
//  extern char *optarg;                              /* for use with getopt() */
//  extern int optind, optopt;                        /* for use with getopt() */
//  int errflg = 0;                                   /* for use with getopt() */
//
//  while ((c = getopt(argc, argv, "X:n:s:o:vh")) != -1)
//*/

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
void SetPart()
{
  //  char ctmp[3];    // unused
  //  int itemp;    // unused

  // int retval;    // unused
  //new TGMsgBox(gClient->GetRoot(),gClient->GetRoot(),
  //       "Look on the particle kind",
  //       "Please set the particle types according to the track",
  //       kMBIconQuestion,kMBYes|kMBNo,&retval);
  //if (retval == kMBNo) 
  // {
  //  return 0;
  // }
  //open dialog window to get the particle
  new MyTestDialog(gClient->GetRoot(), gClient->GetRoot(), 400, 450);
}
