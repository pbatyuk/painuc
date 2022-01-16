//  Reading of images from data file of image, filling data plot.

#include "TPostScript.h"
#include <dirent.h>
#include "../recon/_ROOT.h"
#include "../recon/_STD_UNIX.h"
#include "TPolyMarker.h"
#include "../recon/vertex.h"

#include <cfortran.h>
#include "minuit.h"

typedef unsigned* PUINT;
typedef PUINT* PPUINT;

PPUINT pCDBuffer[2];

void usage(char* argv[]);
int end_input();

PPUINT Alloc_Buffer(int hsize, int vsize);
void Free_Buffer(PPUINT buffer[2], int vsize);
   
void GetNewEvent();
void RemoveCreated();
Bool_t CheckOutFile();
void SetPart(void);
void Draw_region(int xpix, int ypix, int w);
    
void LifeFirst(int image, int iter, int ipass);
void LifeHard(int image, int iter, int ipass);
void LifeSoft(int image, int iter, int ipass);
void RenewPoints(int image);
void GroupPoints(int image, int iw, int ipass);
void SelectPoints(int image, int ipass);
void MakeSegments(int image);
void DeletePoints(int image, float aver, int iw);
void MakeMainSegments(int image, int ipass);
void SelectTracks(int version, int image, int ipass);
void find_track_sum_angle(int n, int k);
float D_to_Line(float xp, float yp, float Slope, float Y_Origin);
void EraseNearestToTrackPoints(int image, int Start, int End, int ipass, int what);

void FindVertex(int image, int ipass);
void LookForRestTrack(int image);
int circle_r(int N);
int circle_i(int N, int *Xp, int *Yp, float *Wp, float *A, float *B, float *R,
	     float *VAR);
void XY_Vertex_Net(int NTracks, float *a, float *b, float *r, 
		   float *Xvertex, float *Yvertex, 
		   float *Dvertex);
void FindPairsOnImages(void);
int  FindFirstPion(void);
void FindNewVertex(int image);
void DigitizeEvent(void);
void FindPointsForCircleFit(int image);
int PixelBr(int ima, int ypix, int xpix);

// for MINUIT usage
void fcn(int npar, double grad[3], double * fcnval,
                double xval[2],int iflag, void (*Dummy)());

FCALLSCSUB6(fcn,FCN,fcn,INT,DOUBLEV,PDOUBLE,DOUBLEV,INT,ROUTINE);

void circle_minuit (int N);
//		    float *A, float *B, float *R, float *pchi);

  //double DistMax(float *x, float *y);
#define OUT_filename "dubtoev.dat"
FILE *out_file;   /* output file handler */

TCanvas *c2,*ch;
TPad *pleft, *pright, *ptmp, *pMsg, *ptest;
TPad *hp0;

TH1F *htest,*hbr;
TH2F *left, *right;
gzFile *im_right, *im_left;
TGFileInfo fi;
TPaveLabel *pl, *pl1;
TObjectTable *fObjTbl;
TPolyMarker *markp, *marke;

Float_t *px_mark, *py_mark;

TStopwatch timer;

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
char iniDir[100] = {"/home/koval/DubTo/automaton/data/"};
//char iniDir[100] = {"/data/dubto/data99_1/cdrom/nov8/"};
//char iniDir[100] = {"/mnt/cdrom/nov12"};
//char iniDir[100] = {"/data/dubto/data99_1/cdrom/nov9"};
//char iniDir[100] = {"/home/gilpont/DubTo/Vera/031100/nov14"};
int maxl=0, minl=4096, maxr=0, minr=4096;
int xlmin=108, xlmax=500, ylmin=68, ylmax=443;
int xrmin=128, xrmax=520, yrmin=52, yrmax=414;
char WorkDir[100];
char eventname[100];
char ofilename[100], gifname[100];
char test_outfile[100];

//block for segments & tracks
#define NAT_LOG 2.718281828
#define  MAXNUMBEROFPOINTS 10000
int tailx[2][2], taily[2][2], iw;
float rad_tail[2][3][4];
short int sostn[8];
short int sosts[658][517], rsost[658][517];
short int flag[517][658];
float aver[2], main_aver[2];

float avermax;
int xrest[2][MAXNUMBEROFPOINTS], yrest[2][MAXNUMBEROFPOINTS];
int xmain[2][MAXNUMBEROFPOINTS], ymain[2][MAXNUMBEROFPOINTS];
int Npoint[2], Npoint_pass, Npoint_main, Not_Used;
float eff;
#define MAXNUMBEROFSEGMENTS 50000
#define MAXNUMBEROFTRACKS 20
int Xmain[2][MAXNUMBEROFTRACKS][100][500], Ymain[2][MAXNUMBEROFTRACKS][100][500];
int Nmain[2][MAXNUMBEROFTRACKS][100];
int NumberOfSegments;
short int StartOfSegment[MAXNUMBEROFSEGMENTS], EndOfSegment[MAXNUMBEROFSEGMENTS];
double cos_a[MAXNUMBEROFSEGMENTS], cos_b[MAXNUMBEROFSEGMENTS], CosTheta[100];
double MaxCosPhi;
int NumberOfTracks, Level, NumberOfTrackPoints;
int LengthOfBestTrack[MAXNUMBEROFTRACKS];
int CurrentTrack[MAXNUMBEROFSEGMENTS];
int BestTrack[MAXNUMBEROFTRACKS][1000];
short int LevelOfSegment[MAXNUMBEROFSEGMENTS];
double CosPhiMax;// = 0.9948;//0.9948; //10 degrees
//const double SinPhiMax = 0.0;/1745;//1 degrees
short int TrackPoint[MAXNUMBEROFPOINTS];
int X[2][MAXNUMBEROFTRACKS][1000], Y[2][MAXNUMBEROFTRACKS][1000];
int ntr[2], np[2][MAXNUMBEROFTRACKS];

float elx[1000], ely[1000];

float xcD[2][MAXNUMBEROFTRACKS], ycD[2][MAXNUMBEROFTRACKS], rcD[2][MAXNUMBEROFTRACKS];
float chi2_circle[2][MAXNUMBEROFTRACKS];
float X_Vert[2], Y_Vert[2], Min_Vertex_Error[2];
int mark[2][MAXNUMBEROFTRACKS];

int nomer_tr[2][MAXNUMBEROFTRACKS];
float angle[2][MAXNUMBEROFTRACKS];
float Angle[2][MAXNUMBEROFTRACKS];
int GoodEvent, NGoodTracks[2];

float xtrack[2][6][10000], ytrack[2][6][10000];
float xt[2][6][10000], yt[2][6][10000];
float wt[2][6][10000], xp[10000], yp[10000];;
float Wpoints[10000];
float Xcircle[10000], Ycircle[10000], Wpoint[10000];
float XC, YC, RC, chi2;

int NTracks[2], pair[MAXNUMBEROFTRACKS], Track[2][MAXNUMBEROFTRACKS];
int PionTrack[2];

// globals for MINUIT
int Npoints_fit;
float fcn_print;
Float_t x_mark[20000], y_mark[20000];

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
  labelfont = fClient->GetFontByName(gEnv->GetValue("Gui.BoldFont",
						    "-adobe-helvetica-bold-r-*-*-14-*-*-*-*-*-iso8859-1"));
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
  labelfont = fClient->GetFontByName(gEnv->GetValue("Gui.ProportionalFont",
						    "-adobe-courier-bold-r-*-*-16-*-*-*-*-*-iso8859-1"));
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

/*
Int_t mb_button_id[9] = { kMBYes, kMBNo, kMBOk, kMBApply,
                          kMBRetry, kMBIgnore, kMBCancel,
                          kMBClose, kMBDismiss };

EMsgBoxIcon mb_icon[4] = { kMBIconStop, kMBIconQuestion,
                           kMBIconExclamation, kMBIconAsterisk };
*/
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

   virtual void CloseWindow();
   virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
};

void SetNewPalette()
{
  TColor *newcolor = new TColor;
  int itmp;
  int ncolors = 37;
  Float_t light;
  Int_t Apalet[50];


  TStyle *palett = new TStyle; //->SetPalette(500);

  for ( itmp=0; itmp<50; itmp++)
    Apalet[itmp] = palett->GetColorPalette(itmp+1);

  for ( itmp=1; itmp<(ncolors+1); itmp++)
    {
      newcolor->SetNumber(itmp);
      light = (Float_t)itmp/ncolors;
      newcolor->SetRGB(light,light,light);
      Apalet[itmp-1]=itmp;
    }

  palett->SetPalette(ncolors+1,Apalet);
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
  palett->SetOptStat(0);
  //  stat->SaveStyle();


  //  left->SetFillColor(20);
  //  TStyle *mystyle=new Tstyle();
  palett->SetPadTopMargin(0.005);
  palett->SetPadBottomMargin(0.005);
  palett->SetPadLeftMargin(0.005);
  palett->SetPadRightMargin(0.005);

  palett->SetHistFillColor(38);
  palett->SetHistFillStyle(1001);

  palett->cd();
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

  event_kind = 1;
  num_tracks=1;
  start_image = 0;
  ipoint[1] = 0; ipoint[2] = 0;

  //  gObjectTable->Print();
}

WaitFrame *WaitMsg;


void GetNewEvent()
{
  char lfile[100], rfile[100];
  char *cptmp1, *cptmp2;
  char ctmp[10000];
  char fout_file[100];
  int maxl=0, minl=4096, maxr=0, minr=4096;

  fi.fFileTypes = (char **)filetypes;
  gSystem->ChangeDirectory(iniDir);
  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen,&fi);
  if (fi.fFilename == NULL) exit(0);
  //  strcpy(iniDir,gSystem->WorkingDirectory());
  strcpy(iniDir,gSystem->DirName(fi.fFilename));
  //  cout<<"just after"<<endl;

  strcpy(lfile,gSystem->BaseName(fi.fFilename));
  strcpy(rfile,lfile);

  if (isupper(rfile[0]))
    rfile[0] = 'R';
  else
    rfile[0] = 'r';

  //   cout<< "lfile = "<<lfile<<endl;
  //   cout<< "rfile = "<<rfile<<endl;

  cptmp1 = gSystem->ConcatFileName(gSystem->DirName(fi.fFilename),(char *)lfile);
  cptmp2 = gSystem->ConcatFileName(gSystem->DirName(fi.fFilename),(char *)rfile);
  //  cptmp1 = gSystem->ConcatFileName(gSystem->WorkingDirectory(),(char *)lfile);
  //  cptmp2 = gSystem->ConcatFileName(gSystem->WorkingDirectory(),(char *)rfile);

  WaitMsg = new WaitFrame(gClient->GetRoot(),300,70);

  cout << "start reading"<<endl;
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
    iitmp = ((int)cptmp2 - (int)cptmp1)/sizeof(char) + 1;
  
  iitmp += 1;     /*  removing first char from name  */
  iitmp *=sizeof(char);

  strcpy(eventname,fout_file + iitmp);
  //  fprintf(stderr,"\n event : %s \n",eventname);
  char dirout[10] = {"dig/"};
  strcpy(ofilename, dirout);
  int i;
  int k = 0;
  for (i = (int) strlen(dirout); 
       i<(int) strlen(dirout) + (int) strlen(eventname); i++)
    {ofilename[i] = eventname[k]; k++;}
  for (i=(int) strlen(dirout) + (int) strlen(eventname);i < 100;i++)
    ofilename[i] = '\0';
  strcat(ofilename,".dig");

  //  fprintf(stderr,"\n event filename : %s \n",ofilename);

  /*  put out Dummy Message Box to waiting  */
  //  cout<<"just before"<<endl;
  //  c2->Update();
  /*  Window_t waitid;
      waitid=WaitMsg->GetId();
      gVirtualX->MapRaised(waitid);
  */
  //  cout <<"just after msg"<<endl;
  

  
  /*  allocate memory buffer for images */

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
  unsigned maxBright = 0;
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
	  
	  istrstream Inbuffer(ctmp,10000);
	  /*    j - x coordinate ;  i - y coordinate   */ 
	  for ( j=0; j<hsize; j++ )
	    {
	      Inbuffer >> Bright;
	      //	      printf(" i,j = %i \n",Bright);
	      pCDBuffer[image][i][j] = Bright;
	      if ( Bright > maxBright ) 
		maxBright = Bright;
	      //	      getc(stdin);
	    }
	}
      if (image == 0) gzclose(im_left);
      else gzclose(im_right);
    }
  cout<< "end of read"<<endl; // ^test
  //  left = new TH2F("left",lfile, 
  left = new TH2F("left","bitmap for left image", 
           (Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);

  //  right = new TH2F("right",rfile,
  right = new TH2F("right","bitmap for right image",
	   (Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);

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
		      htest->Fill(ww);
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
  //  getc(stdin);
  /*
  free(pCDBuffer[0]);
  free(pCDBuffer[1]);
  */
  Free_Buffer(pCDBuffer,vsize);

  pleft->cd();
  left->Draw("col");
  pright->cd();
  right->Draw("col");

  /*
  left->SetMaximum(left->GetMaximum());
  left->SetMinimum(left->GetMinimum());
  right->SetMaximum(right->GetMaximum());
  right->SetMinimum(right->GetMinimum());
  */
  left->SetMaximum(maxl);
  left->SetMinimum(minl);
  right->SetMaximum(maxr);
  right->SetMinimum(minr);

  pMsg->cd();
  pl1 = new TPaveLabel(0.07,0.15,0.93,0.45,eventname,"br");
  pl1->SetFillColor(45);
  pl1->SetTextSize(0.5);
  pl1->Draw();

  c2->Update();

  //  cout << "before close"<<endl;
  WaitMsg->CloseWaitFrame();

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
  
}



int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);

  int length = 0;
  char ctmp[10000];
  char lfile[100], rfile[100];
  char fout_file[100];
  char *cptmp1, *cptmp2;

  SetNewPalette();

  tailx[0][0] = 30;//left image left margin
  tailx[0][1] = 100;//left image right margin
  tailx[1][0] = 50;//right image left margin
  tailx[1][1] = 80;//right image right margin
  taily[0][0] = 105;//left image bottom margin
  taily[0][1] = 55;//left image top margin
  taily[1][0] = 50;//right image bottom margin
  taily[1][1] = 115;//right image top margin
  
  //left image
  //left top corner x, y, r
  rad_tail[0][0][0] = 370.64; rad_tail[0][1][0] = 150.22; rad_tail[0][2][0] = pow(399.54 - 20., 2);
  // right top corner
  rad_tail[0][0][1] = 337.67; rad_tail[0][1][1] = 247.34; rad_tail[0][2][1] = pow(257.52 - 20., 2);
  // right bottom corner
  rad_tail[0][0][2] = 287.39; rad_tail[0][1][2] = 298.12; rad_tail[0][2][2] = pow(309.70 - 20., 2);
  // left bottom corner
  rad_tail[0][0][3] = 381.53; rad_tail[0][1][3] = 391.64; rad_tail[0][2][3] = pow(422.80 - 20., 2);
  //right image
  //left top corner x, y, r
  rad_tail[1][0][0] = 535.84; rad_tail[1][1][0] = 153.91; rad_tail[1][2][0] = pow(502.97 - 20., 2);
  // right top corner
  rad_tail[1][0][1] = 185.84; rad_tail[1][1][1] = 92.71; rad_tail[1][2][1] = pow(468.80 - 20., 2);
  // right bottom corner
  rad_tail[1][0][2] =-121.51; rad_tail[1][1][2] = 622.66; rad_tail[1][2][2] = pow(867.50 - 20., 2);
  // left bottom corner
  rad_tail[1][0][3] = 318.02; rad_tail[1][1][3] = 268.96; rad_tail[1][2][3] = pow(293.52 - 20., 2);

  fObjTbl = new TObjectTable(10000);
  fi.fFileTypes = (char **)filetypes;

  strcpy(WorkDir,gSystem->WorkingDirectory());

  gSystem->ChangeDirectory(iniDir);
       /* print files in current directory in reverse order */
  struct dirent **namelist;
  int n;
  
  n = scandir(".", &namelist, 0, alphasort);
  /*
  if (n < 0)
    perror("scandir");
  else
    while(n--) 
      printf("%s\n", namelist[n]->d_name);
  */
  //new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDOpen,&fi);
  for(int iname = 0; iname <n; iname++)
    {
      length = 0;
      timer.Reset();
      timer.Start(kTRUE);
      fi.fFilename = namelist[iname]->d_name;
      printf("%s \n", fi.fFilename);
      strcpy(lfile,gSystem->BaseName(fi.fFilename));
      if(lfile[0] != 'L' && lfile[0] != 'l') continue;

      
      /*    forming whole filenames islower  */
      strcpy(lfile,gSystem->BaseName(fi.fFilename));
      strcpy(rfile,lfile);
      
      if (isupper(rfile[0]))
	rfile[0] = 'R';
      else
	rfile[0] = 'r';

      cout<< "lfile = "<<lfile<<endl;
      cout<< "rfile = "<<rfile<<endl;

      cptmp1 = gSystem->ConcatFileName(iniDir,(char *)lfile);
      cptmp2 = gSystem->ConcatFileName(iniDir,(char *)rfile);
      //  cptmp1 = lfile;
      //  cptmp2 = rfile;
      cout<< "clfile = "<<cptmp1<<endl;
      cout<< "crfile = "<<cptmp2<<endl;
      //  exit(1);
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
      
      int iitmp = 0; 
      cptmp1 = strcpy(fout_file,cptmp1);
  
      if ((cptmp2 = strrchr(fout_file,'.')) != NULL) // looking for "."
	*cptmp2 = '\0';
      
      if ((cptmp2 = strrchr(fout_file,'/')) != NULL) 
	iitmp = ((int)cptmp2 - (int)cptmp1)/sizeof(char) + 1;
      
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
      istrstream Inbuffer(ctmp,10000);
      
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
	      
	      istrstream Inbuffer(ctmp,10000);
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
      
      left = new TH2F("left","bitmap for left image", 
		      (Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);
      
      right = new TH2F("right","bitmap for right image",
		       (Int_t) hsize,0,(Axis_t) hsize,(Int_t) vsize,0,(Axis_t) vsize);
      
      Axis_t xx,yy;
      Stat_t ww;
      float sumw;
      int inum;
      for (image=0; image<2; image++)
	{
	  sumw = 0.;
	  inum = 0;
	  for (i=taily[image][1]; i<vsize-taily[image][0]-1; i++)
	    {
	      for (j=tailx[image][0]; j<hsize-tailx[image][1]; j++)
		{
		  sumw += pCDBuffer[image][i][j];
		  inum++;
		}
	    }
	  main_aver[image] = 0.;
	  if(inum !=0 )main_aver[image] = sumw/(float) inum;
	  for (i=0; i<vsize; i++)
	    {
	      for (j=0; j<hsize; j++)
		{
		  ww = (Stat_t) pCDBuffer[image][i][j];
		  xx = (Axis_t) j;//(j + 0.5);
		  yy = (Axis_t) (vsize - i - 1);//0.5);  
		  /* histogram filling in oppozitive*/ 
		  //   cout << "brightness ="<<pCDBuffer[image][i][j]<<" weight = "
		  //		     << ww << endl;
		  //	  getc(stdin);		
		  if (image == 0) 
		    {
		      left->Fill(xx,yy,ww);
		    }
		  else
		    {
		      right->Fill(xx,yy,ww);
		    }
		}
	    }
	}
      
      
      cout << "Real time for reading ="<<timer.RealTime()<< endl;
      cout << "CPU time for reading ="<<timer.CpuTime()<< endl;
      
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
      
      cout <<"x , y ="<<ixWindow<<" "<<iyWindow<<endl;

      c2 = new TCanvas("c2","Event digitizing",10,10,ixWindow,iyWindow);
      pleft = new TPad("pleft","Left image",    xleft, yleft, xright, yright, 21);
      pright = new TPad("pright","Right image", xleft + 0.5, yleft, xright + 0.5, yright, 21);
      ptmp = new TPad("ptmp"," ",xlefttmp, ylefttmp, xrighttmp, yrighttmp, 21);
      
      pMsg = new TPad("pMsg"," ",0.01,0.01,xlefttmp-0.01,yrighttmp,28);
      
      pleft->SetEditable(kFALSE);
      pright->SetEditable(kFALSE);
      ptmp->SetEditable(kFALSE);
      pMsg->SetEditable(kFALSE);
      
      pleft->Draw();
      pright->Draw();
      pleft->SetLogz();
      pright->SetLogz();
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
      
      pleft->cd();
      left->Draw("col");
      pright->cd();
      right->Draw("col");
      
      c2->Modified(); 
      c2->Update();

      //new canvas for GIF file
      strcpy(gifname, eventname);
      strcat(gifname, ".gif");
      c2->Print(gifname,"gif");
      
      int column, row;
      for (image=0; image<2; image++)
	{
	  if(image == 0) 
	    pleft->cd();
	  else
	    pright->cd();
	  int ipass = 0;
	  avermax = 4000.;
	LABEL99:
	  Npoint[image] = 0;
	  ntr[image] = 0;
	  for (column=0; column<hsize; column++)
	    {
	      for (row=0; row<vsize; row++)
		{
		  rsost[column][row] = 0;
		  flag[row][column] = 1;
		}
	    }
	LABEL100:
	  ipass++;
	  
	  //if (ipass <= 1)
	  //LifeSoft(image, niter, ipass+1);
	  //else 
	  LifeFirst(image, 1, ipass);
	  //else if (ipass == 2)
	  //LifeHard(image, 3, ipass);
	  
	  sumw = 0.;
	  inum = 0;
	  int npoint = Npoint[image];
	  if(ipass == 1) iw = 0;
	  if(ipass > 1) iw = 0;
	  int maxbr;
	  
	  for (column=0; column<hsize; column++)
	    {
	      for (row=0; row<vsize; row++)
		{
		  rsost[column][row] = 0;
		  if(sosts[column][row] == 1) flag[row][column] = 0;
		}
	    }
	  
	  if(iw != 0) 
	    {       
	      for (column=tailx[image][0]+iw; column<hsize-tailx[image][1]-iw; column += 2*iw)
		{
		  for (row=taily[image][1]+iw; row<vsize-taily[image][0]-iw-1; row+=2*iw)
		    {
		      if ((row>taily[image][1]+iw && row<vsize-taily[image][0]-iw-1) &&
			  (column>tailx[image][0]+iw && column<hsize-tailx[image][1]-iw))
			{
			  maxbr = 0;
			  for (i = column - iw; i< column + iw; i++)
			    {
			      for (j = row - iw; j < row + iw; j++)
				{
				  if (sosts[i][j] == 1 && ((int) pCDBuffer[image][j][i] > maxbr)
				      && rsost[i][j] == 0)
				    {
				      rsost[i][j] = 1;
				      
				      maxbr = pCDBuffer[image][j][i];
				      xrest[image][npoint] = i;
				      yrest[image][npoint] = vsize - j - 1;  
				    }
				}
			    }
			  if (maxbr > 0) npoint ++;
			}
		    }
		}
	    }
	  else
	    {
	      for (column=tailx[image][0]; column<hsize-tailx[image][1]; column++)
		{
		  for (row=taily[image][1]; row<vsize-taily[image][0]-1; row++)
		    {
		      if (sosts[column][row] == 1)
			{
			  if ((row>taily[image][1] && row<vsize-taily[image][0]-1) &&
			      (column>tailx[image][0] && column<hsize-tailx[image][1]))
			    {
			      int corner;
			      float rho;
			      if (column < hsize/2)
				{
				  if (vsize - row -1 > vsize/2)  // left top corner
				    corner = 0;
				  else //left bottom corner
				    corner = 3;
				}
			      else
				{
				  if (vsize - row - 1 > vsize/2)  // right top corner
				    corner = 1;
				  else //right bottom corner
				    corner = 2;
				}
			      rho = pow((float)column - rad_tail[image][0][corner], 2) +
				pow((float)(vsize - row - 1) - rad_tail[image][1][corner], 2);
			      if(rho < rad_tail[image][2][corner]) 
				{
				  rsost[column][row] = 1;
				  
				  xrest[image][npoint] = column;
				  yrest[image][npoint] = vsize - row - 1; 
				  if(vsize - row - 1 < 0 || vsize - row - 1 > vsize)
				    printf("column %d row %d npoint %d\n", column, row, npoint);
				  npoint ++; 
				}
			    }
			}
		    }
		}
	    }
	  
	  Npoint_pass = npoint;
	  printf(" %d points are on the %d image\n", Npoint_pass-Npoint[image], image);
	  //MakeSegments(image);
	  Npoint[image] = Npoint_pass;
	  //if(ipass == 1)SelectPoints(image, ipass);
	  
	  //draw points and track
	  for(int imark=0; imark<Npoint[image]; imark++)
	    {
	      x_mark[imark] = xrest[image][imark];
	      y_mark[imark] = yrest[image][imark];
	    }
 	  px_mark = &x_mark[0];
	  py_mark = &y_mark[0];
	  markp = new TPolyMarker((Int_t) Npoint[image], px_mark, py_mark, "same");
	  fObjTbl->Add(markp);
	  markp->SetMarkerStyle(20);
	  markp->SetMarkerSize(0.2);
	  markp->SetMarkerColor(4);
	  markp->Draw("same");	     
	  c2->Modified(); 
	  c2->Update();
	  //RenewPoints(image);
	  //Npoint_main = Npoint[image];
	  GroupPoints(image, iw, ipass);
	  MakeMainSegments(image, ipass);
	  SelectTracks(1, image, ipass);
	  printf(" Number of tracks %d\n", NumberOfTracks);
	  for (int itr = 0; itr < NumberOfTracks; itr++)
	    {
	      int jtr = itr;
	      j = 0;
	      for (k = 0; k<LengthOfBestTrack[jtr]; k++)
		{
		  i = BestTrack[jtr][k];
		  X[image][itr][j] = xmain[image][StartOfSegment[i]];
		  Y[image][itr][j] = ymain[image][StartOfSegment[i]];
		  j++;
		}
	      k = LengthOfBestTrack[jtr] - 1;
	      i = BestTrack[jtr][k];
	      X[image][itr][j] = xmain[image][EndOfSegment[i]];
	      Y[image][itr][j] = ymain[image][EndOfSegment[i]];
	      j++;
	      np[image][itr] = j;
	    }

	  for (column=0; column<hsize; column++)
	    for (row=0; row<vsize; row++)
	      rsost[column][row] = 0;
	  FindPointsForCircleFit(image);
	  ntr[image] = NumberOfTracks;
	  //draw marked points
	  i = 0;
	  for (column=0; column<hsize; column++)
	    for (row=0; row<vsize; row++)
	      if(rsost[column][row] == 1)
		{
		  x_mark[i] = column;
		  y_mark[i] = vsize - row - 1;
		  i++;
		}
	  printf(" Original points used %d\n", i);
 	  px_mark = &x_mark[0];
	  py_mark = &y_mark[0];
	  markp = new TPolyMarker((Int_t) i, px_mark, py_mark, "same");
	  fObjTbl->Add(markp);
	  markp->SetMarkerStyle(8);
	  markp->SetMarkerSize(0.2);
	  markp->SetMarkerColor(38);
	  if(image == 0) pleft->cd();
	  if(image == 1) pright ->cd();
	  markp->Draw("same");	     
	  
	  Min_Vertex_Error[image] = 10000.;
	  //if(ntr[image] >= 2) 
	  FindVertex(image, ipass);

	  //if(ipass == 1 && 
	  // ((Min_Vertex_Error[image] > 12. && eff > 0.3) || ntr[image] <= 2)) goto LABEL99;
	  //draw points and track
	  i = 0;
	  for(int imark=0; imark<Npoint_main; imark++)
	    {
	      if(TrackPoint[imark] == 0)
		{
		  x_mark[i] = xmain[image][imark];
		  y_mark[i] = ymain[image][imark];
		  i++;
		}
	    }
 	  px_mark = &x_mark[0];
	  py_mark = &y_mark[0];
	  markp = new TPolyMarker((Int_t) i, px_mark, py_mark, "same");
	  fObjTbl->Add(markp);
	  markp->SetMarkerStyle(8);
	  markp->SetMarkerSize(0.25);
	  markp->SetMarkerColor(43);
	  markp->Draw("same");	     
	  i = 0;
	  for(int imark=0; imark<Npoint_main; imark++)
	    {
	      if(TrackPoint[imark] == 1)
		{
		  x_mark[i] = xmain[image][imark];
		  y_mark[i] = ymain[image][imark];
		  i++;
		}
	    }
 	  px_mark = &x_mark[0];
	  py_mark = &y_mark[0];
	  markp = new TPolyMarker((Int_t) i, px_mark, py_mark, "same");
	  fObjTbl->Add(markp);
	  markp->SetMarkerStyle(8);
	  markp->SetMarkerSize(0.25);
	  markp->SetMarkerColor(50);
	  markp->Draw("same");	     

	  LookForRestTrack(image);

	  //draw vertex position
	  float xx = X_Vert[image];
	  float yy = Y_Vert[image];
	  TArc *myarc = new TArc(xx,yy, Min_Vertex_Error[image], 0.,360.);
	  myarc->SetLineColor(41);
	  myarc->Draw("same");
	  fObjTbl->Add(myarc);
	  TMarker *marke = new TMarker(xx, yy, 5);
	  fObjTbl->Add(marke);
	  marke->SetMarkerColor(41);
	  marke->SetMarkerStyle(20);
	  marke->SetMarkerSize(0.3);
	  marke->DrawMarker(xx,yy);
	  //draw tracks associated to the vertex
	  for (int itr=0; itr < ntr[image]; itr++)
	    {
	      //printf("Track %d Points %d Mark %d\n", itr, np[image][itr], mark[image][itr]);
	      if(np[image][itr] ==0) continue; //check for excluded close tracks
	      if(mark[image][itr] != 1 ) continue;
	      for (i=0; i<np[image][itr]; i++)
		{
		  elx[i] = (float) X[image][itr][i];
		  ely[i] = (float) Y[image][itr][i];
		  //printf("elx[%d] %f ely[%d] %f\n", i, elx[i], i, ely[i]);
		}
	      float *elxp = &elx[0];
	      float *elyp = &ely[0];
	      poly1 = new TPolyLine(np[image][itr], elxp, elyp, 
				    "same");
	      poly1->SetLineColor(50+itr);
	      poly1->SetLineWidth(3.);
	      poly1->Draw("same");
	    }
	  c2->Modified(); 
	  c2->Update();
	  
	}

      FindPairsOnImages();
      NGoodTracks[0] = NGoodTracks[1] = 0;
      for(image=0;image<2; image++)
	for (int itr=0; itr < ntr[image]; itr++)
	  if(mark[image][itr] == 1) NGoodTracks[image]++;
      printf("Final number of tracks found is %d %d \n", NGoodTracks[0], NGoodTracks[1]);
      GoodEvent = 2;
      if(NGoodTracks[0] == NGoodTracks[1] && NGoodTracks[0] > 1) 
	{
	  for(image=0; image<2; image++) FindNewVertex(image);
	  GoodEvent = 3;
	  if(NGoodTracks[0] > 2) GoodEvent = 4;
 	  if(FindFirstPion() == 1) GoodEvent = 5;
	}
      float dx_Vert = fabs((X_Vert[1] - 24.) - (X_Vert[0]));
      float dy_Vert = fabs(Y_Vert[1] + 30. - Y_Vert[0]);
      printf("Difference between Vertexes in X: %f in Y %f\n",
	     dx_Vert, dy_Vert);
      if(dx_Vert > 30. || dy_Vert > 30.) 
	{
	  GoodEvent = 1;
	  for(image=0;image<2; image++)
	    for (int itr=0; itr < ntr[image]; itr++)
	      mark[image][itr] = 0;
	}
      printf("Final Event quality is : %d\n", GoodEvent);
      
      //draw tracks pair by pair
      int jtr;
      for(image=0;image<2; image++)
	{
	  jtr = 0;
	  if(image == 0) 
	    pleft->cd();
	  else
	    pright->cd();
	  //draw new vertex if found
	  float xx = X_Vert[image];
	  float yy = Y_Vert[image];
	  TArc *myarc = new TArc(xx,yy, Min_Vertex_Error[image], 0.,360.);
	  myarc->SetLineColor(42);
	  myarc->Draw("same");
	  fObjTbl->Add(myarc);
	  TMarker *marke = new TMarker(xx, yy, 5);
	  fObjTbl->Add(marke);
	  marke->SetMarkerColor(41);
	  marke->SetMarkerStyle(20);
	  marke->SetMarkerSize(0.3);
	  marke->DrawMarker(xx,yy);
	  //draw tracks associated to the vertex and in good coincedence with another image
	  for (int itr=0; itr < ntr[image]; itr++)
	    {
	      printf("Track %d Points %d Mark %d\n", itr, np[image][itr], mark[image][itr]);
	      for (i=0; i<np[image][itr]; i++)
		{
		  elx[i] = (float) X[image][itr][i];
		  ely[i] = (float) Y[image][itr][i];
		  //printf("elx[%d] %f ely[%d] %f\n", i, elx[i], i, ely[i]);
		}
	      float *elxp = &elx[0];
	      float *elyp = &ely[0];
	      poly1 = new TPolyLine(np[image][itr], elxp, elyp, 
				    "same");
	      if(mark[image][itr] != 1 ) poly1->SetLineColor(5);
	      if(mark[image][itr] == 1 ) {poly1->SetLineColor(50+jtr);jtr++;}
	      poly1->SetLineWidth(3.);
	      poly1->Draw("same");
	    }
	  c2->Modified(); 
	  c2->Update();
	  
	}

      char ctmp[1];
      //final canvas for GIF file
      strcpy(gifname, eventname);
      sprintf(ctmp, "%u", GoodEvent);
      strcat(gifname, ".");
      strcat(gifname, ctmp);
      strcat(gifname, ".gif");
      c2->Print(gifname,"gif");

      if(GoodEvent > 1) DigitizeEvent();
      
      Free_Buffer(pCDBuffer,vsize);
      RemoveCreated();
      //int buf1;
      //cin>>buf1;
    }

  maiWin = new MyChoiseFrame(gClient->GetRoot(), 200, 380);

  cout << "Real time after drawing ="<<timer.RealTime()<< endl;
  cout << "CPU time after drawing ="<<timer.CpuTime()<< endl;
  printf("There was running event number %s\n\a", eventname);

  //gGXW->Bell(100);
  //gClient->Bell(100);
  theApp.Run(kTRUE);
 
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

void Draw_region(int xpix, int ypix, int w)
{
  if(w > 0)
    {
      TMarker *mark = new TMarker(xpix,ypix,7);
      fObjTbl->Add(mark);
      mark->SetMarkerStyle(8);
      mark->SetMarkerSize(0.2);
      mark->SetMarkerColor(41);
      mark->Draw("same");	     
    }
  else
    {
      TMarker *mark = new TMarker(xpix,ypix,7);
      fObjTbl->Add(mark);
      mark->SetMarkerStyle(8);
      mark->SetMarkerSize(0.2);
      mark->SetMarkerColor(40);
      mark->Draw("same");	     
    }
}
void MakeSegments(int image)
{
  int Limit0, Limit1;
  NumberOfSegments = 0;
  Limit0 = Npoint[image]; 
  Limit1 = TMath::Min(MAXNUMBEROFPOINTS, Npoint_pass);
  for (int i = Limit0; i < Limit1; i++)
    {
      TrackPoint[i] = 0;
      for(int j = Limit0; j < Limit1; j++)
	{
	  if(i==j) continue;
	  StartOfSegment[NumberOfSegments] = i;
	  EndOfSegment[NumberOfSegments]   = j;
	  double dx = (float) xrest[image][j]-xrest[image][i];
	  double dy = (float) yrest[image][j]-yrest[image][i];
	  double dist = TMath::Sqrt(dx*dx + dy*dy);
	  if (dist < 150. )//&& grad <0.5)//more than 10 pixels!!
	    {
	      //rotate on 5 degree
	      /*
	      double xrj = xrest[image][j]*cosphi + yrest[image][j]*sinphi;
	      double yrj =-xrest[image][j]*sinphi + yrest[image][j]*cosphi;
	      double xri = xrest[image][i]*cosphi + yrest[image][i]*sinphi;
	      double yri =-xrest[image][i]*sinphi + yrest[image][i]*cosphi;
	      dx = (float) xrj-xri;
	      dy = (float) yrj-yri;
	      dist = TMath::Sqrt(dx*dx + dy*dy);
	      if(dx != 0.) tg[NumberOfSegments] = dy/dx;
	      if(dx == 0.) {
	      	tg[NumberOfSegments] = 0.;
	      	printf("There was  zero dx!\n");
	      }
	      */
	      cos_a[NumberOfSegments] = dx/dist;
	      cos_b[NumberOfSegments] = dy/dist;
	      if(MAXNUMBEROFSEGMENTS>NumberOfSegments) 
		NumberOfSegments = NumberOfSegments + 1;
	    }
	}
    }
  TrackPoint[Limit1-1] = 0;
  printf("Segment number is %d\n", NumberOfSegments);
}
void MakeMainSegments(int image, int ipass)
{
  int Limit0, Limit1;
  Limit0 = 0; 
  Limit1 = TMath::Min(MAXNUMBEROFPOINTS, Npoint_main);
  NumberOfSegments = 0;
  for (int i = Limit0; i < Limit1; i++)
    {
      TrackPoint[i] = 0;
      for(int j = Limit0; j < Limit1; j++)
	{
	  if(i==j) continue;
	  StartOfSegment[NumberOfSegments] = i;
	  EndOfSegment[NumberOfSegments]   = j;
	  double dx = (float) xmain[image][j]-xmain[image][i];
	  double dy = (float) ymain[image][j]-ymain[image][i];
	  double dist = TMath::Sqrt(dx*dx + dy*dy);
	  if (dist>1.1 && 
	      ((dist < 60. && ipass == 1) || 
	       (dist < 80. && ipass ==-1) || 
	       (dist < 80. && ipass == 2)))// && grad <0.5)//more than 10 pixels!!
	    //if (dist < 20 )//more than 10 pixels!!
	    {
	      cos_a[NumberOfSegments] = dx/dist;
	      cos_b[NumberOfSegments] = dy/dist;
	      if(MAXNUMBEROFSEGMENTS>NumberOfSegments) 
		NumberOfSegments = NumberOfSegments + 1;
	    }
	}
    }
  TrackPoint[Limit1-1] = 0;
  printf("Main Segment number is %d\n", NumberOfSegments);
}

void SelectTracks(int version, int image, int ipass)
{
  int MinLevelLong = 6;
  int maxnumberoftracks = 4;
  CosPhiMax = 0.9948;//0.9948; //10 degrees
  if(version == 1)
    {
      if(ipass == 1) {
	CosPhiMax = 0.9950;//48;//48;//0.9948; //10 degrees
	MinLevelLong = 5;
      }
      if(ipass > 1) MinLevelLong = 5;
      maxnumberoftracks = 10;
      if(ipass == -1) {
	maxnumberoftracks = 3;
	MinLevelLong = 3;
	CosPhiMax = 0.9948;//48;//48;
      }

    }

  short int HasNeighbour[MAXNUMBEROFSEGMENTS][2];
  int LivingSegments;
  int i, j, k, n;
  int MinLevel; 

  NumberOfTracks = 0;
  MinLevel = MinLevelLong;

  Level = MinLevel;
  while(Level >= MinLevel)
    {
      for(i = 0; i< NumberOfSegments; i++)
	{
	  for(j = 0; j < NumberOfSegments; j++)
	    {
	      if (EndOfSegment[i] == StartOfSegment[j]) 
		{
                  double CosPhi = cos_a[i]*cos_a[j]+cos_b[i]*cos_b[j];
		  // 		  double SinPhi = TMath::Sqrt(1. - CosPhi*CosPhi);
		  //if(cos_b[j]/cos_a[j] - cos_b[i]/cos_a[i] < 0. )
		  // SinPhi = -SinPhi;
		  //if(SinPhi > -SinPhiMax) //take only greater than -1degree angle 
		  if  (CosPhi > CosPhiMax) 
		    {
		      HasNeighbour[i][0] = 1;
		      HasNeighbour[j][0] = 1;
		      LevelOfSegment[i] = 1;
		      LevelOfSegment[j] = 1;
		    }
		}
	    }
	}

      for(i = 0; i< NumberOfSegments; i++)
	{
	  if (TrackPoint[StartOfSegment[i]]==1 ||
	      TrackPoint[EndOfSegment[i]]==1 )
	    {
	      HasNeighbour[i][0] = 0;
	      LevelOfSegment[i] = 0;
	    }
	}

      Level = 0;

      LivingSegments = 1;
      while (LivingSegments > 0)
	{
	  Level++;
	  LivingSegments = 0;

	  for (i = 0; i < NumberOfSegments; i++)
	    {
	      HasNeighbour[i][1] = 0;
	      if (HasNeighbour[i][0]) 
		{
		  for(j = 0; j < NumberOfSegments; j++)
		    {
		      if (EndOfSegment[i] == StartOfSegment[j]) 
			{
			  if (HasNeighbour[j][0]) 
			    {
			      double CosPhi = cos_a[i]*cos_a[j]+cos_b[i]*cos_b[j];
			      //double SinPhi = TMath::Sqrt(1. - CosPhi*CosPhi);
			      //if(cos_b[j]/cos_a[j] - cos_b[i]/cos_a[i] < 0. )
			      //		SinPhi = -SinPhi;
			      //if(SinPhi > -SinPhiMax) //take only greater than -1degree angle 
			      if (CosPhi > CosPhiMax) 
				{
				  HasNeighbour[i][1] = 1;
				  LevelOfSegment[i] = Level + 1;
				  LivingSegments = LivingSegments + 1;
				}
			    }
                        }
		    }
		}
	    }
	  
	  for (i = 0; i < NumberOfSegments; i++)
	    HasNeighbour[i][0] = HasNeighbour[i][1];
	}
      printf("Level is %d\n", Level);
      if (Level >= MinLevel)
	{
	  if (NumberOfTracks >= maxnumberoftracks) return;

	  MaxCosPhi = 0.0;
	  LengthOfBestTrack[NumberOfTracks] = 0;


	  for (k = 0; k< NumberOfSegments; k++)
	    {
	      if (LevelOfSegment[k] == Level) 
		{
                  n = 1;
                  CurrentTrack[n-1] = k;

                  find_track_sum_angle(n, k);
		}
	    }
	  printf("LengthOfBestTrack[NumberOfTracks] %d\n",  
		 LengthOfBestTrack[NumberOfTracks]);
	  //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  //     "delete" track segments 
	  int j_s;
	  for (j = 0; j < LengthOfBestTrack[NumberOfTracks]; j++)
	    { 
	      j_s = BestTrack[NumberOfTracks][j];
	      TrackPoint[StartOfSegment[j_s]] = 1;
	      //delete all points along segments on track except last and first segment!!!
	      //if(j != 0 && j < LengthOfBestTrack[NumberOfTracks]- 1) 
		EraseNearestToTrackPoints(image, StartOfSegment[j_s], EndOfSegment[j_s], ipass, 1);
	      //if(j == 0) TrackPoint[StartOfSegment[j_s]] = 0;
	    }
	  
	  j_s = BestTrack[NumberOfTracks][LengthOfBestTrack[NumberOfTracks]-1];
	  TrackPoint[EndOfSegment[j_s]] = 1;
	  
	  NumberOfTracks = NumberOfTracks + 1;
	  
	}//if    Level > MinLevel
    }//while Level > MinLeve

  //if (MinLevel == MinLevelLong) 
  //  {
  //   MinLevel = MinLevelShort;
  //  goto LABEL100;
  //}
  return;
}
void find_track_sum_angle(int n, int k)
{
  int m, l;
  int i;
  double  SumCosPhi;

  if (Level == n)
    {
      SumCosPhi = 0.0;
      for( i = 0; i < n; i++)
	SumCosPhi += CosTheta[i];
      if ((NumberOfTrackPoints > LengthOfBestTrack[NumberOfTracks]) ||
	  ((NumberOfTrackPoints==LengthOfBestTrack[NumberOfTracks])&&
	   (SumCosPhi/(double) NumberOfTrackPoints > MaxCosPhi)))
	{
	  LengthOfBestTrack[NumberOfTracks] =
	    NumberOfTrackPoints;
	  MaxCosPhi = SumCosPhi/(double) NumberOfTrackPoints;
	  for(i = 0; i < LengthOfBestTrack[NumberOfTracks]; i++)
	    BestTrack[NumberOfTracks][i] = CurrentTrack[i];
	}
      return;
    }


  if (Level > n)
    {
      m = n + 1;
      for (l = 0; l< NumberOfSegments; l++)
	{
	  if (LevelOfSegment[l]==(Level-n))
	    {
              if(EndOfSegment[k]==StartOfSegment[l])
		{
		  CosTheta[m]=cos_a[k]*cos_a[l]+
		    cos_b[k]*cos_b[l];
 		  //double SinPhi = TMath::Sqrt(1. - CosTheta[m]*CosTheta[m]);
		  //if(cos_b[l]/cos_a[l] - cos_b[k]/cos_a[k] < 0. )
		  //SinPhi = -SinPhi;
		  //if(SinPhi > -SinPhiMax) //take only greater than -1degree angle 
		  if (CosTheta[m] > CosPhiMax)
		    {
		      NumberOfTrackPoints = m;
		      CurrentTrack[m-1] = l;
		      find_track_sum_angle(m, l);
		    }
		}
	    }
	}
    }
  return;
}
//////////////////end of automaton part////////////////////
int PixelBr(int ima, int ypix, int xpix)
{
  int bright;
  Int_t bin_arr, binx, biny;

  TAxis *xaxis;
  TAxis *yaxis;
  if(ima == 0) 
    {
      xaxis = left->GetXaxis();
      yaxis = left->GetYaxis();
    }
  else
    {
      xaxis = right->GetXaxis();
      yaxis = right->GetYaxis();
    }
  binx = xaxis->FindBin((Float_t) xpix + 0.5);
  biny = yaxis->FindBin((Float_t) ypix + 0.5);
  //biny = yaxis->FindBin((Float_t) vsize - ypix - 0.5);
  
  //binx = (Int_t) xpix;
  //biny = (Int_t) ypix;
  if ( ima == 0)
    {
      bin_arr = left->GetBin(binx,biny);
      bright = (int) left->GetBinContent(bin_arr);
    }
  else
    {
      bin_arr = right->GetBin(binx,biny);
      bright = (int) right->GetBinContent(bin_arr);
    }
  return bright;
}
void RenewPoints(int image)
{
  int i, j, ip, jp, ii, jj;
  int off[MAXNUMBEROFPOINTS];
  int imain;
  for(ip=0; ip<Npoint[image]; ip++)
    off[ip] = 0;
  for(ip=0; ip<Npoint[image]-1; ip++)
    {
      i = xrest[image][ip];
      j = vsize - yrest[image][ip] - 1;
      for(jp=ip+1; jp<Npoint[image]; jp++)
	{
	  ii = xrest[image][jp];
	  jj = vsize - yrest[image][jp] - 1;
	  if(i == ii && i == jj) off[jp] = 1;
	}
    }
  imain = 0;
  for(ip=0; ip<Npoint[image]; ip++)
    {
      i = xrest[image][ip];
      j = vsize - yrest[image][ip] - 1;
      if(off[ip] == 1) continue;
      xmain[image][imain] = xrest[image][ip];
      ymain[image][imain] = yrest[image][ip];
      imain++;
    }
  Npoint_main = imain;
  printf(" %d MAIN points are on the %d image\n", Npoint_main, image);
}
void DeletePoints(int image, float meanbr, int iwin)
{
  int j;
  float summax = 0.;
  int nmax = 0;
  if (ntr[image] > 0)
    for (int itr = 0; itr < ntr[image]; itr++)
      {
	if(np[image][itr] == 0) continue;
	for (j = 0; j < np[image][itr]; j++)
	  {
	    int xtr = X[image][itr][j];
	    int ytr = vsize - Y[image][itr][j] - 1;
	    float br = (float) pCDBuffer[image][ytr][xtr];
	    if(flag[ytr][xtr] != 0) {
	      summax += br;
	      nmax++;
	    }
	    flag[ytr][xtr] = 0;
	    float stepr = 1.0;
	    float stepf = 0.1;
	    float sf;
	    float brprev;
	    float sr, bc; int xc, yc;
	    for (sf = 0; sf < 360; sf+=stepf)
	      {
		brprev = br;
		for (sr = stepr; sr < iwin; sr +=stepr)
		  {
		    xc = (int) (xtr + sr*cos(sf/180.*(double)TMath::Pi()));
		    yc = (int) (ytr + sr*sin(sf/180.*(double)TMath::Pi()));
		    if(xc == xtr && yc == ytr) continue;
		    if(xc < 0 || xc >= hsize || yc < 0|| yc >= vsize) continue;
		    //if(pCDBuffer[image][yc][xc] == 0) continue;
		    
		    bc = (float) pCDBuffer[image][yc][xc];
		    if(bc < meanbr) break;
		    //float gradient = (bc - brprev)/brprev;
		    //if(sr > iw && gradient*gradient_prev < 1. && gradient > 0.) break;
		    //if(fabs(gradient) < 0.05 && sr > 10) break;
		    //if(gradient > 0.) break;
		    //if(gradient >= 0 && sr > 10) break;
		    //if(bc < br && br/bc < NAT_LOG)
		    if(flag[yc][xc] != 0) {
		      summax += bc;
		      nmax++;
		    }
		    flag[yc][xc] = 0;//rsost[xc][yc] = 1;
		    brprev = bc;
		    //gradient_prev = gradient;
		  }
	      }
	  }
      }
  for(int imark=0; imark<Npoint[image]; imark++)
    {
      int xtr = xrest[image][imark];
      int ytr = vsize - yrest[image][imark] - 1;
      float br = (float) pCDBuffer[image][ytr][xtr];
      if(flag[ytr][xtr] != 0) {
	summax += br;
	nmax++;
      }
      flag[ytr][xtr] = 0;
      float stepr = 1.0;
      float stepf = 0.5;
      float sf;
      float brprev;
      float sr, bc; int xc, yc;
      for (sf = 0; sf < 360; sf+=stepf)
	{
	  brprev = br;
	  for (sr = stepr; sr < iwin; sr +=stepr)
	    {
	      xc = (int) (xtr + sr*cos(sf/180.*(double)TMath::Pi()));
	      yc = (int) (ytr + sr*sin(sf/180.*(double)TMath::Pi()));
	      if(xc == xtr && yc == ytr) continue;
	      if(xc < 0 || xc >= hsize || yc < 0|| yc >= vsize) break;
	      
	      bc = (float) pCDBuffer[image][yc][xc];
	      if(bc < meanbr) break;
	      //float gradient = (bc - brprev)/brprev;
	      //if(sr > 20 && gradient*gradient_prev < 1. && gradient > 0.) break;
	      //if(fabs(gradient) < 0.1 || bc <= aver[image]) break;
	      //if(gradient >= 0) break;
	      //if(bc < br && br/bc < NAT_LOG)
	      if(flag[yc][xc] != 0) {
		summax += bc;
		nmax++;
	      }
	      flag[yc][xc] = 0;//rsost[xc][yc] = 1;
	      brprev = bc;
	      //gradient_prev = gradient;
	    }
	}
    }
  if(nmax !=0) avermax = summax/(float) nmax;
}
void LifeFirst(int image, int maxiter, int ipass)
{
  int i, j;
 int MaxNumberZeros = 1;
 int NumberZeros;
  int iter, irest, irestprev, ialive, NumberAlife;
  //float MinGrad = 1.00;
  float MinGrad = 0.05;//1.2;
  int MinLife = 0;
  int MaxLife = 8;
  int MaxBig = 8;
  int MinBig = 1;
  int Big;
  int column, row, npositive;
  float grad[8], sumgrad;
  for (column=0; column<hsize; column++)
    {  
      for (row=0; row<vsize; row++)
	{
	  sosts[column][row] = 0;
	  rsost[column][row] = 0;
	  if ((column>=tailx[image][0] && column<=hsize-tailx[image][1]) &&
	      (row>=taily[image][1] && row<=vsize-taily[image][0]-1))
	    sosts[column][row] = 1;
	}
    }
  iw = 10;
  for (iter=0; iter<2; iter++)
    {
      float sumgrad_all = 0.;
      irest = 0;
      irestprev = irest;
      for (column=tailx[image][0]+iw; column<hsize-tailx[image][1]-iw; column++)
	{
	  for (row=taily[image][1]+iw; row<vsize-taily[image][0]-iw-1; row++)
	    {
	      if (sosts[column][row] == 1) 
		{
		  rsost[column][row] = 0;
		  int NAlife = 0;
		  int NBig = 0;
		  float sumgrad = 0.;
		  int br0 = pCDBuffer[image][row][column];
		  for(i=column-iw; i<=column+iw; i++)
		    for(j=row-iw; j<=row+iw; j++)
		      {
			if (sosts[i][j] == 1 && i != column && j != row)
			  {
			    int br = pCDBuffer[image][j][i];
			    float gradient = (float) (br0 - br);
			    gradient /= ((float) br0);
			    sumgrad += gradient;
			    if(pCDBuffer[image][row][column] >= pCDBuffer[image][j][i])
			      NBig++;
			    NAlife++;
			  }
		      }
		  if(NAlife > 0) sumgrad = sumgrad/((float) NAlife);
		  if(sumgrad > MinGrad)
		    if(NBig >= (2-iter)*NAlife/4 && NAlife > 0)
		      {
			rsost[column][row] = 1;
			sumgrad_all += sumgrad;
		      }
		}
	    }
	}
      ialive = 0;
      for (column=tailx[image][0]; column<hsize-tailx[image][1]; column++)
	{
	  for (row=taily[image][1]; row<vsize-taily[image][0]-1; row++)
	    {
	      if(rsost[column][row] == sosts[column][row]) irest++;
	      sosts[column][row] = rsost[column][row];
	      if(rsost[column][row] == 1) ialive++;
	    }
	}
      printf(" Third LOOP: %d points are dead, %d are alive with mean grad %f\n", 
	     irest, ialive, sumgrad_all/((float) ialive));
      iw -= 5;
      MinGrad -= 0.01;//sumgrad_all/((float) ialive)/2.;
    }
  //if(ipass > 2)
  {
      int MaxLife;
      if(ipass > 1) {iw = 2; MaxLife = 2;}
      if(ipass == 1) {iw = 1; MaxLife = 1;}
      irest = 0;
      for (iter=0; iter<1; iter++)
	{
	  irestprev = irest;
	  for (column=tailx[image][0]+iw; column<hsize-tailx[image][1]-iw; column++)
	    {
	      for (row=taily[image][1]+iw; row<vsize-taily[image][0]-iw-1; row++)
		{
		  if (sosts[column][row] == 1) 
		    {
		      rsost[column][row] = 0;
		      int NAlife = 0;
		      for(i=column-iw; i<=column+iw; i++)
			for(j=row-iw; j<=row+iw; j++)
			  if (sosts[i][j] == 1 && i != column && j != row)
			    NAlife++;
		      if(NAlife > MaxLife) 
			rsost[column][row] = 1;
		    }
		}
	    }
	  ialive = 0;
	  for (column=tailx[image][0]; column<hsize-tailx[image][1]; column++)
	    {
	      for (row=taily[image][1]; row<vsize-taily[image][0]-1; row++)
		{
		  if(rsost[column][row] == sosts[column][row]) irest++;
		  sosts[column][row] = rsost[column][row];
		  if(rsost[column][row] == 1) ialive++;
		}
	    }
	  printf(" 4th LOOP: %d points are dead, %d are alive \n", 
		 irest, ialive);
	  iter++;
	}
  }
}
void LifeHard(int image, int maxiter, int ipass)
{
  int i, j;
  int iter, irest, irestprev, ialive, NumberAlife;
  float MinGrad = 1.00;
  //float MinGrad = 1.25;
  int MinLife = 4;
  int column, row;
  float grad[8], sumgrad;
  float av;
  printf("Average br on image before %d pass is %f %f\n", ipass, aver[image], MinGrad);
  for (column=0; column<hsize; column++)
    {  
      for (row=0; row<vsize; row++)
	{
	  sosts[column][row] = 0;
	  (pCDBuffer[image][row][column]>aver[image]/(float)(ipass) && 
	   flag[row][column] == 1)?
	    sosts[column][row] = 1:
	    sosts[column][row] = 0;
	}
    }
  iter = 0;
  irest = 100;
  irestprev = 0;
  //while (irest != irestprev)
  ialive = 1000*ipass;
  while (ialive >= 1000*ipass && irest!= irestprev)
    //for (int iter=0; iter<1; iter++)
    {
      irestprev = irest;
      float sumw = 0.;
      int inum = 0;
      for (i=taily[image][1]; i<vsize-taily[image][0]-1; i++)
	{
	  for (j=tailx[image][0]; j<hsize-tailx[image][1]; j++)
	    {
	      sumw += (sosts[j][i] == 1)?pCDBuffer[image][i][j]:0.;
	      inum = (sosts[j][i] == 1)?inum+1:inum;
	    }
	}
      inum?av = sumw/(float) inum/(float) (iter+1): av = 0; 
  
      for (column=1; column<hsize-1; column++)
	{
	  for (row=1; row<vsize-1; row++)
	    {
	      //if ((row>taily[image][1] && row<vsize-taily[image][0]) &&
	      //  (column>tailx[image][0] && column<hsize-tailx[image][1]))
		{
		  sostn[0] = sosts[column][row-1];
		  sostn[1] = sosts[column  ][row+1];
		  sostn[2] = sosts[column-1][row  ];
		  sostn[3] = sosts[column+1][row  ];
		  sostn[4] = sosts[column-1][row-1];
		  sostn[5] = sosts[column+1][row-1];
		  sostn[6] = sosts[column-1][row+1];
		  sostn[7] = sosts[column+1][row+1];
		  
		  grad[0] = pCDBuffer[image][row-1][column];
		  grad[1] = pCDBuffer[image][row+1][column  ];
		  grad[2] = pCDBuffer[image][row  ][column-1];
		  grad[3] = pCDBuffer[image][row  ][column+1];
		  grad[4] = pCDBuffer[image][row-1][column-1];
		  grad[5] = pCDBuffer[image][row-1][column+1];
		  grad[6] = pCDBuffer[image][row+1][column-1];
		  grad[7] = pCDBuffer[image][row+1][column+1];
		  
		  rsost[column][row] = 0;
		  sumgrad = 0.;
		  NumberAlife = 0;
		  for (i = 0; i< 8; i++)
		    if(sostn[i] == 1) NumberAlife++;
		  //for (i = 0; i< 8; i++)
		  //	if(sostn[i] == 1) 
		  //  sumgrad += grad[i];
		  //float meanbr = (sumgrad+pCDBuffer[image][row][column])/
		  //(float) (NumberAlife + 1);
		  //sumgrad = 0;
		  for (i = 0; i< 8; i++)
		    if(sostn[i] == 1)
		      sumgrad += pCDBuffer[image][row][column]/grad[i];
		  if(NumberAlife != 0) sumgrad = sumgrad/(float) NumberAlife;
		  if (sosts[column][row] == 1 && flag[row][column] == 1) 
		    {
		      if(NumberAlife > MinLife &&
			 (sumgrad > MinGrad &&
			pCDBuffer[image][row][column] >= av))
			{
			  //printf("Number of Alive Cells %d Column %d Row %d\n", NumberAlife, column, row);
			  rsost[column][row] = 1;
			  //br[column][row] = meanbr; 
			}
		      /*
		      if(NumberAlife > 7 && NumberZeros < MaxNumberZeros &&
			  pCDBuffer[image][row][column] >= av)
			{
			  rsost[column][row] = 1;
			}
		      */
		    }
		}
	    }
	}
      irest = 0;
      ialive = 0;
      for (column=0; column<hsize; column++)
	{
	  for (row=0; row<vsize; row++)
	    {
	      //if ((row>taily[image][1] && row<vsize-taily[image][0]) &&
	      //  (column>tailx[image][0] && column<hsize-tailx[image][1]))
		{
		  if(rsost[column][row] == sosts[column][row]) irest++;
		  sosts[column][row] = rsost[column][row];
		  if(rsost[column][row] == 1) ialive++;
		}
	    }
	}
      printf(" %d points are dead, %d are alive %6.1f averagebr in %d iteration\n", 
	     irest, ialive, av, iter);
      iter++;
    }
}
void LifeSoft(int image, int iterhard, int ipass)
{
  int i, j;
 int MaxNumberZeros = 1;
 int NumberZeros;
 int iter, irest, irestprev, ialive, NumberAlife;

 float MinGrad = 0.699;//0.4999;
  int MinLife = 4;
  int column, row;
  float grad[8], sumgrad;
  float av;
  for (column=0; column<hsize; column++)
    {
      for (row=0; row<vsize; row++)
	{
	  sosts[column][row] = 0;
	  (pCDBuffer[image][row][column]>main_aver[image]/(float)(ipass-1) && 
	   flag[row][column]==1)?
	    sosts[column][row] = 1:
	    sosts[column][row] = 0;
	}
    }
  iter = 0;
  irest = 100;
  irestprev = 0;
  //while (irest != irestprev)
  //      ialive = 1000000;
  //while (ialive > 1000/ipass)
  while (irest != irestprev)
    {
      irestprev = irest;
      float sumw = 0.;
      int inum = 0;
      for (i=taily[image][1]; i<vsize-taily[image][0]-1; i++)
	{
	  for (j=tailx[image][0]; j<hsize-tailx[image][1]; j++)
	    {
	      sumw += (sosts[j][i] == 1)?pCDBuffer[image][i][j]:0.;
	      inum = (sosts[j][i] == 1)?inum+1:inum;
	    }
	}
      av = sumw/(float) inum;
      
      for (column=0; column<hsize; column++)
	{
	  for (row=0; row<vsize; row++)
	    {
	      if ((row>taily[image][1] && row<vsize-taily[image][0]-1) &&
		  (column>tailx[image][0] && column<hsize-tailx[image][1]))
		{
		  if (sosts[column][row] == 1 && flag[row][column] == 1) 
		    {
		      sostn[0] = sosts[column][row-1];
		      sostn[1] = sosts[column  ][row+1];
		      sostn[2] = sosts[column-1][row  ];
		      sostn[3] = sosts[column+1][row  ];
		      sostn[4] = sosts[column-1][row-1];
		      sostn[5] = sosts[column+1][row-1];
		      sostn[6] = sosts[column-1][row+1];
		      sostn[7] = sosts[column+1][row+1];
		      
		      grad[0] = pCDBuffer[image][row-1][column];
		      grad[1] = pCDBuffer[image][row+1][column  ];
		      grad[2] = pCDBuffer[image][row  ][column-1];
		      grad[3] = pCDBuffer[image][row  ][column+1];
		      grad[4] = pCDBuffer[image][row-1][column-1];
		      grad[5] = pCDBuffer[image][row-1][column+1];
		      grad[6] = pCDBuffer[image][row+1][column-1];
		      grad[7] = pCDBuffer[image][row+1][column+1];
		      
		      rsost[column][row] = 0;
		      sumgrad = 0.;
		      NumberAlife = 0;
		      NumberZeros = 0;
		      for (i = 0; i< 8; i++)
			if(sostn[i] == 1) NumberAlife++;
		      for (i = 0; i< 8; i++)
			if(flag[row][column] == 0) NumberZeros++;
		      //for (i = 0; i< 8; i++)
		      //	if(sostn[i] == 1) 
		      //  sumgrad += grad[i];
		      //float meanbr = (sumgrad+pCDBuffer[image][row][column])/
		      //(float) (NumberAlife + 1);
		      //sumgrad = 0;
		      for (i = 0; i< 8; i++)
			if(sostn[i] == 1)
			  sumgrad += pCDBuffer[image][row][column]/grad[i];
		      //sumgrad += TMath::Log10(pCDBuffer[image][row][column])/
		      //TMath::Log10(grad[i]);
		      if(NumberAlife != 0) sumgrad = sumgrad/(float) NumberAlife;
		      if(NumberAlife > MinLife && 
			 (sumgrad > MinGrad && NumberZeros < MaxNumberZeros &&
			  pCDBuffer[image][row][column] >= av ))
			{
			  rsost[column][row] = 1;
			  //br[column][row] = meanbr; 
			}
		    }
		}
	    }
	}
      irest = 0;
      ialive = 0;
      for (column=0; column<hsize; column++)
	{
	  for (row=0; row<vsize; row++)
	    {
	      if ((row>taily[image][1] && row<vsize-taily[image][0]-1) &&
		  (column>tailx[image][0] && column<hsize-tailx[image][1]))
		{
		  if(rsost[column][row] == sosts[column][row]) irest++;
		  sosts[column][row] = rsost[column][row];
		  if(rsost[column][row] == 1) ialive++;
		}
	    }
	}
      printf(" %d points are daed, %d are alive %6.1f averagebr in %d iteration\n", 
	     irest, ialive, av, iter);
      iter ++;
    }
}




void GroupPoints(int image, int iw, int ipass)
{
  int npoint = 0;
  //if(ipass == 1) iw = 4;
  //if(ipass > 1) iw = 4;
  if(Npoint[image] < 100) iw =0;
  else if(Npoint[image] < 200) iw = 1;
  else if(Npoint[image] < 400) iw = 2;
  else if(Npoint[image] < 900) iw = 3;
  else if(Npoint[image] < 1500) iw = 4;
  else iw = 5;
  printf(" There will be density window used %d\n", iw);
  int maxbr;
  int i, j, ip, bi;
  maxbr = 0;
  int column, row;
  for(ip=0; ip<Npoint[image]; ip++)
    {
      i = xrest[image][ip];
      j = vsize - yrest[image][ip] - 1;
      sosts[i][j] = 1;
    }
  //while (iw < 3)
  if (iw == 0)
    {
      npoint = 0;
      for (column=tailx[image][0]+iw; column<hsize-tailx[image][1]-iw; column++)
	//column += 2*iw)
	{
	  for (row=taily[image][1]+iw; row<vsize-taily[image][0]-iw-1; row++)
	    //row+=2*iw)
	    {
	      if ((row>taily[image][1]+iw && row<vsize-taily[image][0]-iw-1) &&
		  (column>tailx[image][0]+iw && column<hsize-tailx[image][1]-iw))
		{
		  maxbr = 0;
		  for(ip=0; ip<Npoint[image]; ip++)
		    {
		      i = xrest[image][ip];
		      j = vsize - yrest[image][ip] - 1;
		      bi = pCDBuffer[image][j][i];
		      if(i >= column - iw && i <= column + iw && j >= row - iw && j <= row + iw)
			{
			  if (sosts[i][j] == 1 && 
			      bi > maxbr)
			    {
			      sosts[i][j] = 0;

			      maxbr = bi;
			      xmain[image][npoint] = i;
			      ymain[image][npoint] = j;  
			    }
			}
		    }
		  if (maxbr > 0) {
		    sosts[xmain[image][npoint]][ymain[image][npoint]] = 1;
		    ymain[image][npoint] = vsize - ymain[image][npoint] - 1;  
		    npoint ++;
		  }
		}
	    }
	}
      Npoint_main = npoint;
      printf(" %d MAIN points are on the %d image with %d density\n", 
	     Npoint_main, image, iw);
      iw *= 2;
    }
  else
    {
      npoint = 0;
      for (column=tailx[image][0]+iw; column<hsize-tailx[image][1]-iw; 
	   column += 2*iw)
	{
	  for (row=taily[image][1]+iw; row<vsize-taily[image][0]-iw-1; 
	    row+=2*iw)
	    {
	      if ((row>taily[image][1]+iw && row<vsize-taily[image][0]-iw-1) &&
		  (column>tailx[image][0]+iw && column<hsize-tailx[image][1]-iw))
		{
		  maxbr = 0;
		  for(ip=0; ip<Npoint[image]; ip++)
		    {
		      i = xrest[image][ip];
		      j = vsize - yrest[image][ip] - 1;
		      bi = pCDBuffer[image][j][i];
		      if(i >= column - iw && i < column + iw && j >= row - iw && j < row + iw)
			{
			  if (sosts[i][j] == 1 && 
			      bi > maxbr)
			    {
			      sosts[i][j] = 0;

			      maxbr = bi;
			      xmain[image][npoint] = i;
			      ymain[image][npoint] = j;  
			    }
			}
		    }
		  if (maxbr > 0) {
		    sosts[xmain[image][npoint]][ymain[image][npoint]] = 1;
		    ymain[image][npoint] = vsize - ymain[image][npoint] - 1;  
		    npoint ++;
		  }
		}
	    }
	}
      Npoint_main = npoint;
      printf(" %d MAIN points are on the %d image with %d density\n", 
	     Npoint_main, image, iw);
      iw *= 2;
    }
}
void EraseNearestToTrackPoints(int image, int Start, int End, int ipass, int what)
{
  int i;
  float X_Start = xmain[image][Start];
  float Y_Start = ymain[image][Start];
  float X_End   = xmain[image][End];
  float Y_End   = ymain[image][End];
  float Slope, Y_Origin;
  int MD = abs(5 + (ipass-1)*10);
  if(X_End != X_Start) 
    {
      Slope = (Y_End - Y_Start)/(X_End - X_Start);
      Y_Origin = Y_Start - Slope*X_Start;
      Slope = TMath::ATan2((Y_End - Y_Start),(X_End - X_Start));
    }
  else if(X_End == X_Start)
    {
      Slope = TMath::Pi()/2.;
      Y_Origin = X_End;
    }
  else if(Y_End == Y_Start)
    {
      Slope = 0.;
      Y_Origin = Y_End;
    }
  
  for (i = 0; i < TMath::Min(MAXNUMBEROFPOINTS, Npoint_main); i++)
    {
      if(i == Start || i == End) continue;
      if(TrackPoint[i] == what) continue;
      if(xmain[image][i] >= X_Start-MD && xmain[image][i] <= X_End+MD &&
	 ymain[image][i] >= Y_Start-MD && ymain[image][i] <= Y_End+MD)
	if(D_to_Line(xmain[image][i], ymain[image][i], Slope, Y_Origin) < MD)
	  TrackPoint[i] = what;
      if(xmain[image][i] >= X_End-MD && xmain[image][i] <= X_Start+MD &&
	 ymain[image][i] >= Y_End-MD && ymain[image][i] <= Y_Start+MD)
	if(D_to_Line(xmain[image][i], ymain[image][i], Slope, Y_Origin) < MD)
	  TrackPoint[i] = what;
      if(xmain[image][i] >= X_Start-MD && xmain[image][i] <= X_End+MD &&
	 ymain[image][i] >= Y_End-MD && ymain[image][i] <= Y_Start+MD)
	if(D_to_Line(xmain[image][i], ymain[image][i], Slope, Y_Origin) < MD)
	  TrackPoint[i] = what;
      if(xmain[image][i] >= X_End-MD && xmain[image][i] <= X_Start+MD &&
	 ymain[image][i] >= Y_Start-MD && ymain[image][i] <= Y_End+MD)
	if(D_to_Line(xmain[image][i], ymain[image][i], Slope, Y_Origin) < MD)
	  TrackPoint[i] = what;
    }
  
}
float D_to_Line(float xp, float yp, float Slope, float Y_Origin)
{//dimention to the line from the point xp, yp.
  double sinalpha = TMath::Sin(Slope);
  double cosalpha = TMath::Cos(Slope);
  float d = -xp*sinalpha + (yp-Y_Origin)*cosalpha;
  //short int znak = 0;
  //Y_Origin? znak =-1.:znak = 1; 
  //float mu = znak * 1./sqrt(Slope*Slope + 1.);
  //float d = xp*Slope*mu - yp*mu +  Y_Origin*mu;
  return fabs(d);
}
void SelectPoints(int image, int ipass)
{
  int npoint, npoint_win;
  int iww = 20;
  int maxbr;
  int i, j, ip;
  maxbr = 0;
  int column, row;
  npoint_win = 0;
  npoint = 0;
  if(ipass > 1)
    {
      for(ip=0; ip < Npoint[image]; ip++)
	{
	  i = xrest[image][ip];
	  j = vsize - yrest[image][ip] - 1;
	  //bi = pCDBuffer[image][j][i];
	  npoint_win = 0;
	  for(column = i - iww; column <= i + iww; column++)
	    {
	      for(row = j - iww; row <= j + iww; row++)
		{
		  if (flag[row][column] == 0) 
		    npoint_win++;
		}
	    }
	  if (npoint_win > 10) 
	    {
	      xmain[image][npoint] = xrest[image][ip];
	      ymain[image][npoint] = yrest[image][ip];
	      npoint ++;
	    }
	}
      Npoint_main = npoint;
      printf(" %d SELECTED points are on the %d image with %d density\n", 
	     Npoint_main, image, iww);
      for(ip=0; ip < Npoint_main; ip++)
	{
	  xrest[image][ip] = xmain[image][ip];
	  yrest[image][ip] = ymain[image][ip];
	}
      Npoint[image] = Npoint_main;
    }
  //else
    {
      iww = 2;
      npoint_win = 0;
      npoint = 0;
      
      for(ip=0; ip < Npoint[image]; ip++)
	{
	  i = xrest[image][ip];
	  j = vsize - yrest[image][ip] - 1;
	  npoint_win = 0;
	  for(column = i - iww; column <= i + iww; column++)
	    {
	      for(row = j - iww; row <= j + iww; row++)
		{
		  if (flag[row][column] == 0) 
		    npoint_win++;
		}
	    }
	  if (npoint_win > 1) 
	    {
	      xmain[image][npoint] = xrest[image][ip];
	      ymain[image][npoint] = yrest[image][ip];
	      npoint ++;
	    }
	}
      Npoint_main = npoint;
      printf(" %d SELECTED points are on the %d image with %d density\n", 
	     Npoint_main, image, iww);
      for(ip=0; ip < Npoint_main; ip++)
	{
	  xrest[image][ip] = xmain[image][ip];
	  yrest[image][ip] = ymain[image][ip];
	}
      Npoint[image] = Npoint_main;
    }
}
void FindVertex(int image, int ipass)
{
  float xc[MAXNUMBEROFTRACKS], yc[MAXNUMBEROFTRACKS], rc[MAXNUMBEROFTRACKS];
  float Xvertex, Yvertex, Dvertex;
  int itr, itr1, itr2, itr3, itr4;
  int track[MAXNUMBEROFTRACKS];
  int i, j, ii, jj, jtr;
  float Phi_FirstPoint, Phi_LastPoint, Phi_Vertex, buf1, buf2, Buf;
  int Segment_Between, nmarked;
  int Limit0, Limit1;
  short int not_for_vertex[2][20];
  float xvert_prev, yvert_prev, err_vert_prev;
  int t_prev[MAXNUMBEROFTRACKS];
  Limit0 = 0; 
  Limit1 = TMath::Min(MAXNUMBEROFPOINTS, Npoint_main);
  Not_Used = 0;
  for (int i = Limit0; i < Limit1; i++)
    if(TrackPoint[i] == 0) Not_Used++;
  eff = 100.*(float)Not_Used/(float)Npoint_main;
  //printf(" %d (%f \%) points were not used in tracking\n", 
  // Not_Used, eff);
  nmarked = 0;
  track[0] = track[1] = track[2] = track[3] = 0;
  for (itr = 0; itr < ntr[image]; itr++)
    {
      jj = 0;
      float WpMin = 0.;
      float WpMax = 4096.;
      for(i=0; i<np[image][itr]-1; i++) 
	{
	  for(j=0; j<Nmain[image][itr][i]; j++) 
	    {
	      Xcircle[jj] = (float) Xmain[image][itr][i][j];
	      Ycircle[jj] = (float) Ymain[image][itr][i][j];
	      int iy = vsize - Ymain[image][itr][i][j] - 1;
	      int ix = Xmain[image][itr][i][j];
	      Wpoint[jj] = (float) pCDBuffer[image][iy][ix];
	      if (Wpoint[jj] > WpMax) WpMax = Wpoint[jj];
	      if (Wpoint[jj] < WpMin) WpMin = Wpoint[jj];
	      jj++;
	    }
	}
      for(i = 0; i< jj; i++) 
	{
	  if(WpMin != WpMax) Wpoint[i]  = (Wpoint[i]-WpMin)/(WpMax - WpMin);
	  if(WpMin == WpMax) Wpoint[i]  = Wpoint[i]/WpMax;
	}
      xcD[image][itr] = ycD[image][itr] = rcD[image][itr] = 0.;
      int ierr = circle_r(jj);
      printf("IERR %d Track %d Points %d:  XC %f YC %f RD %f ch2 %f\n", ierr, itr, jj, XC, YC, RC, chi2);
      Npoints_fit = jj;
      circle_minuit(jj);
      xcD[image][itr] = XC;
      ycD[image][itr] = YC;
      rcD[image][itr] = RC;
      chi2_circle[image][itr] = chi2;
      for(i=0;i<jj; i++)
	{
	  double buf = atan2(Ycircle[i] - YC, Xcircle[i] - XC);
	  x_mark[i] = XC + RC*cos(buf);
	  y_mark[i] = YC + RC*sin(buf);
	}
      px_mark = &x_mark[0];
      py_mark = &y_mark[0];
      marke = new TPolyMarker((Int_t) jj, px_mark, py_mark, "same");
      fObjTbl->Add(marke);
      marke->SetMarkerStyle(8);
      marke->SetMarkerSize(0.1);
      marke->SetMarkerColor(4);
      marke->Draw("same");	     
    }
  for (itr = 0; itr < 20; itr++)
    {
      mark[image][itr] = 0;
      not_for_vertex[image][itr] = 0;
    }
 LABELVERTEX:
  nmarked = 0;
  itr = 0;
  Min_Vertex_Error[image] = 100000.;
  xvert_prev = yvert_prev = err_vert_prev = 0.;
  if(ntr[image] > 3) 
    {
      for (itr1 = 0; itr1 < ntr[image] - 3; itr1++)
	if(rcD[image][itr1] != 0. && not_for_vertex[image][itr1] == 0)
	  {
	    xc[0] = xcD[image][itr1];
	    yc[0] = ycD[image][itr1];
	    rc[0] = rcD[image][itr1];
	    for (itr2 = itr1 + 1; itr2 < ntr[image] - 2; itr2++)
	      if(rcD[image][itr2] != 0. && not_for_vertex[image][itr2] == 0)
		{
		  xc[1] = xcD[image][itr2];
		  yc[1] = ycD[image][itr2];
		  rc[1] = rcD[image][itr2];
		  for (itr3 = itr2 + 1; itr3 < ntr[image] - 1; itr3++)
		    if(rcD[image][itr3] != 0. && not_for_vertex[image][itr3] == 0)
		      {
			xc[2] = xcD[image][itr3];
			yc[2] = ycD[image][itr3];
			rc[2] = rcD[image][itr3];
			for (itr4 = itr3 + 1; itr4 < ntr[image]; itr4++)
			  if(rcD[image][itr4] != 0. && not_for_vertex[image][itr4] == 0)
			    {
			      Xvertex = Yvertex = Dvertex = 0.;
			      itr = 4;
			      xc[3] = xcD[image][itr4];
			      yc[3] = ycD[image][itr4];
			      rc[3] = rcD[image][itr4];
			      XY_Vertex_Net(itr, &xc[0], &yc[0], 
					    &rc[0], &Xvertex, &Yvertex, 
					    &Dvertex);
			      printf("Dvertex with 4 tracks is %f\n", Dvertex);
			      if(Dvertex < Min_Vertex_Error[image] &&
				 Xvertex > tailx[image][0]+50 && 
				 Xvertex < hsize - tailx[image][1]-50 &&
				 Yvertex > 220. && 
				 Yvertex < 350.) 
				{
				  Min_Vertex_Error[image] = Dvertex;
				  track[0] = itr1;
				  track[1] = itr2;
				  track[2] = itr3;
				  track[3] = itr4;
				  nmarked = 4;
				  X_Vert[image] = Xvertex;
				  Y_Vert[image] = Yvertex;
				  printf(" Best Vertex position is : X %f Y %f Error is %f\n", 
					 X_Vert[image], Y_Vert[image], Min_Vertex_Error[image]);
				  printf(" 4 Track numbers are : %d %d %d %d\n", 
					 track[0], track[1], track[2], track[3]);
				  if(xvert_prev != 0. && yvert_prev != 0. && err_vert_prev < 15. &&
				     sqrt(pow(xvert_prev-Xvertex,2) + pow(yvert_prev-Yvertex,2)) < 15.)
				    //add previos found track to these ones
				    {
				      ii = nmarked;
				      for(int iii = 0; iii < nmarked; iii++)
					{
					  int est = 0;
					  for(int jjj = 0; jjj < nmarked; jjj++)
					    if(track[jjj] == t_prev[iii]) est = 1;
					  if(est == 0) track[ii] = t_prev[iii];
					  if(est == 0) 
					    {
					      printf(" One previously found track %d has been added\n", 
						     track[ii]);
					      ii++;
					    }
					}
				      nmarked = ii;
				    }
				//remember vertex position for a case if there will be one best vertex
				// found not very far from this one
				  xvert_prev = Xvertex;
				  yvert_prev = Yvertex;
				  err_vert_prev = Dvertex;
				  t_prev[0] = itr1;
				  t_prev[1] = itr2;
				  t_prev[2] = itr3;
				  t_prev[3] = itr4;
				}
			    }
		      }
		}
	  }
    }
  if ((Min_Vertex_Error[image] > 5.) &&(ntr[image] >= 3))
    {
      Min_Vertex_Error[image] = 100000.;
      for (itr1 = 0; itr1 < ntr[image] - 2; itr1++)
	if(rcD[image][itr1] != 0. && not_for_vertex[image][itr1] == 0)
	  {
	    xc[0] = xcD[image][itr1];
	    yc[0] = ycD[image][itr1];
	    rc[0] = rcD[image][itr1];
	    for (itr2 = itr1 + 1; itr2 < ntr[image] - 1; itr2++)
	      if(rcD[image][itr2] != 0. && not_for_vertex[image][itr2] == 0)
		{
		  xc[1] = xcD[image][itr2];
		  yc[1] = ycD[image][itr2];
		  rc[1] = rcD[image][itr2];
		  for (itr3 = itr2 + 1; itr3 < ntr[image]; itr3++)
		    if(rcD[image][itr3] != 0. && not_for_vertex[image][itr3] == 0)
		      {
			xc[2] = xcD[image][itr3];
			yc[2] = ycD[image][itr3];
			rc[2] = rcD[image][itr3];
			Xvertex = Yvertex = Dvertex = 0.;
			itr = 3;
			XY_Vertex_Net(itr, &xc[0], &yc[0], 
				      &rc[0], &Xvertex, &Yvertex, 
				      &Dvertex);
			printf("Dvertex with 3 tracks is %f\n", Dvertex);
			if(Dvertex < Min_Vertex_Error[image] &&
			   Xvertex > tailx[image][0]+50 && 
			   Xvertex < hsize - tailx[image][1]-50 &&
			   Yvertex > 220. && 
			   Yvertex < 350.) 
			  {
			    Min_Vertex_Error[image] = Dvertex;
			    track[0] = itr1;
			    track[1] = itr2;
			    track[2] = itr3;
			    nmarked = 3;
			    track[3] = 0; // there is no track
			    X_Vert[image] = Xvertex;
			    Y_Vert[image] = Yvertex;
			    printf(" Best Vertex position is : X %f Y %f Error is %f\n", 
				   X_Vert[image], Y_Vert[image], Min_Vertex_Error[image]);
			    printf(" 3 Track numbers are : %d %d %d\n", 
				   track[0], track[1], track[2]);
			  }
		      }
		}
	  }
    }
    if((ntr[image] >= 2) && ((Min_Vertex_Error[image] > 10. && itr == 3) ||
			     (Min_Vertex_Error[image] > 5. && itr == 4) ||
			     (nmarked == 0)))
    {
      for (itr1 = 0; itr1 < ntr[image] - 1; itr1++)
	if(rcD[image][itr1] != 0. && not_for_vertex[image][itr1] == 0)
	  {
	    xc[0] = xcD[image][itr1];
	    yc[0] = ycD[image][itr1];
	    rc[0] = rcD[image][itr1];
	    for (itr2 = itr1 + 1; itr2 < ntr[image]; itr2++)
	      if(rcD[image][itr2] != 0. && not_for_vertex[image][itr2] == 0)
		{
		  xc[1] = xcD[image][itr2];
		  yc[1] = ycD[image][itr2];
		  rc[1] = rcD[image][itr2];
		  Xvertex = Yvertex = Dvertex = 0.;
		  itr = 2;
		  XY_Vertex_Net(itr, &xc[0], &yc[0], 
				&rc[0], &Xvertex, &Yvertex, 
				&Dvertex);
		  printf("Dvertex with 2 tracks is %f\n", Dvertex);
		  if(Dvertex < Min_Vertex_Error[image] &&
		     Xvertex > tailx[image][0] && 
		     Xvertex < hsize - tailx[image][1] &&
		     Yvertex > 220. && 
		     Yvertex < 350.) 
		    {
		      Min_Vertex_Error[image] = Dvertex;
		      track[0] = itr1;
		      track[1] = itr2;
		      nmarked = 2;
		      track[2] = 0;
		      track[3] = 0;
		      X_Vert[image] = Xvertex;
		      Y_Vert[image] = Yvertex;
		      printf(" Best Vertex position is : X %f Y %f Error is %f\n", 
			     X_Vert[image], Y_Vert[image], Min_Vertex_Error[image]);
		      printf(" 2 Track numbers are : %d %d\n", 
			     track[0], track[1]);
		    }
		  if(nmarked == 2) break;
		}
	    if(nmarked == 2) break;
	  }
    }
  if(nmarked > 0)
    { //if(ntr[image] == 1) mark[image][0] = 1;
      /*
      mark[image][track[0]] = 1;
      mark[image][track[1]] = 1;
      if(track[2] != 0) mark[image][track[2]] = 1;
      if(track[3] != 0) mark[image][track[3]] = 1;
      */
      if(itr != 0 ) 
	{
	  for(i=0; i< nmarked; i++)
	    mark[image][track[i]] = 1; 
	  itr = nmarked;
	}
    }
  else
    mark[image][0] = 1;
  printf("There are %d marked tracks\n", itr);
  
  //try to devide one track if the Vertex lies between of its segments
  if (nmarked != 0) 
    {
      for (i = 0; i < nmarked; i++)
	{
	  //look if this tracks is not crossing pion
	  int fp, lp;
	  float dist_vert0 = sqrt(pow(((float) Y[image][track[i]][0]) - Y_Vert[image], 2)+
				  pow(((float) X[image][track[i]][0]) - X_Vert[image], 2));
	  float dist_vert1 = sqrt(pow(((float) Y[image][track[i]][np[image][track[i]] - 1]) - Y_Vert[image], 2)+
				  pow(((float) X[image][track[i]][np[image][track[i]] - 1]) - X_Vert[image], 2));
	  if(dist_vert0 < dist_vert1) {fp = 0; lp = np[image][track[i]] - 1;}
	  if(dist_vert1 < dist_vert0) {fp = np[image][track[i]] - 1; lp = 0;}
	  double angle = TMath::ATan2((Double_t) (Y[image][track[i]][lp] - Y[image][track[i]][fp]),
					   (Double_t) (X[image][track[i]][lp] - X[image][track[i]][fp]));
	  dist_vert0 = fabs(X[image][track[i]][lp] - X[image][track[i]][fp]);
	  // vertex must be not in the left part of image
	  //if crossing pion goes from the vertex to the left
	  if (X_Vert[image] > 2.*tailx[image][0] && 
	      cos(angle) < -cos(TMath::Pi()/180.*10.) &&
	      dist_vert0 > 250. && 1000. < rcD[image][track[i]])
	    {
	      not_for_vertex[image][track[i]] = 1; 
	      for (jtr = 0; jtr < 20; jtr++)
		mark[image][jtr] = 0;
	      printf(" Crossing left part pion detected, will go to look NEW vertex\n");
	      goto LABELVERTEX;
	    }
	  //if crossing pion goes from the left to the right through the camera and 
	  //vertex is between of the tails
	  dist_vert0 = (X_Vert[image] - X[image][track[i]][fp]);
	  dist_vert1 = fabs(X[image][track[i]][lp] - X[image][track[i]][fp]);
	  if (X_Vert[image] > 2.*tailx[image][0] && 
	      cos(angle) > cos(TMath::Pi()/180.*2.) &&
	      dist_vert0 > 200. && 900. < rcD[image][track[i]] && dist_vert1 > 300.)
	    {
	      not_for_vertex[image][track[i]] = 1; 
	      for (jtr = 0; jtr < 20; jtr++)
		mark[image][jtr] = 0;
	      printf(" Crossing pion detected with %f angle\n          => will go to look NEW vertex\n", 
		     angle*180/3.14);
	      goto LABELVERTEX;
	    }
	}
      for (i = 0; i < itr; i++)
	{
	  buf1 =  TMath::ATan2((float)Y[image][track[i]][0]-ycD[image][track[i]], 
			       (float)X[image][track[i]][0]-xcD[image][track[i]]);
	  buf2 = TMath::ATan2((float)Y[image][track[i]][np[image][track[i]]-1]-ycD[image][track[i]], 
			      (float)X[image][track[i]][np[image][track[i]]-1]-xcD[image][track[i]]);
	  Phi_FirstPoint = TMath::Min(buf1, buf2);
	  Phi_LastPoint  = TMath::Max(buf1, buf2);
	  Phi_Vertex     = TMath::ATan2((float)Y_Vert[image] - ycD[image][track[i]],
					(float)X_Vert[image] - xcD[image][track[i]]);
	  if (Phi_LastPoint - Phi_FirstPoint > TMath::Pi())
	    {
	      Phi_FirstPoint += TMath::Pi()*2.;
	      Buf = Phi_FirstPoint;
	      Phi_FirstPoint = Phi_LastPoint;
	      Phi_LastPoint = Buf;
	    }
	  if(Phi_LastPoint > TMath::Pi() && Phi_Vertex < 0.) Phi_Vertex += TMath::Pi()*2.;
	  if(Phi_Vertex > Phi_FirstPoint && Phi_Vertex < Phi_LastPoint)
	    // vertex is between segments of one track
	    {
	      float L, L1, L2;
	      buf1 =  TMath::ATan2((float)Y[image][track[i]][0]-ycD[image][track[i]], 
				   (float)X[image][track[i]][0]-xcD[image][track[i]]);
	      if(fabs(Phi_Vertex - buf1) < TMath::Pi()) 
		L1 = fabs(rcD[image][track[i]]*(Phi_Vertex-buf1));
	      else
		{
		  L1 = fabs(rcD[image][track[i]]*(2.*TMath::Pi() - (Phi_Vertex-buf1)));
		  if(Phi_Vertex > buf1) Phi_Vertex -= 2.*TMath::Pi();
		  if(Phi_Vertex < buf1) Phi_Vertex += 2.*TMath::Pi();
		}
	      L  = rcD[image][track[i]]*(Phi_LastPoint-Phi_FirstPoint);
	      L2 = L - L1;
	      Segment_Between = -1;
	      //how many segments are from which side of the track
	      for (j = 1; j < np[image][track[i]]; j++)
		{
		  buf2 =  TMath::ATan2((float)Y[image][track[i]][j]-ycD[image][track[i]], 
				       (float)X[image][track[i]][j]-xcD[image][track[i]]);
		  if(fabs(buf2 - buf1) > TMath::Pi())
		    {
		      if(buf2 > buf1) buf2 -= 2.*TMath::Pi();
		      if(buf2 < buf1) buf2 += 2.*TMath::Pi();
		    }
		  if((Phi_Vertex >= buf1 && Phi_Vertex < buf2) ||
		     (Phi_Vertex >= buf2 && Phi_Vertex < buf1)) //Vertex is between given segment
		    {
		      Segment_Between = j - 1;
		      break;
		    }
		}
	      printf("Segment Number between Vertex %d, np %d L1 %f L2 %f\n", 
		     Segment_Between, np[image][track[i]], L1, L2);
	      if(Segment_Between == -1) {printf("there is some error in track %d on image %d\n",
						track[i], image); continue;}
	      if((Segment_Between < (np[image][track[i]]-1)/2 && 
		  (L1 < 30. && (X_Vert[image] - tailx[image][0] > 100  && 
				hsize - tailx[image][1] - X_Vert[image] > 100))) ||
		 Segment_Between == 0 || Segment_Between == 1)
		{
		  // don't devide track on two parts!!! But delete Segments in between!!
		  //and renaitre erased points again!
		  /*
		  int j_s;
		  for (j = 0; j < Segment_Between; j++)
		    {
		      j_s = BestTrack[track[i]][j];
		      TrackPoint[StartOfSegment[j_s]] = 0;
		      TrackPoint[EndOfSegment[j_s]] = 0;
		      EraseNearestToTrackPoints(image, StartOfSegment[j_s], EndOfSegment[j_s], ipass, 0);
		    }
		  */
		  jj = 0;
		  for (j = Segment_Between+1; j < np[image][track[i]]; j++)
		    {
		      X[image][track[i]][jj] = X[image][track[i]][j];
		      Y[image][track[i]][jj] = Y[image][track[i]][j];
		      if(j != np[image][track[i]] - 1) 
			{
			  for (int j1 = 0; j1 < Nmain[image][track[i]][j]; j1++)
			    {
			      Xmain[image][track[i]][jj][j1] = Xmain[image][track[i]][j][j1];
			      Ymain[image][track[i]][jj][j1] = Ymain[image][track[i]][j][j1];
			    }
			  Nmain[image][track[i]][jj] = Nmain[image][track[i]][j];
			}
		      jj++;
		    }
		  np[image][track[i]] = jj;
		  
		}
	      else if(( Segment_Between > (np[image][track[i]]-1)/2 &&
			(L2 < 30. && (X_Vert[image] - tailx[image][0] > 100  && 
				      hsize - tailx[image][1] - X_Vert[image] > 100))) ||
		      Segment_Between == np[image][track[i]]-2 || 
		      Segment_Between == np[image][track[i]]-3)
		{
		  /*
		  int j_s;
		  for (j = Segment_Between; j < np[image][track[i]]-1; j++)
		    {
		      j_s = BestTrack[track[i]][j];
		      TrackPoint[StartOfSegment[j_s]] = 0;
		      TrackPoint[EndOfSegment[j_s]] = 0;
		      EraseNearestToTrackPoints(image, StartOfSegment[j_s], EndOfSegment[j_s], ipass, 0);
		    }
		  */
		  // don't devide track on two parts!!! But delete Segments in between!!
		  np[image][track[i]] = Segment_Between + 1;
		}
	      else
		//devide it on two parts deleting segment inbetween
		{
		  printf(" There will be one addition track %d from track %d\n", ntr[image], track[i]);
		  jj = 0;
		  for (j = 1; j < np[image][track[i]]; j++)
		    {
		      if ( j > Segment_Between)
			{
			  X[image][ntr[image]][jj] = X[image][track[i]][j];
			  Y[image][ntr[image]][jj] = Y[image][track[i]][j];
			  if(j != np[image][track[i]] - 1) 
			    {
			      for (int j1 = 0; j1 < Nmain[image][track[i]][j]; j1++)
				{
				  Xmain[image][ntr[image]][jj][j1] = Xmain[image][track[i]][j][j1];
				  Ymain[image][ntr[image]][jj][j1] = Ymain[image][track[i]][j][j1];
				}
			      Nmain[image][ntr[image]][jj] = Nmain[image][track[i]][j];
			    }
			  jj++;
			}
		    }
		  np[image][track[i]] = Segment_Between + 1;
		  np[image][ntr[image]] = jj;
		  xcD[image][ntr[image]] = xcD[image][track[i]];
		  ycD[image][ntr[image]] = ycD[image][track[i]];
		  rcD[image][ntr[image]] = rcD[image][track[i]];
		  mark[image][ntr[image]] = 1;
		  nmarked++;
		  if(ntr[image] < MAXNUMBEROFTRACKS) ntr[image]++;
		}
	    }
	}
      //if(nmarked <= 2) 
      //there was not found additional track in case of two tracks=> go to the next pass
      //Min_Vertex_Error[image] = 10000.;
      printf("Image %d: There are %d tracks in all\n", image, ntr[image]);
    }
}
void LookForRestTrack(int image)
{
  int i, j, k, itr, jj;
  int nadd = 0;
  float xc[20], yc[20], rc[20];
  float Xvertex, Yvertex, Dvertex;
  float Phi_FirstPoint, Phi_LastPoint, Phi_Vertex, buf1, buf2, Buf;
  int Segment_Between;
  // LABELVERT:
  if(X_Vert[image] != 0. && Y_Vert[image] != 0. && 
     Min_Vertex_Error[image] < 5.) {
    xmain[image][Npoint_main] = (int) X_Vert[image];
    ymain[image][Npoint_main] = (int) Y_Vert[image];
    TrackPoint[Npoint_main] = 0;
    Npoint_main++;
  }
  int ii = 0;
  for (i=0; i<Npoint_main; i++)
    {
      if(TrackPoint[i] == 0){
	xrest[image][ii] = xmain[image][i];
	yrest[image][ii] = ymain[image][i];
	ii++;
      }
    }
  for (i=0; i<ii; i++)
    {
	xmain[image][i] = xrest[image][i];
	ymain[image][i] = yrest[image][i];
    }
  Npoint_main = ii;
  MakeMainSegments(image, -1);
  SelectTracks(1, image, -1);
  printf(" Number of tracks %d\n", NumberOfTracks);
  nadd = 0;
  for (int itr = ntr[image]; itr < NumberOfTracks+ntr[image]; itr++)
    {
      int jtr = itr - ntr[image];
      j = 0;
      for (k = 0; k<LengthOfBestTrack[jtr]; k++)
	{
	  i = BestTrack[jtr][k];
	  X[image][itr][j] = xmain[image][StartOfSegment[i]];
	  Y[image][itr][j] = ymain[image][StartOfSegment[i]];
	  j++;
	}
      k = LengthOfBestTrack[jtr] - 1;
      i = BestTrack[jtr][k];
      X[image][itr][j] = xmain[image][EndOfSegment[i]];
      Y[image][itr][j] = ymain[image][EndOfSegment[i]];
      j++;
      np[image][itr] = j;
    }
  FindPointsForCircleFit(image);
  for (itr = ntr[image]; itr < ntr[image]+NumberOfTracks; itr++)
    {
      mark[image][itr] = 0;
      jj = 0;
      float WpMin = 0.;
      float WpMax = 4096.;
      for(i=0; i<np[image][itr]-1; i++) 
	{
	  for(j=0; j<Nmain[image][itr][i]; j++) 
	    {
	      Xcircle[jj] = (float) Xmain[image][itr][i][j];
	      Ycircle[jj] = (float) Ymain[image][itr][i][j];
	      int iy = vsize - Ymain[image][itr][i][j] - 1;
	      int ix = Xmain[image][itr][i][j];
	      Wpoint[jj] = (float) pCDBuffer[image][iy][ix];
	      if (Wpoint[jj] > WpMax) WpMax = Wpoint[jj];
	      if (Wpoint[jj] < WpMin) WpMin = Wpoint[jj];
	      jj++;
	    }
	}
      for(i = 0; i< jj; i++) 
	{
	  if(WpMin != WpMax) Wpoint[i]  = (Wpoint[i]-WpMin)/(WpMax - WpMin);
	  if(WpMin == WpMax) Wpoint[i]  = Wpoint[i]/WpMax;
	}
      xcD[image][itr] = ycD[image][itr] = rcD[image][itr] = 0.;
      int ierr = circle_r(jj);
      printf("Track %d:  XC %f YC %f RD %f ch2 %f\n", itr, XC, YC, RC, chi2);
      Npoints_fit = jj;
      circle_minuit(jj);
      xcD[image][itr] = XC;
      ycD[image][itr] = YC;
      rcD[image][itr] = RC;
      chi2_circle[image][itr] = chi2;
      for(i=0;i<jj; i++)
	{
	  double buf = atan2(Ycircle[i] - YC, Xcircle[i] - XC);
	  x_mark[i] = XC + RC*cos(buf);
	  y_mark[i] = YC + RC*sin(buf);
	}
      px_mark = &x_mark[0];
      py_mark = &y_mark[0];
      marke = new TPolyMarker((Int_t) jj, px_mark, py_mark, "same");
      fObjTbl->Add(marke);
      marke->SetMarkerStyle(8);
      marke->SetMarkerSize(0.1);
      marke->SetMarkerColor(4);
      if(image == 0) pleft->cd();
      if(image == 1) pright ->cd();
      marke->Draw("same");	     
      if(rcD[image][itr] != 0.)
	{
	  ii = 0;
	  for (i=0; i<ntr[image]; i++)
	    {
	      if(mark[image][i] == 1)
		{
		  xc[ii] = xcD[image][i];
		  yc[ii] = ycD[image][i];
		  rc[ii] = rcD[image][i];
		  ii++;
		}
	    }
	  xc[ii] = xcD[image][itr];
	  yc[ii] = ycD[image][itr];
	  rc[ii] = rcD[image][itr];
	  ii++;
	  Xvertex = Yvertex = Dvertex = 0.;
	  i = ii;
	  XY_Vertex_Net(i, &xc[0], &yc[0], 
			&rc[0], &Xvertex, &Yvertex, 
			&Dvertex);
	  printf("Dvertex with track %d is %f\n", itr, Dvertex);
	  if(Dvertex < 5.)
	    /*&&
	     Xvertex > tailx[image][0] && 
	     Xvertex < hsize - tailx[image][1] &&
	     Yvertex > taily[image][0]+100 && 
	     Yvertex < vsize - taily[image][1]-100 - 1)
	    */ 
	    {
	      Min_Vertex_Error[image] = Dvertex;
	      X_Vert[image] = Xvertex;
	      Y_Vert[image] = Yvertex;
	      mark[image][itr] = 1;
	      printf(" Best Vertex position is : X %f Y %f Error is %f\n", 
		     X_Vert[image], Y_Vert[image], Dvertex);
	      printf(" Now Track numbers are : %d \n", i);
	    }
	}
      nadd++;// = 1;
    }
  ntr[image] += nadd;
  for (itr = 0; itr < ntr[image]; itr++)
    {
      if(mark[image][itr] == 1) 
	{
	  //cut segments after Vertex from smallest side
	  buf1 =  TMath::ATan2((float)Y[image][itr][0]-ycD[image][itr], 
			       (float)X[image][itr][0]-xcD[image][itr]);
	  buf2 = TMath::ATan2((float)Y[image][itr][np[image][itr]-1]-ycD[image][itr], 
			      (float)X[image][itr][np[image][itr]-1]-xcD[image][itr]);
	  Phi_FirstPoint = TMath::Min(buf1, buf2);
	  Phi_LastPoint  = TMath::Max(buf1, buf2);
	  Phi_Vertex     = TMath::ATan2((float)Y_Vert[image] - ycD[image][itr],
					(float)X_Vert[image] - xcD[image][itr]);
	  if (Phi_LastPoint - Phi_FirstPoint > TMath::Pi())
	    {
	      Phi_FirstPoint += TMath::Pi()*2.;
	      Buf = Phi_FirstPoint;
	      Phi_FirstPoint = Phi_LastPoint;
	      Phi_LastPoint = Buf;
	    }
	  if(Phi_LastPoint > TMath::Pi() && Phi_Vertex < 0.) Phi_Vertex += TMath::Pi()*2.;
	  if(Phi_Vertex > Phi_FirstPoint && Phi_Vertex < Phi_LastPoint)
	    // vertex is between segments of one track
	    {
	      float L, L1, L2;
	      buf1 =  TMath::ATan2((float)Y[image][itr][0]-ycD[image][itr], 
				   (float)X[image][itr][0]-xcD[image][itr]);
	      L1 = fabs(rcD[image][itr]*(Phi_Vertex-buf1));
	      L  = rcD[image][itr]*(Phi_LastPoint-Phi_FirstPoint);
	      L2 = L - L1;
	      Segment_Between = -1;
	      //how many segments are from which side of the track
	      for (j = 1; j < np[image][itr]; j++)
		{
		  buf2 =  TMath::ATan2((float)Y[image][itr][j]-ycD[image][itr], 
				       (float)X[image][itr][j]-xcD[image][itr]);
		  if(fabs(buf2 - buf1) > TMath::Pi())
		    {
		      if(buf2 > buf1) buf2 -= 2.*TMath::Pi();
		      if(buf2 < buf1) buf2 += 2.*TMath::Pi();
		    }
		  if((Phi_Vertex >= buf1 && Phi_Vertex < buf2) ||
		     (Phi_Vertex >= buf2 && Phi_Vertex < buf1)) //Vertex is between given segment
		    {
		      Segment_Between = j - 1;
		      break;
		    }
		}
	      printf("Segment Number between Vertex %d, np %d L1 %f L2 %f\n", 
		     Segment_Between, np[image][itr], L1, L2);
	      if(Segment_Between == -1) {printf("error in LookForRest Tracks: track %d on image %d\n",
						itr, image); continue;}
	      if(Segment_Between < (np[image][itr]-1)/2)
		{
		  // don't devide track on two parts!!! But delete Segments in between!!
		  jj = 0;
		  for (j = Segment_Between+1; j < np[image][itr]; j++)
		    {
		      X[image][itr][jj] = X[image][itr][j];
		      Y[image][itr][jj] = Y[image][itr][j];
		      if(j != np[image][itr] - 1)
			{
			  for (int j1 = 0; j1 < Nmain[image][itr][j]; j1++)
			    {
			      Xmain[image][itr][jj][j1] = Xmain[image][itr][j][j1];
			      Ymain[image][itr][jj][j1] = Ymain[image][itr][j][j1];
			    }
			  Nmain[image][itr][jj] = Nmain[image][itr][j];
			}
		      jj++;
		    }
		  np[image][itr] = jj;
		}
	      else 
		{
		  // don't devide track on two parts!!! But delete Segments in between!!
		  np[image][itr] = Segment_Between + 1;
		}
	    }
	}
    }
  if(ntr[image] < 2) Min_Vertex_Error[image] = 10000.;
  //if(NumberOfTracks != 0 && ntr[image] < MAXNUMBEROFTRACKS  - 1)
  // goto LABELVERT;
}
void FindPairsOnImages()
{
  int image, jtr, i, j;
  int f[MAXNUMBEROFTRACKS];
  for(image=0;image<2; image++)
    {
      i = 0;
      int fp, lp, itr;
      //delete one of parallel track on the same image
      for(itr=0; itr<ntr[image]; itr++)
	{
	  if(mark[image][itr] == 1)
	    {
	      float dist_vert0 = sqrt(pow(((float) Y[image][itr][0]) - Y_Vert[image], 2)+
				      pow(((float) X[image][itr][0]) - X_Vert[image], 2));
	      float dist_vert1 = sqrt(pow(((float) Y[image][itr][np[image][itr] - 1]) - Y_Vert[image], 2)+
				      pow(((float) X[image][itr][np[image][itr] - 1]) - X_Vert[image], 2));
	      if(dist_vert0 < dist_vert1) {fp = 0; lp = np[image][itr] - 1;}
	      if(dist_vert1 < dist_vert0) {fp = np[image][itr] - 1; lp = 0;}
	      angle[image][itr] = TMath::ATan2((Double_t) (Y[image][itr][lp] - Y[image][itr][fp]),
					     (Double_t) (X[image][itr][lp] - X[image][itr][fp]));
	      if(angle[image][itr] < 0.) angle[image][itr] += 2.*TMath::Pi();
	      printf("Image %d Angle of track %d is %f\n", image, itr, angle[image][itr]);
	    }
	}
      double CosPhi_min = cos(TMath::Pi()/60.); //3 degrees!!
      for(itr=0; itr<ntr[image] - 1; itr++)
	{
	  if(mark[image][itr] == 1)
	    for(jtr=itr+1; jtr<ntr[image]; jtr++)
	      {
		if(mark[image][jtr] == 1)
		  {
		    double CosPhi = cos(angle[image][itr])*cos(angle[image][jtr]) + 
		      sin(angle[image][itr])*sin(angle[image][jtr]);
		    if(CosPhi > CosPhi_min) 
		      {
			printf(" There are paraller tracks, remain only the longest\n");
			printf(" Image %d itr %d jtr %d cophi %f %f np %d %d\n",
			       image, itr, jtr, CosPhi, CosPhi_min, np[image][itr], np[image][jtr]);
			if(np[image][itr] >= np[image][jtr]) mark[image][jtr] = 0;
			if(np[image][itr] < np[image][jtr]) mark[image][itr] = 0;
		      }
		  }
	      }
	}
    }

  for(image=0;image<2; image++)
    {
      i = 0;
      int fp, lp, itr;
      for(itr=0; itr<ntr[image]; itr++)
	{
	  f[itr] = 0;
	  if(mark[image][itr] == 1)
	    {
	      float dist_vert0 = sqrt(pow(((float) Y[image][itr][0]) - Y_Vert[image], 2)+
				      pow(((float) X[image][itr][0]) - X_Vert[image], 2));
	      float dist_vert1 = sqrt(pow(((float) Y[image][itr][np[image][itr] - 1]) - Y_Vert[image], 2)+
				      pow(((float) X[image][itr][np[image][itr] - 1]) - X_Vert[image], 2));
	      if(dist_vert0 < dist_vert1) {fp = 0; lp = np[image][itr] - 1;}
	      if(dist_vert1 < dist_vert0) {fp = np[image][itr] - 1; lp = 0;}
	      angle[image][i] = TMath::ATan2((Double_t) (Y[image][itr][lp] - Y[image][itr][fp]),
					     (Double_t) (X[image][itr][lp] - X[image][itr][fp]));
	      if(angle[image][i] < 0.) angle[image][i] += 2.*TMath::Pi();
	      printf("Image %d Angle of track %d is %f\n", image, i, angle[image][i]);
	      Track[image][i] = itr;
	      i++;
	    }
	}
      int itr_min = 0;
      for(int ii = 0; ii < i; ii++)
	{
	  float min_angle = 1000.;
	  for(itr = 0; itr < i; itr++)
	    {
	      if(angle[image][itr] < min_angle && f[itr] == 0) 
		{min_angle = angle[image][itr]; nomer_tr[image][ii] = Track[image][itr]; 
		itr_min = itr;}
	    }
	  f[itr_min] = 1;
	}
      //for(int ii=0;ii<i;ii++)
      //printf("Tracks after sorting: %d is %d\n", ii, nomer_tr[image][ii]);
      NTracks[image] = i;
    }
  for(i=0; i<NTracks[1];i++)
    f[i] = 0;
  for(i=0;i<NTracks[0];i++)
    {
      pair[i] = 100;
      double CosPhi_min = cos(TMath::Pi()/9.); //20 degrees!!
      for(j=0;j<NTracks[1]; j++)
	{
	  if(f[j] == 0)
	    {
	      double CosPhi = cos(angle[0][i])*cos(angle[1][j]) + 
		sin(angle[0][i])*sin(angle[1][j]);
	      if(CosPhi > CosPhi_min) {CosPhi_min = CosPhi; pair[i] = j;}
	    }
	}
      if (pair[i] != 100) f[pair[i]] = 1;
      if (pair[i] == 100) mark[0][Track[0][i]] = 0;
      if (pair[i] != 100) printf(" Track %d(image 0) corresponds to track %d(image 1)\n",
	     Track[0][i], Track[1][pair[i]]); 
    }
  for(j=0;j<NTracks[1]; j++)
    if(f[j] == 0) mark[1][Track[1][j]] = 0;
}
void FindNewVertex(int image)
{
  float xc[20], yc[20], rc[20];
  float Xvertex, Yvertex, Dvertex;
  int itr, i, j, jj, NGoodTracks[2];
  int jtr = 0; //number of good tracks
  for (itr = 0; itr < ntr[image]; itr++)
    {
      if(mark[image][itr] == 0) continue;
      jj = 0;
      float WpMin = 0.;
      float WpMax = 4096.;
      for(i=0; i<np[image][itr]-1; i++) 
	{
	  for(j=0; j<Nmain[image][itr][i]; j++) 
	    {
	      Xcircle[jj] = (float) Xmain[image][itr][i][j];
	      Ycircle[jj] = (float) Ymain[image][itr][i][j];
	      int iy = vsize - Ymain[image][itr][i][j] - 1;
	      int ix = Xmain[image][itr][i][j];
	      Wpoint[jj] = (float) pCDBuffer[image][iy][ix];
	      if (Wpoint[jj] > WpMax) WpMax = Wpoint[jj];
	      if (Wpoint[jj] < WpMin) WpMin = Wpoint[jj];
	      jj++;
	    }
	}
      for(i = 0; i< jj; i++) 
	{
	  if(WpMin != WpMax) Wpoint[i]  = (Wpoint[i]-WpMin)/(WpMax - WpMin);
	  if(WpMin == WpMax) Wpoint[i]  = Wpoint[i]/WpMax;
	}
      //xcD[image][itr] = ycD[image][itr] = rcD[image][itr] = 0.;
      int ierr = circle_r(jj);
      printf("Track %d:  XC %f YC %f RD %f ch2 %f\n", itr, XC, YC, RC, chi2);
      Npoints_fit = jj;
      circle_minuit(jj);
      xc[jtr] = xcD[image][itr] = XC;
      yc[jtr] = ycD[image][itr] = YC;
      rc[jtr] = rcD[image][itr] = RC;
      chi2_circle[image][jtr] = chi2;
      for(i=0;i<jj; i++)
	{
	  double buf = atan2(Ycircle[i] - YC, Xcircle[i] - XC);
	  x_mark[i] = XC + RC*cos(buf);
	  y_mark[i] = YC + RC*sin(buf);
	}
      px_mark = &x_mark[0];
      py_mark = &y_mark[0];
      marke = new TPolyMarker((Int_t) jj, px_mark, py_mark, "same");
      fObjTbl->Add(marke);
      marke->SetMarkerStyle(8);
      marke->SetMarkerSize(0.1);
      marke->SetMarkerColor(40);
      if(image == 0) pleft->cd();
      if(image == 1) pright ->cd();
      marke->Draw("same");	     
      jtr++;
    }
  NGoodTracks[image] = jtr;
  Xvertex = Yvertex = Dvertex = 0.;
  XY_Vertex_Net(NGoodTracks[image], &xc[0], &yc[0], 
		&rc[0], &Xvertex, &Yvertex, 
		&Dvertex);
  Min_Vertex_Error[image] = Dvertex;
  X_Vert[image] = Xvertex;
  Y_Vert[image] = Yvertex;
  printf(" Image %d Final Track Number is %d\n", image, NGoodTracks[image]);
  printf(" Image %d Final Vertex position is : X %f Y %f Error is %f\n", 
	 image, X_Vert[image], Y_Vert[image], Dvertex);
}
int FindFirstPion(void)
{
  int image, itr, jtr;
  int FirstPion[2][20];
  for(image=0;image<2; image++)
    {
      int fp, lp, itr;
      //delete one of parallel track on the same image
      for(itr=0; itr<ntr[image]; itr++)
	{
	  //if(mark[image][itr] == 1)
	    {
	      float dist_vert0 = sqrt(pow(((float) Y[image][itr][0]) - Y_Vert[image], 2)+
				      pow(((float) X[image][itr][0]) - X_Vert[image], 2));
	      float dist_vert1 = sqrt(pow(((float) Y[image][itr][np[image][itr] - 1]) - Y_Vert[image], 2)+
				      pow(((float) X[image][itr][np[image][itr] - 1]) - X_Vert[image], 2));
	      if(dist_vert0 < dist_vert1) {fp = 0; lp = np[image][itr] - 1;}
	      if(dist_vert1 < dist_vert0) {fp = np[image][itr] - 1; lp = 0;}
	      Angle[image][itr] = TMath::ATan2((Double_t) (Y[image][itr][lp] - Y[image][itr][fp]),
					     (Double_t) (X[image][itr][lp] - X[image][itr][fp]));
	      FirstPion[image][itr] = 0;
	      if(fabs(Angle[image][itr]) < TMath::Pi()/10. ) FirstPion[image][itr] = 1;

	      if(Angle[image][itr] < 0.) Angle[image][itr] += 2.*TMath::Pi();
	      printf("Image %d Angle of track %d is %f\n", image, itr, Angle[image][itr]);
	    }
	}
    }
  double CosPhi_min = cos(TMath::Pi()/10.); //18 degrees!!
  PionTrack[0] = PionTrack[1] = -1;
  for(itr=0; itr<ntr[0]; itr++)
    {
      if(mark[0][itr] == 1 && FirstPion[0][itr] == 1)
	for(jtr=0; jtr<ntr[1]; jtr++)
	  {
	    if(mark[1][jtr] == 1 && FirstPion[1][jtr] == 1)
	      {
		double CosPhi = cos(Angle[0][itr])*cos(Angle[1][jtr]) + 
		  sin(Angle[0][itr])*sin(Angle[1][jtr]);
		if(CosPhi > CosPhi_min) 
		  {
		    CosPhi_min = CosPhi;
		    printf("Pion Pairs itr %d jtr %d cophi %f %f Angles %f %f\n",
			   itr, jtr, CosPhi, CosPhi_min, Angle[0][itr], Angle[1][jtr]);
		    PionTrack[0] = itr;
		    PionTrack[1] = jtr;
		  }
	      }
	}
    }
  if(CosPhi_min > cos(TMath::Pi()/10.)) return(1); // pairs of Final Pion has been found
  return(0);
}
void DigitizeEvent()
{
  int i, image, track, itr, jtr;

  int ierr;
  float buf;
  float bufmin, bufmax;
  int imin, imax;
  double x_cur, y_cur;
  double phi1, phi2, range, stangle, endangle;
  //float xcenter, ycenter, rcircle;
  double xc_int, yc_int, rc_int;
  double stx, sty, endx, endy;
  float brrr[2][6];
  int NP[2][6];
  float xcent_vx[10], ycent_vx[10], rcent_vx[10];
  Int_t binx, biny, bin_arr;
  float Xvertax, Yvertax, Dvertax;

  int DigTrack[2][MAXNUMBEROFTRACKS], DigTracks;
  int tookit[MAXNUMBEROFTRACKS], itrmin;

  TAxis *xaxis;
  TAxis *yaxis;
  float sumbright_all;
  int binxmax, binymax;
  float rmax;
  int bright;

  //first reorder tracks in images in accordance with clock arrow
  //(first must be pion track)
  for (image = 0; image < 2; image++)
    for(itr = 0; itr < ntr[image]; itr++)
      DigTrack[image][itr] = -1;
  DigTracks = 0;
  if(PionTrack[0] != -1 && PionTrack[1] != -1) 
    {
      DigTrack[0][0] = PionTrack[0];
      DigTrack[1][0] = PionTrack[1];
      printf(" Pion Tracks %d %d will be %d in digitizing\n", 
	     PionTrack[0], PionTrack[1], DigTracks);
      DigTracks++;
    }
  for(itr = 0; itr < NTracks[0]; itr++) tookit[itr] = 0;
  for(itr = 0; itr < NTracks[0]; itr++)
    {
      float MinAngle = 1000.;
      //if((Track[0][itr] == PionTrack[0])) continue;
      for(jtr = 0; jtr < NTracks[0]; jtr++)
	{
	  if(pair[jtr] == 100) continue;
	  if((Track[0][jtr] == PionTrack[0])) continue;
	  //if(itr == jtr) continue;
	  if(tookit[jtr] == 1) continue;
	  if(Angle[0][Track[0][jtr]] < MinAngle) 
	    {MinAngle = Angle[0][Track[0][jtr]];
	    itrmin = jtr; }
	}
      if(MinAngle == 1000.) continue;
      tookit[itrmin] = 1;
      DigTrack[0][DigTracks] = Track[0][itrmin];
      DigTrack[1][DigTracks] = Track[1][pair[itrmin]];
      printf(" Angle %f Tracks %d %d will be %d in digitizing\n", 
	     MinAngle, Track[0][itrmin], Track[1][pair[itrmin]], DigTracks);
      DigTracks++;
    }
  /*   searching for empty filename   */
  char ctmp[3];
  int itemp = 0;
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
  fprintf(out_file,"%s\n", eventname);
  fprintf(out_file," %d      %d \n", GoodEvent, NGoodTracks[0]);
  for (image = 0; image < 2; image++)
    {
      track = 0;
      if(image == 0)
	{
	  pleft ->cd();
	  xaxis = left -> GetXaxis();
	  yaxis = left -> GetYaxis();
	}
      else
	{
	  pright->cd();
	  xaxis = right -> GetXaxis();
	  yaxis = right -> GetYaxis();
	}
      //for (itr=0; itr<ntr[image]; itr++)
      for (int idig=0; idig < DigTracks; idig++)
	{
	  itr = DigTrack[image][idig];
	  printf( " Image %d track %d itr %d mark %d\n", 
		  image, idig, itr, mark[image][itr]);
	  if(mark[image][itr] == 0) continue;
	  //  printout to the output file 

	  fprintf(out_file," %d \n",track);
	  // suppose that start and end angles are between 180degrees
	  bufmin = 1000.;
	  bufmax =-1000.;
	  for(i=0;i<np[image][itr]; i++)
	    {
	      buf = atan2(((float) Y[image][itr][i]) - ycD[image][itr],
			  ((float) X[image][itr][i]) - xcD[image][itr]);
	      if(buf < bufmin)
		{imin = i; bufmin = buf;}
	      if(buf > bufmax)
		{imax = i; bufmax = buf;}
	    }
	  //check for the -Pi - +Pi bounds
	  if((bufmax - bufmin) > (float)TMath::Pi()) 
	    {
	      bufmin = 1000.;
	      bufmax =-1000.;
	      for(i=0;i<np[image][itr]; i++)
		{
		  buf = atan2(((float) Y[image][itr][i]) - ycD[image][itr],
			      ((float) X[image][itr][i]) - xcD[image][itr]);
		  if (buf < 0.) buf += 2.*(float)TMath::Pi();
		  if(buf < bufmin)
		    {imin = i; bufmin = buf;}
		  if(buf > bufmax)
		    {imax = i; bufmax = buf;}
		}
	    }
	  stx = (double) X[image][itr][imin];
	  sty = (double) Y[image][itr][imin];
	  endx = (double) X[image][itr][imax];     
	  endy = (double) Y[image][itr][imax];
	  printf("imin %d imax %d\n", imin, imax);
	  rc_int = (double) rcD[image][itr]; 
	  xc_int = (double) xcD[image][itr]; 
	  yc_int = (double) ycD[image][itr];

	  range = sqrt((endx-stx)*(endx-stx) + (endy-sty)*(endy-sty));
	  range = range/(2.*rc_int);
	  if(range > 1.) range = 1.;
	  range = 2.*asin(range);
      
	  //  simplest way to get angle in 2Pi
	  stangle = atan2((sty-yc_int),(stx-xc_int));
	  endangle = atan2((endy-yc_int),(endx-xc_int));
	  if (stangle > endangle) endangle += 2.*TMath::Pi();
	  phi1 = stangle*180./TMath::Pi(); 
	  phi2 = endangle*180./TMath::Pi();
	  int epx = (int) endx;
	  int epy = (int) endy;
	  
	  double delta_angle = range/1000.;
	  double cur_angle = stangle;
	  double nat_log = 2.718281828;
      
	  char ctmp1[3],ctmp2[3], name[10];
	  sprintf(ctmp1,"%u",image);
	  sprintf(ctmp2,"%u",track);
	  strcpy(name,"h");
	  strcat(name,ctmp1);
	  strcat(name,"tr");
	  strcat(name,ctmp2);
	  hbr = new TH1F(name,name, 41,-20.,20.);
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
	  for(ii = 0; ii<1000000; ii++)
	    {
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
	      if (binxprev_circle != binx || binyprev_circle != biny)
		{
		  // now go along the radius
		  binxprev_circle = binx; binyprev_circle = biny;
		  //bright = PixelBr(image, ypix, xpix);
		  if ( image == 0)
		    {
		      bin_arr = left->GetBin(binx,biny);
		      bright = (int) left->GetBinContent(bin_arr);
		    }
		  else
		    {
		      bin_arr = right->GetBin(binx,biny);
		      bright = (int) right->GetBinContent(bin_arr);
		    }
		  int ifl = 0;
		  if(ifl == 0 && bright > wpmax) {
		    wpmax = bright; xpmax = xpix; 
		    ypmax = ypix; rmax = rc_int;
		    binxmax = binx; binymax = biny; }
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
		      //bright = PixelBr(image, ypix, xpix);
		      if ( image == 0)
			{
			  bin_arr = left->GetBin(binx,biny);
			  bright = (int) left->GetBinContent(bin_arr);
			}
		      else
			{
			  bin_arr = right->GetBin(binx,biny);
			  bright = (int) right->GetBinContent(bin_arr);
			}
		      int ifl = 0;
		      if(ifl == 0 && bright > wpmax) {
			wpmax = bright; xpmax = xpix; 
			ypmax = ypix; rmax = rcur;
			binxmax = binx; binymax = biny; }
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
		      //bright = PixelBr(image, ypix, xpix);
		      if ( image == 0)
			{
			  bin_arr = left->GetBin(binx,biny);
			  bright = (int) left->GetBinContent(bin_arr);
			}
		      else
			{
			  bin_arr = right->GetBin(binx,biny);
			  bright = (int) right->GetBinContent(bin_arr);
			}
		      int ifl = 0;
		      if(ifl == 0 && bright > wpmax) {
			wpmax = bright; xpmax = xpix; 
			ypmax = ypix; rmax = rcur;
			binxmax = binx; binymax = biny; }
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
		  background_cadr = 0.75*main_aver[image];
		  if( npoints > 10000) break;
		  ip = 0;
		  iwidth = 5;
		  hbr->Fill((Axis_t) ip, (Stat_t) wpmax);
		  wt[image][track][npoints] = 
		    (((float) wpmax) -  background_cadr);
		  wpmax = wp[ip] = (int) wt[image][track][npoints];
		  sumbright_all += (float) wpmax;
		  //xp[npoints_bright] = 
		  xt[image][track][npoints] = (float) xpmax;
		  //yp[npoints_bright] = 
		  yt[image][track][npoints] = (float) ypmax;
		  if(wpmax < 0 ) {wt[image][track][npoints]  = 0.; wp[ip]  = wpmax = 0;}
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
			//bright = PixelBr(image, ypix, xpix);
			if ( image == 0)
			  {
			    bin_arr = left->GetBin(binx,biny);
			    bright = (int) left->GetBinContent(bin_arr);
			  }
			else
			  {
			    bin_arr = right->GetBin(binx,biny);
			    bright = (int) right->GetBinContent(bin_arr);
			  }
			hbr->Fill((Axis_t) ip, (Stat_t) bright);
			binxprev = binx; binyprev = biny;
			int ifl = 0;
			if(ifl == 0) {
			  wp[ip] = (int) (((float) bright) - background_cadr);
			  sumbright_all += wp[ip];;
			  if(wp[ip] < 0) wp[ip] = 0;
			  if(ip < iwidth && npoints < 10000 - 1) {
			    xt[image][track][npoints] = (float) xpix; 
			    yt[image][track][npoints] = (float) ypix;
			    wt[image][track][npoints]  = (float) wp[ip];
			    x_mark[npoints] = (Float_t) xpix;
			    y_mark[npoints] = (Float_t) ypix;
			    //Draw_region(image, xpix, ypix, (int) wt[image][track][npoints]);
			    npoints++;
			  }
			  //else
			    //Draw_region(image, xpix, ypix, (int) wp[ip]);
			  npoints_bright++; ip++;
			  if(wp[ip - 1] == 0) break;
			  if(((float) wpmax)/((float) wp[ip - 1]) > nat_log) break;
			  float gradient = (wp[ip - 1] - wp[ip - 2]);
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
			//bright = PixelBr(image, ypix, xpix);
			if ( image == 0)
			  {
			    bin_arr = left->GetBin(binx,biny);
			    bright = (int) left->GetBinContent(bin_arr);
			  }
			else
			  {
			    bin_arr = right->GetBin(binx,biny);
			    bright = (int) right->GetBinContent(bin_arr);
			  }
			binxprev = binx; binyprev = biny;
			hbr->Fill((Axis_t) -ip, (Stat_t) bright);
			int ifl = 0;
			if(ifl == 0) {
			  wp[ip] = (int) (((float) bright) - background_cadr);
			  sumbright_all += wp[ip];;
			  if(wp[ip] < 0) wp[ip] = 0;
			  if(ip < iwidth && npoints < 10000 - 1) {
			    xt[image][track][npoints] = (float) xpix; 
			    yt[image][track][npoints] = (float) ypix;
			    wt[image][track][npoints]  = (float) wp[ip];
			    x_mark[npoints] = (Float_t) xpix;
			    y_mark[npoints] = (Float_t) ypix;
			    //Draw_region(image, xpix, ypix, (int) wt[image][track][npoints]);
			    npoints++;
			  }
			  //else
			  //Draw_region(image, xpix, ypix, (int) wp[ip]);
			  npoints_bright++; ip++;
			  if(wp[ip - 1] == 0) break;
			  if(((float) wpmax)/((float) wp[ip - 1]) > nat_log) break;
			  float gradient = (wp[ip - 1] - wp[ip - 2]);
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
	  hbr->Scale(1./(float) npp);
	  float bright_av_all;
	  //bright_av_all = sumbright_all/(float) npoints_bright;
	  bright_av_all = sumbright_all/(float) npp;
	  brrr[image][track] = bright_av_all*(float) npoints_bright/(float) npp;
	  printf("average brightness on track %d is %f\n", 
		 track, bright_av_all);
	  NP[image][track] = npoints;// 100;
	  float WpMin = 0.;
	  float WpMax = 4096.;
	  for(i=0; i<NP[image][track]; i++) 
	    {
	      Xcircle[i] = (float) xt[image][track][i];
	      Ycircle[i] = (float) yt[image][track][i];
	      Wpoint[i] = ((float)  wt[image][track][i]);
	      if (Wpoint[i] > WpMax) WpMax = Wpoint[i];
	      if (Wpoint[i] < WpMin) WpMin = Wpoint[i];
	    }
	  for(i = 0; i< NP[image][track]; i++) 
	    {
	      if(WpMin != WpMax) Wpoint[i]  = (Wpoint[i]-WpMin)/(WpMax - WpMin);
	      if(WpMin == WpMax) Wpoint[i]  = Wpoint[i]/WpMax;
	    }
 	  px_mark = &x_mark[0];
	  py_mark = &y_mark[0];
	  marke = new TPolyMarker((Int_t) npoints, px_mark, py_mark, "same");
	  fObjTbl->Add(marke);
	  marke->SetMarkerStyle(8);
	  marke->SetMarkerSize(0.2);
	  marke->SetMarkerColor(41);
	  if(image == 0) pleft->cd();
	  if(image == 1) pright ->cd();
	  marke->Draw("same");	     

	  ierr = circle_r(NP[image][track]);
	  printf("C2 FINAL: Image %d track %d rad %f xc %f yc %f chi2 = %f\n",
		 image, track, RC, XC, YC, chi2);	 
	  Npoints_fit = NP[image][track];
	  circle_minuit(NP[image][track]);
	  {
	    // suppose that start and end angles are between 180degrees
	      bufmin = 1000.;
	      bufmax =-1000.;
	      for(i=0;i<NP[image][track]; i++)
		{
		  buf = atan2(Ycircle[i] - YC, Xcircle[i] - XC);
		  x_mark[i] = XC + RC*cos(buf);
		  y_mark[i] = YC + RC*sin(buf);
		  //draw circle for the corridor
		  if(buf < bufmin)
		    {imin = i; bufmin = buf;}
		  if(buf > bufmax)
		    {imax = i; bufmax = buf;}
		}
	      px_mark = &x_mark[0];
	      py_mark = &y_mark[0];
	      marke = new TPolyMarker((Int_t) NP[image][track], px_mark, py_mark, "same");
	      fObjTbl->Add(marke);
	      marke->SetMarkerStyle(8);
	      marke->SetMarkerSize(0.1);
	      marke->SetMarkerColor(51);
	      marke->Draw("same");	     
	      //check for the -Pi - +Pi bounds
	      if((bufmax - bufmin) > (float)TMath::Pi()) 
		{
		  bufmin = 1000.;
		  bufmax =-1000.;
		  for(i=0;i<NP[image][track]; i++)
		    {
		      buf = atan2(Ycircle[i] - YC, Xcircle[i] - XC);
		      if (buf < 0.) buf += 2.*(float)TMath::Pi();
		      if(buf < bufmin)
			{imin = i; bufmin = buf;}
		      if(buf > bufmax)
			{imax = i; bufmax = buf;}
		    }
		}
	      stx = Xcircle[imin];      sty = Ycircle[imin];
	      endx = Xcircle[imax];     endy = Ycircle[imax];
	    }
	  // for vertex search store
	  xcent_vx[track] = XC;
	  ycent_vx[track] = YC;
	  rcent_vx[track] = RC;
	  //  printout to the output file 	 
	  fprintf(out_file,"%7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f \n",
		  RC, XC, YC, stx, sty, endx, endy);
	  fprintf(out_file,"%7d\n", NP[image][track]);
	  int k = 0;
	  for ( int ic=0; ic<NP[image][track]; ic++)
	    {
	      k++;
	      bright = 0;
	      binx = (int) xt[image][track][ic]; 
	      biny = (int) yt[image][track][ic]; 
	      if ( image == 0)
		{
		  bin_arr = left->GetBin(binx,biny);
		  bright = (int) left->GetBinContent(bin_arr);
		}
	      else
		{
		  bin_arr = right->GetBin(binx,biny);
		  bright = (int) right->GetBinContent(bin_arr);
		}
	      fprintf(out_file," %u %u %d ", binx, biny, (int) wt[image][track][ic]);
	      if(k == 5) { k = 0;    fprintf(out_file,"\n"); }
	    }
	  track++;
	  //end over tracks loop
	}
      num_tracks = track;
      //   vertex finding
      Xvertax = Yvertax = Dvertax = 0.;
      if(num_tracks != 1)
	{
	  XY_Vertex_Net(num_tracks, &xcent_vx[0], &ycent_vx[0], 
			&rcent_vx[0], &Xvertax, &Yvertax, 
			&Dvertax);
	  TMarker *mark = new TMarker(Xvertax,Yvertax,5);
	  mark->SetMarkerColor(50);
	  mark->SetMarkerStyle(5);
	  mark->SetMarkerSize(1.5);
	  mark->DrawMarker(Xvertax, Yvertax);
	  //gPad->Modified();
	  //					  gPad->Update();
	  
	  printf(" Vertex coordinates: x %f, y %f; error %f\n",
		 Xvertax, Yvertax, Dvertax );
	}
      c2->Modified(); 
      c2->Update();
      fprintf(out_file,
	      " Vertex position: X: %f, Y: %f, Error: %f \n", 
	      Xvertax, Yvertax, Dvertax);
    }
  fprintf(out_file,
	  " %f %f \n", 
	  main_aver[0], main_aver[1]);
  for(i = 0; i< num_tracks; i++)
    {
      fprintf(out_file,
	      " %f %f \n", 
	      brrr[0][i], brrr[1][i]);
    }
  fclose(out_file);
  ch = new TCanvas("chisto","Histos of brightness",30,30,500,600);
  ch->Divide(2,1);
  for(i=0;i<2;i++){
    ch->cd(i+1);
    hp0 = (TPad*)gPad;
    gPad->UseCurrentStyle();
    hp0->Divide(1,num_tracks); 
    hp0->UseCurrentStyle();
    for(track = 0;track<num_tracks; track++){
      hp0->cd(track+1);
      gPad->SetTopMargin(0.12);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.12);
      gPad->SetBottomMargin(0.12);
      gPad->SetFillColor(33);
      gPad->SetLogy();
      char ctmp1[3],ctmp2[3], name[10];
      sprintf(ctmp1,"%u",i);
      sprintf(ctmp2,"%u",track);
      strcpy(name,"h");
      strcat(name,ctmp1);
      strcat(name,"tr");
      strcat(name,ctmp2);
      TH1F *h1 = (TH1F*)gDirectory->FindObject(name);
      h1->SetMinimum(100.0);
      h1->SetMaximum(4096.0);
      h1->Draw();
    }
  }
  ch->Update();
  // run the Gil Pontecorvo's programme
  //gil_main(test_outfile, 100);
}


int circle_r (int N)
{
//   A routine for fast fit of a set of measured points by a circle.
//   The points are given via their coordinates X(i),Y(i), i=1,2,...,N.
//   The circle will be defined as an equation
//                2       2    2
//           (X-A) + (Y-B)  = R  ,
//   where  A,B  are the coordinates of the centre and  R  is the radius.
//   The routine uses a fast version of the weighted least-squares fit.
//     ========= Parameters:
//   On entry:
//               N - the number of measured points (an integer variable);
//               X - the real array of x-coordinates of the points;
//               Y - the real array of y-coordinates of the points;
//               W - the real array of weights of the points.
//                   W = 1/(Sigma**2)
//         All these are unchanged on exit.
//   comment: the weigth of a measured point usually corresponds to the
//           statistical error of the measurement (a zero or negative
//           weight causes the correspondent data point to be ignored).
//           If you have no such an array, then simply set W(i)=1 for
//           all i=1,2,...,n in your calling (sub)routine.
//           You may also exclude the array W from the routine according
//           to the comments in the text below.
//   On successful exit:
//             A,B  - the coordinates of the centre (real);
//             R    - the value of radius (real);
//             VAR  - the estimate of variance of the statistical errors
//                     (mean square of deviations of the measured points
//                      about the obtained circle);
//             RES  - the vector of (approximate) residuals (a real array).
//   comment: often one does not need to evaluate the variance and the
//            residual vector. In this case one may set the integer
//            parameter NOSTAT to 1 and then the routine will not
//            calculate the value of VAR and the array RES speeding up
//            its main work. Any other value of NOSTAT causes the value
//            of VAR and the array RES to be evaluated.
//   comment: the routine needs several local real variables to be in
//            high precision form (at least 8 bytes per variable).
//            To this aim they are declared as DOUBLE PRECISION in
//            the text below and the routine uses base external
//            functions DSQRT,DABS. But if you work on a 32-bit computer,
//            then remove the declaration DOUBLE PRECISION and use the
//            base external functions SQRT,ABS.
//   Possible errors:
//            If you obtain IERR=0 on exit, then no errors occured.
//            Otherwise you should correct your data according to the
//            following rules:
//   IERR=1:  the number of points with non-zero weight is less than 3
//            (routine needs at least 3 points to construct a circle!);
//   IERR=2:  all the points lie in a straight line; in this case you
//            should fit them by a straight line, not by a circle!
//            BUT the straight line A*x+B is fitted, assuming the case
//            of the robust ellastic fitting.
//   IERR=3:  the case of divergency of the Newton iterational process;
//            in this case the routine will work, too, but with less
//            accuracy.
//          The routine uses the method of N.I.Chernov and G.A.Ososkov
//          published in Comp.Phys.Comm., v.33 (1984), p.329-333.
//          Written by N.Chernov
  double  x, y, w;
  double XX,YY,XY,F,G,H,HH,P,Q,T,GAMIN,GAMMA,FACTOR;
  double DET,EPS,A0,A1,A2,A2A2,XOLD,YOLD,XNEW,DY;
  double XI, YI, DI;
  int i, ITER;
  int IERR;
  double a, b, r;
  //  double res_array[100];
  double R2, RR, res;
    //*RES;
  //DIMENSION X(N),Y(N),W(N),RES(N);
  int ITMAX = 100;
  EPS = 1.e-12;
  XNEW = 0.;

  IERR = 0;

//         clear all accumulators
  int NPOS = 0;
  double SUMW = 0.;
  double XCEN = 0.;
  double YCEN = 0.;
  double SUMRES = 0.;
  XX = 0.;
  YY = 0.;
  XY = 0.;
  P = 0.;
  Q = 0.;
  T = 0.;
  //       compute the centre of gravity (weighted centre) XCEN,YCEN;
  for(i = 0; i < N; i++)
    {
      x =(double) Xcircle[i];//*(Xcircle+i);
      y =(double) Ycircle[i];//*(Ycircle+i);
      w =(double) Wpoint[i];//*(Wpoint+i);
      if (w > 0.) 
	{
	  XCEN += x*w;
	  YCEN += y*w;
	  SUMW += w;
	  NPOS += 1;
	}
    }
//   if you do not want to use weights, then replace the above
//   block with the following one:
//      DO 1 I = 1,N
//         XCEN = XCEN+X(I)
//         YCEN = YCEN+Y(I)
//   1  CONTINUE
//      SUMW = N
//      NPOS = N
 if (NPOS <= 2) return 1;

 XCEN /= SUMW;
 YCEN /= SUMW;
 //          compute the power series of the point coordinates
 for (i = 0; i<N; i++)
   {
     x =(double) Xcircle[i];//*(Xcircle+i);
     y =(double) Ycircle[i];//*(Ycircle+i);
     w =(double) Wpoint[i];//*(Wpoint+i);
     if (w > 0.) 
       {
	 XI = x-XCEN;
	 YI = y-YCEN;
	 DI = XI*XI+YI*YI;
	 XX += XI*XI*w;
	 YY += YI*YI*w;
	 XY += XI*YI*w;
	 P += XI*DI*w;
	 Q += YI*DI*w;
	 T += DI*DI*w;
       }
   }
 //   if you do not want to use weights, then remove all
 //   the multiplications by W(I) and the IF statement
 //   in the above lines.
 //               normalisation of the power series
 XX /= SUMW;
 YY /= SUMW;
 XY /= SUMW;
 P /= SUMW;
 Q /= SUMW;
 T /= SUMW;
 //            compute the coefficients of the main polynomial
 //              A0 + A1*X + A2*X**2 - 4*X**3 + X**4
 F = 3.0*XX+YY;
 G = 3.0*YY+XX;
 H = XY+XY;
 HH = H*H;
 GAMIN = XX+YY;
 //        check that the points with non-zero weight do not coinside
 if (GAMIN <= 0.) return 2;
 FACTOR = GAMIN*GAMIN;
 A2 = (F*G-HH-T)/FACTOR;
 A2A2 = A2+A2;
 FACTOR = FACTOR*GAMIN;
 A1 = (T*(F+G)-2.0*(P*P+Q*Q))/FACTOR;
 FACTOR = FACTOR*GAMIN;
 A0 = (T*(HH-F*G)+2.0*(P*P*G+Q*Q*F)-4.0*P*Q*H)/FACTOR;
 //       find the root of the main polynomial closest to 1
 //       by the Newton method; XNEW will contain the value of root.
 XOLD = 1.0;
 ITER = 0;
 LABEL3:
 ITER += 1;
 YOLD = A0+XOLD*(A1+XOLD*(A2+XOLD*(XOLD-4.0)));
 DY = A1+XOLD*(A2A2+XOLD*(4.0*XOLD-12.0));
 if (fabs(DY) < EPS) goto LABEL4;
 XNEW -= YOLD/DY;
 // cout << "ITER" << ITER << " XNEW " << XNEW << endl;
 // cout << "fabs(XOLD-XNEW)" << fabs(XOLD-XNEW) << " XOLD " << XOLD << endl;
 if (XNEW < EPS) goto LABEL4;
 if (fabs(XOLD-XNEW) < EPS) goto LABEL5;
 XOLD = XNEW;
 if (ITER < ITMAX) goto LABEL3;
 //      the case of divergency: the root is approximately set to 1.
 LABEL4:
 XNEW = 1.0;
 IERR = 3;
 //          compute the parameters of the obtained circle
 LABEL5:
 GAMMA = GAMIN*XNEW;
 DET = (F-GAMMA)*(G-GAMMA)-HH;
 //        check that the data points are not in line.
 if (fabs(DET) <= EPS) return 2;

 a = (P*(G-GAMMA)-Q*H)/DET;
 b = (Q*(F-GAMMA)-P*H)/DET;
 r = sqrt(GAMMA+a*a+b*b);
 a += XCEN;
 b += YCEN;
 XC = (float) a;
 YC = (float) b;
 RC = (float) r;
 //        check that the statistics (VAR,RES) are to be computed
 // if (NOSTAT==1) return 0;
 //       compute the residual vector
 RR = r*r;
 R2 = r*2.0;
 // RES = &res_array[0];
 SUMRES = 0.;
 // printf("a = %f,b = %f, r = %f, RR = %f, R2 = %f\n",a,b,r,RR,R2);
 for(i = 0; i < N; i++)
   {
     x =(double) Xcircle[i];//*(Xcircle+i);
     y =(double) Ycircle[i];//*(Ycircle+i);
     w =(double) Wpoint[i];//*(Wpoint+i);
     res = sqrt(pow(x-a,2) + pow(y-b,2)) - r;
     //((x-a)*(x-a)+(y-b)*(y-b)-RR)/R2;
     //     printf("x = %f,y = %f, res = %f, w = %f\n",x,y,res,w);
     //     res = *(RES+i);
     if (w > 0.) 
	 SUMRES += res*res*w*w;
   }
 //   if you do not want to use weights, then remove
 //   the multiplication by W(I) in the above line.
 //        compute the estimate of variance of statistical errors
 //         (mean square deviation)
 chi2 = (float) (SUMRES/(double) (NPOS-2));
 // printf("chi2 = %f, SUMRES = %f, NPOS = %d\n", VAR, SUMRES, NPOS);
 // cout << "chi2:"<<VAR<<"SUMRES"<<SUMRES<<endl;
 //      VAR = SUMRES/SUMW
 if(IERR == 3) return 3;
 return 0;
}
void circle_minuit (int N)
  //	      float *A, float *B, float *R, float *pchi)
{
  char name1[10],name2[10],name3[10];
  double val1, val2, val3;
  double eval1, eval2, eval3;
  double bnd1, bnd2;
  
  double *pval1=&val1, *pval2=&val2, *pval3=&val3,
    *peval1=&eval1,*peval2=&eval2,*peval3=&eval3,
    *pbnd1=&bnd1,*pbnd2=&bnd2;

  int error_flag=0;
  int ivarbl;
  int *pivarbl=&ivarbl;
  //  void *fcn();
  double f_null =0.;
  MNINIT(5,16,7);
  MNSETI(" just fit ");
  // set of parameters
  MNPARM(1,"X",0.,10.,f_null,f_null,error_flag);
  MNPARM(2,"Y",0.,10.,f_null,f_null,error_flag);
  MNPARM(3,"R",500.,50.,f_null,f_null,error_flag);
  // fitting algorithm
  MNEXCM(fcn,"MIGRAD",0,0,error_flag,0);
  // readout the result
  MNPOUT(1,name1,val1,eval1,bnd1,bnd2,ivarbl);
  MNPOUT(2,name2,val2,eval2,bnd1,bnd2,ivarbl);
  MNPOUT(3,name3,val3,eval3,bnd1,bnd2,ivarbl);
  // printout
  printf(" MINUIT output: chi2 %f \n", fcn_print/((float) N - 2));
  printf("var - %s = %f, err= %f\n",name1,val1,eval1);
  printf("var - %s = %f, err= %f\n",name2,val2,eval2);
  printf("var - %s = %f, err= %f\n",name3,val3,eval3);
  //getc(stdin);
  if(val3 < 10000. && fcn_print != 0. && N > 2)
    {//take the MINUIT values of fit
      XC = val1;
      YC = val2;
      RC = val3;
    }
  else if(RC > 0.) 
    {
      // remain Ososkov's fit values
    }
  else
    {//take the MINUIT values of fit
      XC = val1;
      YC = val2;
      RC = val3;
    }
}
void fcn(int npar, double grad[3], double * fcnval,
           double xval[3],int iflag, void (*Dummy)())
{
  double x0, y0, rad;
  double x, y, w; // current point ...........
  double curr,sum2;
  
  x0 = xval[0];
  y0 = xval[1];
  rad = xval[2];


  switch(iflag){
  case 1:
    //   Initialaise.
    printf(" fcn_c called to init..\n");
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
	      curr = w*w*pow((sqrt(pow(x-x0,2) + pow(y-y0,2)) - rad), 2.);
	      sum2 = sum2 + curr;
	    }
	}
      *fcnval = sum2;
      fcn_print = sum2;
    }
    break;
  }
}
void FindPointsForCircleFit(int image)
{
  int column, row, i, j, itr, jj;
  float MD = 1.2;
  for (itr=ntr[image]; itr < ntr[image]+NumberOfTracks; itr++)
    {
      jj = 0;
      for(i = 0; i < np[image][itr]-1; i++)
	{
	  j = 0;
	  float X_Start = X[image][itr][i];
	  float Y_Start = Y[image][itr][i];
	  float X_End   = X[image][itr][i+1];
	  float Y_End   = Y[image][itr][i+1];
	  float Slope, Y_Origin;
	  if(X_End != X_Start) 
	    {
	      Slope = (Y_End - Y_Start)/(X_End - X_Start);
	      Y_Origin = Y_Start - Slope*X_Start;
	      Slope = TMath::ATan2((Y_End - Y_Start),(X_End - X_Start));
	    }
	  else if(X_End == X_Start)
	    {
	      Slope = TMath::Pi()/2.;
	      Y_Origin = X_End;
	    }
	  else if(Y_End == Y_Start)
	    {
	      Slope = 0.;
	      Y_Origin = Y_End;
	    }
	  for (column=tailx[image][0]; column<hsize-tailx[image][1]; column++)
	    {
	      for (row=taily[image][1]; row<vsize-taily[image][0]-1; row++)
		{
		  if (sosts[column][row] == 1 && rsost[column][row] == 0)
		    {
		      if ((row>taily[image][1] && row<vsize-taily[image][0]-1) &&
			  (column>tailx[image][0] && column<hsize-tailx[image][1]))
			{
			  if(column >= X_Start-MD && column <= X_End+MD &&
			     vsize-row-1 >= Y_Start-MD && vsize-row-1 <= Y_End+MD) goto MARKPOINT;
			  else if(column >= X_End-MD && column <= X_Start+MD &&
			     vsize-row-1 >= Y_End-MD && vsize-row-1 <= Y_Start+MD) goto MARKPOINT;
			  else if(column >= X_Start-MD && column <= X_End+MD &&
			     vsize-row-1 >= Y_End-MD && vsize-row-1 <= Y_Start+MD) goto MARKPOINT;
			  else if(column >= X_End-MD && column <= X_Start+MD &&
			     vsize-row-1 >= Y_Start-MD && vsize-row-1 <= Y_End+MD) goto MARKPOINT;
			  continue;
			MARKPOINT:
			  if(D_to_Line(column, vsize-row-1, Slope, Y_Origin) < MD)
			    {
			      Xmain[image][itr][i][j] = column;
			      Ymain[image][itr][i][j] = vsize-row-1;
			      rsost[column][row] = 1;
			      j++;
			      jj++;
			    }
			}
		    }
		}
	    }
	  Nmain[image][itr][i] = j;
	}
      printf(" Number of real points along %d track is %d\n", itr, jj);
    }
}
