#include "_ROOT.h"
#include "_STD_UNIX.h"
#include "vertex.h"

extern int hsize, vsize;
extern int event_kind, tracks, ipoint[2], itrack[2];
extern float xcent_vx[10], ycent_vx[10], rcent_vx[10];
//extern float Xvertax, Yvertax, Dvertax;
extern TObjectTable *fObjTbl;

const float RTarget = 300.;
unsigned short int Cross;
float X1, Y1, X2, Y2;
float XvertInTarg, YvertInTarg, DvertInTarg;
float a1[3][3], b1[3][3], r1[3][3], a2[3][3], b2[3][3], r2[3][3];


void XY_Vertex_Net(int NTracks, float *a, float *b, float *r, 
		   float *Xvertex, float *Yvertex, 
		   float *Dvertex)
{
  int Step = 0;
  double XvOld, YvOld, DvOld;
  double Xd, Yd;
  double buf;
  int i;
  double xvertex, yvertex, dvertex;
  double xv, yv;
  xvertex = (double) 0.;//*Xvertex;//0.;//(double)hsize/2.;
  yvertex = (double) 300.;*Yvertex;//100.;//(double)vsize/2.;
  dvertex = 0.0;
  for(i = 0; i < NTracks ; i++)
    {
		xv = (double) *(a+i);
		yv = (double) *(b+i);
		printf("xv %f yv %f r %f\n", xv, yv, *(r+i));
      buf = fabs(sqrt((xv-xvertex)*(xv-xvertex)+
		      (yv-yvertex)*(yv-yvertex))-(*(r+i)));
      if (buf > dvertex)
	dvertex = buf;
    }
  printf(" Before vertex search is in position: %f %f %f\n", xvertex, yvertex, dvertex);
  //  xvertex = 0.;
  //  yvertex = 0.;
 label1:
  XvOld = xvertex;
  YvOld = yvertex;
  DvOld = dvertex;
  
  //cout << "step = "<<Step<<" "<<xvertex<<" "<<yvertex<<" "<<dvertex<<endl;
  Xd = Yd = 0.0;
  for(i=0; i < NTracks; i++)
    {
      buf = (sqrt((*(a+i)-xvertex)*(*(a+i)-xvertex)+
	      (*(b+i)-yvertex)*(*(b+i)-yvertex))-*(r+i))/
	(sqrt((*(a+i)-xvertex)*(*(a+i)-xvertex)+
	      (*(b+i)-yvertex)*(*(b+i)-yvertex)));
      Xd += (*(a+i)-xvertex)*buf;
      Yd += (*(b+i)-yvertex)*buf;
    }

  //  float sum_vtx = (float)(xvertex*xvertex+yvertex*yvertex);

  // cout <<"sumxv = "<<sum_vtx<< "rta = "<<RTarget*RTarget<<endl;

  if ((float)((xvertex-*Xvertex)*(xvertex-*Xvertex)+(yvertex-*Yvertex)*(yvertex-*Yvertex)) < (RTarget*RTarget))
    {
      xvertex += Xd/(double)NTracks;
      yvertex += Yd/(double)NTracks;
    }

  dvertex = 0.0;
  for(i = 0; i < NTracks ; i++)
    {
      buf = fabs(sqrt((*(a+i)-xvertex)*(*(a+i)-xvertex)+
		      (*(b+i)-yvertex)*(*(b+i)-yvertex))-(*(r+i)));
      if (buf > dvertex)
	dvertex = buf;
    }
  TArc *myarc = new TArc(xvertex,yvertex,dvertex,0.,360.);
  fObjTbl->Add(myarc);
  myarc->Draw("same");

  //  Tcl_Sleep(ms);
  Step += 1;
    
  if ((Step > 600) ||
      ((((xvertex-XvOld)*(xvertex-XvOld)+
	 (yvertex-YvOld)*(yvertex-YvOld)) < 0.1) &&
		 //		 dvertex<10. || dvertex-DvOld > 0.))
		        (fabs(dvertex-DvOld)< 0.01))) 
    {
      cout << "step = "<<Step<<" "<<xvertex<<" "<<yvertex<<" "<<dvertex<<endl;
      *Xvertex = (float)xvertex;
      *Yvertex = (float)yvertex;
      *Dvertex = (float)dvertex;
      return;
    }
  else
    goto label1;
}
