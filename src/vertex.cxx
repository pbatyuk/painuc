#include "_ROOT.h"
#include "_STD_UNIX.h"
#include "vertex.h"

extern int hsize, vsize;
//extern int event_kind, tracks, ipoint[2], itrack[2];
//extern float xcent_vx[10], ycent_vx[10], rcent_vx[10];
//extern float Xvertax, Yvertax, Dvertax;
//extern TObjectTable *fObjTbl;

const float RTarget = 600.;
//unsigned short int Cross;
//float X1, Y1, X2, Y2;
//float XvertInTarg, YvertInTarg, DvertInTarg;
//float a1[3][3], b1[3][3], r1[3][3], a2[3][3], b2[3][3], r2[3][3];


void XY_Vertex_Net(int NTracks, float *a, float *b, float *r, 
		   float *Xvertex, float *Yvertex, 
		   float *Dvertex)
{
  double xvertex, yvertex, dvertex;
  double ad[20], bd[20], rd[20];
  xvertex = (double)hsize/2.;
  yvertex = (double)vsize/2.;
  //  xvertex = 0.;
  //  yvertex = 0.;
  dvertex = 0.;
  int Step = 0;
  double XvOld, YvOld, DvOld;
  double Xd, Yd;
  double buf;
  int i;
 label1:
  XvOld = xvertex;
  YvOld = yvertex;
  DvOld = dvertex;
  
  Xd = Yd = 0.0;
  if(NTracks < 2) return;
  for(i=0; i < NTracks; i++)
    {
      ad[i] = (double) a[i];
      bd[i] = (double) b[i];
      rd[i] = (double) r[i];
    }
  for(i=0; i < NTracks; i++)
    {
      buf = (sqrt((ad[i]-xvertex)*(ad[i]-xvertex)+
	      (bd[i]-yvertex)*(bd[i]-yvertex))-rd[i])/
	(sqrt((ad[i]-xvertex)*(ad[i]-xvertex)+
	      (bd[i]-yvertex)*(bd[i]-yvertex)));
      Xd += (ad[i]-xvertex)*buf;
      Yd += (bd[i]-yvertex)*buf;
    }

  //  float sum_vtx = (float)(xvertex*xvertex+yvertex*yvertex);

  // cout <<"sumxv = "<<sum_vtx<< "rta = "<<RTarget*RTarget<<endl;

  if ((float)(xvertex*xvertex+yvertex*yvertex) < (RTarget*RTarget))
    {
      xvertex += Xd/(double)NTracks;
      yvertex += Yd/(double)NTracks;
    }

  dvertex = 0.0;
  for(i = 0; i < NTracks ; i++)
    {
      buf = fabs(sqrt((ad[i]-xvertex)*(ad[i]-xvertex)+
		      (bd[i]-yvertex)*(bd[i]-yvertex))-(rd[i]));
      if (buf > dvertex)
	dvertex = buf;
    }
  //TArc *myarc = new TArc(xvertex,yvertex,dvertex,0.,360.);
  //fObjTbl->Add(myarc);
  //myarc->Draw("same");

  //  Tcl_Sleep(ms);
  Step += 1;
    
  if ((Step > 300) ||
      ((((xvertex-XvOld)*(xvertex-XvOld)+
	 (yvertex-YvOld)*(yvertex-YvOld)) < 0.05) &&
       (fabs(dvertex-DvOld)< 0.05))) 
    {
      //cout << "step = "<<Step<<" "<<xvertex<<" "<<yvertex<<" "<<dvertex<<endl;
      *Xvertex = (float)xvertex;
      *Yvertex = (float)yvertex;
      *Dvertex = (float)dvertex;
      return;
    }
  else
    goto label1;
}
