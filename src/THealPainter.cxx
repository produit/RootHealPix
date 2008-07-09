// $Id: THealPainter.cxx,v 1.2 2008/07/09 11:50:10 oxon Exp $
// Author: Akira Okumura 2008/07/07

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include <float.h>

#include "Hoption.h"
#include "Hparam.h"
#include "TClass.h"
#include "TColor.h"
#include "TCrown.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TGraphDelaunay.h"
#include "TF1.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TPainter3dAlgorithms.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TVirtualPadEditor.h"

#include "THealPix.h"
#include "THealPainter.h"
#include "THealPaletteAxis.h"

static THealPix* gCurrentHeal = 0;

static Hoption_t Healoption;
static Hparam_t  Healparam;

static const Int_t kNMAX = 2000;

ClassImp(THealPainter)

//______________________________________________________________________________
THealPainter::THealPainter()
{
  fHeal  = 0;
  fXaxis = 0;
  fYaxis = 0;
  fZaxis = 0;
  fFunctions = 0;
}

//______________________________________________________________________________
THealPainter::~THealPainter()
{
}

//______________________________________________________________________________
Int_t THealPainter::DistancetoPrimitive(Int_t px, Int_t py)
{
  return 0;
}

//______________________________________________________________________________
void THealPainter::DrawPanel()
{
  gCurrentHeal = fHeal;
  if(!gPad){
    Error("DrawPanel", "need to draw HEALPix first");
    return;
  } // if
  TVirtualPadEditor* editor = TVirtualPadEditor::GetPadEditor();
  editor->Show();
  gROOT->ProcessLine(Form("((TCanvas*)0x%lx)->Selected((TVirtualPad*)0x%lx, (TObject*)0x%lx, 1)", gPad->GetCanvas(), gPad, fHeal));
}

//______________________________________________________________________________
void THealPainter::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
}

//______________________________________________________________________________
TList* THealPainter::GetContourList(Double_t contour) const
{
  TGraphDelaunay* dt;
  
  // Check if fHeal contains a TGraphDelaunay
  TList* hl = fHeal->GetListOfFunctions();
  dt = (TGraphDelaunay*)hl->FindObject("TGraphDelaunay");
  if(!dt) return 0;
  /*
  gCurrentHeal = fHeal;
  
  if(!fGraph2DPainter) ((THistPainter*)this)->fGraph2DPainter = new TGraph2DPainter(dt);
  
  return fGraph2DPainter->GetContourList(contour);
  */
  return 0;
}

//______________________________________________________________________________
Bool_t THealPainter::IsInside(Int_t ix, Int_t iy)
{
  return kTRUE;
}

//______________________________________________________________________________
Bool_t THealPainter::IsInside(Double_t x, Double_t y)
{
  return kTRUE;
}

//______________________________________________________________________________
Int_t THealPainter::MakeChopt(Option_t* choptin)
{
  char *l;
  char chopt[128];
  Int_t nch = strlen(choptin);
  strcpy(chopt,choptin);
  
  Healoption.Axis = Healoption.Bar   = Healoption.Curve   = Healoption.Error = 0;
  Healoption.Hist = Healoption.Line  = Healoption.Mark    = Healoption.Fill  = 0;
  Healoption.Same = Healoption.Func  = Healoption.Plus    = Healoption.Scat  = 0;
  Healoption.Star = Healoption.Arrow = Healoption.Box     = Healoption.Text  = 0;
  Healoption.Char = Healoption.Color = Healoption.Contour = Healoption.Logx  = 0;
  Healoption.Logy = Healoption.Logz  = Healoption.Lego    = Healoption.Surf  = 0;
  Healoption.Off  = Healoption.Tri   = Healoption.Proj    = Healoption.AxisPos = 0;
  Healoption.Spec = Healoption.Pie   = 0;
  
  //    special 2D options
  Healoption.List     = 0;
  Healoption.Zscale   = 0;
  Healoption.FrontBox = 1;
  Healoption.BackBox  = 1;
  Healoption.System   = kCARTESIAN;
  
  Healoption.HighRes  = 0;
  
  Healoption.Zero     = 0;
  
  //check for graphical cuts
  MakeCuts(chopt);
  
  for (Int_t i = 0;i < nch; i++) chopt[i] = toupper(chopt[i]);
  Healoption.Scat = 1;
  if(!nch) Healoption.Hist = 1;
  if(fFunctions->First()) Healoption.Func = 2;
  
  l = strstr(chopt,"SPEC");
  if (l) {
    Healoption.Scat = 0;
    Healoption.Spec = 1; strncpy(l,"    ",4);
    return 1;
  }
  
  l = strstr(chopt,"GL");
  if (l) {
    strncpy(l,"  ",2);
  }
  l = strstr(chopt,"X+");
  if (l) {
    Healoption.AxisPos = 10;
    strncpy(l,"  ",2);
  }
  l = strstr(chopt,"Y+");
  if (l) {
    Healoption.AxisPos += 1;
    strncpy(l,"  ",2);
  }
  if((Healoption.AxisPos == 10 || Healoption.AxisPos == 1) && (nch == 2)) Healoption.Hist = 1;
  if(Healoption.AxisPos == 11 && nch == 4) Healoption.Hist = 1;
  
  l = strstr(chopt,"SAMES");
  if (l) {
    if (nch == 5) Healoption.Hist = 1;
    Healoption.Same = 2;
    strncpy(l,"     ",5);
  }
  l = strstr(chopt,"SAME");
  if (l) {
    if (nch == 4) Healoption.Hist = 1;
    Healoption.Same = 1;
    strncpy(l,"    ",4);
  }
  
  l = strstr(chopt,"PIE");
  if (l) {
    Healoption.Pie = 1;
    strncpy(l,"   ",3);
   }
  
  l = strstr(chopt,"LEGO");
  if (l) {
    Healoption.Scat = 0;
    Healoption.Lego = 1; strncpy(l,"    ",4);
    if (l[4] == '1') { Healoption.Lego = 11; l[4] = ' '; }
    if (l[4] == '2') { Healoption.Lego = 12; l[4] = ' '; }
    l = strstr(chopt,"FB"); if (l) { Healoption.FrontBox = 0; strncpy(l,"  ",2); }
    l = strstr(chopt,"BB"); if (l) { Healoption.BackBox = 0;  strncpy(l,"  ",2); }
    l = strstr(chopt,"0");  if (l) { Healoption.Zero = 1;  strncpy(l," ",1); }
  }
  
  l = strstr(chopt,"SURF");
  if (l) {
    Healoption.Scat = 0;
    Healoption.Surf = 1; strncpy(l,"    ",4);
    if (l[4] == '1') { Healoption.Surf = 11; l[4] = ' '; }
    if (l[4] == '2') { Healoption.Surf = 12; l[4] = ' '; }
    if (l[4] == '3') { Healoption.Surf = 13; l[4] = ' '; }
    if (l[4] == '4') { Healoption.Surf = 14; l[4] = ' '; }
    if (l[4] == '5') { Healoption.Surf = 15; l[4] = ' '; }
    if (l[4] == '6') { Healoption.Surf = 16; l[4] = ' '; }
    l = strstr(chopt,"FB");   if (l) { Healoption.FrontBox = 0; strncpy(l,"  ",2); }
    l = strstr(chopt,"BB");   if (l) { Healoption.BackBox = 0;  strncpy(l,"  ",2); }
  }
  
  l = strstr(chopt,"TF3");
  if (l) {
    l = strstr(chopt,"FB");   if (l) { Healoption.FrontBox = 0; strncpy(l,"  ",2); }
    l = strstr(chopt,"BB");   if (l) { Healoption.BackBox = 0;  strncpy(l,"  ",2); }
  }
  
  l = strstr(chopt,"ISO");
  if (l) {
    l = strstr(chopt,"FB");   if (l) { Healoption.FrontBox = 0; strncpy(l,"  ",2); }
    l = strstr(chopt,"BB");   if (l) { Healoption.BackBox = 0;  strncpy(l,"  ",2); }
  }
  
  l = strstr(chopt,"LIST");    if (l) { Healoption.List = 1;  strncpy(l,"    ",4);}
  
  l = strstr(chopt,"CONT");
  if (l) {
    Healoption.Scat = 0;
    Healoption.Contour = 1; strncpy(l,"    ",4);
    if (l[4] == '1') { Healoption.Contour = 11; l[4] = ' '; }
    if (l[4] == '2') { Healoption.Contour = 12; l[4] = ' '; }
    if (l[4] == '3') { Healoption.Contour = 13; l[4] = ' '; }
    if (l[4] == '4') { Healoption.Contour = 14; l[4] = ' '; }
    if (l[4] == '5') { Healoption.Contour = 15; l[4] = ' '; }
  }
  l = strstr(chopt,"HBAR");
  if (l) {
    Healoption.Hist = 0;
    Healoption.Bar = 20; strncpy(l,"    ",4);
    if (l[4] == '1') { Healoption.Bar = 21; l[4] = ' '; }
    if (l[4] == '2') { Healoption.Bar = 22; l[4] = ' '; }
    if (l[4] == '3') { Healoption.Bar = 23; l[4] = ' '; }
    if (l[4] == '4') { Healoption.Bar = 24; l[4] = ' '; }
  }
  l = strstr(chopt,"BAR");
  if (l) {
    Healoption.Hist = 0;
    Healoption.Bar = 10; strncpy(l,"   ",3);
    if (l[3] == '1') { Healoption.Bar = 11; l[3] = ' '; }
    if (l[3] == '2') { Healoption.Bar = 12; l[3] = ' '; }
    if (l[3] == '3') { Healoption.Bar = 13; l[3] = ' '; }
    if (l[3] == '4') { Healoption.Bar = 14; l[3] = ' '; }
  }
  l = strstr(chopt,"+-");   if (l) { Healoption.Plus = 2; strncpy(l,"  ",2); }
  l = strstr(chopt,"-+");   if (l) { Healoption.Plus = 2; strncpy(l,"  ",2); }
  
  l = strstr(chopt,"ARR" ); if (l) { Healoption.Arrow  = 1; strncpy(l,"   ", 3); Healoption.Scat = 0; }
  l = strstr(chopt,"BOX" );
  if (l) {
    Healoption.Scat = 0;
    Healoption.Box  = 1; strncpy(l,"   ", 3);
    if (l[3] == '1') { Healoption.Box = 11; l[3] = ' '; }
  }
  l = strstr(chopt,"COLZ"); if (l) { Healoption.Color  = 2; strncpy(l,"    ",4); Healoption.Scat = 0; Healoption.Zscale = 1;}
  l = strstr(chopt,"COL" ); if (l) { Healoption.Color  = 1; strncpy(l,"   ", 3); Healoption.Scat = 0; }
  l = strstr(chopt,"CHAR"); if (l) { Healoption.Char   = 1; strncpy(l,"    ",4); Healoption.Scat = 0; }
  l = strstr(chopt,"FUNC"); if (l) { Healoption.Func   = 2; strncpy(l,"    ",4); Healoption.Hist = 0; }
  l = strstr(chopt,"HIST"); if (l) { Healoption.Hist   = 2; strncpy(l,"    ",4); Healoption.Func = 0; Healoption.Error = 0;}
  l = strstr(chopt,"AXIS"); if (l) { Healoption.Axis   = 1; strncpy(l,"    ",4); }
  l = strstr(chopt,"AXIG"); if (l) { Healoption.Axis   = 2; strncpy(l,"    ",4); }
  l = strstr(chopt,"SCAT"); if (l) { Healoption.Scat   = 1; strncpy(l,"    ",4); }
  l = strstr(chopt,"TEXT");
  if (l) {
    Int_t angle;
    if (sscanf(&l[4],"%d",&angle) > 0) {
      if (angle < 0)  angle=0;
      if (angle > 90) angle=90;
      Healoption.Text = 1000+angle;
    } else {
      Healoption.Text = 1;
    }
    strncpy(l,"    ",4);
    Healoption.Scat = 0;
  }
  l = strstr(chopt,"POL");  if (l) { Healoption.System = kPOLAR;       strncpy(l,"   ",3); }
  l = strstr(chopt,"CYL");  if (l) { Healoption.System = kCYLINDRICAL; strncpy(l,"   ",3); }
  l = strstr(chopt,"SPH");  if (l) { Healoption.System = kSPHERICAL;   strncpy(l,"   ",3); }
  l = strstr(chopt,"PSR");  if (l) { Healoption.System = kRAPIDITY;    strncpy(l,"   ",3); }
  
  l = strstr(chopt,"TRI");
  if (l) {
    Healoption.Scat = 0;
    Healoption.Color  = 0;
    Healoption.Tri = 1; strncpy(l,"   ",3);
    l = strstr(chopt,"FB");   if (l) { Healoption.FrontBox = 0; strncpy(l,"  ",2); }
    l = strstr(chopt,"BB");   if (l) { Healoption.BackBox = 0;  strncpy(l,"  ",2); }
  }
  
  l = strstr(chopt,"AITOFF");
  if (l) {
    Healoption.Proj = 1; strncpy(l,"     ",6);       //Aitoff projection
  }
  l = strstr(chopt,"MERCATOR");
  if (l) {
    Healoption.Proj = 2; strncpy(l,"       ",8);     //Mercator projection
  }
  l = strstr(chopt,"SINUSOIDAL");
  if (l) {
    Healoption.Proj = 3; strncpy(l,"         ",10);  //Sinusoidal projection
  }
  l = strstr(chopt,"PARABOLIC");
  if (l) {
    Healoption.Proj = 4; strncpy(l,"        ",9);    //Parabolic projection
  }
  if (Healoption.Proj > 0) {
    Healoption.Scat = 0;
    Healoption.Contour = 14;
  }
  
  if (strstr(chopt,"A"))   Healoption.Axis = -1;
  if (strstr(chopt,"B"))   Healoption.Bar  = 1;
  if (strstr(chopt,"C")) { Healoption.Curve =1; Healoption.Hist = -1;}
  if (strstr(chopt,"F"))   Healoption.Fill =1;
  if (strstr(chopt,"][")) {Healoption.Off  =1; Healoption.Hist =1;}
  if (strstr(chopt,"F2"))  Healoption.Fill =2;
  if (strstr(chopt,"L")) { Healoption.Line =1; Healoption.Hist = -1;}
  if (strstr(chopt,"P")) { Healoption.Mark =1; Healoption.Hist = -1;}
  if (strstr(chopt,"Z"))   Healoption.Zscale =1;
  if (strstr(chopt,"*"))   Healoption.Star =1;
  if (strstr(chopt,"+"))   Healoption.Plus =1;
  if (strstr(chopt,"-"))   Healoption.Plus =-1;
  if (strstr(chopt,"H"))   Healoption.Hist =2;
  if (strstr(chopt,"P0"))  Healoption.Mark =10;
  if (strstr(chopt,"E")) {
    if (Healoption.Error == 0) {
      Healoption.Error = 100;
      Healoption.Scat  = 0;
    }
    if (Healoption.Text) {
      Healoption.Text += 2000;
      Healoption.Error = 0;
    }
  }
  
  if (strstr(chopt,"9"))  Healoption.HighRes = 1;
  
  if (Healoption.Surf == 15) {
    if (Healoption.System == kPOLAR || Healoption.System == kCARTESIAN) {
      Healoption.Surf = 13;
      Warning("MakeChopt","option SURF5 is not supported in Cartesian and Polar modes");
    }
  }
  
  //      Copy options from current style
  Healoption.Logx = gPad->GetLogx();
  Healoption.Logy = gPad->GetLogy();
  Healoption.Logz = gPad->GetLogz();
  
  //       Check options incompatibilities
  if (Healoption.Bar  == 1) Healoption.Hist = -1;
  if (Healoption.Same && Healoption.Plus) {
    Error("MakeChopt", "select only one of the options S,+");
    return 0;
  }
  if (Healoption.Plus) {
    if (Healoption.Line || Healoption.Curve || Healoption.Text || Healoption.Mark) {
      Error("MakeChopt", "options L,C,T,P are incompatible with options U and K");
      if (Healoption.Hist && Healoption.Bar) return 0;
    }
  }
  if (Healoption.Error || Healoption.Func || Healoption.Star) {
    if (Healoption.Plus) {
      Error("MakeChopt", "U, + options incompatible with errors/function");
      return 0;
    }
  }
  return 1;
}


//______________________________________________________________________________
Int_t THealPainter::MakeCuts(char* choptin)
{
  fNcuts = 0;
  char *left = (char*)strchr(choptin,'[');
  if (!left) return 0;
  char *right = (char*)strchr(choptin,']');
  if (!right) return 0;
  Int_t nch = right-left;
  if (nch < 2) return 0;
  char *cuts = left+1;
  *right = 0;
  char *comma, *minus;
  Int_t i;
  while(1) {
    comma = strchr(cuts,',');
    if (comma) *comma = 0;
    minus = strchr(cuts,'-');
    if (minus) cuts = minus+1;
    while (*cuts == ' ') cuts++;
    Int_t nc = strlen(cuts);
    while (cuts[nc-1] == ' ') {cuts[nc-1] = 0; nc--;}
    TIter next(gROOT->GetListOfSpecials());
    TCutG *cut=0;
    TObject *obj;
    while ((obj = next())) {
      if (!obj->InheritsFrom(TCutG::Class())) continue;
      if (strcmp(obj->GetName(),cuts)) continue;
      cut = (TCutG*)obj;
      break;
    }
    if (cut) {
      fCuts[fNcuts] = cut;
      fCutsOpt[fNcuts] = 1;
      if (minus) fCutsOpt[fNcuts] = -1;
      fNcuts++;
    }
    if (!comma) break;
    cuts = comma+1;
  }
  for (i=0;i<=nch;i++) left[i] = ' ';
  return fNcuts;
}

//______________________________________________________________________________
void THealPainter::Paint(Option_t* option)
{
  printf("*****Paint*****\n");
  THealPix* oldheal = gCurrentHeal;
  gCurrentHeal = fHeal;
  THealPix* healsave = fHeal;

  Double_t minsav = fHeal->GetMinimumStored();

  if(!MakeChopt(option)){
    return;
  } // if

  fXbuf = new Double_t[kNMAX];
  fYbuf = new Double_t[kNMAX];

  TView* view = gPad->GetView();
  if(view){
    if(!Healoption.Lego && !Healoption.Surf && !Healoption.Tri){
      delete view;
      gPad->SetView(0);
    } // if
  } // if

  PaintTable(option);
  fHeal->SetMinimum(minsav);
  gCurrentHeal = oldheal;

  delete [] fXbuf;
  delete [] fYbuf;
}

//______________________________________________________________________________
void THealPainter::PaintAxis(Bool_t drawGridOnly)
{
  printf("*****PaintAxis*****\n");
  if(Healoption.Axis == -1){
    return;
  } // if
  if(Healoption.Same && Healoption.Axis <= 0){
    return;
  } // if
  
  static char chopt[10] = "";
  Double_t gridl = 0;
  Int_t ndiv, ndivsave;
  Int_t useHealparam = 0;
  Double_t umin, umax, uminsave, umaxsave;
  Short_t xAxisPos = Healoption.AxisPos/10;
  Short_t yAxisPos = Healoption.AxisPos - 10*xAxisPos;
  
  Double_t axmin = gPad->GetUxmin();
  Double_t axmax = gPad->GetUxmax();
  Double_t aymin = gPad->GetUymin();
  Double_t aymax = gPad->GetUymax();
  char* cw = 0;
  TGaxis axis;
  
  // In case of option 'cont4' or in case of option 'same' over a 'cont4 plot'
  // Healparam must be use for the axis limits.
  if(Healoption.Contour == 14){
    useHealparam = 1;
  } // if
  if(Healoption.Same){
    TObject *obj;
    TIter next(gPad->GetListOfPrimitives());
    while((obj = next())){
      if(strstr(obj->GetDrawOption(), "cont4")){
	useHealparam = 1;
 	break;
      } // if
    } // while
  } // if
  
  // Paint X axis
  Int_t ndivx = fXaxis->GetNdivisions();
  if (ndivx > 1000) {
    Int_t nx2   = ndivx/100;
    Int_t nx1   = TMath::Max(1, ndivx%100);
    ndivx = 100*nx2 + Int_t(Float_t(nx1)*gPad->GetAbsWNDC());
  }
  axis.SetTextAngle(0);
  axis.ImportAxisAttributes(fXaxis);
  
  chopt[0] = 0;
  strcat(chopt, "SDH");
  if(ndivx < 0){
    strcat(chopt, "N");
  } // if
  if(gPad->GetGridx()){
    gridl = (aymax - aymin)/(gPad->GetY2() - gPad->GetY1());
    strcat(chopt, "W");
  } // if
  
  // Define X-Axis limits
  if(Healoption.Logx){
    strcat(chopt, "G");
    ndiv = TMath::Abs(ndivx);
    if(useHealparam){
      umin = TMath::Power(10, Healparam.xmin);
      umax = TMath::Power(10, Healparam.xmax);
    } else {
      umin = TMath::Power(10, axmin);
      umax = TMath::Power(10, axmax);
    } // if
  } else {
    ndiv = TMath::Abs(ndivx);
    if(useHealparam){
      umin = Healparam.xmin;
      umax = Healparam.xmax;
    } else {
      umin = axmin;
      umax = axmax;
    } // if
  } // if
  
  // The main X axis can be on the bottom or on the top of the pad
  Double_t xAxisYPos1, xAxisYPos2;
  if(xAxisPos == 1){
    // Main X axis top
    xAxisYPos1 = aymax;
    xAxisYPos2 = aymin;
  } else {
    // Main X axis bottom
    xAxisYPos1 = aymin;
    xAxisYPos2 = aymax;
  } // if
  
  // Paint the main X axis (always)
  uminsave = umin;
  umaxsave = umax;
  ndivsave = ndiv;
  axis.SetOption(chopt);
  if(xAxisPos){
    strcat(chopt, "-");
    gridl = -gridl;
  } // if
  axis.PaintAxis(axmin, xAxisYPos1, axmax, xAxisYPos1,
		 umin, umax,  ndiv, chopt, gridl, drawGridOnly);
  
  // Paint additional X axis (if needed)
  if(gPad->GetTickx()){
    if(xAxisPos){
      cw = strstr(chopt,"-");
      *cw = 'z';
    } else {
      strcat(chopt, "-");
    } // if
    if(gPad->GetTickx() < 2){
      strcat(chopt, "U");
    } // if
    if((cw=strstr(chopt,"W"))){
      *cw='z';
    } // if
    axis.SetTitle("");
    axis.PaintAxis(axmin, xAxisYPos2, axmax, xAxisYPos2,
		   uminsave, umaxsave, ndivsave, chopt, gridl, drawGridOnly);
  }
  
  // Paint Y axis
  Int_t ndivy = fYaxis->GetNdivisions();
  axis.ImportAxisAttributes(fYaxis);
  
  chopt[0] = 0;
  strcat(chopt, "SDH");
  if(ndivy < 0){
    strcat(chopt, "N");
  } // if
  if(gPad->GetGridy()){
    gridl = (axmax - axmin)/(gPad->GetX2() - gPad->GetX1());
    strcat(chopt, "W");
  } // if
  
  // Define Y-Axis limits
  if(Healoption.Logy){
    strcat(chopt, "G");
    ndiv = TMath::Abs(ndivy);
    if(useHealparam){
      umin = TMath::Power(10, Healparam.ymin);
      umax = TMath::Power(10, Healparam.ymax);
    } else {
      umin = TMath::Power(10, aymin);
      umax = TMath::Power(10, aymax);
    } // if
  } else {
    ndiv = TMath::Abs(ndivy);
    if(useHealparam){
      umin = Healparam.ymin;
      umax = Healparam.ymax;
    } else {
      umin = aymin;
      umax = aymax;
    } // if
  } // if
  
  // Display axis as time
  if(fYaxis->GetTimeDisplay()){
    strcat(chopt, "t");
    if(strlen(fYaxis->GetTimeFormatOnly()) == 0){
      axis.SetTimeFormat(fYaxis->ChooseTimeFormat(Healparam.ymax - Healparam.ymin));
    } // if
  } // if
  
  // The main Y axis can be on the left or on the right of the pad
  Double_t yAxisXPos1, yAxisXPos2;
  if(yAxisPos == 1){
    // Main Y axis left
    yAxisXPos1 = axmax;
    yAxisXPos2 = axmin;
  } else {
    // Main Y axis right
    yAxisXPos1 = axmin;
    yAxisXPos2 = axmax;
  } // if
  
  // Paint the main Y axis (always)
  uminsave = umin;
  umaxsave = umax;
  ndivsave = ndiv;
  axis.SetOption(chopt);
  if(yAxisPos){
    strcat(chopt, "+L");
    gridl = -gridl;
  } // if
  axis.PaintAxis(yAxisXPos1, aymin, yAxisXPos1, aymax,
		 umin, umax,  ndiv, chopt, gridl, drawGridOnly);
  
  // Paint the additional Y axis (if needed)
  if(gPad->GetTicky()){
    if(gPad->GetTicky() < 2){
      strcat(chopt, "U");
      axis.SetTickSize(-fYaxis->GetTickLength());
    } else {
      strcat(chopt, "+L");
    } // if
    if((cw=strstr(chopt,"W"))){
      *cw='z';
    } // if
    axis.SetTitle("");
    axis.PaintAxis(yAxisXPos2, aymin, yAxisXPos2, aymax,
		   uminsave, umaxsave,  ndivsave, chopt, gridl, drawGridOnly);
  } // if
}

//______________________________________________________________________________
void THealPainter::PaintBarH(Option_t* option)
{
  printf("*****PaintBarH*****\n");
  gPad->SetVertical(kFALSE);
  
  PaintInitH();
  
  TAxis* xaxis = fXaxis;
  TAxis* yaxis = fYaxis;
  if (!strcmp(xaxis->GetName(), "xaxis")) {
    fXaxis = yaxis;
    fYaxis = xaxis;
  }

  PaintFrame();
  
  Int_t bar = Healoption.Bar - 20;
  Double_t xmin, xmax, ymin, ymax, umin, umax, w;
  Double_t offset = fHeal->GetBarOffset();
  Double_t width  = fHeal->GetBarWidth();
  TBox box;
  Int_t hcolor = fHeal->GetFillColor();
  Int_t hstyle = fHeal->GetFillStyle();
  box.SetFillColor(hcolor);
  box.SetFillStyle(hstyle);
  for(Int_t bin = fYaxis->GetFirst(); bin <= fYaxis->GetLast(); bin++) {
    ymin = gPad->YtoPad(fYaxis->GetBinLowEdge(bin));
    ymax = gPad->YtoPad(fYaxis->GetBinUpEdge(bin));
    xmin = gPad->GetUxmin();
    xmax = gPad->XtoPad(fHeal->GetBinContent(bin));
    if(xmax < gPad->GetUxmin()) continue;
    if(xmax > gPad->GetUxmax()) xmax = gPad->GetUxmax();
    if(xmin < gPad->GetUxmin()) xmin = gPad->GetUxmin();
    w = (ymax - ymin)*width;
    ymin += offset*(ymax - ymin);
    ymax = ymin + w;
    if(bar < 1){
      box.PaintBox(xmin, ymin, xmax, ymax);
    } else {
      umin = ymin + bar*(ymax - ymin)/10.;
      umax = ymax - bar*(ymax - ymin)/10.;
      box.SetFillColor(TColor::GetColorDark(hcolor)); //dark
      box.PaintBox(xmin, ymin, xmax, umin);
      box.SetFillColor(hcolor);
      box.PaintBox(xmin, umin, xmax, umax);
      box.SetFillColor(TColor::GetColorBright(hcolor)); //bright
      box.PaintBox(xmin, umax, xmax, ymax);
    } // if
  } // bin

  PaintTitle();
 //    Draw box with histogram statistics and/or fit parameters
  if(Healoption.Same != 1 && !fHeal->TestBit(THealPix::kNoStats)){  // bit set via TH1::SetStats
    TIter next(fFunctions);
    TObject* obj = 0;
    while ((obj = next())) {
      if (obj->InheritsFrom(TF1::Class())) break;
      obj = 0;
    } // while
    PaintStat(gStyle->GetOptStat(), (TF1*)obj);
  } // if

  PaintAxis(kFALSE);
  fXaxis = xaxis;
  fYaxis = yaxis;  
}

//______________________________________________________________________________
void THealPainter::PaintColorLevels(Option_t* option)
{
  printf("*****PaintColorLevels*****\n");

  Double_t z, zc, xk, xstep, yk, ystep, xlow, xup, ylow, yup;
  
  Double_t zmin = fHeal->GetMinimum();
  Double_t zmax = fHeal->GetMaximum();
  if(Healoption.Logz){
    if(zmin > 0){
      zmin = TMath::Log10(zmin);
      zmax = TMath::Log10(zmax);
    } else {
      return;
    } // if
  } // if
  
  Double_t dz = zmax - zmin;
  
  if(dz <= 0){
    return;
  } // if
  
  Style_t fillsav = fHeal->GetFillStyle();
  Style_t colsav  = fHeal->GetFillColor();
  fHeal->SetFillStyle(1001);
  fHeal->TAttFill::Modify();
  
  // Initialize the levels on the Z axis
  Int_t ncolors = gStyle->GetNumberOfColors();
  Int_t ndiv = fHeal->GetContour();
  if(ndiv == 0){
    ndiv = gStyle->GetNumberContours();
    fHeal->SetContour(ndiv);
  } // if
  Int_t ndivz = TMath::Abs(ndiv);
  if(fHeal->TestBit(THealPix::kUserContour) == 0){
    fHeal->SetContour(ndiv);
  } // if
  Double_t scale = ndivz/dz;
  
  Int_t color;
  for(Int_t j = Healparam.yfirst; j <= Healparam.ylast; j++){
    yk    = fYaxis->GetBinLowEdge(j);
    ystep = fYaxis->GetBinWidth(j);
    if(Healoption.System == kPOLAR && yk < 0){
      yk= 2*TMath::Pi() + yk;
    } // if
    for(Int_t i = Healparam.xfirst; i <= Healparam.xlast; i++){
      Int_t bin = j*(fXaxis->GetNbins() + 2) + i;
      xk    = fXaxis->GetBinLowEdge(i);
      xstep = fXaxis->GetBinWidth(i);
      if(!IsInside(xk + 0.5*xstep, yk + 0.5*ystep)){
	continue;
      } // if
      z = fHeal->GetBinContent(bin);
      if(z == 0 && (zmin >= 0 || Healoption.Logz)){
	continue; // don't draw the empty bins for histograms with positive content
      } // if
      if(Healoption.Logz){
	z = z > 0 ? TMath::Log10(z) : zmin;
      } // if
      if(z < zmin){
	continue;
      } // if
      xup  = xk + xstep;
      xlow = xk;
      if(Healoption.Logx){
	if (xup > 0)  xup  = TMath::Log10(xup);
	else continue;
	if (xlow > 0) xlow = TMath::Log10(xlow);
	else continue;
      } // if
      yup  = yk + ystep;
      ylow = yk;
      if(Healoption.System != kPOLAR){
	if(Healoption.Logy){
	  if (yup > 0)  yup  = TMath::Log10(yup);
	  else continue;
	  if (ylow > 0) ylow = TMath::Log10(ylow);
	  else continue;
	} // if
	if(xup  < gPad->GetUxmin()) continue;
	if(yup  < gPad->GetUymin()) continue;
	if(xlow > gPad->GetUxmax()) continue;
	if(ylow > gPad->GetUymax()) continue;
	if(xlow < gPad->GetUxmin()) xlow = gPad->GetUxmin();
	if(ylow < gPad->GetUymin()) ylow = gPad->GetUymin();
	if(xup  > gPad->GetUxmax()) xup  = gPad->GetUxmax();
	if(yup  > gPad->GetUymax()) yup  = gPad->GetUymax();
      } // if
      
      if(fHeal->TestBit(THealPix::kUserContour)){
	zc = fHeal->GetContourLevelPad(0);
	if (z < zc) continue;
	color = -1;
	for(Int_t k = 0; k < ndiv; k++){
	  zc = fHeal->GetContourLevelPad(k);
	  if(z < zc){
	    continue;
	  } else {
	    color++;
	  } // if
	} // k
      } else {
	color = Int_t(0.01 + (z - zmin)*scale);
      } // if
      
      Int_t theColor = Int_t((color + 0.99)*Float_t(ncolors)/Float_t(ndivz));
      if(theColor > ncolors - 1){
	theColor = ncolors - 1;
      } // if
      fHeal->SetFillColor(gStyle->GetColorPalette(theColor));
      fHeal->TAttFill::Modify();
      if(Healoption.System != kPOLAR){
	gPad->PaintBox(xlow, ylow, xup, yup);
      } else {
	TCrown crown(0,0,xlow,xup,ylow*TMath::RadToDeg(),yup*TMath::RadToDeg());
	crown.SetFillColor(gStyle->GetColorPalette(theColor));
	crown.Paint();
      } // if
    } // i
  } // j
  
  if(Healoption.Zscale){
    PaintPalette();
  } // if
  
  fHeal->SetFillStyle(fillsav);
  fHeal->SetFillColor(colsav);
  fHeal->TAttFill::Modify();
}

//______________________________________________________________________________
void THealPainter::PaintContour(Option_t* option)
{
}

//______________________________________________________________________________
Int_t THealPainter::PaintInit()
{
  printf("*****PaintInit*****\n");
  return 1;
}

//______________________________________________________________________________
Int_t THealPainter::PaintInitH()
{
  printf("*****PaintInitH*****\n");
  /*
  static const char* where = "PaintInitH";
  Double_t yMARGIN = gStyle->GetHistTopMargin();
  Int_t maximum = 0;
  Int_t minimum = 0;
  if(fHeal->GetMaximumStored() != -1111) maximum = 1;
  if(fHeal->GetMinimumStored() != -1111) minimum = 1;
  
  //     Compute X axis parameters
  
  Int_t last      = fXaxis->GetLast();
  Int_t first     = fXaxis->GetFirst();
  Healparam.xlowedge = fXaxis->GetBinLowEdge(first);
  Healparam.xbinsize = fXaxis->GetBinWidth(first);
  Healparam.xlast    = last;
  Healparam.xfirst   = first;
  Healparam.ymin     = Healparam.xlowedge;
  Healparam.ymax     = fXaxis->GetBinLowEdge(last) + fXaxis->GetBinWidth(last);
  
  //       if log scale in Y, replace ymin,max by the log
  if(Healoption.Logy){
    if(Healparam.xlowedge <= 0 ){
      Healparam.xlowedge = 0.1*Healparam.xbinsize;
      Healparam.ymin  = Healparam.xlowedge;
    } // if
    if(Healparam.ymin <= 0 || Healparam.ymax <= 0){
      Error(where, "cannot set Y axis to log scale");
      return 0;
    } // if
    Healparam.xfirst= fXaxis->FindFixBin(Healparam.ymin);
    Healparam.xlast = fXaxis->FindFixBin(Healparam.ymax);
    Healparam.ymin  = TMath::Log10(Healparam.ymin);
    Healparam.ymax  = TMath::Log10(Healparam.ymax);
    if(Healparam.xlast > last) Healparam.xlast = last;
  } // if
  
  //     Compute Y axis parameters
  Double_t bigp = TMath::Power(10, 32);
  Double_t xmax = -bigp;
  Double_t xmin = bigp;
  Double_t allchan = 0;

  for(Int_t i = first; i <= last; i++){
    Dobule_t c1 = fHeal->GetBinContent(i);
    xmax = TMath::Max(xmax, c1);
    xmin = TMath::Min(xmin, c1);
    if(Healoption.Error){
      Double_t e1 = fHeal->GetBinError(i);
      xmax = TMath::Max(xmax, c1 + e1);
      xmin = TMath::Min(xmin, c1 - e1);
    } // if
    allchan += c1;
  }
  
  //     Take into account maximum , minimum
  
  if (Healoption.Logx && xmin <= 0) {
    if (xmax >= 1) xmin = TMath::Max(.5,xmax*1e-10);
    else           xmin = 0.001*xmax;
  }
  Double_t xm = xmin;
  if (maximum) xmax = fHeal->GetMaximumStored();
  if (minimum) xm   = fHeal->GetMinimumStored();
  if (Healoption.Logx && xm <= 0) {
    Error(where, "log scale requested with zero or negative argument (%f)", xm);
    return 0;
  }
  else xmin = xm;
  if (xmin >= xmax && !Healoption.Plus) {
    if (Healoption.Logx) {
      if (xmax > 0) xmin = 0.001*xmax;
      else {
	if (!Healoption.Same) Error(where, "log scale is requested but maximum is less or equal 0 (%f)", xmax);
	return 0;
      }
    }
    else {
      if (xmin > 0) {
	xmin = 0;
	xmax *= 2;
      } else if (xmin < 0) {
	xmax = 0;
	xmin *= 2;
      } else {
	xmin = -1;
	xmax = 1;
      }
    }
  }
  
  //     take into account normalization factor
  Healparam.allchan = allchan;
  Double_t factor = allchan;
  if (fHeal->GetNormFactor() > 0) factor = fHeal->GetNormFactor();
  if (allchan) factor /= allchan;
  if (factor == 0) factor = 1;
  Healparam.factor = factor;
  xmax = factor*xmax;
  xmin = factor*xmin;
  
  //         For log scales, histogram coordinates are LOG10(ymin) and
  //         LOG10(ymax). Final adjustment (if not option "Same"
  //         or "+" for ymax) of ymax and ymin for logarithmic scale, if
  //         Maximum and Minimum are not defined.
  if (Healoption.Logx) {
    if (xmin <=0 || xmax <=0) {
      Error(where, "Cannot set Y axis to log scale");
      return 0;
    }
    xmin = TMath::Log10(xmin);
    if (!minimum) xmin += TMath::Log10(0.5);
    xmax = TMath::Log10(xmax);
    if (!maximum && !Healoption.Plus) xmax += TMath::Log10(2*(0.9/0.95));
    if (!Healoption.Same) {
      Healparam.xmin = xmin;
      Healparam.xmax = xmax;
    }
    return 1;
  }
  
  //         final adjustment of ymin for linear scale.
  //         if minimum is not set , then ymin is set to zero if >0
  //         or to ymin - margin if <0.
  if (!minimum) {
    if (xmin >= 0) xmin = 0;
    else           xmin -= yMARGIN*(xmax-xmin);
  }
  
  //         final adjustment of YMAXI for linear scale (if not option "Same"):
  //         decrease histogram height to MAX% of allowed height if HMAXIM
  //         has not been called.
  if (!maximum && !Healoption.Plus) {
    xmax += yMARGIN*(xmax-xmin);
  }
  Healparam.xmin = xmin;
  Healparam.xmax = xmax;
  */
  return 1;
}

//______________________________________________________________________________
void THealPainter::PaintFrame()
{
  printf("*****PaintFrame*****\n");
  if(Healoption.Same){
    return;
  } // if

  RecalculateRange();

  if(Healoption.Lego || Healoption.Surf || Healoption.Tri ||
     Healoption.Contour == 14 || Healoption.Error >= 100) {
    TObject* frame = gPad->FindObject("TFrame");
    if(frame){
      gPad->GetListOfPrimitives()->Remove(frame);
    } // if
    return;
  } // if
  gPad->PaintPadFrame(Healparam.xmin, Healparam.ymin, Healparam.xmax, Healparam.ymax);
}

//______________________________________________________________________________
void THealPainter::PaintLego(Option_t* option)
{
}

//______________________________________________________________________________
void THealPainter::PaintPalette()
{
  THealPaletteAxis* palette = (THealPaletteAxis*)fFunctions->FindObject("palette");
  TView* view = gPad->GetView();
  if(palette){
    if(view){
      if(!palette->TestBit(THealPaletteAxis::kHasView)) {
	delete palette;
	palette = 0;
      } // if
    } else {
      if(palette->TestBit(THealPaletteAxis::kHasView)){
	delete palette; palette = 0;
      } // if
    } // if
  } // if
  
  if(!palette){
    Double_t xup  = gPad->GetUxmax();
    Double_t x2   = gPad->PadtoX(gPad->GetX2());
    Double_t ymin = gPad->PadtoY(gPad->GetUymin());
    Double_t ymax = gPad->PadtoY(gPad->GetUymax());
    Double_t xr   = 0.05*(gPad->GetX2() - gPad->GetX1());
    Double_t xmin = gPad->PadtoX(xup + 0.1*xr);
    Double_t xmax = gPad->PadtoX(xup + xr);
    if(xmax > x2){
      xmax = gPad->PadtoX(gPad->GetX2() - 0.01*xr);
    } // if
    palette = new THealPaletteAxis(xmin, ymin, xmax, ymax, fHeal);
    palette->Paint();
  } // if
}

//______________________________________________________________________________
void THealPainter::PaintStat(Int_t dostat, TF1* fit)
{
}

//______________________________________________________________________________
void THealPainter::PaintSurface(Option_t* option)
{
}

//______________________________________________________________________________
void THealPainter::PaintTable(Option_t* option)
{
  printf("*****PaintTable*****\n");
  if(!TableInit()){
    return;
  } // if

  PaintFrame();

  if(fHeal->GetEntries() != 0 && Healoption.Axis <= 0){
    //    if(Healoption.Scat)    PaintScatterPlot(option);
    //    if(Healoption.Arrow)   PaintArrows(option);
    //    if(Healoption.Box)     PaintBoxes(option);
    if(Healoption.Color)   PaintColorLevels(option);
    //    if(Healoption.Contour) PaintContour(option);
    //    if(Healoption.Text)    PaintText(option);
    //    if(Healoption.Error >= 100)   Paint2DErrors(option);
  } // if

  if(Healoption.Lego) PaintLego(option);
  if(Healoption.Surf && !Healoption.Contour) PaintSurface(option);
  if(Healoption.Tri) PaintTriangles(option);
  
  if(!Healoption.Lego && !Healoption.Surf &&
     !Healoption.Tri  && !(Healoption.Error >= 100)) PaintAxis(kFALSE);

  PaintTitle();
}

//______________________________________________________________________________
void THealPainter::PaintTitle()
{
  printf("*****PaintTitle*****\n");
  if(Healoption.Same){
    return;
  } // if

  if(fHeal->TestBit(THealPix::kNoTitle)){
    return;
  } // if

  Int_t nt = strlen(fHeal->GetTitle());
  TPaveText* title = 0;
  TObject* obj;
  TIter next(gPad->GetListOfPrimitives());
  while((obj = next())) {
    if(!obj->InheritsFrom(TPaveText::Class())){
      continue;
    } // if
    title = (TPaveText*)obj;
    if(strcmp(title->GetName(), "title")){
      title = 0;
      continue;
    } // if
    break;
  } // while

  if(nt == 0 || gStyle->GetOptTitle() <= 0){
    if(title){
      delete title;
    } // if
    return;
  } // if

  Double_t ht = gStyle->GetTitleH();
  Double_t wt = gStyle->GetTitleW();
  if(ht <= 0){
    ht = 1.1*gStyle->GetTitleFontSize();
  } // if
  if(ht <= 0){
    ht = 0.05;
  } // if
  if(wt <= 0){
    TLatex l;
    l.SetTextSize(ht);
    l.SetTitle(fHeal->GetTitle());
    // adjustment in case the title has several lines (#splitline)
    ht = TMath::Max(ht, 1.2*l.GetYsize()/(gPad->GetY2() - gPad->GetY1()));
    Double_t wndc = l.GetXsize()/(gPad->GetX2() - gPad->GetX1());
    wt = TMath::Min(0.7, 0.02+wndc);
  } // if

  if(title){
    TText *t0 = (TText*)title->GetLine(0);
    if(t0){
      if(!strcmp(t0->GetTitle(), fHeal->GetTitle())){
	return;
      } // if
      t0->SetTitle(fHeal->GetTitle());
      if(wt > 0){
	title->SetX2NDC(title->GetX1NDC() + wt);
      } // if
    }
    return;
  }

  Int_t talh = gStyle->GetTitleAlign()/10;
  if(talh < 1){
    talh = 1;
    if(talh > 3){
      talh = 3;
    } // if
  } // if
  Int_t talv = gStyle->GetTitleAlign()%10;
  if(talv < 1){
    talv = 1;
    if(talv > 3){
      talv = 3;
    } // if
  } // if
  Double_t xpos, ypos;
  xpos = gStyle->GetTitleX();
  ypos = gStyle->GetTitleY();
  if(talh == 2) xpos = xpos-wt/2.;
  if(talh == 3) xpos = xpos-wt;
  if(talv == 2) ypos = ypos+ht/2.;
  if(talv == 1) ypos = ypos+ht;

  TPaveText *ptitle = new TPaveText(xpos, ypos - ht, xpos + wt, ypos, "blNDC");

  //     box with the histogram title
  ptitle->SetFillColor(gStyle->GetTitleFillColor());
  ptitle->SetFillStyle(gStyle->GetTitleStyle());
  ptitle->SetName("title");
  ptitle->SetBorderSize(gStyle->GetTitleBorderSize());
  ptitle->SetTextColor(gStyle->GetTitleTextColor());
  ptitle->SetTextFont(gStyle->GetTitleFont(""));
  if(gStyle->GetTitleFont("")%10 > 2){
    ptitle->SetTextSize(gStyle->GetTitleFontSize());
  } // if
  ptitle->AddText(fHeal->GetTitle());
  ptitle->SetBit(kCanDelete);
  ptitle->Draw();
  ptitle->Paint();
}

//______________________________________________________________________________
void THealPainter::PaintTriangles(Option_t *option)
{
}

//______________________________________________________________________________
void THealPainter::ProcessMessage(const char *mess, const TObject *obj)
{
}

//______________________________________________________________________________
void THealPainter::RecalculateRange()
{
  if(Healoption.Same){
    return;
  } // if

  Double_t xmin = Healparam.xmin;
  Double_t xmax = Healparam.xmax;
  Double_t ymin = Healparam.ymin;
  Double_t ymax = Healparam.ymax;

  Double_t xmin_aid, ymin_aid, xmax_aid, ymax_aid;
  if(Healoption.Proj == 1){
  } else if(Healoption.Proj == 2){
  } else if(Healoption.Proj == 3){
  } else if(Healoption.Proj == 4){
  } // if

  Healparam.xmin = xmin;
  Healparam.xmax = xmax;
  Healparam.ymin = ymin;
  Healparam.ymax = ymax;

  Double_t dx = xmax - xmin;
  Double_t dy = ymax - ymin;
  Double_t dxr  = dx/(1 - gPad->GetLeftMargin()   - gPad->GetRightMargin());
  Double_t dyr  = dy/(1 - gPad->GetBottomMargin() - gPad->GetTopMargin());

  gPad->Range(xmin - dxr*gPad->GetLeftMargin(),
	      ymin - dyr*gPad->GetBottomMargin(),
	      xmax + dxr*gPad->GetRightMargin(),
	      ymax + dyr*gPad->GetTopMargin());
  gPad->RangeAxis(xmin, ymin, xmax, ymax);
}

//______________________________________________________________________________
void THealPainter::SetHealPix(THealPix* heal)
{
  if(heal == 0) return;
  fHeal = heal;
  fXaxis = heal->GetXaxis();
  fYaxis = heal->GetYaxis();
  fZaxis = heal->GetZaxis();
  fFunctions = fHeal->GetListOfFunctions();
}

//______________________________________________________________________________
Int_t THealPainter::TableInit()
{
  printf("*****TableInit*****\n");  

  static const char* where = "TableInit";
  
  Int_t first, last;
  Double_t yMARGIN = gStyle->GetHistTopMargin();
  Double_t zmin, zmax;
  Int_t maximum = 0;
  Int_t minimum = 0;
  if(fHeal->GetMaximumStored() != -1111) maximum = 1;
  if(fHeal->GetMinimumStored() != -1111) minimum = 1;
  
  //    -----------------  Compute X axis parameters
  first           = fXaxis->GetFirst();
  last            = fXaxis->GetLast();
  Healparam.xlast    = last;
  Healparam.xfirst   = first;
  Healparam.xlowedge = fXaxis->GetBinLowEdge(first);
  Healparam.xbinsize = fXaxis->GetBinWidth(first);
  Healparam.xmin     = Healparam.xlowedge;
  Healparam.xmax     = fXaxis->GetBinLowEdge(last) + fXaxis->GetBinWidth(last);
  
  //       if log scale in X, replace xmin,max by the log
  if(Healoption.Logx){
    //   find the first edge of a bin that is > 0
    if(Healparam.xlowedge <= 0){
      Healparam.xlowedge = fXaxis->GetBinUpEdge(fXaxis->FindFixBin(0.01*Healparam.xbinsize));
      Healparam.xmin  = Healparam.xlowedge;
    } // if
    if(Healparam.xmin <= 0 || Healparam.xmax <= 0){
      Error(where, "cannot set X axis to log scale");
      return 0;
    } // if
    Healparam.xfirst= fXaxis->FindFixBin(Healparam.xmin);
    if(Healparam.xfirst < first){
      Healparam.xfirst = first;
    } // if
    Healparam.xlast = fXaxis->FindFixBin(Healparam.xmax);
    if(Healparam.xlast > last){
      Healparam.xlast = last;
    } // if
    Healparam.xmin  = TMath::Log10(Healparam.xmin);
    Healparam.xmax  = TMath::Log10(Healparam.xmax);
  }
  
  //    -----------------  Compute Y axis parameters
  first           = fYaxis->GetFirst();
  last            = fYaxis->GetLast();
  Healparam.ylast    = last;
  Healparam.yfirst   = first;
  Healparam.ylowedge = fYaxis->GetBinLowEdge(first);
  Healparam.ybinsize = fYaxis->GetBinWidth(first);
  if(!Healparam.ybinsize){
    Healparam.ybinsize = 1;
  } // if
  Healparam.ymin     = Healparam.ylowedge;
  Healparam.ymax     = fYaxis->GetBinLowEdge(last) + fYaxis->GetBinWidth(last);
  
  //       if log scale in Y, replace ymin,max by the log
  if(Healoption.Logy){
    if(Healparam.ylowedge <=0){
      Healparam.ylowedge = fYaxis->GetBinUpEdge(fYaxis->FindFixBin(0.01*Healparam.ybinsize));
      Healparam.ymin  = Healparam.ylowedge;
    } // if
    if(Healparam.ymin <=0 || Healparam.ymax <=0){
      Error(where, "cannot set Y axis to log scale");
      return 0;
    } // if
    Healparam.yfirst = fYaxis->FindFixBin(Healparam.ymin);
    if(Healparam.yfirst < first){
      Healparam.yfirst = first;
    } // if
    Healparam.ylast = fYaxis->FindFixBin(Healparam.ymax);
    if(Healparam.ylast > last){
      Healparam.ylast = last;
    } // if
    Healparam.ymin  = TMath::Log10(Healparam.ymin);
    Healparam.ymax  = TMath::Log10(Healparam.ymax);
  }
  
  //    -----------------  Compute Z axis parameters
  Double_t bigp = TMath::Power(10, 32);
  zmax = -bigp;
  zmin = bigp;
  Double_t allchan = 0;
  for(Int_t i = 0; i <= fHeal->GetNpix(); i++){
    Double_t c1 = fHeal->GetBinContent(i);
    zmax = TMath::Max(zmax,c1);
    if(Healoption.Error){
      Double_t e1 = fHeal->GetBinError(i);
      zmax = TMath::Max(zmax, c1 + e1);
    } // if
    zmin = TMath::Min(zmin, c1);
    allchan += c1;
  } // i
  
  //     Take into account maximum , minimum
  
  if(maximum) zmax = fHeal->GetMaximumStored();
  if(minimum) zmin = fHeal->GetMinimumStored();
  if(Healoption.Logz && zmax <= 0) {
    if(!Healoption.Same){
      Error(where, "log scale is requested but maximum is less or equal 0 (%f)", zmax);
    } // if
    return 0;
  }
  if(zmin >= zmax && !Healoption.Plus){
    if(Healoption.Logz){
      if(zmax > 0){
	zmin = 0.001*zmax;
      } else {
	if(!Healoption.Same){
	  Error(where, "log scale is requested but maximum is less or equal 0 (%f)", zmax);
	} // if
	return 0;
      } // if
    } // if
  } // if
  
  //     take into account normalization factor
  Healparam.allchan = allchan;
  Double_t factor = allchan;
  if(fHeal->GetNormFactor() > 0) factor = fHeal->GetNormFactor();
  if(allchan) factor /= allchan;
  if(factor == 0) factor = 1;
  Healparam.factor = factor;
  zmax = factor*zmax;
  zmin = factor*zmin;
  Double_t c1 = zmax;
  if(TMath::Abs(zmin) > TMath::Abs(c1)){
    c1 = zmin;
  } // if
  
  //         For log scales, histogram coordinates are log10(ymin) and
  //         log10(ymax). Final adjustment (if not option "Same")
  //         or "+" for ymax) of ymax and ymin for logarithmic scale, if
  //         Maximum and Minimum are not defined.
  if(Healoption.Logz){
    if(zmin <= 0){
      zmin = TMath::Min((Double_t)1, (Double_t)0.001*zmax);
      fHeal->SetMinimum(zmin);
    } // if
    zmin = TMath::Log10(zmin);
    if(!minimum){
      zmin += TMath::Log10(0.5);
    } // if
    zmax = TMath::Log10(zmax);
    if(!maximum && !Healoption.Plus){
      zmax += TMath::Log10(2*(0.9/0.95));
    } // if
    goto LZMIN;
  }
  
  //         final adjustment of YMAXI for linear scale (if not option "Same"):
  //         decrease histogram height to MAX% of allowed height if HMAXIM
  //         has not been called.
  //         MAX% is the value in percent which has been set in HPLSET
   //         (default is 90%).
  if(!maximum && !Healoption.Plus){
    zmax += yMARGIN*(zmax - zmin);
  } // if
  
  //         final adjustment of ymin for linear scale.
  //         if minimum is not set , then ymin is set to zero if >0
  //         or to ymin - yMARGIN if <0.
  if(!minimum){
    if(gStyle->GetHistMinimumZero()){
      if(zmin >= 0){
	zmin = 0;
      } else {
	zmin -= yMARGIN*(zmax - zmin);
      } // if
    } else {
      Double_t dzmin = yMARGIN*(zmax - zmin);
      if(zmin >= 0 && (zmin - dzmin <= 0)){
	zmin  = 0;
      } else {
	zmin -= dzmin;
      } // if
    } // if
  } // if
  
 LZMIN:
  Healparam.zmin = zmin;
  Healparam.zmax = zmax;
  
  //     Set bar offset and width
  Healparam.baroffset = fHeal->GetBarOffset();
  Healparam.barwidth  = fHeal->GetBarWidth();
  
  return 1;
}
