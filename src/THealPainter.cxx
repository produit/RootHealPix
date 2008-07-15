// $Id: THealPainter.cxx,v 1.4 2008/07/15 14:48:58 oxon Exp $
// Author: Akira Okumura 2008/07/07

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include <float.h>

#include "Healoption.h"
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

static Healoption_t Healoption;
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
Bool_t THealPainter::IsInside(Int_t bin /*ix*/, Int_t /*iy*/)
{
  for (Int_t i=0;i<fNcuts;i++) {
    Double_t x, y;
    fHeal->GetBinCenter(bin, x, y);
    if(!fHeal->IsDegree()){
      x *= TMath::RadToDeg();
      y *= TMath::RadToDeg();
    } // if
    if(Healoption.System == kGalactic || Healoption.System == kLatLong){
      x -= 180;
      y = 90 - y;
    } else if(Healoption.System == kCelestial){
      y = 90 - y;
    } // if

    if (fCutsOpt[i] > 0) {
      if (!fCuts[i]->IsInside(x,y)) return kFALSE;
    } else {
      if (fCuts[i]->IsInside(x,y))  return kFALSE;
    }
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t THealPainter::IsInside(Double_t x, Double_t y)
{
  for(Int_t i = 0; i < fNcuts; i++){
    if(fCutsOpt[i] > 0){
      if(!fCuts[i]->IsInside(x,y)) return kFALSE;
    } else {
      if(fCuts[i]->IsInside(x,y))  return kFALSE;
    } // if
  } // if
  return kTRUE;
}

//______________________________________________________________________________
Int_t THealPainter::MakeChopt(Option_t* choptin)
{
  if(gPad->GetLogx() || gPad->GetLogy()){
    Error("MakeChopt", "logarithmic X/Y axis is not supported.");
  } // if

  char *l;
  char chopt[128];
  Int_t nch = strlen(choptin);
  strcpy(chopt,choptin);
  
  Healoption.Axis = Healoption.Heal    = Healoption.Same    = Healoption.Func =
  Healoption.Scat = Healoption.Color   = Healoption.Contour = Healoption.Logz =
  Healoption.Lego = Healoption.Surf    = Healoption.Off     = Healoption.Tri  =
  Healoption.AxisPos = 0;
  
  //    special 2D options
  Healoption.List     = 0;
  Healoption.Zscale   = 0;
  Healoption.FrontBox = 1;
  Healoption.BackBox  = 1;
  Healoption.System   = kThetaPhi;
  Healoption.Proj     = kEquirect;
  
  Healoption.HighRes  = 0;
  
  Healoption.Zero     = 0;
  
  //check for graphical cuts
  MakeCuts(chopt);
  
  for (Int_t i = 0;i < nch; i++) chopt[i] = toupper(chopt[i]);
  Healoption.Scat = 1;
  if(!nch) Healoption.Heal = 1;
  if(fFunctions->First()) Healoption.Func = 2;
  
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
  if((Healoption.AxisPos == 10 || Healoption.AxisPos == 1) && (nch == 2)) Healoption.Heal = 1;
  if(Healoption.AxisPos == 11 && nch == 4) Healoption.Heal = 1;
  
  // Coordinate options
  l = strstr(chopt, "THETAPHI");
  if(l){
    if(nch == 8) Healoption.Heal = 1;
    Healoption.System = kThetaPhi;
    strncpy(l, "        ", 8);
  } // if
  l = strstr(chopt, "GALACTIC");
  if(l){
    if(nch == 8) Healoption.Heal = 1;
    Healoption.System = kGalactic;
    strncpy(l, "        ", 8);
  } // if
  l = strstr(chopt, "CELESTIAL");
  if(l){
    if(nch == 9) Healoption.Heal = 1;
    Healoption.System = kCelestial;
    strncpy(l, "         ", 9);
  } // if
  l = strstr(chopt, "LATLONG");
  if(l){
    if(nch == 7) Healoption.Heal = 1;
    Healoption.System = kLatLong;
    strncpy(l, "       ", 7);
  } // if

  // Projection options
  l = strstr(chopt, "EQUIRECT");
  if(l){
    if(nch == 8) Healoption.Heal = 1;
    Healoption.Proj = kEquirect;
    strncpy(l, "        ", 8);
  } // if
  l = strstr(chopt, "HAMMER");
  if(l){
    if(nch == 6) Healoption.Heal = 1;
    Healoption.Proj = kHammer;
    strncpy(l, "      ", 6);
  } // if

  l = strstr(chopt,"SAMES");
  if (l) {
    if (nch == 5) Healoption.Heal = 1;
    Healoption.Same = 2;
    strncpy(l,"     ",5);
  }
  l = strstr(chopt,"SAME");
  if (l) {
    if (nch == 4) Healoption.Heal = 1;
    Healoption.Same = 1;
    strncpy(l,"    ",4);
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
  l = strstr(chopt,"COLZ"); if (l) { Healoption.Color  = 2; strncpy(l,"    ",4); Healoption.Scat = 0; Healoption.Zscale = 1;}
  l = strstr(chopt,"COL" ); if (l) { Healoption.Color  = 1; strncpy(l,"   ", 3); Healoption.Scat = 0; }
  l = strstr(chopt,"FUNC"); if (l) { Healoption.Func   = 2; strncpy(l,"    ",4); Healoption.Heal = 0; }
  l = strstr(chopt,"HEAL"); if (l) { Healoption.Heal   = 2; strncpy(l,"    ",4); Healoption.Func = 0;}
  l = strstr(chopt,"AXIS"); if (l) { Healoption.Axis   = 1; strncpy(l,"    ",4); }
  l = strstr(chopt,"AXIG"); if (l) { Healoption.Axis   = 2; strncpy(l,"    ",4); }
  l = strstr(chopt,"SCAT"); if (l) { Healoption.Scat   = 1; strncpy(l,"    ",4); }

  l = strstr(chopt, "AITOFF");
  if(l){
    Healoption.Proj = kEquirect;
    strncpy(l, "      ", 6);
  } // if
  l = strstr(chopt, "LAMBERT");
  if(l){
    Healoption.Proj = kLambert;
    strncpy(l, "       ", 7);
  } // if
  l = strstr(chopt, "HAMMER");
  if (l) {
    Healoption.Proj = 2;
    strncpy(l, "      ", 6);
  } // if
  l = strstr(chopt, "EQUIRECT");
  if (l) {
    Healoption.Proj = kEquirect;
    strncpy(l, "       ", 7);
  } // if
  
  if (strstr(chopt,"A"))   Healoption.Axis = -1;
  if (strstr(chopt,"][")) {Healoption.Off  =1; Healoption.Heal =1;}
  if (strstr(chopt,"Z"))   Healoption.Zscale =1;
  if (strstr(chopt,"H"))   Healoption.Heal =2;
  if (strstr(chopt,"9"))  Healoption.HighRes = 1;
  
  if (Healoption.Surf == 15) {
    if (Healoption.System == kPOLAR || Healoption.System == kCARTESIAN) {
      Healoption.Surf = 13;
      Warning("MakeChopt","option SURF5 is not supported in Cartesian and Polar modes");
    }
  }
  
  //      Copy options from current style
  Healoption.Logz = gPad->GetLogz(); // Z only
  
  //       Check options incompatibilities

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

  if(Healoption.Proj == kLambert){
    printf("lambert\n");
    fXaxis->Set(1, -1., 1.);
    fYaxis->Set(1, -1., 1.);
  } else {
    if(Healoption.System == kThetaPhi){
      fXaxis->Set(1, 0., 360.);
      fYaxis->Set(1, 0., 180.);
    } else if(Healoption.System == kGalactic){
      fXaxis->Set(1, -180., 180.);
      fYaxis->Set(1, -90., 90.);
    } else if(Healoption.System == kCelestial){
      fXaxis->Set(1, 0., 360.);
      fYaxis->Set(1, -90., 90.);
    } else if(Healoption.System == kLatLong){
      fXaxis->Set(1, -180., 180.);
      fYaxis->Set(1, -90., 90.);
    } // if
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
  if(Healoption.Func){
    Healoption_t hoptsave = Healoption;
    Hparam_t  hparsave = Healparam;
    PaintFunction(option);
    SetHealPix(healsave);
    Healoption = hoptsave;
    Healparam  = hparsave;
  }
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
  
  if(Healoption.System == kGalactic || Healoption.System == kLatLong){
    fXaxis->Set(1, -180., 180.);
    fYaxis->Set(1, -90., 90.);
  } else if(Healoption.System == kCelestial){
    fXaxis->Set(1, 0., 360.);
    fYaxis->Set(1, -90., 90.);
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
   if (ndivx < 0) strcat(chopt, "N");
   if (gPad->GetGridx()) {
      gridl = (aymax-aymin)/(gPad->GetY2() - gPad->GetY1());
      strcat(chopt, "W");
   }

   // Define X-Axis limits
   ndiv = TMath::Abs(ndivx);
   if (useHealparam) {
     umin = Healparam.xmin;
     umax = Healparam.xmax;
   } else {
     umin = axmin;
     umax = axmax;
   }

   // Display axis as time
   if (fXaxis->GetTimeDisplay()) {
      strcat(chopt,"t");
      if (strlen(fXaxis->GetTimeFormatOnly()) == 0) {
         axis.SetTimeFormat(fXaxis->ChooseTimeFormat(Healparam.xmax-Healparam.xmin));
      }
   }

   // The main X axis can be on the bottom or on the top of the pad
   Double_t xAxisYPos1, xAxisYPos2;
   if (xAxisPos == 1) {
      // Main X axis top
      xAxisYPos1 = aymax;
      xAxisYPos2 = aymin;
   } else {
      // Main X axis bottom
      xAxisYPos1 = aymin;
      xAxisYPos2 = aymax;
   }

   // Paint the main X axis (always)
   uminsave = umin;
   umaxsave = umax;
   ndivsave = ndiv;
   axis.SetOption(chopt);
   if (xAxisPos) {
      strcat(chopt, "-");
      gridl = -gridl;
   }
   axis.PaintAxis(axmin, xAxisYPos1,
                  axmax, xAxisYPos1,
                  umin, umax,  ndiv, chopt, gridl, drawGridOnly);

   // Paint additional X axis (if needed)
   if (gPad->GetTickx()) {
      if (xAxisPos) {
         cw=strstr(chopt,"-");
         *cw='z';
      } else {
         strcat(chopt, "-");
      }
      if (gPad->GetTickx() < 2) strcat(chopt, "U");
      if ((cw=strstr(chopt,"W"))) *cw='z';
      axis.SetTitle("");
      axis.PaintAxis(axmin, xAxisYPos2,
                     axmax, xAxisYPos2,
                     uminsave, umaxsave,  ndivsave, chopt, gridl, drawGridOnly);
   }

   // Paint Y axis
   Int_t ndivy = fYaxis->GetNdivisions();
   axis.ImportAxisAttributes(fYaxis);

   chopt[0] = 0;
   strcat(chopt, "SDH");
   if (ndivy < 0) strcat(chopt, "N");
   if (gPad->GetGridy()) {
      gridl = (axmax-axmin)/(gPad->GetX2() - gPad->GetX1());
      strcat(chopt, "W");
   }

   // Define Y-Axis limits
   ndiv = TMath::Abs(ndivy);
   if (useHealparam) {
     umin = Healparam.ymin;
     umax = Healparam.ymax;
   } else {
     umin = aymin;
     umax = aymax;
   }

   // Display axis as time
   if (fYaxis->GetTimeDisplay()) {
      strcat(chopt,"t");
      if (strlen(fYaxis->GetTimeFormatOnly()) == 0) {
         axis.SetTimeFormat(fYaxis->ChooseTimeFormat(Healparam.ymax-Healparam.ymin));
      }
   }

   // The main Y axis can be on the left or on the right of the pad
   Double_t yAxisXPos1, yAxisXPos2;
   if (yAxisPos == 1) {
      // Main Y axis left
      yAxisXPos1 = axmax;
      yAxisXPos2 = axmin;
   } else {
      // Main Y axis right
      yAxisXPos1 = axmin;
      yAxisXPos2 = axmax;
   }

   // Paint the main Y axis (always)
   uminsave = umin;
   umaxsave = umax;
   ndivsave = ndiv;
   axis.SetOption(chopt);
   if (yAxisPos) {
      strcat(chopt, "+L");
      gridl = -gridl;
   }
   axis.PaintAxis(yAxisXPos1, aymin,
                  yAxisXPos1, aymax,
                  umin, umax,  ndiv, chopt, gridl, drawGridOnly);

   // Paint the additional Y axis (if needed)
   if (gPad->GetTicky()) {
      if (gPad->GetTicky() < 2) {
         strcat(chopt, "U");
         axis.SetTickSize(-fYaxis->GetTickLength());
      } else {
         strcat(chopt, "+L");
      }
      if ((cw=strstr(chopt,"W"))) *cw='z';
      axis.SetTitle("");
      axis.PaintAxis(yAxisXPos2, aymin,
                     yAxisXPos2, aymax,
                     uminsave, umaxsave,  ndivsave, chopt, gridl, drawGridOnly);
   }

  /*
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
  if(Healoption.System == kThetaPhi || Healoption.System == kLatLong){
    strcat(chopt, "SDHR+");
  } else {
    strcat(chopt, "SDHL-");
    axis.SetLabelOffset(-TMath::Abs(axis.GetLabelOffset()));
  } // if
  if(ndivx < 0){
    strcat(chopt, "N");
  } // if
  if(gPad->GetGridx()){
    if(Healoption.System == kThetaPhi || Healoption.System == kLatLong){
      gridl = (aymax - aymin)/(gPad->GetY2() - gPad->GetY1());
    } else {
      gridl = -(aymax - aymin)/(gPad->GetY2() - gPad->GetY1());
    } // if
    strcat(chopt, "W");
  } // if
  
  // Define X-Axis limits
  ndiv = TMath::Abs(ndivx);
  if(useHealparam){
    umin = Healparam.xmin;
    umax = Healparam.xmax;
  } else {
    umin = axmin;
    umax = axmax;
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
    if(Healoption.System == kThetaPhi || Healoption.System == kLatLong){
      strcat(chopt, "R+");
    } else {
      char* c = strstr(chopt, "R+");
      if(c) strncpy(c, "  ", 2);
      strcat(chopt, "L-");
      axis.SetLabelOffset(-TMath::Abs(axis.GetLabelOffset()));
    } // if
    gridl = -gridl;
  } // if
  if(Healoption.System == kThetaPhi || Healoption.System == kLatLong){
    axis.PaintAxis(axmin, xAxisYPos1, axmax, xAxisYPos1,
		   umin, umax,  ndiv, chopt, gridl, drawGridOnly);
  } else {
    axis.PaintAxis(axmax, xAxisYPos1, axmin, xAxisYPos1,
		   umin, umax,  ndiv, chopt, gridl, drawGridOnly);
  } // if
  
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
  printf("%s\n", chopt);
  // Paint Y axis
  Int_t ndivy = fYaxis->GetNdivisions();
  axis.ImportAxisAttributes(fYaxis);
  axis.SetLabelOffset(-axis.GetLabelOffset());

  chopt[0] = 0;
  strcat(chopt, "SDHR+");
  if(ndivy < 0){
    strcat(chopt, "N");
  } // if
  if(gPad->GetGridy()){
    gridl = (axmax - axmin)/(gPad->GetX2() - gPad->GetX1());
    strcat(chopt, "W");
  } // if
  
  // Define Y-Axis limits
  ndiv = TMath::Abs(ndivy);
  if(useHealparam){
    umin = Healparam.ymin;
    umax = Healparam.ymax;
  } else {
    umin = aymin;
    umax = aymax;
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
    char* c = strstr(chopt, "R+");
    if(c) strncpy(c, "  ", 2);
    strcat(chopt, "L-");
    axis.SetLabelOffset(-axis.GetLabelOffset());
    gridl = -gridl;
  } // if

  axis.PaintAxis(yAxisXPos1, aymax, yAxisXPos1 - 1e-3, aymin,
		 umin, umax,  ndiv, chopt, gridl, drawGridOnly);
  
  // Paint the additional Y axis (if needed)
  if(gPad->GetTicky()){
    if(gPad->GetTicky() < 2){
      strcat(chopt, "U");
      axis.SetTickSize(-fYaxis->GetTickLength());
    } else {
      strcat(chopt, "-L");
    } // if
    if((cw=strstr(chopt,"W"))){
      *cw='z';
    } // if
    axis.SetTitle("");
    axis.PaintAxis(yAxisXPos2, aymax, yAxisXPos2 - 1e-3, aymin,
		   uminsave, umaxsave,  ndivsave, chopt, gridl, drawGridOnly);
  } // if
  printf("%s\n", chopt);
  */
}

//______________________________________________________________________________
void THealPainter::PaintColorLevels(Option_t* option)
{
  printf("*****PaintColorLevels*****\n");

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
  Double_t xmin = gPad->GetUxmin();
  Double_t xmax = gPad->GetUxmax();
  Double_t ymin = gPad->GetUymin();
  Double_t ymax = gPad->GetUymax();

  for(Int_t bin = 0; bin < fHeal->GetNpix(); bin++){
    Double_t z = fHeal->GetBinContent(bin);
    if(z == 0 && (zmin >= 0 || Healoption.Logz)){
      continue; // don't draw the empty bins for histograms with positive content
    } // if
    if(Healoption.Logz){
      z = z > 0 ? TMath::Log10(z) : zmin;
    } // if
    if(z < zmin){
      continue;
    } // if

    Double_t theta, phi; // x/y center
    fHeal->GetBinCenter(bin, theta, phi);
    if(!fHeal->IsDegree()){
      theta *= TMath::RadToDeg();
      phi *= TMath::RadToDeg();
    } // if
    if(Healoption.System == kGalactic || Healoption.System == kLatLong){
      theta = 90 - theta;
      phi -= 180;
    } else if(Healoption.System == kCelestial){
      theta = 90 - theta;
    } // if

    if(!IsInside(phi, theta)){
      continue;
    } // if

    Double_t x[5], y[5];
    Int_t n = fHeal->GetBinVertices(bin, x, y);
    if(Healoption.System == kGalactic || Healoption.System == kLatLong){
      if(x[0] >= 180. && x[1] >= 180. && x[2] >= 180. && x[3] >= 180. &&
	 (n == 4 || (n == 5 && x[4] >= 180.))){
	for(Int_t j = 0; j < n; j++){
	  x[j] -= 360.;
	} // j
      } // if
      for(Int_t j = 0; j < n; j++){
	y[j] = 90. - y[j];
      } // j
    } else if(Healoption.System == kCelestial){
      for(Int_t j = 0; j < n; j++){
	y[j] = 90. - y[j];
      } // j
    } // if

    if(Healoption.Proj == kHammer){
      if(Healoption.System == kGalactic || Healoption.System == kLatLong){
	for(Int_t j = 0; j < n; j++){
	  Double_t lng = x[j]*TMath::DegToRad();
	  Double_t lat = y[j]*TMath::DegToRad();
	  x[j] = 180*TMath::Cos(lat)*TMath::Sin(lng/2.)
	    /TMath::Sqrt(1. + TMath::Cos(lat)*TMath::Cos(lng/2.));
	  y[j] = 90.*TMath::Sin(lat)
	    /TMath::Sqrt(1. + TMath::Cos(lat)*TMath::Cos(lng/2.));
	} // j
      } else if(Healoption.System == kCelestial){
	for(Int_t j = 0; j < n; j++){
	  Double_t lng = (x[j] - 180.)*TMath::DegToRad();
	  Double_t lat = y[j]*TMath::DegToRad();
	  x[j] = 180. + 180*TMath::Cos(lat)*TMath::Sin(lng/2.)
	    /TMath::Sqrt(1. + TMath::Cos(lat)*TMath::Cos(lng/2.));
	  y[j] = 90.*TMath::Sin(lat)
	    /TMath::Sqrt(1. + TMath::Cos(lat)*TMath::Cos(lng/2.));
	} // j
      } else {
	for(Int_t j = 0; j < n; j++){
	  Double_t lng = (x[j] - 180.)*TMath::DegToRad();
	  Double_t lat = (90. - y[j])*TMath::DegToRad();
	  x[j] = 180. + 180*TMath::Cos(lat)*TMath::Sin(lng/2.)
	    /TMath::Sqrt(1. + TMath::Cos(lat)*TMath::Cos(lng/2.));
	  y[j] = 90. - 90.*TMath::Sin(lat)
	    /TMath::Sqrt(1. + TMath::Cos(lat)*TMath::Cos(lng/2.));
	} // j
      } // if
    } // if

    if(fHeal->TestBit(THealPix::kUserContour)){
      Double_t zc = fHeal->GetContourLevelPad(0);
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
    gPad->PaintFillArea(n, x, y);
  } // bin
  
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
  printf("*****PaintContour*****\n");  
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
     Healoption.Contour == 14) {
    TObject* frame = gPad->FindObject("TFrame");
    if(frame){
      gPad->GetListOfPrimitives()->Remove(frame);
    } // if
    return;
  } // if
  gPad->PaintPadFrame(Healparam.xmin, Healparam.ymin, Healparam.xmax, Healparam.ymax);
}

//______________________________________________________________________________
void THealPainter::PaintFunction(Option_t*)
{
  TObjOptLink* lnk = (TObjOptLink*)fFunctions->FirstLink();
  TObject* obj;
  
  while(lnk){
    obj = lnk->GetObject();
    TVirtualPad *padsave = gPad;
    if(obj->InheritsFrom(TF1::Class())){
      if (obj->TestBit(TF1::kNotDraw) == 0) obj->Paint("lsame");
    } else  {
      obj->Paint(lnk->GetOption());
    } // if
    lnk = (TObjOptLink*)lnk->Next();
    padsave->cd();
  } // while
}

//______________________________________________________________________________
void THealPainter::PaintLego(Option_t* option)
{
  printf("*****PaintLego*****\n");  
}

//______________________________________________________________________________
void THealPainter::PaintPalette()
{
  printf("*****PaintPalette*****\n");  
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
void THealPainter::PaintScatterPlot(Option_t* option)
{
  printf("*****PaintScatterPlot*****\n");  
}

//______________________________________________________________________________
void THealPainter::PaintStat(Int_t dostat, TF1* fit)
{
  printf("*****PaintStat*****\n");  
}

//______________________________________________________________________________
void THealPainter::PaintSurface(Option_t* option)
{
  printf("*****PaintSurface*****\n");  
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
    if(Healoption.Scat)    PaintScatterPlot(option);
    if(Healoption.Color)   PaintColorLevels(option);
    if(Healoption.Contour) PaintContour(option);
  } // if

  if(Healoption.Lego) PaintLego(option);
  if(Healoption.Surf && !Healoption.Contour) PaintSurface(option);
  if(Healoption.Tri) PaintTriangles(option);
  
  if(!Healoption.Lego && !Healoption.Surf && !Healoption.Tri) PaintAxis(kFALSE);

  PaintTitle();

  TF1* fit  = 0;
  TIter next(fFunctions);
  TObject *obj;
  while((obj = next())){
    if (obj->InheritsFrom(TF1::Class())) {
      fit = (TF1*)obj;
      break;
    } // if
  } // while
  if(Healoption.Same != 1){
    if (!fHeal->TestBit(THealPix::kNoStats)) {  // bit set via TH1::SetStats
      PaintStat(gStyle->GetOptStat(), fit);
    } // if
  } // if
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
void THealPainter::PaintTriangles(Option_t* option)
{
  printf("*****PaintTriangles*****\n");  
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
  
  //    -----------------  Compute Z axis parameters
  Double_t bigp = TMath::Power(10, 32);
  zmax = -bigp;
  zmin = bigp;
  Double_t allchan = 0;
  for(Int_t i = 0; i <= fHeal->GetNpix(); i++){
    Double_t c1 = fHeal->GetBinContent(i);
    zmax = TMath::Max(zmax,c1);
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
  if(zmin >= zmax){
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
    if(!maximum){
      zmax += TMath::Log10(2*(0.9/0.95));
    } // if
    goto LZMIN;
  }
  
  //         final adjustment of YMAXI for linear scale (if not option "Same"):
  //         decrease histogram height to MAX% of allowed height if HMAXIM
  //         has not been called.
  //         MAX% is the value in percent which has been set in HPLSET
   //         (default is 90%).
  if(!maximum){
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
