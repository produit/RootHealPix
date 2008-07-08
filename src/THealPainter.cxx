// $Id: THealPainter.cxx,v 1.1 2008/07/08 08:35:49 oxon Exp $
// Author: Akira Okumura 2008/07/07

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include "TClass.h"
#include "THealPainter.h"
#include "TPad.h"
#include "TROOT.h"
#include "TVirtualPadEditor.h"

THealPix* gCurrentHeal = 0;

ClassImp(THealPainter)

//______________________________________________________________________________
THealPainter::THealPainter()
{
  fHealPix = 0;
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
  gCurrentHeal = fHealPix;
  if(!gPad){
    Error("DrawPanel", "need to draw HEALPix first");
    return;
  } // if
  TVirtualPadEditor* editor = TVirtualPadEditor::GetPadEditor();
  editor->Show();
  gROOT->ProcessLine(Form("((TCanvas*)0x%lx)->Selected((TVirtualPad*)0x%lx, (TObject*)0x%lx, 1)", gPad->GetCanvas(), gPad, fHealPix));
}

//______________________________________________________________________________
void THealPainter::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
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
void THealPainter::Paint(Option_t* option)
{
}

//______________________________________________________________________________
void THealPainter::ProcessMessage(const char *mess, const TObject *obj)
{
}

//______________________________________________________________________________
void THealPainter::SetHealPix(THealPix* heal)
{
  if(heal == 0) return;
  fHealPix = heal;
}
