// $Id: THealPainter.h,v 1.2 2008/07/09 11:50:10 oxon Exp $
// Author: Akira Okumura 2008/07/07

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_PAINTER
#define T_HEAL_PAINTER

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THealPainter                                                         //
//                                                                      //
// helper class to draw HEALPixs                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "TVirtualHealPainter.h"

class TAxis;
class TCutG;
class THealPix;
const Int_t kMaxCuts = 16;

class THealPainter : public TVirtualHealPainter {

protected:
   THealPix*      fHeal;              //pointer to HEALPix to paint
   TAxis*         fXaxis;             //pointer to X axis
   TAxis*         fYaxis;             //pointer to Y axis
   TAxis*         fZaxis;             //pointer to Z axis
   TList*         fFunctions;         //pointer to histogram list of functions
   Double_t*      fXbuf;              //X buffer coordinates
   Double_t*      fYbuf;              //Y buffer coordinates
   Int_t          fNcuts;             //Number of graphical cuts
   Int_t          fCutsOpt[kMaxCuts]; //sign of each cut
   TCutG*         fCuts[kMaxCuts];    //Pointers to graphical cuts

public:
   THealPainter();
   virtual ~THealPainter();
   virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
   virtual void   DrawPanel();
   virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual TList* GetContourList(Double_t contour) const;
   virtual Bool_t IsInside(Int_t x, Int_t y);
   virtual Bool_t IsInside(Double_t x, Double_t y);
   virtual Int_t  MakeChopt(Option_t* option);
   virtual Int_t  MakeCuts(char* cutsopt);
   virtual void   Paint(Option_t* option = "");
   virtual void   PaintAxis(Bool_t drawGridOnly = kFALSE);
   virtual void   PaintBarH(Option_t* option);
   virtual void   PaintColorLevels(Option_t* option);
   virtual void   PaintContour(Option_t* option);
   virtual void   PaintFrame();
   virtual Int_t  PaintInit();
   virtual Int_t  PaintInitH();
   virtual void   PaintLego(Option_t* option);
   virtual void   PaintPalette();
   virtual void   PaintStat(Int_t dostat, TF1* fit);
   virtual void   PaintSurface(Option_t* option);
   virtual void   PaintTable(Option_t* option);
   virtual void   PaintTitle();
   virtual void   PaintTriangles(Option_t *option);
   virtual void   ProcessMessage(const char* mess, const TObject* obj);
   virtual void   RecalculateRange();
   virtual void   SetHealPix(THealPix *hp);
   virtual Int_t  TableInit();

   ClassDef(THealPainter, 0)  //Helper class to draw HEALPix
};

#endif // T_HEAL_PAINTER
