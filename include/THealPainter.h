// $Id: THealPainter.h,v 1.1 2008/07/08 08:35:48 oxon Exp $
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

class THealPix;
class TAxis;

class THealPainter : public TVirtualHealPainter {

protected:
   THealPix*      fHealPix; //pointer to HEALPix to paint

public:
   THealPainter();
   virtual ~THealPainter();
   virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py);
   virtual void   DrawPanel();
   virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual Bool_t IsInside(Int_t x, Int_t y);
   virtual Bool_t IsInside(Double_t x, Double_t y);
   virtual void   Paint(Option_t *option="");
   virtual void   ProcessMessage(const char* mess, const TObject* obj);

   virtual void   SetHealPix(THealPix *hp);

   ClassDef(THealPainter, 0)  //Helper class to draw HEALPix
};

#endif // T_HEAL_PAINTER
