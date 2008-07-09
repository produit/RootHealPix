// $Id: TVirtualHealPainter.h,v 1.2 2008/07/09 11:50:10 oxon Exp $
// Author: Akira Okumura 2008/07/07

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_VIRTUAL_HEAL_PAINTER
#define T_VIRTUAL_HEAL_PAINTER

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualHealPainter                                                  //
//                                                                      //
// Abstract base class for HEALPix painters                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TF1;
class TClass;
class THealPix;

class TVirtualHealPainter : public TObject {
private:
  static TClass* fgPainter; //Pointer to class painter

public:
  TVirtualHealPainter() {}
  ~TVirtualHealPainter() {}

  virtual Int_t  DistancetoPrimitive(Int_t px, Int_t py) = 0;
  virtual void   DrawPanel() = 0;
  virtual void   ExecuteEvent(Int_t event, Int_t px, Int_t py) = 0;
  virtual TList* GetContourList(Double_t contour) const = 0;
  virtual Bool_t IsInside(Int_t x, Int_t y) = 0;
  virtual Bool_t IsInside(Double_t x, Double_t y) = 0;
  virtual void   Paint(Option_t* option = "") = 0;
  virtual void   PaintStat(Int_t dostat, TF1* fit) = 0;
  virtual void   ProcessMessage(const char* mess, const TObject* obj) = 0;
  virtual void   SetHealPix(THealPix* hp) = 0;
  virtual Int_t  MakeCuts(char* cutsopt) = 0;

  static TVirtualHealPainter* HealPainter(THealPix* obj);
  static void                 SetPainter(const char* painter);

  ClassDef(TVirtualHealPainter, 0);
};

#endif // T_VIRTUAL_HEAL_PAINTER
