// $Id: TVirtualHealPainter.cxx,v 1.1 2008/07/08 08:12:06 oxon Exp $
// Author: Akira Okumura 2008/07/07

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include "TClass.h"
#include "TROOT.h"
#include "TVirtualHealPainter.h"

TClass* TVirtualHealPainter::fgPainter = 0;

ClassImp(TVirtualHealPainter)

//______________________________________________________________________________
TVirtualHealPainter* TVirtualHealPainter::HealPainter(THealPix* obj)
{
  if(!fgPainter){
    TVirtualHealPainter::SetPainter("THealPainter");
    if(!fgPainter){
      return 0;
    } // if
  } // if
  
  //create an instance of the histogram painter
  TVirtualHealPainter* p = (TVirtualHealPainter*)fgPainter->New();
  if(p){
    p->SetHealPix(obj);
  } // if
  
  return p;
}

//______________________________________________________________________________
void TVirtualHealPainter::SetPainter(const char *painter)
{
  // Static function to set an alternative histogram painter.
  fgPainter = TClass::GetClass(painter);
}
