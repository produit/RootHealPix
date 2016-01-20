// $Id: THealPaletteAxis.h,v 1.1 2008/07/09 11:50:10 oxon Exp $
// Author: Akira Okumura 2008/07/09

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_PALETTE_AXIS
#define T_HEAL_PALETTE_AXIS


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THealPaletteAxis                                                     //
//                                                                      //
// class used to display a color palette axis for 2-d plots             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPave.h"
#include "TGaxis.h"

class THealPix;

class THealPaletteAxis : public TPave {

protected:
   TGaxis    fAxis; //palette axis
   THealPix* fHeal; //pointer to parent histogram
   TString   fName; //Pave name

public:
   // THealPaletteAxis status bits
   enum { kHasView   = BIT(11)};

   THealPaletteAxis();
   THealPaletteAxis(Double_t x1, Double_t y1, Double_t x2 ,Double_t y2, THealPix* hp);
   THealPaletteAxis(const THealPaletteAxis& palette);
   virtual ~THealPaletteAxis();
   void          Copy(TObject &palette) const;
   virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);
   virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
   TGaxis*       GetAxis() {return &fAxis;}
   Option_t*     GetName() const {return fName.Data();}
   virtual char* GetObjectInfo(Int_t px, Int_t py) const;
   virtual void  Paint(Option_t* option = "");
   virtual void  SavePrimitive(std::ostream& out, Option_t* option = "");
   virtual void  SetName(const char* name = "") {fName = name;} // *MENU*
   virtual void  SetLabelColor(Int_t labelcolor) {fAxis.SetLabelColor(labelcolor);} // *MENU*
   virtual void  SetLabelFont(Int_t labelfont) {fAxis.SetLabelFont(labelfont);} // *MENU*
   virtual void  SetLabelOffset(Float_t labeloffset) {fAxis.SetLabelOffset(labeloffset);} // *MENU*
   virtual void  SetLabelSize(Float_t labelsize) {fAxis.SetLabelSize(labelsize);} // *MENU*
   virtual void  SetTitleOffset(Float_t titleoffset = 1) {fAxis.SetTitleOffset(titleoffset);} // *MENU*
   virtual void  SetTitleSize(Float_t titlesize) {fAxis.SetTitleSize(titlesize);} // *MENU*
   virtual void  SetLineColor(Color_t linecolor) {fAxis.SetLineColor(linecolor);} // *MENU*
   virtual void  SetLineWidth(Width_t linewidth) {fAxis.SetLineWidth(linewidth);} // *MENU*
   virtual void  UnZoom();  // *MENU*

   ClassDef(THealPaletteAxis, 1)  //class used to display a color palette axis for HEALPix
};

#endif // T_HEAL_PALETTE_AXIS
