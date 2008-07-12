// $Id: THealPixCube.cxx,v 1.1 2008/07/12 09:14:50 oxon Exp $
// Author: Akira Okumura 2008/07/11

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THealPixCube                                                         //
//                                                                      //
// Cubic HEALPix Image class                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "THealPixCube.h"

ClassImp(THealPixCube)

//______________________________________________________________________________
THealPixCube::THealPixCube(): TNamed()
{
  fN       = 0;
  fHeals   = 0;
  fMaximum = -1111;
  fMinimum = -1111;
  fWaxis.SetName("waxis");
  fWaxis.SetParent(this);
}

//______________________________________________________________________________
THealPixCube::THealPixCube(const char* name, const char* title,
			   Int_t nbins, Double_t wlow, Double_t wup)
  : TNamed(name, title)
{
  fN = (nbins <= 0) ? 1 : nbins;
  fWaxis.Set(fN, wlow, wup);
  fHeals = new THealPix*[fN];
  for(Int_t i = 0; i < fN; i++){
    fHeals[i] = 0;
  } // i
  fMaximum = -1111;
  fMinimum = -1111;
}

//______________________________________________________________________________
THealPixCube::~THealPixCube()
{
  for(Int_t i = 0; i < fN; i++){
    SafeDelete(fHeals[i]);
  } // i
  delete [] fHeals;
  fHeals = 0;
}

//______________________________________________________________________________
THealPix* THealPixCube::operator[](Int_t n) const
{
  if(n < 0 || n >= fN){
    return 0;
  } // if

  return fHeals[n];
}

//______________________________________________________________________________
THealPixCubeF::THealPixCubeF(): THealPixCube()
{
}

//______________________________________________________________________________
THealPixCubeF::THealPixCubeF(const char* name, const char* title, Int_t order,
			     Int_t nbins, Double_t wlow, Double_t wup,
			     Bool_t nested)
  : THealPixCube(name, title, nbins, wlow, wup)
{
  for(Int_t i = 0; i < fN; i++){
    fHeals[i] = new THealPixF(Form("%s%d", name, i), title, order, nested);
  } // i
}

//______________________________________________________________________________
THealPixCubeF::~THealPixCubeF()
{
}

//______________________________________________________________________________
THealPixCubeD::THealPixCubeD(): THealPixCube()
{
}

//______________________________________________________________________________
THealPixCubeD::~THealPixCubeD()
{
}

//______________________________________________________________________________
THealPixCubeD::THealPixCubeD(const char* name, const char* title, Int_t order,
			     Int_t nbins, Double_t wlow, Double_t wup,
			     Bool_t nested)
  : THealPixCube(name, title, nbins, wlow, wup)
{
  for(Int_t i = 0; i < fN; i++){
    fHeals[i] = new THealPixD(Form("%s%d", name, i), title, order, nested);
  } // i
}
