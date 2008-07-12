// $Id: LinkDef.h,v 1.16 2008/07/12 09:14:49 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class THealAlm<Float_t>-;
#pragma link C++ class THealAlm<Double_t>-;

#pragma link C++ function THealAlm<Float_t>::Alm2Map(THealPix&);
#pragma link C++ function THealAlm<Double_t>::Alm2Map(THealPix&);

#pragma link C++ function operator*(Double_t, THealAlm<Float_t>&);
#pragma link C++ function operator*(Double_t, THealAlm<Double_t>&);
#pragma link C++ function operator*(THealAlm<Float_t>&,  Double_t);
#pragma link C++ function operator*(THealAlm<Double_t>&, Double_t);
#pragma link C++ function operator+(THealAlm<Float_t>&,  THealAlm<Float_t>&);
#pragma link C++ function operator+(THealAlm<Double_t>&, THealAlm<Double_t>&);
#pragma link C++ function operator-(THealAlm<Float_t>&,  THealAlm<Float_t>&);
#pragma link C++ function operator-(THealAlm<Double_t>&, THealAlm<Double_t>&);
#pragma link C++ function operator*(THealAlm<Float_t>&,  THealAlm<Float_t>&);
#pragma link C++ function operator*(THealAlm<Double_t>&, THealAlm<Double_t>&);
#pragma link C++ function operator/(THealAlm<Float_t>&,  THealAlm<Float_t>&);
#pragma link C++ function operator/(THealAlm<Double_t>&, THealAlm<Double_t>&);

#pragma link C++ class THealPix-;
#pragma link C++ class std::vector<THealPix*>+;
#pragma link C++ class THealPixF+;
#pragma link C++ class THealPixD+;

#pragma link C++ function THealPix::Map2Alm(THealAlm<Float_t>&, Bool_t);
#pragma link C++ function THealPix::Map2Alm(THealAlm<Double_t>&, Bool_t);

void THealPix::Map2Alm(THealAlm<Float_t>& alm, const
		       std::vector<Double_t>& weight, Bool_t add = kFALSE);
void THealPix::Map2Alm(THealAlm<Double_t>& alm, const
		       std::vector<Double_t>& weight, Bool_t add = kFALSE);

#pragma link C++ function THealPix::Map2Alm(THealAlm<Float_t>&, std::vector<Double_t>&, Bool_t);
#pragma link C++ function THealPix::Map2Alm(THealAlm<Double_t>&, std::vector<Double_t>&, Bool_t);

#pragma link C++ function operator*(Double_t, THealPixF&);
#pragma link C++ function operator*(THealPixF&, Double_t);
#pragma link C++ function operator+(THealPixF&, THealPixF&);
#pragma link C++ function operator-(THealPixF&, THealPixF&);
#pragma link C++ function operator*(THealPixF&, THealPixF&);
#pragma link C++ function operator/(THealPixF&, THealPixF&);

#pragma link C++ function operator*(Double_t, THealPixD&);
#pragma link C++ function operator*(THealPixD&, Double_t);
#pragma link C++ function operator+(THealPixD&, THealPixD&);
#pragma link C++ function operator-(THealPixD&, THealPixD&);
#pragma link C++ function operator*(THealPixD&, THealPixD&);
#pragma link C++ function operator/(THealPixD&, THealPixD&);

#pragma link C++ class THealPixCube+;
#pragma link C++ class THealPixCubeF+;
#pragma link C++ class THealPixCubeD+;

#pragma link C++ namespace THealUtil;
#pragma link C++ function THealUtil::FitsReportError(Int_t);
#pragma link C++ function THealUtil::SaveToFits(const char*, const std::vector<THealPix*>&);

#pragma link C++ class THealPainter+;
#pragma link C++ class THealPaletteAxis+;
#pragma link C++ class TVirtualHealPainter+;

#endif
