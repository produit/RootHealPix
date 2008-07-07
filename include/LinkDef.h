// $Id: LinkDef.h,v 1.12 2008/07/07 07:04:03 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class THealAlm<Double_t>-;
#pragma link C++ class THealAlm<Float_t>-;

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

namespace THealUtil {
  // Work around a CINT instantiation problem ...
  void Map2Alm(const THealPix& map, THealAlm<Float_t>& alm, Bool_t add_alm = kFALSE);
  void Map2Alm(const THealPix& map, THealAlm<Double_t>& alm, Bool_t add_alm = kFALSE);
  void Map2Alm(const THealPix& map, THealAlm<Float_t>& alm, const
	       std::vector<Double_t>& weight, Bool_t add_alm = kFALSE);
  void Map2Alm(const THealPix& map, THealAlm<Double_t>& alm, const
	       std::vector<Double_t>& weight, Bool_t add_alm = kFALSE);
}
#pragma link C++ namespace THealUtil;
#pragma link C++ function THealUtil::Alm2Map(const THealAlm<Float_t>&, THealPix&);
#pragma link C++ function THealUtil::Alm2Map(const THealAlm<Double_t>&, THealPix&);
#pragma link C++ function THealUtil::Map2Alm(const THealPix&, THealAlm<Float_t>&, Bool_t);
#pragma link C++ function THealUtil::Map2Alm(const THealPix&, THealAlm<Double_t>&, Bool_t);
#pragma link C++ function THealUtil::Map2Alm(const THealPix&, THealAlm<Float_t>&, std::vector<Double_t>&, Bool_t);
#pragma link C++ function THealUtil::Map2Alm(const THealPix&, THealAlm<Double_t>&, std::vector<Double_t>&, Bool_t);

#pragma link C++ function THealUtil::FitsReportError(Int_t);
#pragma link C++ function THealUtil::Nside2Npix(Int_t);
#pragma link C++ function THealUtil::Order2Nside(Int_t);
#pragma link C++ function THealUtil::Isqrt(UInt_t);
#pragma link C++ function THealUtil::Modulo(Double_t, Double_t);
#pragma link C++ function THealUtil::Modulo(Int_t, Int_t);
#pragma link C++ function THealUtil::Modulo(Long_t, Long_t);
#pragma link C++ function THealUtil::SaveToFits(const char*, const std::vector<THealPix*>&);
#endif
