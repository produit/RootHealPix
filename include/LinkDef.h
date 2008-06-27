// $Id: LinkDef.h,v 1.5 2008/06/27 18:35:55 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class THealAlm<Double_t>;
#pragma link C++ class THealAlm<Float_t>;

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

#pragma link C++ class std::complex<Float_t>;
#pragma link C++ class std::complex<Double_t>;

#pragma link C++ class THealPix+;
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

#pragma link C++ namespace THealUtil;
#pragma link C++ function THealUtil::FitsReportError(Int_t);
#pragma link C++ function THealUtil::Nside2Npix(Int_t);
#pragma link C++ function THealUtil::Order2Nside(Int_t);
#pragma link C++ function THealUtil::Isqrt(UInt_t);
#pragma link C++ function THealUtil::Modulo(Double_t, Double_t);
#pragma link C++ function THealUtil::Modulo(Int_t, Int_t);
#pragma link C++ function THealUtil::Modulo(Long_t, Long_t);
#pragma link C++ function THealUtil::SaveToFits(const char*, const std::vector<THealPix*>&);
#endif
