// $Id: LinkDef.h,v 1.3 2008/06/25 06:03:18 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class THealPix+;
#pragma link C++ class std::vector<THealPix*>+;
#pragma link C++ class THealPixF+;
#pragma link C++ class THealPixD+;

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
