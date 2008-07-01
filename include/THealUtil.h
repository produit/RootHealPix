// $Id: THealUtil.h,v 1.4 2008/07/01 23:03:35 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_UTIL
#define T_HEAL_UTIL

#include <vector>
#include <cmath>

#include "Rtypes.h"
#include "TMath.h"

namespace THealUtil {
  inline Int_t    Nside2Npix(Int_t nside) { return 12*nside*nside;}
  inline Int_t    Order2Nside(Int_t order) { return 1<<order;}

  Bool_t          FitsReportError(Int_t status);
  inline UInt_t   Isqrt(UInt_t v);
  void            GetChunkInfo(Int_t nrings, Int_t& nchunks, Int_t& chunksize);
  inline Double_t Modulo(Double_t v1, Double_t v2);
  inline Int_t    Modulo(Int_t v1, Int_t v2);
  inline Long_t   Modulo(Long_t v1, Long_t v2);
  Bool_t          SaveToFits(const char* fname, const THealPix* hp);
  Bool_t          SaveToFits(const char* fname, const std::vector<THealPix*>& hp);
  Bool_t          SaveToFits(const char* fname, const std::vector<const THealPix*>& hp);
};

//_____________________________________________________________________________
UInt_t THealUtil::Isqrt(UInt_t v)
{
  // Original code is isqrt of HEALPix C++
  return UInt_t(TMath::Sqrt(v + .5));
}

//_____________________________________________________________________________
Double_t THealUtil::Modulo(Double_t v1, Double_t v2)
{
  // Original code is modulo of HEALPix C++
  return (v1 >= 0) ? ((v1 < v2) ? v1 : fmod(v1, v2)) : (fmod(v1, v2) + v2);
}

//_____________________________________________________________________________
Int_t THealUtil::Modulo(Int_t v1, Int_t v2)
{
  // Original code is modulo of HEALPix C++
  return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
}

//_____________________________________________________________________________
Long_t THealUtil::Modulo(Long_t v1, Long_t v2)
{
  // Original code is modulo of HEALPix C++
  return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
}

#endif // T_HEAL_UTIL
