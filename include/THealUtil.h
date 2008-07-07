// $Id: THealUtil.h,v 1.6 2008/07/07 07:19:05 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.

   This is a port of HEALPix C++ package to ROOT system.
   Original code is available at <http://healpix.jpl.nasa.gov> under GPL.
******************************************************************************/

#ifndef T_HEAL_UTIL
#define T_HEAL_UTIL

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THealUtil                                                            //
//                                                                      //
// Encapsulate HEALPix utility functions.                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>

#include "Rtypes.h"
#include "TMath.h"

template<typename T> class THealAlm;
class THealPix;

namespace THealUtil {
  template<typename T> void Alm2Map(const THealAlm<T>& alm, THealPix& map);
  template<typename T> void Map2Alm(const THealPix& map, THealAlm<T>& alm, Bool_t add_alm = kFALSE);
  template<typename T> void Map2Alm(const THealPix& map, THealAlm<T>& alm, const std::vector<double>& weight, Bool_t add_alm = kFALSE);

  inline Int_t    Nside2Npix(Int_t nside) { return 12*nside*nside;}
  inline Int_t    Order2Nside(Int_t order) { return 1<<order;}

  Bool_t          FitsReportError(Int_t status);
  inline UInt_t   Isqrt(UInt_t v);
  inline Double_t Modulo(Double_t v1, Double_t v2);
  inline Int_t    Modulo(Int_t v1, Int_t v2);
  inline Long_t   Modulo(Long_t v1, Long_t v2);
  Bool_t          SaveToFits(const char* fname, const THealPix* hp);
  Bool_t          SaveToFits(const char* fname, const std::vector<THealPix*>& hp);
  Bool_t          SaveToFits(const char* fname, const std::vector<const THealPix*>& hp);
};

//______________________________________________________________________________
UInt_t THealUtil::Isqrt(UInt_t v)
{
  // Calculate square root of an integer
  // Original code is isqrt of HEALPix C++
  return UInt_t(TMath::Sqrt(v + .5));
}

//______________________________________________________________________________
Double_t THealUtil::Modulo(Double_t v1, Double_t v2)
{
  // Calculate modulo of two doubles
  // Original code is modulo of HEALPix C++
  return (v1 >= 0) ? ((v1 < v2) ? v1 : fmod(v1, v2)) : (fmod(v1, v2) + v2);
}

//______________________________________________________________________________
Int_t THealUtil::Modulo(Int_t v1, Int_t v2)
{
  // Calculate modulo of two integers
  // Original code is modulo of HEALPix C++
  return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
}

//______________________________________________________________________________
Long_t THealUtil::Modulo(Long_t v1, Long_t v2)
{
  // Calculate modulo of two long integers
  // Original code is modulo of HEALPix C++
  return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
}

#endif // T_HEAL_UTIL
