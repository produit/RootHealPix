// $Id: THealUtil.h,v 1.1 2008/06/24 08:16:43 oxon Exp $
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

namespace THealUtil {
  inline Int_t NsideToNpix(Int_t nside) { return 12*nside*nside;}
  inline Int_t OrderToNside(Int_t order) { return 1<<order;}

  Bool_t       FitsReportError(Int_t status);
  inline Double_t Modulo(Double_t v1, Double_t v2);
  inline Int_t    Modulo(Int_t v1, Int_t v2);
  inline Long_t   Modulo(Long_t v1, Long_t v2);
  Bool_t       SaveToFits(const char* fname, const std::vector<THealPix*>& hp);
};

//_____________________________________________________________________________
Double_t THealUtil::Modulo(Double_t v1, Double_t v2)
{
  return (v1 >= 0) ? ((v1 < v2) ? v1 : fmod(v1, v2)) : (fmod(v1, v2) + v2);
}

//_____________________________________________________________________________
Int_t THealUtil::Modulo(Int_t v1, Int_t v2)
{
  return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
}

//_____________________________________________________________________________
Long_t THealUtil::Modulo(Long_t v1, Long_t v2)
{
  return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
}

#endif // T_HEAL_UTIL
