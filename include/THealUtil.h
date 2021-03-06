// $Id: THealUtil.h,v 1.11 2009/01/12 20:51:24 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
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

class THealPix;
class THealPixCube;

namespace THealUtil {
  Bool_t          FitsReportError(Int_t status);
  Bool_t          SaveToFits(const char* fname, const THealPix* hp);
  Bool_t          SaveToFits(const char* fname, const THealPixCube* hpc);
  Bool_t          SaveToFits(const char* fname, const std::vector<THealPix*>& hp);
  Bool_t          SaveToFits(const char* fname, const std::vector<const THealPix*>& hp);
};

#endif // T_HEAL_UTIL
