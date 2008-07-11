// $Id: THealUtil.h,v 1.8 2008/07/11 23:57:48 oxon Exp $
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

class THealPix;

namespace THealUtil {
  Bool_t          FitsReportError(Int_t status);
  Bool_t          SaveToFits(const char* fname, const THealPix* hp);
  Bool_t          SaveToFits(const char* fname, const std::vector<THealPix*>& hp);
  Bool_t          SaveToFits(const char* fname, const std::vector<const THealPix*>& hp);
};

#endif // T_HEAL_UTIL
