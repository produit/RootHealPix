// $Id: THealUtil.cxx,v 1.5 2008/07/01 23:03:35 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include <cmath>
#include <vector>

#include <fitsio.h>

#include "TMath.h"

#include "THealPix.h"
#include "THealUtil.h"

namespace THealUtil {

Bool_t FitsReportError(Int_t status)
{
  if(status != 0){
      fits_report_error(stderr, status);
      return kFALSE;
  } // if
  
  return kTRUE;
}
  
//______________________________________________________________________________
void GetChunkInfo(Int_t nrings, Int_t& nchunks, Int_t& chunksize)
{
  // Original code is get_chunk_info of HEALPix C++
  nchunks = nrings/TMath::Max(100, nrings/10) + 1;
  chunksize = (nrings + nchunks - 1)/nchunks;
}

//______________________________________________________________________________
Bool_t SaveToFits(const char* fname, const THealPix* hp)
{
  std::vector<THealPix*> vec;
  vec.push_back((THealPix*)hp);

  return SaveToFits(fname, vec);
}

//______________________________________________________________________________
Bool_t SaveToFits(const char* fname, const std::vector<THealPix*>& hp)
{
  std::vector<const THealPix*> vec;
  for(UInt_t i = 0; i < vec.size(); i++){
    vec.push_back(hp[i]);
  } // i

  return SaveToFits(fname, vec);
}

//______________________________________________________________________________
Bool_t SaveToFits(const char* fname, const std::vector<const THealPix*>& hp)
{
  fitsfile* fptr;
  Int_t status = 0;
  
  fits_create_file(&fptr, fname, &status);
  if(!THealUtil::FitsReportError(status)){
    return kFALSE;
  } // if

  char** type = new char*[hp.size()];
  char** unit = new char*[hp.size()];
  char** form = new char*[hp.size()];
  for(UInt_t i = 0; i < hp.size(); i++){
    type[i] = new char[81];
    unit[i] = new char[81];
    form[i] = new char[81];
  } // i

  for(UInt_t i = 0; i < hp.size(); i++){
    strcpy(type[i], hp[i]->GetTitle());
    strcpy(unit[i], hp[i]->GetUnit().c_str());
    strcpy(form[i], Form("%d%s", hp[i]->GetNrows(), hp[i]->GetTypeString().c_str()));
  } // i

  fits_insert_btbl(fptr, 0, hp.size(), (char**)type, (char**)form, (char**)unit, "xtension", 0, &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(fptr, &status);
    return kFALSE;
  } // if

  fits_write_key(fptr, TSTRING, (char*)"PIXTYPE", (char*)"HEALPIX", (char*)"HEALPIX pixelisation", &status);
  fits_write_key(fptr, TSTRING, (char*)"ORDERING", (char*)hp[0]->GetSchemeString().c_str(), (char*)"Pixel ordering scheme, either RING or NESTED", &status);
  Long_t nside = hp[0]->GetNside();
  fits_write_key(fptr, TINT, (char*)"NSIDE", &nside, (char*)"Resolution parameter for HEALPIX", &status);
  Long_t first = 0;
  fits_write_key(fptr, TLONG, (char*)"FIRSTPIX", &first, (char*)"First pixel # (0 based)", &status);
  Long_t last = hp[0]->GetNpix() - 1;
  fits_write_key(fptr, TLONG, (char*)"LASTPIX", &last, (char*)"Last pixel # (0 based)", &status);
  fits_write_key(fptr, TSTRING, (char*)"INDXSCHM", (char*)"IMPLICIT", (char*)"Indexing: IMPLICIT or EXPLICIT", &status);
  Long_t grain = 0;
  fits_write_key(fptr, TLONG, (char*)"GRAIN", &grain, (char*)"Grain of pixel indexing", &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(fptr, &status);
    return kFALSE;
  } // if

  for(UInt_t i = 0; i < hp.size(); i++){
    if(hp[i]->GetType() == TDOUBLE){
      fits_write_col(fptr, hp[i]->GetType(), i + 1, 1, 1, hp[i]->GetNpix(),
		     ((THealPixD*)hp[i])->GetArray(), &status);
      Long_t nrows;
      fits_get_num_rows(fptr, &nrows, &status);
    } // if
    if(!THealUtil::FitsReportError(status)){
      fits_close_file(fptr, &status);
      return kFALSE;
    } // if
  } // i

  fits_close_file(fptr, &status);

  return kTRUE;
}

} // namespace THealUtil
