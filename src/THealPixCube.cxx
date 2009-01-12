// $Id: THealPixCube.cxx,v 1.5 2009/01/12 20:51:25 oxon Exp $
// Author: Akira Okumura 2008/07/11

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THealPixCube                                                         //
//                                                                      //
// Cubic HEALPix Image class                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "THealPixCube.h"
#include "THealUtil.h"

ClassImp(THealPixCube)

//______________________________________________________________________________
THealPixCube::THealPixCube(): TNamed()
{
  fN       = 0;
  fHeals   = 0;
  fMaximum = -1111;
  fMinimum = -1111;
  fWaxis.SetName("waxis");
  fWaxis.SetParent(this);
}

//______________________________________________________________________________
THealPixCube::THealPixCube(const char* name, const char* title,
			   Int_t nbins, Double_t wlow, Double_t wup)
  : TNamed(name, title)
{
  fN = (nbins <= 0) ? 1 : nbins;
  fWaxis.Set(fN, wlow, wup);
  fHeals = new THealPix*[fN];
  for(Int_t i = 0; i < fN; i++){
    fHeals[i] = 0;
  } // i
  fMaximum = -1111;
  fMinimum = -1111;
}

//______________________________________________________________________________
THealPixCube::~THealPixCube()
{
  for(Int_t i = 0; i < fN; i++){
    SafeDelete(fHeals[i]);
  } // i
  delete [] fHeals;
  fHeals = 0;
}

//______________________________________________________________________________
THealPix* THealPixCube::GetSlice(Int_t n) const
{
  if(n < 0 || n >= fN){
    return 0;
  } // if
  
  return fHeals[n];
}

//______________________________________________________________________________
THealPix* THealPixCube::operator[](Int_t n) const
{
  if(n < 0 || n >= fN){
    return 0;
  } // if

  return fHeals[n];
}

ClassImp(THealPixCubeF)

//______________________________________________________________________________
THealPixCubeF::THealPixCubeF(): THealPixCube()
{
}

//______________________________________________________________________________
THealPixCubeF::~THealPixCubeF()
{
}

//______________________________________________________________________________
THealPixCubeF::THealPixCubeF(const char* name, const char* title, Int_t order,
			     Int_t nbins, Double_t wlow, Double_t wup,
			     Bool_t nested)
  : THealPixCube(name, title, nbins, wlow, wup)
{
  for(Int_t i = 0; i < fN; i++){
    fHeals[i] = new THealPixF(Form("%s%d", name, i), title, order, nested);
  } // i
}

//______________________________________________________________________________
Int_t THealPixCubeF::GetType() const
{
  // Return type of content as CFITSIO type 'TFLOAT'
  return TFLOAT;
}

//_____________________________________________________________________________
std::string THealPixCubeF::GetTypeString() const
{
  // Return type of content as FITS type string 'E'
  return "E";
}

//______________________________________________________________________________
THealPixCubeF* THealPixCubeF::ReadFits(const char* fname, const char* colname)
{
  THealPix::HealHeader_t head;
  fitsfile* fptr = 0;

  if(!THealPix::ReadFitsHeader(&fptr, fname, colname, head)){
    return 0;
  } // if

  THealPixCubeF* hpcf;
  if(head.nrows*head.repeat == head.npix){ // CMB-like FITS
    hpcf = new THealPixCubeF(fname, colname, head.order, 1, 0, 0, head.isNested);
  } else if(head.nrows == head.npix){ // HEALPix Cube
    hpcf = new THealPixCubeF(fname, colname, head.order, head.repeat, 0, 0, head.isNested);
  } else {
    return 0;
  } // if

  for(Int_t i = 0; i < hpcf->GetN(); i++){
    (*hpcf)[i]->SetUnit(head.tunit);
  } // i

  Int_t status = 0;

  if(head.nrows*head.repeat == head.npix){ // CMB-like FITS
    for(Int_t i = 0; i < head.nrows; i++){
      fits_read_col(fptr, TFLOAT, head.colnum, i + 1, 1, head.repeat, 0,
		    &(((THealPixF*)((*hpcf)[0]))->GetArray()[i*head.repeat]),
		    0, &status);
      if(!THealUtil::FitsReportError(status)){
	delete hpcf;
	return 0;
      } // if
    } // i
  } else if(head.nrows == head.npix){ // HEALPix Cube
    for(Int_t i = 0; i < head.nrows; i++){
      for(Int_t j = 0; j < hpcf->GetN(); j++){
	fits_read_col(fptr, TFLOAT, head.colnum, i + 1, j + 1, 1, 0,
		      &(((THealPixF*)((*hpcf)[j]))->GetArray()[i]), 0, &status);
	if(!THealUtil::FitsReportError(status)){
	  delete hpcf;
	  return 0;
	} // if
      } // j
    } // i
  } // if

  for(Int_t i = 0; i < hpcf->GetN(); i++){
    Double_t entries = 0;
    for(Int_t j = 0; j < (*hpcf)[i]->GetNpix(); j++){
      if((*hpcf)[i]->GetBinContent(j) != 0){
	entries++;
      } // if
    } // j
    (*hpcf)[i]->SetEntries(entries);
  } // i

  return hpcf;
}

ClassImp(THealPixCubeD)

//______________________________________________________________________________
THealPixCubeD::THealPixCubeD(): THealPixCube()
{
}

//______________________________________________________________________________
THealPixCubeD::~THealPixCubeD()
{
}

//______________________________________________________________________________
THealPixCubeD::THealPixCubeD(const char* name, const char* title, Int_t order,
			     Int_t nbins, Double_t wlow, Double_t wup,
			     Bool_t nested)
  : THealPixCube(name, title, nbins, wlow, wup)
{
  for(Int_t i = 0; i < fN; i++){
    fHeals[i] = new THealPixD(Form("%s%d", name, i), title, order, nested);
  } // i
}

//______________________________________________________________________________
Int_t THealPixCubeD::GetType() const
{
  // Return type of content as CFITSIO type 'TDOUBLE'
  return TDOUBLE;
}

//_____________________________________________________________________________
std::string THealPixCubeD::GetTypeString() const
{
  // Return type of content as FITS type string 'D'
  return "D";
}

//______________________________________________________________________________
THealPixCubeD* THealPixCubeD::ReadFits(const char* fname, const char* colname)
{
  THealPix::HealHeader_t head;
  fitsfile* fptr = 0;

  if(!THealPix::ReadFitsHeader(&fptr, fname, colname, head)){
    return 0;
  } // if

  THealPixCubeD* hpcd;
  if(head.nrows*head.repeat == head.npix){ // CMB-like FITS
    hpcd = new THealPixCubeD(fname, colname, head.order, 1, 0, 0, head.isNested);
  } else if(head.nrows == head.npix){ // HEALPix Cube
    hpcd = new THealPixCubeD(fname, colname, head.order, head.repeat, 0, 0, head.isNested);
  } else {
    return 0;
  } // if

  for(Int_t i = 0; i < hpcd->GetN(); i++){
    (*hpcd)[i]->SetUnit(head.tunit);
  } // i

  Int_t status = 0;

  if(head.nrows*head.repeat == head.npix){ // CMB-like FITS
    for(Int_t i = 0; i < head.nrows; i++){
      fits_read_col(fptr, TDOUBLE, head.colnum, i + 1, 1, head.repeat, 0,
		    &(((THealPixD*)((*hpcd)[0]))->GetArray()[i*head.repeat]),
		    0, &status);
      if(!THealUtil::FitsReportError(status)){
	delete hpcd;
	return 0;
      } // if
    } // i
  } else if(head.nrows == head.npix){ // HEALPix Cube
    for(Int_t i = 0; i < head.nrows; i++){
      for(Int_t j = 0; j < hpcd->GetN(); j++){
	fits_read_col(fptr, TDOUBLE, head.colnum, i + 1, j + 1, 1, 0,
		      &(((THealPixD*)((*hpcd)[j]))->GetArray()[i]), 0, &status);
	if(!THealUtil::FitsReportError(status)){
	  delete hpcd;
	  return 0;
	} // if
      } // j
    } // i
  } // if

  for(Int_t i = 0; i < hpcd->GetN(); i++){
    Double_t entries = 0;
    for(Int_t j = 0; j < (*hpcd)[i]->GetNpix(); j++){
      if((*hpcd)[i]->GetBinContent(j) != 0){
	entries++;
      } // if
    } // j
    (*hpcd)[i]->SetEntries(entries);
  } // i

  return hpcd;
}
