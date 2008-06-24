// $Id: THealPix.cxx,v 1.1 2008/06/24 08:16:44 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include <cstring>
#include <iostream>
#include "fitsio.h"

#include "TArrayD.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"

#include "THealPix.h"
#include "THealUtil.h"

Bool_t THealPix::fgAddDirectory = kTRUE;
THealPix::THealTable THealPix::fgTable = THealPix::THealTable();

THealPix::THealTable::THealTable()
{
  for(Int_t i = 0; i < 0x100; i++){
    fCtab[i] = (i&0x1) | ((i&0x2 ) << 7) | ((i&0x4 ) >> 1) | ((i&0x8 ) << 6)
      | ((i&0x10) >> 2) | ((i&0x20) << 5) | ((i&0x40) >> 3) | ((i&0x80) << 4);
    fUtab[i] = (i&0x1) | ((i&0x2 ) << 1) | ((i&0x4 ) << 2) | ((i&0x8 ) << 3)
      | ((i&0x10) << 4) | ((i&0x20) << 5) | ((i&0x40) << 6) | ((i&0x80) << 7);
  } // i
}

ClassImp(THealPix)

//_____________________________________________________________________________
THealPix::THealPix()
: TNamed(),
  fOrder(0),
  fNside(1),
  fNpix(12),
  fIsNested(kFALSE),
  fUnit(""),
  fDirectory(0)
{
}

//_____________________________________________________________________________
THealPix::~THealPix()
{
  if(!TestBit(kNotDeleted)){
    return;
  } // if

  if(fDirectory){
    fDirectory->Remove(this);
    fDirectory = 0;
  } // if
}

//_____________________________________________________________________________
THealPix::THealPix(const char* name, const char* title, Int_t order,
		   Bool_t isNested)
  : TNamed(name, title), fIsNested(isNested), fUnit("")
{
  Build();
  fOrder = order < 0 ? 0 : order;
  fNside = THealUtil::OrderToNside(fOrder);
  fNpix  = THealUtil::NsideToNpix(fNside);
}

//_____________________________________________________________________________
THealPix::THealPix(const THealPix& hp) : TNamed()
{
  ((THealPix&)hp).Copy(*this);
}

//_____________________________________________________________________________
void THealPix::AddDirectory(Bool_t add)
{
  fgAddDirectory = add;
}

//______________________________________________________________________________
void THealPix::AddBinContent(Int_t)
{
  AbstractMethod("AddBinContent");
}

//______________________________________________________________________________
void THealPix::AddBinContent(Int_t, Double_t)
{
  AbstractMethod("AddBinContent");
}

//_____________________________________________________________________________
Bool_t THealPix::AddDirectoryStatus()
{
  return fgAddDirectory;
}

//_____________________________________________________________________________
Int_t THealPix::AngToPix(Double_t theta, Double_t phi) const
{
  Double_t z = TMath::Cos(theta);
  Double_t zabs = TMath::Abs(z);
  Double_t z0 = 2./3.;
  // Conversion [0, 2pi) => [0, 4)
  Double_t tt = THealUtil::Modulo(phi, TMath::TwoPi())/TMath::PiOver2();
  Int_t ncap = 2*fNside*(fNside - 1);

  if(!fIsNested){ // RING
    if(zabs <= z0){ // Equatorial region
      Double_t temp1 = fNside*(0.5 + tt);
      Double_t temp2 = fNside*z*0.75;
      Int_t jp = Int_t(temp1 - temp2); // index of  ascending edge line
      Int_t jm = Int_t(temp1 + temp2); // index of descending edge line
      
      // ring number counted from z = 2/3
      Int_t ir = fNside + 1 + jp - jm; // in {1,2n+1}
      Int_t kshift = 1 - (ir&1); // kshift=1 if ir even, 0 otherwise
      
      Int_t ip = (jp + jm -fNside + kshift + 1)/2; // in {0,4n-1}
      ip = THealUtil::Modulo(ip, 4*fNside);
      
      return ncap + (ir - 1)*4*fNside + ip;
    } else { // North & South polar caps
      Double_t tp = tt - Int_t(tt);
      Double_t tmp = fNside*TMath::Sqrt(3*(1 - zabs));
      
      Int_t jp = Int_t(tp*tmp); // increasing edge line index
      Int_t jm = Int_t((1. - tp)*tmp); // decreasing edge line index
      
      Int_t ir = jp + jm + 1; // ring number counted from the closest pole
      Int_t ip = Int_t(tt*ir); // in {0,4*ir-1}
      ip = THealUtil::Modulo(ip, 4*ir);
      
      if(z>0){
	return 2*ir*(ir - 1) + ip;
      } else {
	return fNpix - 2*ir*(ir + 1) + ip;
      } // if
    } // if
  } else { // scheme_ == NEST
    const Int_t order_max = 13;
    const Int_t ns_max = 1<<order_max;
    Int_t face_num, ix, iy;
    
    if(zabs <= z0){ // Equatorial region
      Double_t temp1 = ns_max*(0.5 + tt);
      Double_t temp2 = ns_max*z*0.75;
      Int_t jp = Int_t(temp1 - temp2); // index of  ascending edge line
      Int_t jm = Int_t(temp1 + temp2); // index of descending edge line
      Int_t ifp = jp >> order_max;  // in {0,4}
      Int_t ifm = jm >> order_max;
      if(ifp == ifm){           // faces 4 to 7
	face_num = (ifp == 4) ? 4: ifp + 4;
      } else if(ifp < ifm){       // (half-)faces 0 to 3
	face_num = ifp;
      } else {                     // (half-)faces 8 to 11
	face_num = ifm + 8;
      } // if
      
      ix = jm & (ns_max-1);
      iy = ns_max - (jp & (ns_max-1)) - 1;
    } else { // polar region, zabs > 2/3
      Int_t ntt = Int_t(tt);
      Double_t tp = tt - ntt;
      Double_t tmp = ns_max*TMath::Sqrt(3.*(1 - zabs));
      
      Int_t jp = Int_t(tp*tmp); // increasing edge line index
      Int_t jm = Int_t((1. - tp)*tmp); // decreasing edge line index
      if(jp >= ns_max){
	jp = ns_max - 1; // for points too close to the boundary
      } // if
      if(jm >= ns_max){
	jm = ns_max - 1;
      } // if
      if(z >= 0){
	face_num = ntt;  // in {0,3}
	ix = ns_max - jm - 1;
	iy = ns_max - jp - 1;
      } else {
	face_num = ntt + 8; // in {8,11}
	ix =  jp;
	iy =  jm;
      } // if
    } // if
    
    Int_t ipf = XYToPix(ix, iy);
    
    ipf >>= (2*(order_max - fOrder));  // in {0, fNside**2 - 1}
    
    return ipf + (face_num << (2*fOrder));    // in {0, 12*fNside**2 - 1}
  } // if
}

//_____________________________________________________________________________
void THealPix::Build()
{
  fDirectory  = 0;

  SetTitle(fTitle.Data());
  
  if(THealPix::AddDirectoryStatus()){
    fDirectory = gDirectory;
    if(fDirectory){
      fDirectory->Append(this, kTRUE);
    } // if
  } // if
}

//_____________________________________________________________________________
void THealPix::Copy(TObject& newhp) const
{
  THealPix::Copy(newhp);
}

//______________________________________________________________________________
Int_t THealPix::Fill(Double_t theta, Double_t phi)
{
  /*
  if(fBuffer){
    return BufferFill(x, y, 1);
  } // if
  */

  if(fIsDegree){
    theta *= TMath::DegToRad();
    phi   *= TMath::DegToRad();
  } // if

  Int_t bin = AngToPix(theta, phi);
  AddBinContent(bin);

  return bin;
}

//______________________________________________________________________________
Int_t THealPix::Fill(Double_t theta, Double_t phi, Double_t w)
{
  /*
  if(fBuffer){
    return BufferFill(x, y, w);
  } // if
  */
  if(fIsDegree){
    theta *= TMath::DegToRad();
    phi   *= TMath::DegToRad();
  } // if

  Int_t bin = AngToPix(theta, phi);
  AddBinContent(bin, w);

  return bin;
}

//______________________________________________________________________________
Double_t THealPix::GetBinContent(Int_t) const
{
  AbstractMethod("GetBinContent");
  return 0;
}

//_____________________________________________________________________________
void THealPix::Draw(Option_t* option)
{
  TString opt = option;
  opt.ToLower();
  if(gPad){
    if(!gPad->IsEditable()){
      gROOT->MakeDefCanvas();
    } // if
  } // if

  AppendPad(option);
}

//_____________________________________________________________________________
void THealPix::Paint(Option_t* /*option*/)
{
  // Not implemented yet
}

//_____________________________________________________________________________
Bool_t THealPix::ReadFitsHeader(fitsfile** fptr, const char* fname,
				const char* colname,
				HealHeader_t& head)
{
  Int_t status = 0;

  fits_open_file(fptr, fname, READONLY, &status);
  if(!THealUtil::FitsReportError(status)){
    return kFALSE;
  } // if

  Int_t hdutype;
  fits_movabs_hdu(*fptr, 2, &hdutype, &status);
  if(!THealUtil::FitsReportError(status) or hdutype != BINARY_TBL){
    return kFALSE;
  } // if

  Long_t nrows;
  Int_t ncols;
  fits_get_num_rows(*fptr, &nrows, &status);
  fits_get_num_cols(*fptr, &ncols, &status);

  char comment[FLEN_COMMENT];
  char pixtype[FLEN_VALUE];
  char ordering[FLEN_VALUE];

  fits_read_key(*fptr, TSTRING, "PIXTYPE", pixtype, comment, &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(*fptr, &status);
    return kFALSE;
  } // if
  
  fits_read_key(*fptr, TSTRING, "ORDERING", ordering, comment, &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(*fptr, &status);
    return kFALSE;
  } // if

  if(strcmp(pixtype, "HEALPIX") != 0){
    fits_close_file(*fptr, &status);
    return kFALSE;
  } // if
  
  Bool_t isNested;
  if(strcmp(ordering, "RING") == 0){
    isNested = kFALSE;
  } else if(strcmp(ordering, "NESTED") == 0){
    isNested = kTRUE;
  } else {
    fits_close_file(*fptr, &status);
    std::cerr << Form("Unknown ordering keyword %s found.\n", ordering);
    return kFALSE;
  } // if
  
  Int_t colnum = 0;
  fits_get_colnum(*fptr, CASEINSEN, (char*)colname, &colnum, &status);
  
  if(colnum == 0){
    std::cerr << Form("Column %s not found.\n", colname);
    return kFALSE;
  } // if

  char ttype[FLEN_VALUE], typechar[FLEN_VALUE],
    tdisp[FLEN_VALUE];
  Long_t repeat, nulval;
  Double_t scale, zero;
    
  fits_get_bcolparms(*fptr, colnum, ttype, head.tunit, typechar, &repeat,
		     &scale, &zero, &nulval, tdisp, &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(*fptr, &status);
    return kFALSE;
  } // if

  Long_t nside;
  fits_read_key_lng(*fptr, "NSIDE", &nside, comment, &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(*fptr, &status);
    return kFALSE;
  } // if

  // Check if nside == 2^0 or 2^2 ... or 2^13 (8192)
  Int_t order = 0;
  for(; order < 15; order++){
    if(nside == 1<<order){
      break;
    } // if
    if(order == 14){
      fits_close_file(*fptr, &status);
      std::cerr << Form("Invalid Nside %d was given.\n", nside);
      return kFALSE;
    } // if
  } // i

  Int_t npix = THealUtil::NsideToNpix(nside);

  if(npix%nrows != 0){
    fits_close_file(*fptr, &status);
    std::cerr << Form("Npix(%d) is not dividable by nrows(%d).\n", npix, nrows);
    return kFALSE;
  } // if

  head.isNested = isNested;
  head.order    = order;
  head.nside    = nside;
  head.npix     = npix;
  head.nrows    = nrows;
  head.colnum   = colnum;
  
  return kTRUE;
}

//_____________________________________________________________________________
void THealPix::SetUnit(const char* unit)
{
  fUnit = std::string(unit);
}

//_____________________________________________________________________________
Int_t THealPix::GetNrows() const
{
  return fOrder < 4 ? 1 : 1024;
}

std::string THealPix::GetOrderingTypeString() const
{
  if(fIsNested){
    return "NESTED";
  } else {
    return "RING";
  } // if
}

//_____________________________________________________________________________
std::string THealPix::GetTypeString() const
{
  if(fType == TDOUBLE){
    return "D";
  } // if

  return "D";
}

//_____________________________________________________________________________
Int_t THealPix::XYToPix(Int_t x, Int_t y) const
{
  return fgTable.U(x&0xff) | (fgTable.U(x>>8)<<16) | (fgTable.U(y&0xff)<<1)
    | (fgTable.U(y>>8)<<17);
}

ClassImp(THealPixD)

//_____________________________________________________________________________
// THealPixD methods
//_____________________________________________________________________________
THealPixD::THealPixD(): THealPix(), TArrayD()
{
}

//_____________________________________________________________________________
THealPixD::THealPixD(const char* name, const char* title, Int_t order,
		     Bool_t isNested)
: THealPix(name, title, order, isNested)
{
  fType = TDOUBLE;
  TArrayD::Set(fNpix);
}

//_____________________________________________________________________________
THealPixD::THealPixD(const THealPixD& hpd) : THealPix(), TArrayD()
{
  ((THealPixD&)hpd).Copy(*this);
}

//_____________________________________________________________________________
THealPixD::~THealPixD()
{
}

//_____________________________________________________________________________
void THealPixD::Copy(TObject& newhp) const
{
  THealPix::Copy(newhp);
}

//______________________________________________________________________________
Double_t THealPixD::GetBinContent(Int_t bin) const
{
  if(!fArray){
    return 0;
  } // if

  bin = (bin < 0 ? 0 : bin) < fNpix ? bin : fNpix - 1;

  return Double_t(fArray[bin]);
}

//______________________________________________________________________________
THealPixD* THealPixD::ReadFits(const char* fname, const char* colname)
{
  THealPix::HealHeader_t head;
  fitsfile* fptr = 0;

  if(!THealPix::ReadFitsHeader(&fptr, fname, colname, head)){
    return 0;
  } // if

  THealPixD* hpd = new THealPixD(fname, colname, head.order, head.isNested);
  hpd->SetUnit(head.tunit);

  Long_t npercol = head.npix/head.nrows;
  Int_t status = 0;

  for(Int_t i = 0; i < head.nrows; i++){
    fits_read_col(fptr, TDOUBLE, head.colnum, i+1, 1, npercol, 0,
		  &(hpd->GetArray()[i*npercol]), 0, &status);
    if(!THealUtil::FitsReportError(status)){
      delete hpd;
      return 0;
    } // if
  } // i

  return hpd;
}
