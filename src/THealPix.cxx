// $Id: THealPix.cxx,v 1.15 2008/06/30 19:28:10 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include <cstring>
#include <iostream>
#include "fitsio.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TVector3.h"

#include "THealPix.h"
#include "THealUtil.h"

Bool_t THealPix::fgAddDirectory = kTRUE;
Bool_t THealPix::fgDefaultSumw2 = kFALSE;
THealPix::THealTable THealPix::fgTable = THealPix::THealTable();
const Int_t THealPix::fgJrll[] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
const Int_t THealPix::fgJpll[] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};

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
  fIsNested(kFALSE),
  fUnit(""),
  fEntries(0),
  fTsumw(0),		 
  fTsumw2(0),						 
  fDirectory(0)
{
  SetOrder(0);
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

  if(order < 0)  order = 0;
  if(order > 13) order = 13;

  SetOrder(order);
  if(fgDefaultSumw2){
    Sumw2();
  } // if
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
void THealPix::Add(const THealPix* hp1, Double_t c1)
{
  if(!hp1){
    Error("Add", "Attempt to add a non-existing HEALPix");
    return;
  } // if
  
  // Check HEALPix compatibility
  if(fOrder != hp1->GetOrder() || fIsNested != hp1->IsNested()){
    Error("Add", "Attempt to add HEALPixs with different number of order or scheme");
    return;
  } // if

  // Create Sumw2 if hp1 has Sumw2 set
  if(fSumw2.fN == 0 && hp1->GetSumw2N() != 0){
    Sumw2();
  } // if

  // Add statistics
  fEntries += c1*hp1->GetEntries();
 
  Double_t factor = 1;
  //if (hp1->GetNormFactor() != 0) factor = hp1->GetNormFactor()/hp1->GetSumOfWeights();;
  for(Int_t bin = 0; bin < fNpix; bin++){
    //special case where histograms have the kIsAverage bit set
    if(this->TestBit(kIsAverage) && hp1->TestBit(kIsAverage)){
      Double_t v1 = hp1->GetBinContent(bin);
      Double_t v2 = this->GetBinContent(bin);
      Double_t e1 = hp1->GetBinError(bin);
      Double_t e2 = this->GetBinError(bin);
      Double_t w1 = 1., w2 = 1.;
      if(e1 > 0) w1 = 1./(e1*e1);
      if(e2 > 0) w2 = 1./(e2*e2);
      SetBinContent(bin, (w1*v1 + w2*v2)/(w1 + w2));
      if(fSumw2.fN){
	fSumw2.fArray[bin] = 1./(w1 + w2);
      } // if
    } else {
      Double_t cu = c1*factor*hp1->GetBinContent(bin);
      AddBinContent(bin, cu);
      if(fSumw2.fN){
	Double_t e1 = factor*hp1->GetBinError(bin);
	fSumw2.fArray[bin] += c1*c1*e1*e1;
      } // if
    } // if
  } // i

}

//______________________________________________________________________________
void THealPix::Add(const THealPix* hp1, const THealPix* hp2, Double_t c1, Double_t c2)
{
  if(!hp1 || !hp2){
    Error("Add", "Attempt to add a non-existing HEALPix");
    return;
  } // if

  Bool_t normWidth = kFALSE;
  if(hp1 == hp2 && c2 < 0){
    c2 = 0;
    normWidth = kTRUE;
  } // if
  // Check HEALPix compatibility
  if(fOrder != hp1->GetOrder() || fIsNested != hp1->IsNested()
  || fOrder != hp2->GetOrder() || fIsNested != hp2->IsNested()){
    Error("Add", "Attempt to add HEALPixs with different number of order or scheme");
    return;
  } // if

  // Create Sumw2 if hp1 or hp2 have Sumw2 set
  if(fSumw2.fN == 0 && (hp1->GetSumw2N() !=0 || hp2->GetSumw2N()  != 0)){
    Sumw2();
  } // if

  // Add statistics
  Double_t nEntries = c1*hp1->GetEntries() + c2*hp2->GetEntries();
  
  for(Int_t bin = 0; bin < fNpix; bin++){
    if(hp1->TestBit(kIsAverage) && hp2->TestBit(kIsAverage)){
      Double_t v1 = hp1->GetBinContent(bin);
      Double_t v2 = hp2->GetBinContent(bin);
      Double_t e1 = hp1->GetBinError(bin);
      Double_t e2 = hp2->GetBinError(bin);
      Double_t w1 = 1., w2 = 1.;
      if (e1 > 0) w1 = 1./(e1*e1);
      if (e2 > 0) w2 = 1./(e2*e2);
      SetBinContent(bin, (w1*v1 + w2*v2)/(w1 + w2));
      if(fSumw2.fN){
	fSumw2.fArray[bin] = 1./(w1 + w2);
      } // if
    } else {
      if(normWidth){
	Double_t w = 4.*TMath::Pi()/fNpix;
	Double_t cu = c1*hp1->GetBinContent(bin)/w;
	SetBinContent(bin, cu);
	if(fSumw2.fN){
	  Double_t e1 = hp1->GetBinError(bin)/w;
	  fSumw2.fArray[bin] = c1*c1*e1*e1;
	} // if
      } else {
	Double_t cu = c1*hp1->GetBinContent(bin) + c2*hp2->GetBinContent(bin);
	SetBinContent(bin, cu);
	if(fSumw2.fN){
	  Double_t e1 = hp1->GetBinError(bin);
	  Double_t e2 = hp2->GetBinError(bin);
	  fSumw2.fArray[bin] = c1*c1*e1*e1 + c2*c2*e2*e2;
	} // if
      } // if
    } // if
  } // i

  SetEntries(nEntries);
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
void THealPix::Build()
{
  fEntries    = 0;
  fTsumw      = 0;
  fTsumw2     = 0;
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
void THealPix::Copy(TObject& obj) const
{
  if(((THealPix&)obj).fDirectory){
    // We are likely to change the hash value of this object
    // with TNamed::Copy, to keep things correct, we need to
    // clean up its existing entries.
    ((THealPix&)obj).fDirectory->Remove(&obj);
    ((THealPix&)obj).fDirectory = 0;
  } // if
  TNamed::Copy(obj);
  ((THealPix&)obj).fOrder    = fOrder;
  ((THealPix&)obj).fNside    = fNside;
  ((THealPix&)obj).fNpix     = fNpix;
  ((THealPix&)obj).fNpFace   = fNpFace;
  ((THealPix&)obj).fNcap     = fNcap;
  ((THealPix&)obj).fIsDegree = fIsDegree;
  ((THealPix&)obj).fIsNested = fIsNested;
  ((THealPix&)obj).fType     = fType;
  ((THealPix&)obj).fUnit     = fUnit;
  ((THealPix&)obj).fTsumw    = fTsumw;
  ((THealPix&)obj).fTsumw2   = fTsumw2;

  TArray* a = dynamic_cast<TArray*>(&obj);
  if(a){
    a->Set(fNpix);
  } // if
  for(Int_t i = 0; i < fNpix; i++){
    ((THealPix&)obj).SetBinContent(i, this->GetBinContent(i));
  } // i
  ((THealPix&)obj).fEntries  = fEntries;

  fSumw2.Copy(((THealPix&)obj).fSumw2);

  if(fgAddDirectory && gDirectory){
    gDirectory->Append(&obj);
    ((THealPix&)obj).fDirectory = gDirectory;
  } // if
}

//______________________________________________________________________________
void THealPix::Divide(const THealPix* hp1)
{
  if(!hp1){
    Error("Divide", "Attempt to divide by a non-existing HEALPix");
    return;
  } // if
  
  // Check HEALPix compatibility
  if(fOrder != hp1->GetOrder() || fIsNested != hp1->IsNested()){
    Error("Divide", "Attempt to devide HEALPixs with different number of order or scheme");
    return;
  } // if

  // Create Sumw2 if hp1 has Sumw2 set
  if(fSumw2.fN == 0 && hp1->GetSumw2N() != 0){
    Sumw2();
  } // if

  // Reset statistics
  Double_t nEntries = fEntries;
  fEntries = fTsumw = 0;

  for(Int_t i = 0; i < fNpix; i++){
    Double_t c0 = GetBinContent(i);
    Double_t c1 = hp1->GetBinContent(i);
    Double_t w;
    if(c1) w = c0/c1;
    else   w = 0;
    SetBinContent(i, w);
    fEntries++;
    if(fSumw2.fN){
      Double_t e0 = GetBinError(i);
      Double_t e1 = hp1->GetBinError(i);
      Double_t c12 = c1*c1;
      if(!c1){
	fSumw2.fArray[i] = 0;
	continue;
      } // if
      fSumw2.fArray[i] = (e0*e0*c1*c1 + e1*e1*c0*c0)/(c12*c12);
    } // if
  } // i

  SetEntries(nEntries);
}

//______________________________________________________________________________
Int_t THealPix::Fill(Double_t theta, Double_t phi)
{
  /*
  if(fBuffer){
    return BufferFill(x, y, 1);
  } // if
  */

  Int_t bin = FindBin(theta, phi);
  fEntries++;
  AddBinContent(bin);
  if(fSumw2.fN){
    fSumw2.fArray[bin]++;
  } // if
  fTsumw++;
  fTsumw2++;
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

  Int_t bin = FindBin(theta, phi);
  fEntries++;
  AddBinContent(bin, w);
  if(fSumw2.fN){
    fSumw2.fArray[bin] += w*w;
  } // if
  Double_t z = (w > 0 ? w : -w);
  fTsumw  += z;
  fTsumw2 += z*z;
  return bin;
}

//______________________________________________________________________________
Int_t THealPix::Fill(const Double_t* x)
{
  TVector3 vec(x);

  if(fIsDegree){
    return Fill(vec.Theta()*TMath::RadToDeg(), vec.Phi()*TMath::RadToDeg());
  } else {
    return Fill(vec.Theta(), vec.Phi());
  } // if
}

//______________________________________________________________________________
Int_t THealPix::Fill(const Double_t* x, Double_t w)
{
  TVector3 vec(x);

  if(fIsDegree){
    return Fill(vec.Theta()*TMath::RadToDeg(), vec.Phi()*TMath::RadToDeg(), w);
  } else {
    return Fill(vec.Theta(), vec.Phi(), w);
  } // if
}

//______________________________________________________________________________
Int_t THealPix::Fill(const Float_t* x)
{
  TVector3 vec(x);

  if(fIsDegree){
    return Fill(vec.Theta()*TMath::RadToDeg(), vec.Phi()*TMath::RadToDeg());
  } else {
    return Fill(vec.Theta(), vec.Phi());
  } // if
}

//______________________________________________________________________________
Int_t THealPix::Fill(const Float_t* x, Double_t w)
{
  TVector3 vec(x);

  if(fIsDegree){
    return Fill(vec.Theta()*TMath::RadToDeg(), vec.Phi()*TMath::RadToDeg(), w);
  } else {
    return Fill(vec.Theta(), vec.Phi(), w);
  } // if
}

//_____________________________________________________________________________
Int_t THealPix::FindBin(Double_t theta, Double_t phi) const
{
  // Oribinal code is Healpix_Base::ang2pix_z_phi of HEALPix C++
  if(fIsDegree){
    theta *= TMath::DegToRad();
    phi   *= TMath::DegToRad();
  } // if

  // Theta must be in range [0, pi] # NOT [0, pi)
  while(theta < 0)               theta += TMath::TwoPi();
  while(theta >= TMath::TwoPi()) theta -= TMath::TwoPi();
  if(theta > TMath::Pi()){
    theta -= TMath::Pi();
    phi   += TMath::Pi();
  } // if

  // Phi must be in range [0, 2pi)
  while(phi   < 0)               phi   += TMath::TwoPi();
  while(phi   >= TMath::TwoPi()) phi   -= TMath::TwoPi();

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
    
    Int_t ipf = XY2Pix(ix, iy);
    
    ipf >>= (2*(order_max - fOrder));  // in {0, fNside**2 - 1}
    
    return ipf + (face_num << (2*fOrder));    // in {0, 12*fNside**2 - 1}
  } // if
}

//______________________________________________________________________________
Double_t THealPix::GetBinContent(Int_t) const
{
  AbstractMethod("GetBinContent");
  return 0;
}

//______________________________________________________________________________
Bool_t THealPix::GetDefaultSumw2()
{
  // static function
  // return kTRUE if THealPix::Sumw2 must be called when creating new HEALPix.
  // see THealPix::SetDefaultSumw2.

  return fgDefaultSumw2;
}

//______________________________________________________________________________
Double_t THealPix::GetEntries() const
{
  // return the current number of entries
  return fEntries;
}

//______________________________________________________________________________
void THealPix::GetRingInfo(Int_t ring, Int_t& startpix, Int_t& ringpix,
			   Double_t &costheta, Double_t& sintheta,
			   Bool_t& shifted) const
{
  if(fIsNested){
    Error("GetRingInfo", "Only RING scheme is supported");
    return;
  } // if

  Int_t northring = (ring > 2*fNside) ? 4*fNside - ring : ring;

  if(northring < fNside){
    costheta = 1 - northring*northring*4./fNpix;
    sintheta = TMath::Sin(2*TMath::ASin(northring/(TMath::Sqrt(6.)*fNside)));
    ringpix = 4*northring;
    shifted = kTRUE;
    startpix = 2*northring*(northring - 1);
  } else {
    costheta = (2*fNside - northring)*(8.*fNside/fNpix);
    sintheta = TMath::Sqrt(1 - costheta*costheta);
    ringpix = 4*fNside;
    shifted = ((northring - fNside) & 1) == 0;
    startpix = fNcap + (northring - fNside)*ringpix;
  } // if

  if(northring != ring){ // southern hemisphere
    costheta = -costheta;
    startpix = fNpix - startpix - ringpix;
  } // if
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
void THealPix::Multiply(const THealPix* hp1)
{
  if(!hp1){
    Error("Multiply", "Attempt to multiply a non-existing HEALPix");
    return;
  } // if
  
  // Check HEALPix compatibility
  if(fOrder != hp1->GetOrder() || fIsNested != hp1->IsNested()){
    Error("Multiply", "Attempt to multiply HEALPixs with different number of order or scheme");
    return;
  } // if

  // Create Sumw2 if hp1 has Sumw2 set
  if(fSumw2.fN == 0 && hp1->GetSumw2N() != 0){
    Sumw2();
  } // if

  // Reset statistics
  Double_t nEntries = fEntries;
  fEntries = fTsumw = 0;
  
  for(Int_t bin = 0; bin < fNpix; bin++){
    Double_t c0 = GetBinContent(bin);
    Double_t c1 = hp1->GetBinContent(bin);
    Double_t w  = c0*c1;
    SetBinContent(bin, w);
    fEntries++;
    if(fSumw2.fN){
      Double_t e0 = GetBinError(bin);
      Double_t e1 = hp1->GetBinError(bin);
      fSumw2.fArray[bin] = (e0*e0*c1*c1 + e1*e1*c0*c0);
    } // if
  } // i

  SetEntries(nEntries);
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

  Int_t npix = THealUtil::Nside2Npix(nside);

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
THealPix* THealPix::Rebin(Int_t neworder, const char* newname)
{
  // Rebin this HEALPix
  // -case 1 ngroup>0
  // -case 2 ngroup<0
  if(neworder > 13 || neworder < 0){
    Error("Rebin", "Illegal value of neworder=%d",neworder);
    return 0;
  } // if

  THealPix* hpnew = this;
  if(newname && strlen(newname) > 0){
    hpnew = (THealPix*)Clone(newname);
  } // if
  if(neworder == fOrder){
    return hpnew;
  } // if

  // Save old bin contents into a new array
  Double_t entries = fEntries;
  Double_t* oldBins = new Double_t[fNpix];
  for(Int_t i = 0; i < fNpix; i++){
    oldBins[i] = GetBinContent(i);
  } // i
  Double_t *oldErrors = 0;
  if (fSumw2.fN != 0) {
    oldErrors = new Double_t[fNpix];
    for(Int_t bin = 0; bin < fNpix; bin++){
      oldErrors[bin] = GetBinError(bin);
    } // if
  } // if

  hpnew->SetOrder(neworder);
  hpnew->SetBinsLength(hpnew->GetNpix());

  Int_t newnpix  = hpnew->GetNpix();
  if(IsNested()){ // NESTED
    if(neworder < fOrder){ // rebin to lower level
      Int_t div = fNpix/newnpix;
      for(Int_t i = 0; i < newnpix; i++){
	Double_t content = 0;
	Double_t error   = 0;
	for(Int_t j = i*div; j < (i + 1)*div; j++){
	  content += oldBins[j];
	  if(oldErrors){
	    error += oldErrors[j]*oldErrors[j];
	  } // if
	} // j
	hpnew->SetBinContent(i, content);
	if(oldErrors){
	  hpnew->SetBinError(i, TMath::Sqrt(error));
	} // if
      } // i
    } else { // rebin to higher level
      Int_t div = newnpix/fNpix;
      for(Int_t i = 0; i < newnpix; i++){
	hpnew->SetBinContent(i, oldBins[i/div]/div);
	if(oldErrors){
	  hpnew->SetBinError(i, oldErrors[i/div]/TMath::Sqrt(div));
	} // if
      } // i
    } // if
  } else { // RING
    if(neworder < fOrder){ // rebin to lower level
      Int_t div = fNpix/newnpix;
      for(Int_t i = 0; i < newnpix; i++){
	Double_t content = 0;
	Double_t error   = 0;
	Int_t inest = hpnew->Ring2Nest(i);
	for(Int_t jnest = inest*div; jnest < (inest + 1)*div; jnest++){
	  Int_t j = Nest2Ring(jnest);
	  content += oldBins[j];
	  if(oldErrors){
	    error += oldErrors[j]*oldErrors[j];
	  } // if
	} // j
	hpnew->SetBinContent(i, content);
	if(oldErrors){
	  hpnew->SetBinError(i, TMath::Sqrt(error));
	} // if
      } // i
    } else { // rebin to higher level
      Int_t div = newnpix/fNpix;
      for(Int_t i = 0; i < newnpix; i++){
	Int_t inest = hpnew->Ring2Nest(i);
	Int_t p = Nest2Ring(inest/div);
	hpnew->SetBinContent(i, oldBins[p]/div);
	if(oldErrors){
	  hpnew->SetBinError(i, oldErrors[p]/TMath::Sqrt(div));
	} // if
      } // i
    } // if
  } // if

  hpnew->SetEntries(entries); //was modified by SetBinContent
  delete [] oldBins;
  if(oldErrors){
    delete [] oldErrors;
  } // if
  return hpnew;
}

//_____________________________________________________________________________
void THealPix::SetBinContent(Int_t, Double_t)
{
  AbstractMethod("SetBinContent");
}

//______________________________________________________________________________
void THealPix::SetDefaultSumw2(Bool_t sumw2)
{
  // static function.
  // When this static function is called with sumw2=kTRUE, all new
  // HEALPix will automatically activate the storage of
  // the sum of squares of errors, ie THealPix::Sumw2 is automatically called.
  
  fgDefaultSumw2 = sumw2;
}

//_____________________________________________________________________________
void THealPix::SetOrder(Int_t order)
{
  if(order < 0)  order = 0;
  if(order > 13) order = 13;

  fOrder = order;
  fNside = THealUtil::Order2Nside(fOrder);
  fNpix  = THealUtil::Nside2Npix(fNside);
  fNpFace= fNside*fNside;
  fNcap  = 2*(fNpFace - fNside);
}

//_____________________________________________________________________________
void THealPix::SetUnit(const char* unit)
{
  fUnit = std::string(unit);
}

//______________________________________________________________________________
void THealPix::Sumw2()
{
  // Create structure to store sum of squares of weights*-*-*-*-*-*-*-*
  //
  //     if HEALPix is already filled, the sum of squares of weights
  //     is filled with the existing bin contents
  //
  //     The error per bin will be computed as sqrt(sum of squares of weight)
  //     for each bin.
  //
  //  This function is automatically called when the HEALPix is created
  //  if the static function THealPix::SetDefaultSumw2 has been called before.

   if(!fgDefaultSumw2 && fSumw2.fN){
      Warning("Sumw2", "Sum of squares of weights structure already created");
      return;
   } // if

   fSumw2.Set(fNpix);

   for(Int_t bin = 0; bin < fNpix; bin++) {
     fSumw2.fArray[bin] = GetBinContent(bin);
   } // bin
}

//______________________________________________________________________________
Double_t THealPix::GetBinError(Int_t bin) const
{
  //   -*-*-*-*-*Return value of error associated to bin number bin*-*-*-*-*
  //             ==================================================
  //
  //    if the sum of squares of weights has been defined (via Sumw2),
  //    this function returns the sqrt(sum of w2).
  //    otherwise it returns the sqrt(contents) for this bin.
  //
  //   -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  if(bin < 0){
    bin = 0;
  } // if
  if(bin >= fNpix){
    bin = fNpix - 1;
  } // if

  if(fSumw2.fN){
    Double_t err2 = fSumw2.fArray[bin];
    return TMath::Sqrt(err2);
  } // if
  Double_t error2 = TMath::Abs(GetBinContent(bin));
  return TMath::Sqrt(error2);
}

//_____________________________________________________________________________
Int_t THealPix::GetNrows() const
{
  return fOrder < 4 ? 1 : 1024;
}

//_____________________________________________________________________________
Double_t THealPix::GetPixelArea(Bool_t degree2) const
{
  Double_t area = TMath::Pi()*4/fNpix;
  if(degree2){
    area *= TMath::RadToDeg()*TMath::RadToDeg();
  } // if

  return area;
}

//_____________________________________________________________________________
std::string THealPix::GetSchemeString() const
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
  } else if(fType == TFLOAT){
    return "F";
  } // if

  return "D";
}

//_____________________________________________________________________________
void THealPix::Scale(Double_t c1, Option_t* option)
{
  TString opt = option;
  opt.ToLower();
  Double_t ent = fEntries;
  if(opt.Contains("width")){
    Add(this, this, c1, -1);
  } else {
    Add(this, this, c1, 0);
  } // if
  fEntries = ent;
}

//______________________________________________________________________________
void THealPix::SetBinError(Int_t bin, Double_t error)
{
  // see convention for numbering bins in THealPix::GetBin
  if(!fSumw2.fN){
    Sumw2();
  } // if
  if(bin <0 || bin >= fSumw2.fN){
    return;
  } // if
  fSumw2.fArray[bin] = error*error;
}

//______________________________________________________________________________
void THealPix::SetDirectory(TDirectory *dir)
{
  // By default when a HEALPix is created, it is added to the list
  // of HEALPix objects in the current directory in memory.
  // Remove reference to this HEALPix from current directory and add
  // reference to new directory dir. dir can be 0 in which case the
  // HEALPix does not belong to any directory.

  if(fDirectory == dir){
    return;
  } // if
  if(fDirectory){
    fDirectory->Remove(this);
  } // if

  fDirectory = dir;

  if(fDirectory){
    fDirectory->Append(this);
  } // if
}

//_____________________________________________________________________________
void THealPix::Nest2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const
{
  // Original code is Healpix_Base::nest2xyf of HEALPix C++
  face = pix>>(2*fOrder);
  Pix2XY(pix & (fNpFace - 1), x, y);
}

//_____________________________________________________________________________
Int_t THealPix::XYF2Nest(Int_t x, Int_t y, Int_t face) const
{
  // Original code is Healpix_Base::xyf2nest of HEALPix C++
  return (face<<(2*fOrder)) + XY2Pix(x, y);
}

//_____________________________________________________________________________
void THealPix::Ring2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const
{
  // Original code is Healpix_Base::ring2xyf of HEALPix C++
  Int_t iring, iphi, kshift, nr;
  
  Int_t nl2 = 2*fNside;
  
  if(pix < fNcap){ // North Polar cap
    iring = Int_t(0.5*(1 + THealUtil::Isqrt(1 + 2*pix))); //counted from North pole
    iphi  = (pix + 1) - 2*iring*(iring - 1);
    kshift = 0;
    nr = iring;
    face = 0;
    Int_t tmp = iphi - 1;
    if(tmp >= 2*iring){
      face = 2;
      tmp -= 2*iring;
    } // if
    if(tmp >= iring){
      face++;
    } // if
  } else if(pix < fNpix - fNcap){ // Equatorial region
    Int_t ip = pix - fNcap;
    if(fOrder >= 0){
      iring = (ip>>(fOrder + 2)) + fNside; // counted from North pole
      iphi  = (ip&(4*fNside - 1)) + 1;
    } else {
      iring = (ip/(4*fNside)) + fNside; // counted from North pole
      iphi = (ip%(4*fNside)) + 1;
    } // if
    kshift = (iring + fNside)&1;
    nr = fNside;
    UInt_t ire = iring - fNside + 1;
    UInt_t irm = nl2 + 2 - ire;
    Int_t ifm, ifp;
    if(fOrder >= 0){
      ifm = (iphi - ire/2 + fNside -1) >> fOrder;
      ifp = (iphi - irm/2 + fNside -1) >> fOrder;
    } else {
      ifm = (iphi - ire/2 + fNside -1) / fNside;
      ifp = (iphi - irm/2 + fNside -1) / fNside;
    } // if
    if(ifp == ifm){ // faces 4 to 7
      face = (ifp==4) ? 4 : ifp+4;
    } else if(ifp<ifm){ // (half-)faces 0 to 3
      face = ifp;
    } else { // (half-)faces 8 to 11
      face = ifm + 8;
    } // if
  } else { // South Polar cap
    Int_t ip = fNpix - pix;
    iring = Int_t(0.5*(1 + THealUtil::Isqrt(2*ip - 1))); //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring - 1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2 - iring;
    face = 8;
    Int_t tmp = iphi - 1;
    if(tmp >= 2*nr){
      face = 10;
      tmp -= 2*nr;
    } // if
    if(tmp >= nr){
      ++face;
    } // if
  } // if

  Int_t irt = iring - (fgJrll[face]*fNside) + 1;
  Int_t ipt = 2*iphi- fgJpll[face]*nr - kshift -1;
  if(ipt >= nl2){
    ipt -= 8*fNside;
  } // if

  x =  (ipt-irt) >>1;
  y =(-(ipt+irt))>>1;
}

//_____________________________________________________________________________
Int_t THealPix::XYF2Ring(Int_t x, Int_t y, Int_t face) const
{
  // Original code is Healpix_Base::xyf2ring of HEALPix C++
  Int_t nl4 = 4*fNside;
  Int_t jr = (fgJrll[face]*fNside) - x - y  - 1;

  Int_t nr, kshift, n_before;
  if(jr < fNside){
    nr = jr;
    n_before = 2*nr*(nr - 1);
    kshift = 0;
  } else if (jr > 3*fNside){
    nr = nl4 - jr;
    n_before = fNpix - 2*(nr + 1)*nr;
    kshift = 0;
  } else {
    nr = fNside;
    n_before = fNcap + (jr - fNside)*nl4;
    kshift = (jr - fNside)&1;
  } // if

  Int_t jp = (fgJpll[face]*nr + x - y + 1 + kshift)/2;
  if(jp > nl4){
    jp -= nl4;
  } else {
    if(jp < 1){
      jp += nl4;
    } // if
  } // if

  return n_before + jp - 1;
}

//_____________________________________________________________________________
Int_t THealPix::Nest2Ring(Int_t pix) const
{
  // Original code is Healpix_Base::nest2ring of HEALPix C++
  Int_t x, y, face;
  Nest2XYF(pix, x, y, face);
  return XYF2Ring(x, y, face);
}

//_____________________________________________________________________________
Int_t THealPix::Ring2Nest(Int_t pix) const
{
  // Original code is Healpix_Base::ring2nest of HEALPix C++
  Int_t x, y, face;
  Ring2XYF(pix, x, y, face);
  return XYF2Nest(x, y, face);
}

//_____________________________________________________________________________
void THealPix::Pix2XY(Int_t pix, Int_t& x, Int_t& y) const
{
  // Original code is Healpix_Base::pix2xy of HEALPix C++
  Int_t raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  x = fgTable.C(raw&0xff) | (fgTable.C(raw>>8)<<4);
  raw = ((pix&0xaaaa)>>1) | ((pix&0xaaaa0000)>>16);
  y = fgTable.C(raw&0xff) | (fgTable.C(raw>>8)<<4);
}

//_____________________________________________________________________________
Int_t THealPix::XY2Pix(Int_t x, Int_t y) const
{
  // Original code is Healpix_Base::xy2pix of HEALPix C++
  return fgTable.U(x&0xff) | (fgTable.U(x>>8)<<16) | (fgTable.U(y&0xff)<<1)
    | (fgTable.U(y>>8)<<17);
}

ClassImp(THealPixF)

//_____________________________________________________________________________
// THealPixF methods
//_____________________________________________________________________________
THealPixF::THealPixF(): THealPix(), TArrayF()
{
  if(fgDefaultSumw2){
    Sumw2();
  } // if
}

//_____________________________________________________________________________
THealPixF::THealPixF(const char* name, const char* title, Int_t order,
		     Bool_t isNested)
: THealPix(name, title, order, isNested)
{
  fType = TFLOAT;
  TArrayF::Set(fNpix);
  if(fgDefaultSumw2){
    Sumw2();
  } // if
}

//_____________________________________________________________________________
THealPixF::THealPixF(const THealPixF& hpd) : THealPix(), TArrayF()
{
  ((THealPixF&)hpd).Copy(*this);
}

//_____________________________________________________________________________
THealPixF::~THealPixF()
{
}

//_____________________________________________________________________________
void THealPixF::Copy(TObject& newhp) const
{
  THealPix::Copy(newhp);
}

//______________________________________________________________________________
Double_t THealPixF::GetBinContent(Int_t bin) const
{
  if(!fArray){
    return 0;
  } // if

  bin = bin < 0 ? 0 : bin;
  bin = bin < fNpix ? bin : fNpix - 1;

  return Double_t(fArray[bin]);
}

//______________________________________________________________________________
THealPixF* THealPixF::ReadFits(const char* fname, const char* colname)
{
  THealPix::HealHeader_t head;
  fitsfile* fptr = 0;

  if(!THealPix::ReadFitsHeader(&fptr, fname, colname, head)){
    return 0;
  } // if

  THealPixF* hpf = new THealPixF(fname, colname, head.order, head.isNested);
  hpf->SetUnit(head.tunit);

  Long_t npercol = head.npix/head.nrows;
  Int_t status = 0;

  for(Int_t i = 0; i < head.nrows; i++){
    fits_read_col(fptr, TFLOAT, head.colnum, i+1, 1, npercol, 0,
		  &(hpf->GetArray()[i*npercol]), 0, &status);
    if(!THealUtil::FitsReportError(status)){
      delete hpf;
      return 0;
    } // if
  } // i

  return hpf;
}

//_____________________________________________________________________________
void THealPixF::SetBinContent(Int_t bin, Double_t content)
{
  if(bin < 0 || fNpix <= 0){
    return;
  } // if
  fArray[bin] = content;
  fEntries++;
  fTsumw = 0;
}

//_____________________________________________________________________________
void THealPixF::SetBinsLength(Int_t n)
{
  // Set total number of bins
  // Reallocate bin contents array
  if(n < 0){
    n = fNpix;
  } // if
  TArrayF::Set(n);
}

//______________________________________________________________________________
THealPixF& THealPixF::operator=(const THealPixF& hp1)
{
   // Operator =

  if(this != &hp1){
    ((THealPixF&)hp1).Copy(*this);
  } // if
  return *this;
}

//______________________________________________________________________________
THealPixF THealPixF::__add__(const THealPixF& hp1) const
{
  // Python operator +
  return *this + hp1;
}

//______________________________________________________________________________
THealPixF THealPixF::__div__(const THealPixF& hp1) const
{
  // Python operator /
  return *this / hp1;
}

//______________________________________________________________________________
THealPixF THealPixF::__mul__(const THealPixF& hp1) const
{
  // Python operator *
  return *this * hp1;
}

//______________________________________________________________________________
THealPixF THealPixF::__mul__(Double_t c1) const
{
  // Python operator *
  return *this*c1;
}

//______________________________________________________________________________
THealPixF THealPixF::__rmul__(Double_t c1) const
{
  // Python operator *
  return *this*c1;
}

//______________________________________________________________________________
THealPixF THealPixF::__sub__(const THealPixF& hp1) const
{
  // Python operator -
  return *this - hp1;
}

//______________________________________________________________________________
THealPixF operator*(Double_t c1, const THealPixF& hp1)
{
   // Operator *

   THealPixF hpnew = hp1;
   hpnew.Scale(c1);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixF operator+(const THealPixF& hp1, const THealPixF& hp2)
{
   // Operator +

   THealPixF hpnew = hp1;
   hpnew.Add(&hp2, 1);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixF operator-(const THealPixF& hp1, const THealPixF& hp2)
{
   // Operator -

   THealPixF hpnew = hp1;
   hpnew.Add(&hp2, -1);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixF operator*(const THealPixF& hp1, const THealPixF& hp2)
{
   // Operator *

   THealPixF hpnew = hp1;
   hpnew.Multiply(&hp2);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixF operator/(const THealPixF& hp1, const THealPixF& hp2)
{
   // Operator /

   THealPixF hpnew = hp1;
   hpnew.Divide(&hp2);
   hpnew.SetDirectory(0);
   return hpnew;
}

ClassImp(THealPixD)

//_____________________________________________________________________________
// THealPixD methods
//_____________________________________________________________________________
THealPixD::THealPixD(): THealPix(), TArrayD()
{
  if(fgDefaultSumw2){
    Sumw2();
  } // if
}

//_____________________________________________________________________________
THealPixD::THealPixD(const char* name, const char* title, Int_t order,
		     Bool_t isNested)
: THealPix(name, title, order, isNested)
{
  fType = TDOUBLE;
  TArrayD::Set(fNpix);
  if(fgDefaultSumw2){
    Sumw2();
  } // if
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

  bin = bin < 0 ? 0 : bin;
  bin = bin < fNpix ? bin : fNpix - 1;

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

//_____________________________________________________________________________
void THealPixD::SetBinContent(Int_t bin, Double_t content)
{
  if(bin < 0 || fNpix <= 0){
    return;
  } // if
  fArray[bin] = content;
  fEntries++;
  fTsumw = 0;
}

//_____________________________________________________________________________
void THealPixD::SetBinsLength(Int_t n)
{
  // Set total number of bins
  // Reallocate bin contents array
  if(n < 0){
    n = fNpix;
  } // if
  TArrayD::Set(n);
}

//______________________________________________________________________________
THealPixD& THealPixD::operator=(const THealPixD& hp1)
{
   // Operator =

  if(this != &hp1){
    ((THealPixD&)hp1).Copy(*this);
  } // if
  return *this;
}

//______________________________________________________________________________
THealPixD THealPixD::__add__(const THealPixD& hp1) const
{
  // Python operator +
  return *this + hp1;
}

//______________________________________________________________________________
THealPixD THealPixD::__div__(const THealPixD& hp1) const
{
  // Python operator /
  return *this / hp1;
}

//______________________________________________________________________________
THealPixD THealPixD::__mul__(const THealPixD& hp1) const
{
  // Python operator *
  return *this * hp1;
}

//______________________________________________________________________________
THealPixD THealPixD::__mul__(Double_t c1) const
{
  // Python operator *
  return *this*c1;
}

//______________________________________________________________________________
THealPixD THealPixD::__rmul__(Double_t c1) const
{
  // Python operator *
  return *this*c1;
}

//______________________________________________________________________________
THealPixD THealPixD::__sub__(const THealPixD& hp1) const
{
  // Python operator -
  return *this - hp1;
}

//______________________________________________________________________________
THealPixD operator*(Double_t c1, const THealPixD& hp1)
{
   // Operator *

   THealPixD hpnew = hp1;
   hpnew.Scale(c1);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixD operator+(const THealPixD& hp1, const THealPixD& hp2)
{
   // Operator +

   THealPixD hpnew = hp1;
   hpnew.Add(&hp2, 1);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixD operator-(const THealPixD& hp1, const THealPixD& hp2)
{
   // Operator -

   THealPixD hpnew = hp1;
   hpnew.Add(&hp2, -1);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixD operator*(const THealPixD& hp1, const THealPixD& hp2)
{
   // Operator *

   THealPixD hpnew = hp1;
   hpnew.Multiply(&hp2);
   hpnew.SetDirectory(0);
   return hpnew;
}

//______________________________________________________________________________
THealPixD operator/(const THealPixD& hp1, const THealPixD& hp2)
{
   // Operator /

   THealPixD hpnew = hp1;
   hpnew.Divide(&hp2);
   hpnew.SetDirectory(0);
   return hpnew;
}
