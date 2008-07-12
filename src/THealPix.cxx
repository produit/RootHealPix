// $Id: THealPix.cxx,v 1.31 2008/07/12 22:06:48 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.

   This is a port of HEALPix C++ package to ROOT system.
   Original code is available at <http://healpix.jpl.nasa.gov> under GPL.
******************************************************************************/

#include <complex>
#include <cstring>
#include <iostream>
#include "fitsio.h"

#include "TDirectory.h"
#include "THashList.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVector3.h"

#include "THealFFT.h"
#include "THealPix.h"
#include "THealUtil.h"
#include "TVirtualHealPainter.h"

namespace {
  inline UInt_t Isqrt(UInt_t v)
  {
    // Calculate square root of an integer
    // Original code is isqrt of HEALPix C++
    return UInt_t(TMath::Sqrt(v + .5));
  }
  
//______________________________________________________________________________
  inline Double_t Modulo(Double_t v1, Double_t v2)
  {
    // Calculate modulo of two doubles
    // Original code is modulo of HEALPix C++
    return (v1 >= 0) ? ((v1 < v2) ? v1 : fmod(v1, v2)) : (fmod(v1, v2) + v2);
  }
  
//______________________________________________________________________________
  inline Int_t Modulo(Int_t v1, Int_t v2)
  {
    // Calculate modulo of two integers
    // Original code is modulo of HEALPix C++
    return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
  }
  
//______________________________________________________________________________
  inline Long_t Modulo(Long_t v1, Long_t v2)
  {
    // Calculate modulo of two long integers
    // Original code is modulo of HEALPix C++
    return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
  }

//______________________________________________________________________________
  void get_chunk_info (int nrings, int &nchunks, int &chunksize)
  {
    nchunks = nrings/TMath::Max(100,nrings/10) + 1;
    chunksize = (nrings+nchunks-1)/nchunks;
  }

//______________________________________________________________________________
  void fill_work (const std::complex<double> *datain, int nph, int mmax,
		  bool shifted, const std::vector<std::complex<double> > &shiftarr,
		  std::vector<std::complex<double> > &work)
  {
    for (int m=1; m<nph; ++m) work[m]=0;
    work[0]=datain[0];
    
    int cnt1=0, cnt2=nph;
    for(int m=1; m<=mmax; ++m){
      if (++cnt1==nph) cnt1=0;
      if (--cnt2==-1) cnt2=nph-1;
      std::complex<double> tmp = shifted ? (datain[m]*shiftarr[m]) : datain[m];
      work[cnt1] += tmp;
      work[cnt2] += conj(tmp);
    }
  }

//______________________________________________________________________________
  void read_work (const std::vector<std::complex<double> >& work, int nph,
		  int mmax, bool shifted,
		  const std::vector<std::complex<double> > &shiftarr,
		  std::complex<double> *dataout)
  {
    int cnt2=0;
    for(int m=0; m<=mmax; ++m){
      dataout[m] = work[cnt2];
      if (++cnt2==nph) cnt2=0;
    }
    if (shifted)
      for (int m=0; m<=mmax; ++m) dataout[m] *= shiftarr[m];
  }

//______________________________________________________________________________
  template<typename T>
  void fft_map2alm (int nph, int mmax, bool shifted, double weight,
		    THealFFT::rfft &plan, const T *mapN, const T *mapS,
		    std::complex<double> *phas_n, std::complex<double> *phas_s,
		    const std::vector<std::complex<double> > &shiftarr,
		    std::vector<std::complex<double> > &work)
  {
    for (int m=0; m<nph; ++m) work[m] = mapN[m]*weight;
    plan.forward_c(work);
    read_work (work, nph, mmax, shifted, shiftarr, phas_n);
    if (mapN!=mapS){
      for (int m=0; m<nph; ++m) work[m] = mapS[m]*weight;
      plan.forward_c(work);
      read_work (work, nph, mmax, shifted, shiftarr, phas_s);
    } else {
      for (int m=0; m<=mmax; ++m) phas_s[m]=0;
    } // if
  }
  
//______________________________________________________________________________
  void recalc_map2alm (int nph, int mmax, THealFFT::rfft &plan,
		       std::vector<std::complex<double> > &shiftarr)
  {
    if (plan.size() == nph) return;
    plan.Set (nph);
    double f1 = TMath::Pi()/nph;
    for (int m=0; m<=mmax; ++m){
      if (m<nph)
	shiftarr[m] = std::complex<double>(cos(m*f1),-sin(m*f1));
      else
	shiftarr[m]=-shiftarr[m-nph];
    } // m
  }

//______________________________________________________________________________
  void recalc_alm2map (int nph, int mmax, THealFFT::rfft &plan,
		       std::vector<std::complex<double> > &shiftarr)
  {
    if (plan.size() == nph) return;
    plan.Set (nph);
    double f1 = TMath::Pi()/nph;
    for (int m=0; m<=mmax; ++m){
      if (m<nph)
	shiftarr[m] = std::complex<double>(cos(m*f1),sin(m*f1));
      else
	shiftarr[m]=-shiftarr[m-nph];
    }
  }

//______________________________________________________________________________
  template<typename T>
  void fft_alm2map (int nph, int mmax, bool shifted, THealFFT::rfft &plan,
		    T *mapN, T *mapS, std::complex<double> *b_north,
		    const std::complex<double> *b_south,
		    const std::vector<std::complex<double> > &shiftarr,
		    std::vector<std::complex<double> > &work)
  {
    fill_work (b_north, nph, mmax, shifted, shiftarr, work);
    plan.backward_c(work);
    for (int m=0; m<nph; ++m) mapN[m] = work[m].real();
    if (mapN==mapS) return;
    fill_work (b_south, nph, mmax, shifted, shiftarr, work);
    plan.backward_c(work);
    for (int m=0; m<nph; ++m) mapS[m] = work[m].real();
  }
}

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
THealPix::THealPix() : TNamed(), TAttLine(), TAttFill(), TAttMarker()
{
  SetOrder(0);
  fDirectory = 0;
  fFunctions = new TList;
  fPainter   = 0;
  fEntries   = 0;
  fNormFactor= 0;
  fTsumw     = fTsumw2 = 0;
  fMaximum   = -1111;
  fMinimum   = -1111;
  fUnit      = "";
  fIsDegree  = kFALSE;
  fIsNested  = kFALSE;
  fXaxis.SetName("xaxis");
  fYaxis.SetName("yaxis");
  fZaxis.SetName("zaxis");
  fXaxis.SetParent(this);
  fYaxis.SetParent(this);
  fZaxis.SetParent(this);
  UseCurrentStyle();
}

//_____________________________________________________________________________
THealPix::~THealPix()
{
  if(!TestBit(kNotDeleted)){
    return;
  } // if

  if(fFunctions){
    fFunctions->SetBit(kInvalidObject);
    TObject* obj = 0;
    //special logic to support the case where the same object is
    //added multiple times in fFunctions.
    //This case happens when the same object is added with different
    //drawing modes
    //In the loop below we must be careful with objects (eg TCutG) that may
    // have been added to the list of functions of several histograms
    //and may have been already deleted.
    while((obj = fFunctions->First())){
      while(fFunctions->Remove(obj)){ }
      if(!obj->TestBit(kNotDeleted)){
	break;
      } // if
      delete obj;
      obj = 0;
    } // while
    delete fFunctions;
    fFunctions = 0;
  } // if

  if(fDirectory){
    fDirectory->Remove(this);
    fDirectory = 0;
  } // if

  delete fPainter;
  fPainter = 0;
}

//_____________________________________________________________________________
THealPix::THealPix(const char* name, const char* title, Int_t order,
		   Bool_t isNested)
  : TNamed(name, title), TAttLine(), TAttFill(), TAttMarker()
{
  fIsNested = isNested;
  fUnit     = "";
  Build();

  if(order < 0)  order = 0;
  if(order > 13) order = 13;

  SetOrder(order);
  if(fgDefaultSumw2){
    Sumw2();
  } // if
}

//_____________________________________________________________________________
THealPix::THealPix(const THealPix& hp) : TNamed(), TAttLine(), TAttFill(), TAttMarker()
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
  if(hp1->GetNormFactor() != 0){
    factor = hp1->GetNormFactor()/hp1->GetSumOfWeights();;
  } // if
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
  fDirectory  = 0;
  fPainter    = 0;
  fEntries    = 0;
  fNormFactor = 0;
  fIsDegree   = kFALSE;
  fTsumw      = 0;
  fTsumw2     = 0;
  fMaximum    = -1111;
  fMinimum    = -1111;
  fUnit       = "";
  fXaxis.SetName("xaxis");
  fYaxis.SetName("yaxis");
  fZaxis.SetName("zaxis");
  fXaxis.Set(1, 0., 360.);
  fYaxis.Set(1, 0., 180.);
  fZaxis.Set(1, 0., 1.);
  fXaxis.SetNdivisions(408, kFALSE);
  fYaxis.SetNdivisions(408, kFALSE);
  fXaxis.SetParent(this);
  fYaxis.SetParent(this);
  fZaxis.SetParent(this);

  SetTitle(fTitle.Data());
  
  fFunctions = new TList;

  UseCurrentStyle();

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
  ((THealPix&)obj).fNormFactor= fNormFactor;
  ((THealPix&)obj).fBarOffset = fBarOffset;
  ((THealPix&)obj).fBarWidth  = fBarWidth;
  ((THealPix&)obj).fOrder     = fOrder;
  ((THealPix&)obj).fNside     = fNside;
  ((THealPix&)obj).fNpix      = fNpix;
  ((THealPix&)obj).fNside2    = fNside2;
  ((THealPix&)obj).f3Nside2   = f3Nside2;
  ((THealPix&)obj).f2over3Nside= f2over3Nside;
  ((THealPix&)obj).f4Nside    = f4Nside;
  ((THealPix&)obj).f2Nside    = f2Nside;
  ((THealPix&)obj).fNcap      = fNcap;
  ((THealPix&)obj).fIsDegree  = fIsDegree;
  ((THealPix&)obj).fIsNested  = fIsNested;
  ((THealPix&)obj).fUnit      = fUnit;
  ((THealPix&)obj).fTsumw     = fTsumw;
  ((THealPix&)obj).fTsumw2    = fTsumw2;
  ((THealPix&)obj).fMaximum   = fMaximum;
  ((THealPix&)obj).fMinimum   = fMinimum;
  ((THealPix&)obj).fOption    = fOption;

  TArray* a = dynamic_cast<TArray*>(&obj);
  if(a){
    a->Set(fNpix);
  } // if
  for(Int_t i = 0; i < fNpix; i++){
    ((THealPix&)obj).SetBinContent(i, this->GetBinContent(i));
  } // i
  ((THealPix&)obj).fEntries  = fEntries;
  
  TAttLine::Copy(((THealPix&)obj));
  TAttFill::Copy(((THealPix&)obj));
  TAttMarker::Copy(((THealPix&)obj));
  fXaxis.Copy(((THealPix&)obj).fXaxis);
  fYaxis.Copy(((THealPix&)obj).fYaxis);
  fZaxis.Copy(((THealPix&)obj).fZaxis);
  ((THealPix&)obj).fXaxis.SetParent(&obj);
  ((THealPix&)obj).fYaxis.SetParent(&obj);
  ((THealPix&)obj).fZaxis.SetParent(&obj);
  fContour.Copy(((THealPix&)obj).fContour);
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
  Double_t tt = ::Modulo(phi, TMath::TwoPi())/TMath::PiOver2();

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
      ip = ::Modulo(ip, f4Nside);
      
      return fNcap + (ir - 1)*f4Nside + ip;
    } else { // North & South polar caps
      Double_t tp = tt - Int_t(tt);
      Double_t tmp = fNside*TMath::Sqrt(3*(1 - zabs));
      
      Int_t jp = Int_t(tp*tmp); // increasing edge line index
      Int_t jm = Int_t((1. - tp)*tmp); // decreasing edge line index
      
      Int_t ir = jp + jm + 1; // ring number counted from the closest pole
      Int_t ip = Int_t(tt*ir); // in {0,4*ir-1}
      ip = ::Modulo(ip, 4*ir);
      
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
Double_t THealPix::GetAverage() const
{
  Double_t total = 0;
  for(Int_t i = 0; i < fNpix; i++){
    total += GetBinContent(i);
  } // i
  
  return total/fNpix;
}

//_____________________________________________________________________________
Double_t THealPix::GetBinArea(Bool_t degree2) const
{
  Double_t area = TMath::Pi()*4/fNpix;
  if(degree2){
    area *= TMath::RadToDeg()*TMath::RadToDeg();
  } // if

  return area;
}

//______________________________________________________________________________
void THealPix::GetBinCenter(Int_t bin, Double_t& theta, Double_t& phi) const
{
  GetBinCenter(bin, &theta, &phi);
}

//______________________________________________________________________________
void THealPix::GetBinCenter(Int_t bin, Double_t* theta, Double_t* phi) const
{
  // Oribinal code is Healpix_Base::pix2ang of HEALPix C++
  // See Gorski et al. ApJ 622 (2005) for equation numbering

  if(fIsNested){
    bin = Nest2Ring(bin);
  } // if

  if(bin < fNcap){ // North Polar cap
    Int_t i = Int_t(.5*(1 + ::Isqrt(1 + 2*bin))); //counted from North pole
    Int_t j = (bin + 1) - 2*i*(i - 1); // (3)
    *theta = TMath::ACos(1.0 - (i*i)/f3Nside2); // (4)
    *phi   = (j - .5)*TMath::Pi()/(2.*i); // (5)
  } else if(bin < (fNpix - fNcap)){ // Equatorial region
    Int_t ip = bin - fNcap;
    Int_t i = ip/f4Nside + fNside; // (6) counted from North pole
    Int_t j = ip%f4Nside + 1;      // (7)
    // 1 if i+Nside is odd, 1/2 otherwise
    Double_t fodd = ((i + fNside)&1) ? 1 : 0.5;
    
    *theta = TMath::ACos(4./3. - i*f2over3Nside); // (8)
    *phi   = (j - fodd) * TMath::Pi()/f2Nside;    // (9)
  } else { // South Polar cap
    Int_t ip = fNpix - bin;
    Int_t i = Int_t(0.5*(1 + ::Isqrt(2*ip - 1))); //counted from South pole
    Int_t j = 4*i + 1 - (ip - 2*i*(i - 1)); //(3)
    
    *theta = TMath::ACos(-1. + (i*i)/f3Nside2); // (4)
    *phi   = (j - .5)*TMath::PiOver2()/i;       // (5)
  } // if

  if(fIsDegree){
    *theta *= TMath::RadToDeg();
    *phi   *= TMath::RadToDeg();
  } // if
}

//______________________________________________________________________________
Double_t THealPix::GetBinContent(Int_t) const
{
  AbstractMethod("GetBinContent");
  return 0;
}

//______________________________________________________________________________
Int_t THealPix::GetBinVertices(Int_t bin, Double_t* x, Double_t* y) const
{
  if(fIsNested){
    bin = Nest2Ring(bin);
  } // if

  if(bin < fNcap){ // North polar cap
    Int_t i = Int_t(.5*(1 + ::Isqrt(1 + 2*bin)));
    Int_t j = (bin + 1) - 2*i*(i - 1); // (3)
    Double_t theta0 = TMath::RadToDeg()*TMath::ACos(1 - (i*i)/f3Nside2); // (4)
    Double_t theta1 = TMath::RadToDeg()*TMath::ACos(1 - ((i+1)*(i+1))/f3Nside2);
    Double_t theta2 = TMath::RadToDeg()*TMath::ACos(1 - ((i-1)*(i-1))/f3Nside2);
    //    Double_t phi0 = (j -  .5)*90./i; // (5)
    Double_t phi1 = (j     )*90./i; // eastward
    Double_t phi2 = (j - 1.)*90./i; // westward
    Double_t phi3 = (j + (j - 1)/i)*90./(i + 1); // southward
    if(i != 1){
      Double_t phi4 = (j - (j - 1)/i - 1. )*90./(i - 1); // northward
      x[0] = phi3; y[0] = theta1;
      x[1] = phi1; y[1] = theta0;
      x[2] = phi4; y[2] = theta2;
      x[3] = phi2; y[3] = theta0;
      return 4;
    } else {
      x[0] = phi3; y[0] = theta1;
      x[1] = phi1; y[1] = theta0;
      x[2] = phi1; y[2] = theta2;
      x[3] = phi2; y[3] = theta2;
      x[4] = phi2; y[4] = theta0;
      return 5;
    } // if
  } else if(bin < fNcap + f4Nside){
    Int_t i = Int_t(.5*(1 + ::Isqrt(1 + 2*bin)));
    Int_t j = (bin + 1) - 2*i*(i - 1); // (3)
    Double_t theta0 = TMath::RadToDeg()*TMath::ACos(1 - (i*i)/f3Nside2); // (4)
    Double_t theta1 = TMath::RadToDeg()*TMath::ACos(1 - ((i+1)*(i+1))/f3Nside2);
    Double_t theta2 = TMath::RadToDeg()*TMath::ACos(1 - ((i-1)*(i-1))/f3Nside2);
    Double_t phi0 = (j -  .5)*90./i; // (5)
    Double_t phi1 = (j     )*90./i; // eastward
    Double_t phi2 = (j - 1.)*90./i; // westward
    Double_t phi3 = (j - (j - 1)/i - 1. )*90./(i - 1); // northward
    x[0] = phi0; y[0] = theta1;
    x[1] = phi1; y[1] = theta0;
    x[2] = phi3; y[2] = theta2;
    x[3] = phi2; y[3] = theta0;
    return 4;
  } else if(bin < fNpix - fNcap - f4Nside){ // Eqatorial region
    Int_t pdash = bin - fNcap;
    Int_t i = pdash/f4Nside + fNside; // (6)
    Int_t j = pdash%f4Nside + 1;      // (7)
    Double_t z = 4./3. - i*f2over3Nside; // (8)
    Double_t theta0 = TMath::RadToDeg()*TMath::ACos(z);
    Double_t theta1 = TMath::RadToDeg()*TMath::ACos(z + f2over3Nside);
    Double_t theta2 = TMath::RadToDeg()*TMath::ACos(z - f2over3Nside);

    Double_t fodd = ((i + fNside)&1) ? 1 : 0.5;
    Double_t phi0 = (j - fodd)*180./f2Nside; // (9)
    Double_t phi1 = phi0 + 180./f4Nside;
    Double_t phi2 = phi0 - 180./f4Nside;
    x[0] = phi0; y[0] = theta1; // counter clock wise from south
    x[1] = phi1; y[1] = theta0;
    x[2] = phi0; y[2] = theta2;
    x[3] = phi2; y[3] = theta0;
    return 4;
  } else if(bin < fNpix - fNcap){
    Int_t ip = fNpix - bin;
    Int_t i = Int_t(.5*(1 + ::Isqrt(2*ip - 1)));
    Int_t j = 4*i + 1 - (ip - 2*i*(i - 1)); // (3)
    Double_t theta0 = TMath::RadToDeg()*TMath::ACos(-1 + (i*i)/f3Nside2); // (4)
    Double_t theta1 = TMath::RadToDeg()*TMath::ACos(-1 + ((i+1)*(i+1))/f3Nside2);
    Double_t theta2 = TMath::RadToDeg()*TMath::ACos(-1 + ((i-1)*(i-1))/f3Nside2);
    Double_t phi0 = (j -  .5)*90./i; // (5)
    Double_t phi1 = (j     )*90./i; // eastward
    Double_t phi2 = (j - 1.)*90./i; // westward
    Double_t phi3 = (j - (j - 1)/i - 1. )*90./(i - 1); // northward
    x[0] = phi0; y[0] = theta1;
    x[1] = phi1; y[1] = theta0;
    x[2] = phi3; y[2] = theta2;
    x[3] = phi2; y[3] = theta0;
    return 4;
  } else { // South polar cap
    Int_t ip = fNpix - bin;
    Int_t i = Int_t(.5*(1 + ::Isqrt(2*ip - 1)));
    Int_t j = 4*i + 1 - (ip - 2*i*(i - 1)); // (3)
    Double_t theta0 = TMath::RadToDeg()*TMath::ACos(-1 + (i*i)/f3Nside2); // (4)
    Double_t theta1 = TMath::RadToDeg()*TMath::ACos(-1 + ((i+1)*(i+1))/f3Nside2);
    Double_t theta2 = TMath::RadToDeg()*TMath::ACos(-1 + ((i-1)*(i-1))/f3Nside2);
    //    Double_t phi0 = (j -  .5)*90./i; // (5)
    Double_t phi1 = (j     )*90./i; // eastward
    Double_t phi2 = (j - 1.)*90./i; // westward
    Double_t phi3 = (j + (j - 1)/i)*90./(i + 1); // southward
    if(i != 1){
      Double_t phi4 = (j - (j - 1)/i - 1. )*90./(i - 1); // northward
      x[0] = phi3; y[0] = theta1;
      x[1] = phi1; y[1] = theta0;
      x[2] = phi4; y[2] = theta2;
      x[3] = phi2; y[3] = theta0;
      return 4;
    } else {
      x[0] = phi3; y[0] = theta1;
      x[1] = phi1; y[1] = theta0;
      x[2] = phi1; y[2] = theta2;
      x[3] = phi2; y[3] = theta2;
      x[4] = phi2; y[4] = theta0;
      return 5;
    } // if
  } // if
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

  Int_t northring = (ring > f2Nside) ? f4Nside - ring : ring;

  if(northring < fNside){
    costheta = 1 - northring*northring*4./fNpix;
    sintheta = TMath::Sin(2*TMath::ASin(northring/(TMath::Sqrt(6.)*fNside)));
    ringpix = 4*northring;
    shifted = kTRUE;
    startpix = 2*northring*(northring - 1);
  } else {
    costheta = (f2Nside - northring)*(8.*fNside/fNpix);
    sintheta = TMath::Sqrt(1 - costheta*costheta);
    ringpix = f4Nside;
    shifted = ((northring - fNside) & 1) == 0;
    startpix = fNcap + (northring - fNside)*ringpix;
  } // if

  if(northring != ring){ // southern hemisphere
    costheta = -costheta;
    startpix = fNpix - startpix - ringpix;
  } // if
}

//______________________________________________________________________________
Double_t THealPix::GetSumOfWeights() const
{
  //   -*-*-*-*-*-*Return the sum of weights excluding under/overflows*-*-*-*-*
  //               ===================================================
  Double_t sum =0;
  for(Int_t i = 0; i <= fNpix; i++) {
    sum += GetBinContent(i);
  } // i
  return sum;
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
    if(opt.Contains("same")){
      if(gPad->GetX1() == 0 && gPad->GetX2() == 1 &&
	 gPad->GetY1() == 0 && gPad->GetY2() == 1 &&
	 gPad->GetListOfPrimitives()->GetSize() == 0){
	opt.ReplaceAll("same", "");
      } // if
    } else {
      //the following statement is necessary in case one attempts to draw
      //a temporary histogram already in the current pad
      if(TestBit(kCanDelete)){
	gPad->GetListOfPrimitives()->Remove(this);
      } // if
      gPad->Clear();
    } // if
  } else {
    if(opt.Contains("same")){
      opt.ReplaceAll("same", "");
    } // if
  }
  AppendPad(option);
}

//_____________________________________________________________________________
template<typename T>
void THealPix::Map2Alm(THealAlm<T>& alm, Bool_t add) const
{
  std::vector<Double_t> weight(2*fNside, 1);
  Map2Alm(alm, weight, add);
}

//_____________________________________________________________________________
template<typename T>
void THealPix::Map2Alm(THealAlm<T>& alm, const std::vector<double>& weight,
		       Bool_t add) const
{
  if(fIsNested){
    // to be modified
    return;
  } // if

  if((Int_t)weight.size() < f2Nside){
    // to be modified
    return;
  } // if

  int lmax = alm.GetLmax(), mmax = alm.GetMmax();

  int nchunks, chunksize;
  get_chunk_info(f2Nside,nchunks,chunksize);

  std::complex<double> *phas_n[chunksize], *phas_s[chunksize];
  for(Int_t i = 0; i < chunksize; i++){
    phas_n[i] = new std::complex<double>[mmax+1];
    phas_s[i] = new std::complex<double>[mmax+1];
  } // i

  std::vector<double> cth(chunksize), sth(chunksize);
  double normfact = TMath::Pi()/(3*fNside2);

  if (!add) alm.SetToZero();

  for(int chunk=0; chunk<nchunks; ++chunk){
    int llim=chunk*chunksize, ulim=TMath::Min(llim+chunksize,f2Nside);
    std::vector<std::complex<double> > shiftarr(mmax+1), work(f4Nside);
    THealFFT::rfft plan;

    for (int ith=llim; ith<ulim; ++ith){
      int istart_north, istart_south, nph;
      bool shifted;
      GetRingInfo (ith+1,istart_north,nph,cth[ith-llim],sth[ith-llim],shifted);
      istart_south = fNpix - istart_north - nph;
      recalc_map2alm (nph, mmax, plan, shiftarr);
      if(GetTypeString() == "D"){
	const THealPixD* tmp = dynamic_cast<THealPixD*>(this);
	const Double_t* p1 = &(tmp->GetArray()[istart_north]);
	const Double_t* p2 = &(tmp->GetArray()[istart_south]);
	fft_map2alm (nph, mmax, shifted, weight[ith]*normfact, plan,
		     p1, p2,	phas_n[ith-llim], phas_s[ith-llim], shiftarr, work);
      } else if(GetTypeString() == "F"){
	const THealPixF* tmp = dynamic_cast<THealPixF*>(this);
	const Float_t* p1 = &(tmp->GetArray()[istart_north]);
	const Float_t* p2 = &(tmp->GetArray()[istart_south]);
	fft_map2alm (nph, mmax, shifted, weight[ith]*normfact, plan,
		     p1, p2,	phas_n[ith-llim], phas_s[ith-llim], shiftarr, work);
      } // if
    } // ith
    
    THealFFT::Ylmgen generator(lmax,mmax,1e-30);
    std::vector<double> Ylm;
    std::vector<std::complex<double> > alm_tmp(lmax+1);
    for (int m=0; m<=mmax; ++m)
      {
      for (int l=m; l<=lmax; ++l) alm_tmp[l] = std::complex<double>(0.,0.);
      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if (l<=lmax)
          {
	   std::complex<double> p1 = phas_n[ith][m]+phas_s[ith][m],
                                p2 = phas_n[ith][m]-phas_s[ith][m];
          if ((l-m)&1) goto middle;
start:    alm_tmp[l].real() += p1.real()*Ylm[l]; alm_tmp[l].imag() += p1.imag()*Ylm[l];
          if (++l>lmax) goto end;
middle:   alm_tmp[l].real() += p2.real()*Ylm[l]; alm_tmp[l].imag() += p2.imag()*Ylm[l];
          if (++l<=lmax) goto start;
end:      ;
          }
        }
      std::complex<T> *palm = alm.GetMstart(m);
      for (int l=m; l<=lmax; ++l)
        { palm[l].real() += alm_tmp[l].real(); palm[l].imag() += alm_tmp[l].imag(); }
      }
    }

  for(Int_t i = 0; i < chunksize; i++){
    delete phas_n[i];
    delete phas_s[i];
  } // i
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
void THealPix::Paint(Option_t* option)
{
  // Control routine to paint any kind of HEALPixs
  GetPainter(option);

  if(fPainter){
    if(strlen(option) > 0){
      fPainter->Paint(option);
    } else {
      fPainter->Paint(fOption.Data());
    } // if
  } // if
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

  char ttype[FLEN_VALUE], tform[FLEN_VALUE], tdisp[FLEN_VALUE];
  Long_t repeat, nulval;
  Double_t scale, zero;
    
  // For example,
  // TTYPE1  = 'diffuse'
  // TFORM1  = '1024D'
  // gives ttype = 'diffuse', repeat = 1024 and tform = 'D'
  fits_get_bcolparms(*fptr, colnum, ttype, head.tunit, tform, &repeat,
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

  Int_t npix = Nside2Npix(nside);

  if(!(nrows*repeat == npix || nrows == npix || repeat >= 1)){
    fits_close_file(*fptr, &status);
    std::cerr << Form("Incompatible type of header (Npix = %d, Nrows = %d, Repeat = %d).\n", npix, nrows, repeat);
    return kFALSE;
  } // if

  head.isNested = isNested;
  head.order    = order;
  head.nrows    = nrows;
  head.repeat   = repeat;
  head.colnum   = colnum;
  head.npix     = npix;
  
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

  // If overwrite this, clone this to keep original information
  THealPix* hpold = (hpnew == this) ? (THealPix*)Clone() : this;

  hpnew->SetOrder(neworder);
  hpnew->SetBins(hpnew->GetNpix());

  Int_t newnpix  = hpnew->GetNpix();
  if(hpnew->IsNested()){ // NESTED
    if(neworder < hpold->GetOrder()){ // rebin to lower level
      Int_t div = hpold->GetNpix()/newnpix;
      for(Int_t i = 0; i < newnpix; i++){
	Double_t content = 0;
	Double_t error   = 0;
	for(Int_t j = i*div; j < (i + 1)*div; j++){
	  content += hpold->GetBinContent(j);
	  if(hpold->GetSumw2N()){
	    error += hpold->GetBinError(j)*hpold->GetBinError(j);
	  } // if
	} // j
	hpnew->SetBinContent(i, content);
	if(hpold->GetSumw2N()){
	  hpnew->SetBinError(i, TMath::Sqrt(error));
	} // if
      } // i
    } else { // rebin to higher level
      Int_t div = newnpix/hpold->GetNpix();
      for(Int_t i = 0; i < newnpix; i++){
	hpnew->SetBinContent(i, hpold->GetBinContent(i/div)/div);
	if(hpold->GetSumw2N()){
	  hpnew->SetBinError(i, hpold->GetBinError(i/div)/TMath::Sqrt(div));
	} // if
      } // i
    } // if
  } else { // RING
    if(neworder < hpold->GetOrder()){ // rebin to lower level
      Int_t div = hpold->GetNpix()/newnpix;
      for(Int_t i = 0; i < newnpix; i++){
	Double_t content = 0;
	Double_t error   = 0;
	Int_t inest = hpnew->Ring2Nest(i);
	for(Int_t jnest = inest*div; jnest < (inest + 1)*div; jnest++){
	  Int_t j = hpold->Nest2Ring(jnest);
	  content += hpold->GetBinContent(j);
	  if(hpold->GetSumw2N()){
	    error += hpold->GetBinError(j)*hpold->GetBinError(j);
	  } // if
	} // j
	hpnew->SetBinContent(i, content);
	if(hpold->GetSumw2N()){
	  hpnew->SetBinError(i, TMath::Sqrt(error));
	} // if
      } // i
    } else { // rebin to higher level
      Int_t div = newnpix/hpold->GetNpix();
      for(Int_t i = 0; i < newnpix; i++){
	Int_t inest = hpnew->Ring2Nest(i);
	Int_t p = hpold->Nest2Ring(inest/div);
	hpnew->SetBinContent(i, hpold->GetBinContent(p)/div);
	if(hpold->GetSumw2N()){
	  hpnew->SetBinError(i, hpold->GetBinError(p)/TMath::Sqrt(div));
	} // if
      } // i
    } // if
  } // if

  hpnew->SetEntries(hpold->GetEntries()); //was modified by SetBinContent

  if(hpold != this){
    delete hpold;
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

  fOrder   = order;
  fNside   = Order2Nside(fOrder);
  fNpix    = Nside2Npix(fNside);
  fNside2  = fNside*fNside;
  f2Nside  = 2*fNside;
  fNcap    = f2Nside*(fNside - 1);
  f3Nside2 = 3*fNside2;
  f2over3Nside = 2./(3.*fNside);
  f4Nside  = 4*fNside;

}

//_____________________________________________________________________________
void THealPix::SetUnit(const char* unit)
{
  fUnit = std::string(unit);
}

//_____________________________________________________________________________
void THealPix::Streamer(TBuffer& b)
{
  if(b.IsReading()){
    UInt_t R__s, R__c;
    Version_t R__v = b.ReadVersion(&R__s, &R__c);
    fDirectory = 0;
    b.ReadClassBuffer(THealPix::Class(), this, R__v, R__s, R__c);
    SetOrder(fOrder);
  } else {
    b.WriteClassBuffer(THealPix::Class(), this);
  }
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
void THealPix::UseCurrentStyle()
{
  if(gStyle->IsReading()) {
    fXaxis.ResetAttAxis("X");
    fYaxis.ResetAttAxis("Y");
    fZaxis.ResetAttAxis("Z");
    SetBarOffset(gStyle->GetBarOffset());
    SetBarWidth(gStyle->GetBarWidth());
    SetFillColor(gStyle->GetHistFillColor());
    SetFillStyle(gStyle->GetHistFillStyle());
    SetLineColor(gStyle->GetHistLineColor());
    SetLineStyle(gStyle->GetHistLineStyle());
    SetLineWidth(gStyle->GetHistLineWidth());
    SetMarkerColor(gStyle->GetMarkerColor());
    SetMarkerStyle(gStyle->GetMarkerStyle());
    SetMarkerSize(gStyle->GetMarkerSize());
    Int_t dostat = gStyle->GetOptStat();
    if (gStyle->GetOptFit() && !dostat) dostat = 1000000001;
    //    SetStats(dostat);
  } else {
    gStyle->SetBarOffset(fBarOffset);
    gStyle->SetBarWidth(fBarWidth);
    gStyle->SetHistFillColor(GetFillColor());
    gStyle->SetHistFillStyle(GetFillStyle());
    gStyle->SetHistLineColor(GetLineColor());
    gStyle->SetHistLineStyle(GetLineStyle());
    gStyle->SetHistLineWidth(GetLineWidth());
    gStyle->SetMarkerColor(GetMarkerColor());
    gStyle->SetMarkerStyle(GetMarkerStyle());
    gStyle->SetMarkerSize(GetMarkerSize());
    gStyle->SetOptStat(TestBit(kNoStats));
  } // if

  TIter next(GetListOfFunctions());
  TObject* obj;
  
  while((obj = next())){
    obj->UseCurrentStyle();
  } // while
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

//______________________________________________________________________________
Int_t THealPix::GetContour(Double_t* levels)
{
  //  Return contour values into array levels if pointer levels is non zero
  //
  //  The function returns the number of contour levels.
  //  see GetContourLevel to return one contour only
  //
  
  Int_t nlevels = fContour.fN;
  if(levels){
    if(nlevels == 0){
      nlevels = 20;
      SetContour(nlevels);
    } else {
      if(TestBit(kUserContour) == 0) SetContour(nlevels);
    } // if
    for(Int_t i = 0; i < nlevels; i++){
      levels[i] = fContour.fArray[i];
    } // i
  } // if
  return nlevels;
}

//______________________________________________________________________________
Double_t THealPix::GetContourLevel(Int_t level) const
{
  // Return value of contour number level
  // see GetContour to return the array of all contour levels
  
  if(level <0 || level >= fContour.fN) return 0;
  Double_t zlevel = fContour.fArray[level];
  return zlevel;
}

//______________________________________________________________________________
Double_t THealPix::GetContourLevelPad(Int_t level) const
{
  // Return the value of contour number "level" in Pad coordinates ie: if the Pad
  // is in log scale along Z it returns le log of the contour level value.
  // see GetContour to return the array of all contour levels
  
  if(level <0 || level >= fContour.fN) return 0;
  Double_t zlevel = fContour.fArray[level];
  
  // In case of user defined contours and Pad in log scale along Z,
  // fContour.fArray doesn't contain the log of the contour whereas it does
  // in case of equidistant contours.
  if(gPad && gPad->GetLogz() && TestBit(kUserContour)){
    if(zlevel <= 0) return 0;
    zlevel = TMath::Log10(zlevel);
  } // if
  return zlevel;
}

//_____________________________________________________________________________
Double_t THealPix::GetMaximum(Double_t maxval) const
{
  if(fMaximum != -1111){
    return fMaximum;
  } // if

  Double_t maximum = -FLT_MAX;
  for(Int_t i = 0; i < fNpix; i++){
    Double_t value = GetBinContent(i);
    if(value > maximum && value < maxval){
      maximum = value;
    } // if
  } // i

  return maximum;
}

//_____________________________________________________________________________
Int_t THealPix::GetMaximumBin() const
{
  Int_t bin = 0;
  Double_t maximum = -FLT_MAX;
  for(Int_t i = 0; i < fNpix; i++){
    Double_t value = GetBinContent(i);
    if(value > maximum){
      maximum = value;
      bin = i;
    } // if
  } // i

  return bin;
}

//_____________________________________________________________________________
Double_t THealPix::GetMinimum(Double_t minval) const
{
  if(fMinimum != -1111){
    return fMinimum;
  } // if

  Double_t minimum = FLT_MAX;
  for(Int_t i = 0; i < fNpix; i++){
    Double_t value = GetBinContent(i);
    if(value < minimum && value > minval){
      minimum = value;
    } // if
  } // i

  return minimum;
}

//_____________________________________________________________________________
Int_t THealPix::GetMinimumBin() const
{
  Int_t bin = 0;
  Double_t minimum = FLT_MAX;
  for(Int_t i = 0; i < fNpix; i++){
    Double_t value = GetBinContent(i);
    if(value < minimum){
      minimum = value;
      bin = i;
    } // if
  } // i

  return bin;
}

//_____________________________________________________________________________
Int_t THealPix::GetNrows() const
{
  return fOrder < 4 ? 1 : 1024;
}

//_____________________________________________________________________________
TVirtualHealPainter* THealPix::GetPainter(Option_t* option)
{
  // return pointer to patiner
  // if painter does not exist, it is created
  if(!fPainter){
    TString opt = option;
    opt.ToLower();
    if(opt.Contains("gl") || gStyle->GetCanvasPreferGL()){
      fPainter = 0; // to be modified
    } // if
  } // if

  if(!fPainter){
    fPainter = TVirtualHealPainter::HealPainter(this);
  } // if

  return fPainter;
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

//______________________________________________________________________________
TAxis* THealPix::GetXaxis() const
{
  // return a pointer to the X axis object
  
  return &((THealPix*)this)->fXaxis;
}

//______________________________________________________________________________
TAxis* THealPix::GetYaxis() const
{
  // return a pointer to the Y axis object
  
  return &((THealPix*)this)->fYaxis;
}

//______________________________________________________________________________
TAxis* THealPix::GetZaxis() const
{
  // return a pointer to the Z axis object
  
  return &((THealPix*)this)->fZaxis;
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

  //if contours set, must also scale contours
  Int_t ncontours = GetContour();
  if(ncontours == 0) return;
  Double_t* levels = fContour.GetArray();
  for(Int_t i = 0; i < ncontours; i++){
    levels[i] *= c1;
  } // i
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
void THealPix::SetBins(Int_t n)
{
  SetBinsLength(n);
  if(fSumw2.fN){
    fSumw2.Set(n);
  } // if
}

//______________________________________________________________________________
void THealPix::SetContour(Int_t  nlevels, const Double_t* levels)
{
  //   -*-*-*-*-*-*Set the number and values of contour levels*-*-*-*-*-*-*-*-*
  //               ===========================================
  //
  //  By default the number of contour levels is set to 20.
  //
  //  if argument levels = 0 or missing, equidistant contours are computed
  //
  
  ResetBit(kUserContour);
  if(nlevels <= 0){
    fContour.Set(0);
    return;
  } // if
  fContour.Set(nlevels);
  
  //   -  Contour levels are specified
  if(levels){
    SetBit(kUserContour);
    for(Int_t i=0; i < nlevels; i++){
      fContour.fArray[i] = levels[i];
    } // i
  } else {
    //   - contour levels are computed automatically as equidistant contours
    Double_t zmin = GetMinimum();
    Double_t zmax = GetMaximum();
    if((zmin == zmax) && (zmin != 0)){
      zmax += 0.01*TMath::Abs(zmax);
      zmin -= 0.01*TMath::Abs(zmin);
    } // if
    Double_t dz = (zmax - zmin)/Double_t(nlevels);
    if(gPad && gPad->GetLogz()){
      if (zmax <= 0) return;
      if (zmin <= 0) zmin = 0.001*zmax;
      zmin = TMath::Log10(zmin);
      zmax = TMath::Log10(zmax);
      dz   = (zmax-zmin)/Double_t(nlevels);
    } // if
    for(Int_t i=0; i < nlevels; i++) {
      fContour.fArray[i] = zmin + dz*Double_t(i);
    } // i
  } // if
}

//______________________________________________________________________________
void THealPix::SetContourLevel(Int_t level, Double_t value)
{
  //   -*-*-*-*-*-*-*-*-*Set value for one contour level*-*-*-*-*-*-*-*-*-*-*-*
  //                     ===============================
  if(level <0 || level >= fContour.fN) return;
  SetBit(kUserContour);
  fContour.fArray[level] = value;
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

//______________________________________________________________________________
void THealPix::SetMaximum(Double_t maximum)
{
  fMaximum = maximum;
}

//______________________________________________________________________________
void THealPix::SetMinimum(Double_t minimum)
{
  fMinimum = minimum;
}

//______________________________________________________________________________
void THealPix::SetName(const char* name)
{
  if(fDirectory){
    fDirectory->Remove(this);
  } // if
  fName = name;
  if(fDirectory){
    fDirectory->Append(this);
  } // if
}

//______________________________________________________________________________
void THealPix::SetNameTitle(const char* name, const char* title)
{
  if(fDirectory){
    fDirectory->Remove(this);
  } // if
  fName  = name;
  SetTitle(title);
  if(fDirectory){
    fDirectory->Append(this);
  } // if
}

//_____________________________________________________________________________
void THealPix::Ring2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const
{
  // Original code is Healpix_Base::ring2xyf of HEALPix C++
  Int_t iring, iphi, kshift, nr;
  
  if(pix < fNcap){ // North Polar cap
    iring = Int_t(0.5*(1 + ::Isqrt(1 + 2*pix))); //counted from North pole
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
      iphi  = (ip&(f4Nside - 1)) + 1;
    } else {
      iring = (ip/(f4Nside)) + fNside; // counted from North pole
      iphi = (ip%(f4Nside)) + 1;
    } // if
    kshift = (iring + fNside)&1;
    nr = fNside;
    UInt_t ire = iring - fNside + 1;
    UInt_t irm = f2Nside + 2 - ire;
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
    iring = Int_t(0.5*(1 + ::Isqrt(2*ip - 1))); //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring - 1));
    kshift = 0;
    nr = iring;
    iring = f4Nside - iring;
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
  if(ipt >= f2Nside){
    ipt -= 8*fNside;
  } // if

  x =  (ipt-irt) >>1;
  y =(-(ipt+irt))>>1;
}

//_____________________________________________________________________________
Int_t THealPix::XYF2Ring(Int_t x, Int_t y, Int_t face) const
{
  // Original code is Healpix_Base::xyf2ring of HEALPix C++
  Int_t jr = (fgJrll[face]*fNside) - x - y  - 1;

  Int_t nr, kshift, n_before;
  if(jr < fNside){
    nr = jr;
    n_before = 2*nr*(nr - 1);
    kshift = 0;
  } else if (jr > 3*fNside){
    nr = f4Nside - jr;
    n_before = fNpix - 2*(nr + 1)*nr;
    kshift = 0;
  } else {
    nr = fNside;
    n_before = fNcap + (jr - fNside)*f4Nside;
    kshift = (jr - fNside)&1;
  } // if

  Int_t jp = (fgJpll[face]*nr + x - y + 1 + kshift)/2;
  if(jp > f4Nside){
    jp -= f4Nside;
  } else {
    if(jp < 1){
      jp += f4Nside;
    } // if
  } // if

  return n_before + jp - 1;
}

ClassImp(THealPixF)

//_____________________________________________________________________________
// THealPixF methods
//_____________________________________________________________________________
THealPixF::THealPixF(): THealPix(), TArrayF()
{
  SetBinsLength(fNpix);
  if(fgDefaultSumw2){
    Sumw2();
  } // if
}

//_____________________________________________________________________________
THealPixF::THealPixF(const char* name, const char* title, Int_t order,
		     Bool_t isNested)
: THealPix(name, title, order, isNested)
{
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

//_____________________________________________________________________________
Int_t THealPixF::GetType() const
{
  return TFLOAT;
}

//_____________________________________________________________________________
std::string THealPixF::GetTypeString() const
{
  return "F";
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

  Int_t status = 0;

  if(head.nrows*head.repeat == hpf->GetNpix()){ // CMB-like FITS
    for(Int_t i = 0; i < head.nrows; i++){
      fits_read_col(fptr, TFLOAT, head.colnum, i + 1, 1, head.repeat, 0,
		    &(hpf->GetArray()[i*head.repeat]), 0, &status);
      if(!THealUtil::FitsReportError(status)){
	delete hpf;
	return 0;
      } // if
    } // i
  } else if(head.nrows == hpf->GetNpix()){ // HEALPix Cube
    // Read only the first layer
    for(Int_t i = 0; i < head.nrows; i++){
      fits_read_col(fptr, TFLOAT, head.colnum, i + 1, 1, 1, 0,
		    &(hpf->GetArray()[i]), 0, &status);
      if(!THealUtil::FitsReportError(status)){
	delete hpf;
	return 0;
      } // if
    } // i
  } // if

  Double_t entries = 0;
  for(Int_t i = 0; i  < hpf->GetNpix(); i++){
    if(hpf->GetBinContent(i) != 0){
      entries++;
    } // if
  } // i
  hpf->SetEntries(entries);

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
  SetBinsLength(fNpix);
  if(fgDefaultSumw2){
    Sumw2();
  } // if
}

//_____________________________________________________________________________
THealPixD::THealPixD(const char* name, const char* title, Int_t order,
		     Bool_t isNested)
: THealPix(name, title, order, isNested)
{
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

//_____________________________________________________________________________
Int_t THealPixD::GetType() const
{
  return TDOUBLE;
}

//_____________________________________________________________________________
std::string THealPixD::GetTypeString() const
{
  return "D";
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

  Int_t status = 0;

  if(head.nrows*head.repeat == hpd->GetNpix()){ // CMB-like FITS
    for(Int_t i = 0; i < head.nrows; i++){
      fits_read_col(fptr, TDOUBLE, head.colnum, i + 1, 1, head.repeat, 0,
		    &(hpd->GetArray()[i*head.repeat]), 0, &status);
      if(!THealUtil::FitsReportError(status)){
	delete hpd;
	return 0;
      } // if
    } // i
  } else if(head.nrows == hpd->GetNpix()){ // HEALPix Cube
    // Read only the first layer
    for(Int_t i = 0; i < head.nrows; i++){
      fits_read_col(fptr, TDOUBLE, head.colnum, i + 1, 1, 1, 0,
		    &(hpd->GetArray()[i]), 0, &status);
      if(!THealUtil::FitsReportError(status)){
	delete hpd;
	return 0;
      } // if
    } // i
  } // if

  Double_t entries = 0;
  for(Int_t i = 0; i  < hpd->GetNpix(); i++){
    if(hpd->GetBinContent(i) != 0){
      entries++;
    } // if
  } // i
  hpd->SetEntries(entries);

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

