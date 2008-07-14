// $Id: THealAlm.cxx,v 1.9 2008/07/14 05:03:16 oxon Exp $
// Author: Akira Okumura 2008/06/26

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.

   This is a port of HEALPix C++ package to ROOT system.
   Original code is available at <http://healpix.jpl.nasa.gov> under GPL.
******************************************************************************/

#include "TMath.h"

#include "THealAlm.h"
#include "THealPix.h"

namespace {
  /*
  void get_chunk_info (int nrings, int &nchunks, int &chunksize)
  {
    nchunks = nrings/TMath::Max(100,nrings/10) + 1;
    chunksize = (nrings+nchunks-1)/nchunks;
  }
  */
//______________________________________________________________________________
void recalc_alm2map (int nph, int mmax, THealFFT::rfft &plan,
  std::vector<std::complex<double> > &shiftarr)
  {
  if (plan.size() == nph) return;
  plan.Set (nph);
  double f1 = TMath::Pi()/nph;
  for (int m=0; m<=mmax; ++m)
    {
    if (m<nph)
      shiftarr[m] = std::complex<double>(cos(m*f1),sin(m*f1));
    else
      shiftarr[m]=-shiftarr[m-nph];
    }
  }
  /*
//______________________________________________________________________________
void fill_work (const std::complex<double> *datain, int nph, int mmax,
  bool shifted, const std::vector<std::complex<double> > &shiftarr,
  std::vector<std::complex<double> > &work)
  {
  for (int m=1; m<nph; ++m) work[m]=0;
  work[0]=datain[0];

  int cnt1=0, cnt2=nph;
  for (int m=1; m<=mmax; ++m)
    {
    if (++cnt1==nph) cnt1=0;
    if (--cnt2==-1) cnt2=nph-1;
    std::complex<double> tmp = shifted ? (datain[m]*shiftarr[m]) : datain[m];
    work[cnt1] += tmp;
    work[cnt2] += conj(tmp);
    }
  }

//______________________________________________________________________________
template<typename T> void fft_alm2map (int nph, int mmax, bool shifted,
  THealFFT::rfft &plan, T *mapN, T *mapS, std::complex<double> *b_north,
  const std::complex<double> *b_south, const std::vector<std::complex<double> > &shiftarr,
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
  */
}

templateClassImp(THealAlm)

template<typename T>
std::complex<T> THealAlm<T>::fgAlmDummy = std::complex<T>(0, 0);
template<typename T>
const std::complex<T> THealAlm<T>::fgAlmDummyConst = std::complex<T>(0, 0);

//_____________________________________________________________________________
template<typename T>
THealAlm<T>::~THealAlm()
{
}

//_____________________________________________________________________________
template<typename T>
THealAlm<T>::THealAlm(Int_t lmax, Int_t mmax)
  : TObject()
{
  ResetLM(lmax, mmax);
}

//_____________________________________________________________________________
template<typename T>
THealAlm<T>::THealAlm(const THealAlm<T>& ha) : TObject()
{
  ((THealAlm<T>&)ha).Copy(*this);
}

//______________________________________________________________________________
template<typename T>
void THealAlm<T>::Add(const THealAlm<T>* ha1, Double_t c1)
{
  if(!ha1){
    Error("Add", "Attempt to add a non-existing Alm");
    return;
  } // if

  if(fLmax != ha1->GetLmax() || fMmax != ha1->GetMmax()){
    Error("Add", "Attempt to add Alms with different number of lmax or mmax");
    return;
  } // if

  for(Int_t l = 0; l <= fLmax; l++){
    for(Int_t m = 0; m <= TMath::Min(l, fMmax); m++){
      std::complex<T> comp = ((*ha1)(l, m));
      comp *= c1;
      (*this)(l, m) += comp;
    } // m
  } // l
}

//______________________________________________________________________________
template<typename T>
void THealAlm<T>::Add(const THealAlm<T>* ha1, const THealAlm<T>* ha2, Double_t c1, Double_t c2)
{
  if(!ha1 || !ha2){
    Error("Add", "Attempt to add a non-existing Alm");
    return;
  } // if

  if(fLmax != ha1->GetLmax() || fMmax != ha1->GetMmax()
  || fLmax != ha2->GetLmax() || fMmax != ha2->GetMmax()){
    Error("Add", "Attempt to add Alms with different number of lmax or mmax");
    return;
  } // if

  for(Int_t l = 0; l <= fLmax; l++){
    for(Int_t m = 0; m <= TMath::Min(l, fMmax); m++){
      std::complex<T> comp1 = ((*ha1)(l, m));
      std::complex<T> comp2 = ((*ha2)(l, m));
      comp1 *= c1;
      comp2 *= c2;
      (*this)(l, m) = comp1 + comp2;
    } // m
  } // l
}

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::Alm2Map(THealPix& map) const
{
  int nside = map.GetNside();
  
  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);
  
  std::complex<double> *b_north[chunksize], *b_south[chunksize];
  for(Int_t i = 0; i < chunksize; i++){
    b_north[i] = new std::complex<double>[fMmax+1];
    b_south[i] = new std::complex<double>[fMmax+1];
  } // i
  std::vector<double> cth(chunksize),sth(chunksize);
  
  for(int chunk=0; chunk<nchunks; ++chunk){
    int llim=chunk*chunksize, ulim=TMath::Min(llim+chunksize,2*nside);
    for(int ith=llim; ith<ulim; ++ith){
      int nph, istart_north;
      bool shifted;
      map.GetRingInfo (ith+1,istart_north,nph, cth[ith-llim],sth[ith-llim],
		       shifted);
    }
    
    THealFFT::Ylmgen generator(fLmax,fMmax,1e-30);
    std::vector<double> Ylm;
    std::vector<std::complex<double> > alm_tmp(fLmax+1);
    for(int m=0; m<=fMmax; ++m){
      for(int l=m; l<=fLmax; ++l){
        alm_tmp[l]=(*this)(l,m);
      } // l
      
      for(int ith=0; ith<ulim-llim; ++ith){
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if(l<=fLmax){
	  std::complex<double> p1(0, 0), p2(0, 0);
	  
          if ((l-m)&1) goto middle;
	start:    p1.real() += alm_tmp[l].real()*Ylm[l]; p1.imag() += alm_tmp[l].imag()*Ylm[l];
          if (++l>fLmax) goto end;
	middle:   p2.real() += alm_tmp[l].real()*Ylm[l]; p2.imag() += alm_tmp[l].imag()*Ylm[l];
          if (++l<=fLmax) goto start;
	end:      b_north[ith][m] = p1+p2; b_south[ith][m] = p1-p2;
	}
        else
          {
	    b_north[ith][m] = b_south[ith][m] = 0;
          }
      }
    }

    std::vector<std::complex<double> > shiftarr(fMmax+1), work(4*nside);
    THealFFT::rfft plan;
    for(int ith=llim; ith<ulim; ++ith){
      int istart_north, istart_south, nph;
      double dum1, dum2;
      bool shifted;
      map.GetRingInfo (ith+1,istart_north,nph,dum1,dum2,shifted);
      istart_south = map.GetNpix()-istart_north-nph;
      recalc_alm2map (nph, fMmax, plan, shiftarr);
      if(map.GetTypeString() == "D"){
	fft_alm2map (nph, fMmax, shifted, plan,
		     &(dynamic_cast<THealPixD&>(map))[istart_north],
		     &(dynamic_cast<THealPixD&>(map))[istart_south],
		     b_north[ith-llim], b_south[ith-llim], shiftarr, work);
      } else if(map.GetTypeString() == "F"){
	fft_alm2map (nph, fMmax, shifted, plan,
		     &(dynamic_cast<THealPixF&>(map))[istart_north],
		     &(dynamic_cast<THealPixF&>(map))[istart_south],
		     b_north[ith-llim], b_south[ith-llim], shiftarr, work);
      } // if
    }
  }
  
  Double_t entries = 0;
  for(Int_t i = 0; i < map.GetNpix(); i++){
    if(map.GetBinContent(i) != 0){
      entries++;
    } // if
  } // i
  map.SetEntries(entries);

  for(Int_t i = 0; i < chunksize; i++){
    delete b_north[i];
    delete b_south[i];
  } // i
}

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::Copy(TObject& obj) const
{
  TObject::Copy(obj);
  ((THealAlm&)obj).fAlm  = fAlm;
  ((THealAlm&)obj).fLmax = fLmax;
  ((THealAlm&)obj).fMmax = fMmax;
  ((THealAlm&)obj).fTval = fTval;
}

//______________________________________________________________________________
template<typename T>
void THealAlm<T>::Divide(const THealAlm<T>* ha1)
{
  if(!ha1){
    Error("Divide", "Attempt to divide by a non-existing Alm");
    return;
  } // if

  if(fLmax != ha1->GetLmax() || fMmax != ha1->GetMmax()){
    Error("Divide", "Attempt to devide Alms with different number of lmax or mmax");
    return;
  } // if

  for(Int_t l = 0; l <= fLmax; l++){
    for(Int_t m = 0; m <= TMath::Min(l, fMmax); m++){
      if(std::norm((*ha1)(l, m)) != 0){
	(*this)(l, m) /= (*ha1)(l, m);
      } else {
	(*this)(l, m) = 0;
      } // if
    } // m
  } // l
}

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::Multiply(const THealAlm<T>* ha1)
{
  if(!ha1){
    Error("Multiply", "Attempt to multiply a non-existing Alm");
    return;
  } // if

  if(fLmax != ha1->GetLmax() || fMmax != ha1->GetMmax()){
    Error("Divide", "Attempt to multiply Alms with different number of lmax or mmax");
    return;
  } // if

  for(Int_t l = 0; l <= fLmax; l++){
    for(Int_t m = 0; m <= TMath::Min(l, fMmax); m++){
      (*this)(l, m) *= (*ha1)(l, m);
    } // m
  } // l
}

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::ResetLM(Int_t lmax, Int_t mmax)
{
  fLmax = (lmax > 0) ? lmax : 0;
  mmax  = (mmax > 0) ? mmax : 0;
  fMmax = (mmax <= fLmax) ? mmax : fLmax;
  fTval = 2*fLmax + 1;
  Int_t n = (fMmax + 1)*(fMmax + 2)/2 + (fMmax + 1)*(fLmax - fMmax);

  fAlm.resize(n, fgAlmDummyConst);
}

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::Scale(Double_t c1)
{
  Add(this, this, c1, 0);
}

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::SetToZero()
{
  for(UInt_t i = 0; i < fAlm.size(); i++){
    fAlm[i] = fgAlmDummyConst;
  } // i
}

//______________________________________________________________________________
template<typename T>
void THealAlm<T>::Streamer(TBuffer& b)
{
  if(b.IsReading()){
    UInt_t R__s, R__c;
    Version_t R__v = b.ReadVersion(&R__s, &R__c);
    b.ReadClassBuffer(THealPix::Class(), this, R__v, R__s, R__c);
    fTval = 2*fLmax + 1;
  } else {
    b.WriteClassBuffer(THealPix::Class(), this);
  }
}

//______________________________________________________________________________
template<typename T>
std::complex<T>& THealAlm<T>::operator()(Int_t l, Int_t m)
{
  if(0 <= l && l <= fLmax && 0 <= m && m <= TMath::Min(fMmax, l)){
    return fAlm[((m*(fTval - m))>>1) + l];
  } // if

  fgAlmDummy = fgAlmDummyConst;
  return fgAlmDummy;
}

//______________________________________________________________________________
template<typename T>
const std::complex<T>& THealAlm<T>::operator()(Int_t l, Int_t m) const
{
  if(0 <= l && l <= fLmax && 0 <= m && m <= TMath::Min(fMmax, l)){
    return fAlm[((m*(fTval - m))>>1) + l];
  } // if

  return fgAlmDummyConst;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T>& THealAlm<T>::operator=(const THealAlm<T>& ha1)
{
   // Operator =

  if(this != &ha1){
    ((THealAlm<T>&)ha1).Copy(*this);
  } // if
  return *this;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> THealAlm<T>::__add__(const THealAlm<T>& ha1) const
{
  // Python operator +
  return *this + ha1;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> THealAlm<T>::__div__(const THealAlm<T>& ha1) const
{
  // Python operator /
  return *this / ha1;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> THealAlm<T>::__mul__(const THealAlm<T>& ha1) const
{
  // Python operator *
  return *this * ha1;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> THealAlm<T>::__mul__(Double_t c1) const
{
  // Python operator *
  return *this*c1;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> THealAlm<T>::__rmul__(Double_t c1) const
{
  // Python operator *
  return *this*c1;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> THealAlm<T>::__sub__(const THealAlm<T>& ha1) const
{
  // Python operator -
  return *this - ha1;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> operator*(Double_t c1, const THealAlm<T>& ha1)
{
   // Operator *
   THealAlm<T> hanew = ha1;
   hanew.Scale(c1);
   return hanew;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> operator+(const THealAlm<T>& ha1, const THealAlm<T>& ha2)
{
   // Operator +
   THealAlm<T> hanew = ha1;
   hanew.Add(&ha2, 1);
   return hanew;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> operator-(const THealAlm<T>& ha1, const THealAlm<T>& ha2)
{
   // Operator -
   THealAlm<T> hanew = ha1;
   hanew.Add(&ha2, -1);
   return hanew;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> operator*(const THealAlm<T>& ha1, const THealAlm<T>& ha2)
{
   // Operator *
   THealAlm<T> hanew = ha1;
   hanew.Multiply(&ha2);
   return hanew;
}

//______________________________________________________________________________
template<typename T>
THealAlm<T> operator/(const THealAlm<T>& ha1, const THealAlm<T>& ha2)
{
   // Operator /
   THealAlm<T> hanew = ha1;
   hanew.Divide(&ha2);
   return hanew;
}

template class THealAlm<Float_t>;
template class THealAlm<Double_t>;
template THealAlm<Float_t> operator*<Float_t>(Double_t c1, const THealAlm<Float_t>& ha1);
template THealAlm<Double_t> operator*<Double_t>(Double_t c1, const THealAlm<Double_t>& ha1);
template THealAlm<Float_t> operator*<Float_t>(const THealAlm<Float_t>& ha1, Double_t c1);
template THealAlm<Double_t> operator*<Double_t>(const THealAlm<Double_t>& ha1, Double_t c1);
template THealAlm<Float_t> operator+<Float_t>(const THealAlm<Float_t>& ha1, const THealAlm<Float_t>& ha2);
template THealAlm<Double_t> operator+<Double_t>(const THealAlm<Double_t>& ha1, const THealAlm<Double_t>& ha2);
template THealAlm<Float_t> operator-<Float_t>(const THealAlm<Float_t>& ha1, const THealAlm<Float_t>& ha2);
template THealAlm<Double_t> operator-<Double_t>(const THealAlm<Double_t>& ha1, const THealAlm<Double_t>& ha2);
template THealAlm<Float_t> operator*<Float_t>(const THealAlm<Float_t>& ha1, const THealAlm<Float_t>& ha2);
template THealAlm<Double_t> operator*<Double_t>(const THealAlm<Double_t>& ha1, const THealAlm<Double_t>& ha2);
template THealAlm<Float_t> operator/<Float_t>(const THealAlm<Float_t>& ha1, const THealAlm<Float_t>& ha2);
template THealAlm<Double_t> operator/<Double_t>(const THealAlm<Double_t>& ha1, const THealAlm<Double_t>& ha2);
