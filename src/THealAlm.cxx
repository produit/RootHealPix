// $Id: THealAlm.cxx,v 1.6 2008/07/04 23:06:09 oxon Exp $
// Author: Akira Okumura 2008/06/26

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include "TMath.h"

#include "THealAlm.h"
#include "THealPix.h"

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
