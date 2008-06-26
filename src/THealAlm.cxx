// $Id: THealAlm.cxx,v 1.1 2008/06/26 23:23:51 oxon Exp $
// Author: Akira Okumura 2008/06/26

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include "THealAlm.h"
#include "THealPix.h"

templateClassImp(THealAlm)

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

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::ResetLM(Int_t lmax, Int_t mmax)
{
  fLmax = (lmax >= 0) ? lmax : 0;
  mmax  = (mmax >= 0) ? mmax : 0;
  fMmax = (mmax <= fLmax) ? mmax : fLmax;
  fTval = 2*fLmax + 1;
  Int_t n = (fMmax + 1)*(fMmax + 2)/2 + (fMmax + 1)*(fLmax - fMmax);
  fAlm.resize(n, std::vector<std::complex<T> >(0, 0));
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
}

//_____________________________________________________________________________
template<typename T>
void THealAlm<T>::Multiply(const THealAlm<T>* ha1)
{
  if(!ha1){
    Error("Multiply", "Attempt to multiply a non-existing Alm");
    return;
  } // if
}

//______________________________________________________________________________
template<typename T>
std::vector<std::complex<T> >& THealAlm<T>::operator()(Int_t l, Int_t m)
{
  if(0 <= l && l <= fLmax && -fMmax <= m && m <= fMmax){
    return fAlm[((m*(fTval - m))>>1) + l];
  } // if

  fAlmDummy = std::complex<T>(0, 0);
  return fAlmDummy;
}

//______________________________________________________________________________
template<typename T>
const std::vector<std::complex<T> >& THealAlm<T>::operator()(Int_t l, Int_t m) const
{
  if(0 <= l && l <= fLmax && -fMmax <= m && m <= fMmax){
    return fAlm[((m*(fTval - m))>>1) + l];
  } // if

  fAlmDummy = std::complex<T>(0, 0);
  return fAlmDummy;
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
