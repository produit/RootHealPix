// $Id: THealAlm.h,v 1.3 2008/07/04 22:11:33 oxon Exp $
// Author: Akira Okumura 2008/06/26

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_ALM
#define T_HEAL_ALM

#include <complex>
#include <vector>

#include "TObject.h"

#include "THealFFT.h"

class THealPix;
template<typename T> class THealAlm;
template<typename T> THealAlm<T> operator*(Double_t c1, const THealAlm<T>& ha1);
template<typename T> THealAlm<T> operator*(const THealAlm<T>& ha1, Double_t c1);
template<typename T> THealAlm<T> operator+(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
template<typename T> THealAlm<T> operator-(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
template<typename T> THealAlm<T> operator*(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
template<typename T> THealAlm<T> operator/(const THealAlm<T>& ha1, const THealAlm<T>& ha2);

template<typename T> class THealAlm : public TObject {
private:
  std::vector<std::complex<T> > fAlm;
  static std::complex<T> fgAlmDummy; //!
  static const std::complex<T> fgAlmDummyConst; //!
  Int_t fLmax; // Maximum of l
  Int_t fMmax; // Maximum of m
  Int_t fTval; //!to be used for faster calculation

protected:
  virtual void Copy(TObject& hanew) const;
  
public:
  THealAlm(Int_t lmax = 0, Int_t mmax = 0);
  THealAlm(const THealAlm<T>& ha);
  virtual ~THealAlm();
  virtual void  Add(const THealAlm<T>* ha1, Double_t c1 = 1);
  virtual void  Add(const THealAlm<T>* ha1, const THealAlm<T>* ha2, Double_t c1, Double_t c2);
  virtual Int_t GetLmax() const {return fLmax;}
  virtual Int_t GetMmax() const {return fMmax;}
  virtual std::complex<T>* GetMstart(Int_t m) {return &(fAlm[(m*(fTval - m))>>1]);}
  virtual void  Divide(const THealAlm<T>* ha1);
  virtual void  Multiply(const THealAlm<T>* ha1);
  virtual void  ResetLM(Int_t lmax, Int_t mmax);
  virtual void  Scale(Double_t c1);
  virtual void  SetAlm(Int_t l, Int_t m, std::complex<T> a) {(*this)(l, m) = a;}
  virtual void  SetToZero();

        std::complex<T>& operator()(Int_t l, Int_t m);
  const std::complex<T>& operator()(Int_t l, Int_t m) const;

  THealAlm<T>& operator=(const THealAlm<T>& ha1);
  friend THealAlm<T>  operator*<>(Double_t c1, const THealAlm<T>& ha1);
  friend THealAlm<T>  operator*<>(const THealAlm<T>& ha1, Double_t c1);
  friend THealAlm<T>  operator+<>(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
  friend THealAlm<T>  operator-<>(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
  friend THealAlm<T>  operator*<>(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
  friend THealAlm<T>  operator/<>(const THealAlm<T>& ha1, const THealAlm<T>& ha2);

  // Operators for PyROOT
  virtual THealAlm<T>  __add__(const THealAlm<T>& ha1) const;
  virtual THealAlm<T>  __div__(const THealAlm<T>& ha1) const;
  virtual THealAlm<T>  __mul__(Double_t c1) const;
  virtual THealAlm<T>  __mul__(const THealAlm<T>& ha1) const;
  virtual THealAlm<T>  __rmul__(Double_t c1) const;
  virtual THealAlm<T>  __sub__(const THealAlm<T>& ha1) const;

  ClassDef(THealAlm, 1);
};

template<typename T> inline THealAlm<T> operator*(const THealAlm<T>& ha1, Double_t c1)
{
  return operator*(c1, ha1);
}

#endif // T_HEAL_ALM
