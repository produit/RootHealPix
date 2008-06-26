// $Id: THealAlm.h,v 1.1 2008/06/26 23:23:50 oxon Exp $
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

class THealPix;
template<typename T> class THealAlm : public TObject {
private:
  std::vector<std::complex<T> > fAlm;
  std::complex<T> fAlmDummy;
  Int_t fLmax;
  Int_t fMmax;
  Int_t fTval;
  
protected:
  virtual void Copy(TObject& hanew) const;
  
public:
  THealAlm(Int_t lmax = 0, Int_t mmax = 0);
  THealAlm(const THealAlm<T>& ha);
  virtual ~THealAlm();
  virtual void  Add(const THealAlm<T>* ha1, Double_t c1 = 1);
  virtual Int_t GetLmax() const {return fLmax;}
  virtual Int_t GetMmax() const {return fMmax;}
  virtual void  Divide(const THealAlm<T>* ha1);
  virtual void  Multiply(const THealAlm<T>* ha1);
  virtual void  ResetLM(Int_t lmax, Int_t mmax);
  
        std::vector<std::complex<T> >& operator()(Int_t l, Int_t m);
  const std::vector<std::complex<T> >& operator()(Int_t l, Int_t m) const;

          THealAlm<T>& operator=(const THealAlm<T>& ha1);
  friend  THealAlm<T>  operator*(Double_t c1, const THealAlm<T>& ha1);
  friend  THealAlm<T>  operator*(const THealAlm<T>& ha1, Double_t c1);
  friend  THealAlm<T>  operator+(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
  friend  THealAlm<T>  operator-(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
  friend  THealAlm<T>  operator*(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
  friend  THealAlm<T>  operator/(const THealAlm<T>& ha1, const THealAlm<T>& ha2);

  // Operators for PyROOT
  virtual THealAlm<T>  __add__(const THealAlm<T>& ha1) const;
  virtual THealAlm<T>  __div__(const THealAlm<T>& ha1) const;
  virtual THealAlm<T>  __mul__(Double_t c1) const;
  virtual THealAlm<T>  __mul__(const THealAlm<T>& ha1) const;
  virtual THealAlm<T>  __rmul__(Double_t c1) const;
  virtual THealAlm<T>  __sub__(const THealAlm<T>& ha1) const;

  ClassDef(THealAlm, 0);
};

template<typename T> THealAlm<T> operator*(Double_t c1, const THealAlm<T>& ha1);
template<typename T> inline THealAlm<T> operator*(const THealAlm<T>& ha1, Double_t c1)
{
  return operator*(c1, ha1);
}
template<typename T> THealAlm<T> operator+(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
template<typename T> THealAlm<T> operator-(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
template<typename T> THealAlm<T> operator*(const THealAlm<T>& ha1, const THealAlm<T>& ha2);
template<typename T> THealAlm<T> operator/(const THealAlm<T>& ha1, const THealAlm<T>& ha2);

#endif // T_HEAL_ALM
