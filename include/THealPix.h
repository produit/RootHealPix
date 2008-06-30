// $Id: THealPix.h,v 1.13 2008/06/30 19:26:09 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_PIX
#define T_HEAL_PIX

#include <string>

#include "TArrayD.h"
#include "TArrayF.h"
#include "TNamed.h"

#include <fitsio.h>

class TDirectory;

class THealPix : public TNamed {
protected:
  class THealTable {
  private:
    Short_t fCtab[0x100];
    Short_t fUtab[0x100];
  public:
    THealTable();
    inline Short_t C(Int_t i) const {return fCtab[i];}
    inline Short_t U(Int_t i) const {return fUtab[i];}
  };

  struct HealHeader_t {
    Int_t order, colnum;
    Int_t nside, npix, nrows;
    Bool_t isNested;
    char tunit[FLEN_VALUE];;
  };

protected:
  Int_t         fOrder;         //Order of resolution
  Int_t         fNside;         //= 2^fOrder
  Int_t         fNpix;          //= 12*fNside^2
  Int_t         fNpFace;        //= fNside*fNside (used for faster calculation)
  Int_t         fNcap;          //= fNpFace*fNside (used for faster calculation)
  Bool_t        fIsDegree;      //deg = true, rad = false
  Bool_t        fIsNested;      //RING = false, NESTED = true
  Int_t         fType;          //Data type : TDOUBLE, TFLOAT, TINT
  std::string   fUnit;          //Unit of data (used in FITS header)
  Double_t      fEntries;       //Number of entries
  Double_t      fTsumw;         //Total Sum of weights
  Double_t      fTsumw2;        //Total Sum of squares of weights
  TArrayD       fSumw2;         //Total Sum of squares of weights
  TDirectory*   fDirectory;     //!Pointer to directory holding this HEALPixxs
  static Bool_t fgAddDirectory; //!flag to add HEALPixs to the directory
  static Bool_t fgDefaultSumw2; //!flag to call THealPix::Sumw2 automatically at histogram creation time
  static THealTable fgTable;    //!ctab and utab holder
  static const Int_t fgJrll[];  //!
  static const Int_t fgJpll[];  //!

private:
  void Build();
  void FillTable();

  THealPix& operator=(const THealPix&); // Not implemented

protected:
  THealPix();
  THealPix(const char* name, const char* title, Int_t order,
	   Bool_t isNested = kFALSE);
  virtual void  Copy(TObject& hpnew) const;
  static Bool_t ReadFitsHeader(fitsfile** fptr, const char* fname,
			       const char* colname, HealHeader_t& head);

public:
  // THealPix status bits
  enum {
    kCanRebin  = BIT(11), // can rebin axis
    kIsAverage = BIT(18)  // Bin contents are average (used by Add)
  };

  THealPix(const THealPix&);
  virtual ~THealPix();
  
  virtual void     Add(const THealPix *hp1, Double_t c1 = 1);
  virtual void     Add(const THealPix* hp1, const THealPix* hp2, Double_t c1 = 1, Double_t c2 = 1);
  virtual void     AddBinContent(Int_t bin);
  virtual void     AddBinContent(Int_t bin, Double_t w);
  static  void     AddDirectory(Bool_t add = kTRUE);
  static  Bool_t   AddDirectoryStatus();
  virtual void     Divide(const THealPix* hp1);
  virtual void     Draw(Option_t* option = "");
  virtual Int_t    Fill(Double_t theta, Double_t phi);
  virtual Int_t    Fill(Double_t theta, Double_t phi, Double_t w);
  virtual Int_t    Fill(const Double_t* x);
  virtual Int_t    Fill(const Double_t* x, Double_t w);
  virtual Int_t    FindBin(Double_t theta, Double_t phi) const;
  virtual Double_t GetBinContent(Int_t bin) const;
  virtual Double_t GetBinError(Int_t bin) const;
  TDirectory*      GetDirectory() const {return fDirectory;}
  static  Bool_t   GetDefaultSumw2();
  virtual Double_t GetEntries() const;
  virtual Int_t    GetNside() const { return fNside;}
  virtual Int_t    GetNpix() const { return fNpix;}
  virtual Int_t    GetNrows() const;
  virtual Int_t    GetOrder() const { return fOrder;}
  virtual Double_t GetPixelArea(Bool_t degree2 = kFALSE) const;
  virtual std::string GetSchemeString() const;
  virtual TArrayD* GetSumw2() {return &fSumw2;}
  virtual const TArrayD* GetSumw2() const {return &fSumw2;}
  virtual void     GetRingInfo(Int_t ring, Int_t& startpix, Int_t& ringpix, Double_t &costheta, Double_t& sintheta, Bool_t& shifted) const;
  virtual Int_t    GetSumw2N() const {return fSumw2.fN;}
  virtual Int_t    GetType() const { return fType;}
  virtual std::string GetTypeString() const;
  virtual std::string GetUnit() const { return fUnit;}
  virtual Bool_t   IsNested() const {return fIsNested;}
  virtual void     Multiply(const THealPix* hp1);
  virtual void     Paint(Option_t* option = "");
  virtual THealPix* Rebin(Int_t neworder, const char* newname = "");
  virtual void     Scale(Double_t c1 = 1, Option_t* option = "");
  virtual void     SetBinContent(Int_t bin, Double_t content);
  virtual void     SetBinError(Int_t bin, Double_t error);
  virtual void     SetBinsLength(Int_t = -1) {} // redefined in derived cplasses
  static  void     SetDefaultSumw2(Bool_t sumw2=kTRUE);
  virtual void     SetDegree(Bool_t isDegree = kTRUE) {fIsDegree = isDegree;}
  virtual void     SetDirectory(TDirectory *dir);
  virtual void     SetEntries(Double_t n) {fEntries = n;}
  virtual void     SetOrder(Int_t order);
  virtual void     SetUnit(const char* unit);
  virtual void     Sumw2();
  virtual void     Nest2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const;
  virtual Int_t    XYF2Nest(Int_t x, Int_t y, Int_t face) const;
  virtual void     Ring2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const;
  virtual Int_t    XYF2Ring(Int_t x, Int_t y, Int_t face) const;
  virtual Int_t    Nest2Ring(Int_t pix) const;
  virtual Int_t    Ring2Nest(Int_t pix) const;
  virtual void     Pix2XY(Int_t pix, Int_t& x, Int_t& y) const;
  virtual Int_t    XY2Pix(Int_t x, Int_t y) const;

  ClassDef(THealPix, 0);
};

//_____________________________________________________________________________
class THealPixF : public THealPix, public TArrayF {
public:
  THealPixF();
  THealPixF(const char* name, const char* title, Int_t order,
	    Bool_t isNested = kFALSE);
  THealPixF(const THealPixF& hpd);
  virtual ~THealPixF();
  
  virtual void       AddBinContent(Int_t bin) {++fArray[bin];}
  virtual void       AddBinContent(Int_t bin, Double_t w)
                                  {fArray[bin] += Float_t(w);}
  virtual Double_t   GetBinContent(Int_t bin) const;
  
  virtual void       Copy(TObject& newhp) const;
  static  THealPixF* ReadFits(const char* fname, const char* colname);
  virtual void       SetBinContent(Int_t bin, Double_t content);
  virtual void       SetBinsLength(Int_t n=-1);
          THealPixF& operator=(const THealPixF &hp1);
  friend  THealPixF  operator*(Double_t c1, const THealPixF &hp1);
  friend  THealPixF  operator*(const THealPixF &hp1, Double_t c1);
  friend  THealPixF  operator+(const THealPixF &hp1, const THealPixF &hp2);
  friend  THealPixF  operator-(const THealPixF &hp1, const THealPixF &hp2);
  friend  THealPixF  operator*(const THealPixF &hp1, const THealPixF &hp2);
  friend  THealPixF  operator/(const THealPixF &hp1, const THealPixF &hp2);

  // Operators for PyROOT
  virtual THealPixF  __add__(const THealPixF& hp1) const;
  virtual THealPixF  __div__(const THealPixF& hp1) const;
  virtual THealPixF  __mul__(Double_t c1) const;
  virtual THealPixF  __mul__(const THealPixF& hp1) const;
  virtual THealPixF  __rmul__(Double_t c1) const;
  virtual THealPixF  __sub__(const THealPixF& hp1) const;

  ClassDef(THealPixF, 0);
};

THealPixF operator*(Double_t c1, const THealPixF &hp1);
inline
THealPixF operator*(const THealPixF &hp1, Double_t c1) {return operator*(c1, hp1);}
THealPixF operator+(const THealPixF &hp1, const THealPixF &hp2);
THealPixF operator-(const THealPixF &hp1, const THealPixF &hp2);
THealPixF operator*(const THealPixF &hp1, const THealPixF &hp2);
THealPixF operator/(const THealPixF &hp1, const THealPixF &hp2);

//_____________________________________________________________________________
class THealPixD : public THealPix, public TArrayD {
public:
  THealPixD();
  THealPixD(const char* name, const char* title, Int_t order,
	    Bool_t isNested = kFALSE);
  THealPixD(const THealPixD& hpd);
  virtual ~THealPixD();
  
  virtual void       AddBinContent(Int_t bin) {++fArray[bin];}
  virtual void       AddBinContent(Int_t bin, Double_t w)
                                  {fArray[bin] += Double_t(w);}
  virtual Double_t   GetBinContent(Int_t bin) const;
  
  virtual void       Copy(TObject& newhp) const;
  static  THealPixD* ReadFits(const char* fname, const char* colname);
  virtual void       SetBinContent(Int_t bin, Double_t content);
  virtual void       SetBinsLength(Int_t n=-1);
          THealPixD& operator=(const THealPixD& hp1);
  friend  THealPixD  operator*(Double_t c1, const THealPixD& hp1);
  friend  THealPixD  operator*(const THealPixD& hp1, Double_t c1);
  friend  THealPixD  operator+(const THealPixD& hp1, const THealPixD& hp2);
  friend  THealPixD  operator-(const THealPixD& hp1, const THealPixD& hp2);
  friend  THealPixD  operator*(const THealPixD& hp1, const THealPixD& hp2);
  friend  THealPixD  operator/(const THealPixD& hp1, const THealPixD& hp2);

  // Operators for PyROOT
  virtual THealPixD  __add__(const THealPixD& hp1) const;
  virtual THealPixD  __div__(const THealPixD& hp1) const;
  virtual THealPixD  __mul__(Double_t c1) const;
  virtual THealPixD  __mul__(const THealPixD& hp1) const;
  virtual THealPixD  __rmul__(Double_t c1) const;
  virtual THealPixD  __sub__(const THealPixD& hp1) const;

  ClassDef(THealPixD, 0);
};

THealPixD operator*(Double_t c1, const THealPixD& hp1);
inline
THealPixD operator*(const THealPixD& hp1, Double_t c1) {return operator*(c1, hp1);}
THealPixD operator+(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator-(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator*(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator/(const THealPixD& hp1, const THealPixD& hp2);

#endif // T_HEAL_PIX
