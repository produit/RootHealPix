// $Id: THealPix.h,v 1.4 2008/06/25 05:39:26 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_PIX
#define T_HEAL_PIX

#include <string>

#include "TArrayD.h"
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
  TDirectory*   fDirectory;     //!Pointer to directory holding this HEALPixxs
  static Bool_t fgAddDirectory; //!flag to add HEALPixs to the directory
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
  THealPix(const THealPix&);
  virtual ~THealPix();
  
  virtual void     AddBinContent(Int_t bin);
  virtual void     AddBinContent(Int_t bin, Double_t w);
  static  void     AddDirectory(Bool_t add = kTRUE);
  static  Bool_t   AddDirectoryStatus();
  virtual void     Draw(Option_t* option = "");
  virtual Int_t    Fill(Double_t theta, Double_t phi);
  virtual Int_t    Fill(Double_t theta, Double_t phi, Double_t w);
  virtual Int_t    FindBin(Double_t theta, Double_t phi) const;
  virtual Double_t GetBinContent(Int_t bin) const;
  virtual Double_t GetEntries() const;
  virtual Int_t    GetNside() const { return fNside;}
  virtual Int_t    GetNpix() const { return fNpix;}
  virtual Int_t    GetNrows() const;
  virtual Int_t    GetOrder() const { return fOrder;}
  virtual std::string GetSchemeString() const;
  virtual Int_t    GetType() const { return fType;}
  virtual std::string GetTypeString() const;
  virtual std::string GetUnit() const { return fUnit;}
  virtual Bool_t   IsNested() const {return fIsNested;}
  virtual void     Paint(Option_t *option="");
  virtual THealPix* Rebin(Int_t neworder, const char* newname = "");
  virtual void     SetBinContent(Int_t bin, Double_t content);
  virtual void     SetBinsLength(Int_t = -1) {} // redefined in derived cplasses
  virtual void     SetEntries(Double_t n) {fEntries = n;}
  virtual void     SetOrder(Int_t order);
  virtual void     SetUnit(const char* unit);
  virtual void     SetDegree(Bool_t isDegree = kTRUE) {fIsDegree = isDegree;}

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

  ClassDef(THealPixD, 0);
};

#endif // T_HEAL_PIX
