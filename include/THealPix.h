// $Id: THealPix.h,v 1.27 2008/07/11 23:57:48 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.

   This is a port of HEALPix C++ package to ROOT system.
   Original code is available at <http://healpix.jpl.nasa.gov> under GPL.
******************************************************************************/

#ifndef T_HEAL_PIX
#define T_HEAL_PIX

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THealPix                                                             //
//                                                                      //
// HEALPix-like 2D histogram base class.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

#include "TArrayD.h"
#include "TArrayF.h"
#include "TAttFill.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TAxis.h"
#include "TNamed.h"

#include <fitsio.h>

class TDirectory;
template<typename T> class THealAlm;
class TList;
class TVirtualHealPainter;

class THealPix : public TNamed, public TAttLine, public TAttFill, public TAttMarker {
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
  TAxis         fXaxis;         //X axis descriptor
  TAxis         fYaxis;         //Y axis descriptor
  TAxis         fZaxis;         //Z axis descriptor
  TList*        fFunctions;     //->Pointer to list of functions (fits and user)
  Short_t       fBarOffset;     //(1000*offset) for bar charts or legos
  Short_t       fBarWidth;      //(1000*width) for bar charts or legos
  Int_t         fOrder;         //Order of resolution
  Int_t         fNside;         //!= 2^fOrder
  Int_t         fNpix;          //!= 12*fNside^2
  Int_t         fNside2;        //!= fNside*fNside (used for faster calculation)
  Int_t         fNcap;          //!= 2*fNside*(fNside-1) (used for faster calculation)
  Double_t      f3Nside2;       //!= 3*fNside*fNside (used for faster calculation)
  Double_t      f2over3Nside;   //!= 2/(3*fNside) (used for faster calculation)
  Int_t         f4Nside;        //!= 4*fNside
  Int_t         f2Nside;        //!= 2*fNside
  Bool_t        fIsDegree;      //deg = true, rad = false
  Bool_t        fIsNested;      //RING = false, NESTED = true
  std::string   fUnit;          //Unit of data (used in FITS header)
  Double_t      fEntries;       //Number of entries
  Double_t      fMaximum;       //Maximum value for plotting
  Double_t      fMinimum;       //Minimum value for plotting
  Double_t      fNormFactor;    //Normalization factor
  TArrayD       fContour;       //Array to display contour levels
  Double_t      fTsumw;         //Total Sum of weights
  Double_t      fTsumw2;        //Total Sum of squares of weights
  TArrayD       fSumw2;         //Total Sum of squares of weights
  TString       fOption;        //histogram options
  TDirectory*   fDirectory;     //!Pointer to directory holding this HEALPixxs
  TVirtualHealPainter* fPainter;//!pointer to HEALPix painter
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
    kNoStats     = BIT(9),  // don't draw stats box
    kUserContour = BIT(10), // user specified contour levels
    kCanRebin    = BIT(11), // can rebin axis
    kLogX        = BIT(15), // X-axis in log scale
    kIsZoomed    = BIT(16), // bit set when zooming on Y axis
    kNoTitle     = BIT(17), // don't draw the histogram title
    kIsAverage   = BIT(18)  // Bin contents are average (used by Add)
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
  virtual Int_t    Fill(const Float_t* x);
  virtual Int_t    Fill(const Float_t* x, Double_t w);
  virtual Int_t    FindBin(Double_t theta, Double_t phi) const;
  virtual Double_t GetAverage() const;
  virtual Float_t  GetBarOffset() const {return Float_t(0.001*Float_t(fBarOffset));}
  virtual Float_t  GetBarWidth() const  {return Float_t(0.001*Float_t(fBarWidth));}
  virtual Double_t GetBinArea(Bool_t degree2 = kFALSE) const;
  virtual void     GetBinCenter(Int_t bin, Double_t& theta, Double_t& phi) const;
  virtual void     GetBinCenter(Int_t bin, Double_t* theta, Double_t* phi) const;
  virtual Double_t GetBinContent(Int_t bin) const;
  virtual Double_t GetBinError(Int_t bin) const;
  virtual Int_t    GetBinVertices(Int_t bin, Double_t* x, Double_t* y) const;
  virtual Int_t    GetContour(Double_t* levels = 0);
  virtual Double_t GetContourLevel(Int_t level) const;
  virtual Double_t GetContourLevelPad(Int_t level) const;
  TDirectory*      GetDirectory() const {return fDirectory;}
  static  Bool_t   GetDefaultSumw2();
  virtual Double_t GetEntries() const;
  TList*           GetListOfFunctions() const {return fFunctions;}
  virtual Double_t GetMaximum(Double_t maxval = FLT_MAX) const;
  virtual Int_t    GetMaximumBin() const;
  virtual Int_t    GetMaximumStored() const {return fMaximum;}
  virtual Double_t GetMinimum(Double_t minval = -FLT_MAX) const;
  virtual Int_t    GetMinimumBin() const;
  virtual Int_t    GetMinimumStored() const {return fMinimum;}
  virtual Int_t    GetNside() const { return fNside;}
  virtual Double_t GetNormFactor() const {return fNormFactor;}
  virtual Int_t    GetNpix() const { return fNpix;}
  virtual Int_t    GetNrows() const;
  Option_t*        GetOption() const {return fOption.Data();}
  virtual Int_t    GetOrder() const { return fOrder;}
  TVirtualHealPainter* GetPainter(Option_t* option = "");
  virtual std::string GetSchemeString() const;
  virtual Double_t GetSumOfWeights() const;
  virtual TArrayD* GetSumw2() {return &fSumw2;}
          TAxis*   GetXaxis() const;
          TAxis*   GetYaxis() const;
          TAxis*   GetZaxis() const;
  virtual const TArrayD* GetSumw2() const {return &fSumw2;}
  virtual void     GetRingInfo(Int_t ring, Int_t& startpix, Int_t& ringpix, Double_t &costheta, Double_t& sintheta, Bool_t& shifted) const;
  virtual Int_t    GetSumw2N() const {return fSumw2.fN;}
  virtual Int_t    GetType() const = 0;
  virtual std::string GetTypeString() const = 0;
  virtual std::string GetUnit() const { return fUnit;}
  virtual Bool_t   IsDegree() const {return fIsDegree;}
  virtual Bool_t   IsNested() const {return fIsNested;}
  template<typename T>
          void     Map2Alm(THealAlm<T>& alm, Bool_t add = kFALSE) const;
  template<typename T>
          void     Map2Alm(THealAlm<T>& alm, const std::vector<double>& weight, Bool_t add = kFALSE) const;
  virtual void     Multiply(const THealPix* hp1);
  virtual void     Paint(Option_t* option = "");
  virtual THealPix* Rebin(Int_t neworder, const char* newname = "");
  virtual void     Scale(Double_t c1 = 1, Option_t* option = "");
  virtual void     SetBarOffset(Float_t offset = 0.25) {fBarOffset = Short_t(1000*offset);}
  virtual void     SetBarWidth(Float_t width = 0.5) {fBarWidth = Short_t(1000*width);}
  virtual void     SetBinContent(Int_t bin, Double_t content);
  virtual void     SetBinError(Int_t bin, Double_t error);
  virtual void     SetBins(Int_t n);
  virtual void     SetBinsLength(Int_t = -1) {} // redefined in derived cplasses
  virtual void     SetContour(Int_t nlevels, const Double_t* levels = 0);
  virtual void     SetContourLevel(Int_t level, Double_t value);
  static  void     SetDefaultSumw2(Bool_t sumw2=kTRUE);
  virtual void     SetDegree(Bool_t isDegree = kTRUE) {fIsDegree = isDegree;}
  virtual void     SetDirectory(TDirectory *dir);
  virtual void     SetEntries(Double_t n) {fEntries = n;}
  virtual void     SetMaximum(Double_t maximum = -1111);
  virtual void     SetMinimum(Double_t minimum = -1111);
  virtual void     SetName(const char* name);
  virtual void     SetNameTitle(const char* name, const char* title);
  virtual void     SetNormFactor(Double_t factor = 1) {fNormFactor = factor;}
  virtual void     SetOption(Option_t* option = "") {fOption = option;}
  virtual void     SetOrder(Int_t order);
  virtual void     SetUnit(const char* unit);
  virtual void     SetXTitle(const char* title) {fXaxis.SetTitle(title);}
  virtual void     SetYTitle(const char* title) {fYaxis.SetTitle(title);}
  virtual void     SetZTitle(const char* title) {fZaxis.SetTitle(title);}
  virtual void     Sumw2();
  void             UseCurrentStyle();
  virtual void     Nest2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const;
  virtual Int_t    XYF2Nest(Int_t x, Int_t y, Int_t face) const;
  virtual void     Ring2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const;
  virtual Int_t    XYF2Ring(Int_t x, Int_t y, Int_t face) const;
  virtual Int_t    Nest2Ring(Int_t pix) const;
  virtual Int_t    Ring2Nest(Int_t pix) const;
  virtual void     Pix2XY(Int_t pix, Int_t& x, Int_t& y) const;
  virtual Int_t    XY2Pix(Int_t x, Int_t y) const;

  static  Int_t    Nside2Npix(Int_t nside) {return 12*nside*nside;}
  static  Int_t    Order2Nside(Int_t order) {return 1<<order;}

  ClassDef(THealPix, 1);
};

//_____________________________________________________________________________
inline void THealPix::Nest2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const
{
  // Original code is Healpix_Base::nest2xyf of HEALPix C++
  face = pix>>(2*fOrder);
  Pix2XY(pix & (fNside2 - 1), x, y);
}

//_____________________________________________________________________________
inline Int_t THealPix::XYF2Nest(Int_t x, Int_t y, Int_t face) const
{
  // Original code is Healpix_Base::xyf2nest of HEALPix C++
  return (face<<(2*fOrder)) + XY2Pix(x, y);
}

//_____________________________________________________________________________
inline Int_t THealPix::Nest2Ring(Int_t pix) const
{
  // Original code is Healpix_Base::nest2ring of HEALPix C++
  Int_t x, y, face;
  Nest2XYF(pix, x, y, face);
  return XYF2Ring(x, y, face);
}

//_____________________________________________________________________________
inline Int_t THealPix::Ring2Nest(Int_t pix) const
{
  // Original code is Healpix_Base::ring2nest of HEALPix C++
  Int_t x, y, face;
  Ring2XYF(pix, x, y, face);
  return XYF2Nest(x, y, face);
}

//_____________________________________________________________________________
inline void THealPix::Pix2XY(Int_t pix, Int_t& x, Int_t& y) const
{
  // Original code is Healpix_Base::pix2xy of HEALPix C++
  Int_t raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  x = fgTable.C(raw&0xff) | (fgTable.C(raw>>8)<<4);
  raw = ((pix&0xaaaa)>>1) | ((pix&0xaaaa0000)>>16);
  y = fgTable.C(raw&0xff) | (fgTable.C(raw>>8)<<4);
}

//_____________________________________________________________________________
inline Int_t THealPix::XY2Pix(Int_t x, Int_t y) const
{
  // Original code is Healpix_Base::xy2pix of HEALPix C++
  return fgTable.U(x&0xff) | (fgTable.U(x>>8)<<16) | (fgTable.U(y&0xff)<<1)
    | (fgTable.U(y>>8)<<17);
}

//______________________________________________________________________________
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
  virtual void       Copy(TObject& newhp) const;
  virtual Double_t   GetBinContent(Int_t bin) const;
  virtual Int_t      GetType() const;
  virtual std::string GetTypeString() const;
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

  ClassDef(THealPixF, 1);
};

THealPixF operator*(Double_t c1, const THealPixF &hp1);
inline
THealPixF operator*(const THealPixF &hp1, Double_t c1) {return operator*(c1, hp1);}
THealPixF operator+(const THealPixF &hp1, const THealPixF &hp2);
THealPixF operator-(const THealPixF &hp1, const THealPixF &hp2);
THealPixF operator*(const THealPixF &hp1, const THealPixF &hp2);
THealPixF operator/(const THealPixF &hp1, const THealPixF &hp2);

//______________________________________________________________________________
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
  virtual void       Copy(TObject& newhp) const;
  virtual Double_t   GetBinContent(Int_t bin) const;
  virtual Int_t      GetType() const;
  virtual std::string GetTypeString() const;
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

  ClassDef(THealPixD, 1);
};

THealPixD operator*(Double_t c1, const THealPixD& hp1);
inline
THealPixD operator*(const THealPixD& hp1, Double_t c1) {return operator*(c1, hp1);}
THealPixD operator+(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator-(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator*(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator/(const THealPixD& hp1, const THealPixD& hp2);

#endif // T_HEAL_PIX
