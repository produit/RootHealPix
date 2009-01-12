// $Id: THealPix.h,v 1.36 2009/01/12 20:51:24 oxon Exp $
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

#include <complex>
#include <float.h>
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

#include "THealFFT.h"

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

  struct HealHeader_t {
    Int_t order, colnum;
    Int_t nside, npix, nrows, repeat;
    Bool_t isNested;
    char tunit[FLEN_VALUE];
    char ttype[FLEN_VALUE];
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
  virtual THealPix* DrawCopy(Option_t* option = "") const;
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
  virtual Double_t GetMaximumStored() const {return fMaximum;}
  virtual Double_t GetMinimum(Double_t minval = -FLT_MAX) const;
  virtual Int_t    GetMinimumBin() const;
  virtual Double_t GetMinimumStored() const {return fMinimum;}
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
  static  Bool_t   ReadFitsHeader(fitsfile** fptr, const char* fname,
				  const char* colname, HealHeader_t& head);
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

  ClassDef(THealPix, 1); // 2D HEALPix histogram
};

//______________________________________________________________________________
namespace {
  inline UInt_t Isqrt(UInt_t v);
  inline Double_t Modulo(Double_t v1, Double_t v2);
  inline Int_t Modulo(Int_t v1, Int_t v2);
  inline Long_t Modulo(Long_t v1, Long_t v2);
  void get_chunk_info (int nrings, int &nchunks, int &chunksize);
  void fill_work (const std::complex<double> *datain, int nph, int mmax,
		  bool shifted, const std::vector<std::complex<double> > &shiftarr,
		  std::vector<std::complex<double> > &work);
  void read_work (const std::vector<std::complex<double> >& work, int nph,
		  int mmax, bool shifted,
		  const std::vector<std::complex<double> > &shiftarr,
		  std::complex<double> *dataout);
  template<typename T>
  void fft_map2alm (int nph, int mmax, bool shifted, double weight,
		    THealFFT::rfft &plan, const T *mapN, const T *mapS,
		    std::complex<double> *phas_n, std::complex<double> *phas_s,
		    const std::vector<std::complex<double> > &shiftarr,
		    std::vector<std::complex<double> > &work);
  void recalc_map2alm (int nph, int mmax, THealFFT::rfft &plan,
		       std::vector<std::complex<double> > &shiftarr);
  void recalc_alm2map (int nph, int mmax, THealFFT::rfft &plan,
		       std::vector<std::complex<double> > &shiftarr);
  template<typename T>
  void fft_alm2map (int nph, int mmax, bool shifted, THealFFT::rfft &plan,
		    T *mapN, T *mapS, std::complex<double> *b_north,
		    const std::complex<double> *b_south,
		    const std::vector<std::complex<double> > &shiftarr,
		    std::vector<std::complex<double> > &work);
} // namespace

//______________________________________________________________________________
template<typename T>
void THealPix::Map2Alm(THealAlm<T>& alm, Bool_t add) const
{
  std::vector<Double_t> weight(2*fNside, 1);
  Map2Alm(alm, weight, add);
}

//______________________________________________________________________________
inline void THealPix::Nest2XYF(Int_t pix, Int_t& x, Int_t& y, Int_t& face) const
{
  // Original code is Healpix_Base::nest2xyf of HEALPix C++
  face = pix>>(2*fOrder);
  Pix2XY(pix & (fNside2 - 1), x, y);
}

//______________________________________________________________________________
inline Int_t THealPix::XYF2Nest(Int_t x, Int_t y, Int_t face) const
{
  // Original code is Healpix_Base::xyf2nest of HEALPix C++
  return (face<<(2*fOrder)) + XY2Pix(x, y);
}

//______________________________________________________________________________
inline Int_t THealPix::Nest2Ring(Int_t pix) const
{
  // Original code is Healpix_Base::nest2ring of HEALPix C++
  Int_t x, y, face;
  Nest2XYF(pix, x, y, face);
  return XYF2Ring(x, y, face);
}

//______________________________________________________________________________
inline Int_t THealPix::Ring2Nest(Int_t pix) const
{
  // Original code is Healpix_Base::ring2nest of HEALPix C++
  Int_t x, y, face;
  Ring2XYF(pix, x, y, face);
  return XYF2Nest(x, y, face);
}

//______________________________________________________________________________
inline void THealPix::Pix2XY(Int_t pix, Int_t& x, Int_t& y) const
{
  // Original code is Healpix_Base::pix2xy of HEALPix C++
  Int_t raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  x = fgTable.C(raw&0xff) | (fgTable.C(raw>>8)<<4);
  raw = ((pix&0xaaaa)>>1) | ((pix&0xaaaa0000)>>16);
  y = fgTable.C(raw&0xff) | (fgTable.C(raw>>8)<<4);
}

//______________________________________________________________________________
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
  virtual THealPix*  DrawCopy(Option_t* option = "") const;
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

  ClassDef(THealPixF, 1); // 2D HEALPix histogram (one float per channel)
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
  virtual THealPix*  DrawCopy(Option_t* option = "") const;
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

  ClassDef(THealPixD, 1); // 2D HEALPix histogram (one double per channel)
};

THealPixD operator*(Double_t c1, const THealPixD& hp1);
inline
THealPixD operator*(const THealPixD& hp1, Double_t c1) {return operator*(c1, hp1);}
THealPixD operator+(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator-(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator*(const THealPixD& hp1, const THealPixD& hp2);
THealPixD operator/(const THealPixD& hp1, const THealPixD& hp2);

//______________________________________________________________________________
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
	const THealPixD* tmp = dynamic_cast<const THealPixD*>(this);
	const Double_t* p1 = &(tmp->GetArray()[istart_north]);
	const Double_t* p2 = &(tmp->GetArray()[istart_south]);
	fft_map2alm (nph, mmax, shifted, weight[ith]*normfact, plan,
		     p1, p2,	phas_n[ith-llim], phas_s[ith-llim], shiftarr, work);
      } else if(GetTypeString() == "E"){
	const THealPixF* tmp = dynamic_cast<const THealPixF*>(this);
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

//______________________________________________________________________________
namespace {
  UInt_t Isqrt(UInt_t v)
  {
    // Calculate square root of an integer
    // Original code is isqrt of HEALPix C++
    return UInt_t(TMath::Sqrt(v + .5));
  }
  
//______________________________________________________________________________
  Double_t Modulo(Double_t v1, Double_t v2)
  {
    // Calculate modulo of two doubles
    // Original code is modulo of HEALPix C++
    return (v1 >= 0) ? ((v1 < v2) ? v1 : fmod(v1, v2)) : (fmod(v1, v2) + v2);
  }
  
//______________________________________________________________________________
  Int_t Modulo(Int_t v1, Int_t v2)
  {
    // Calculate modulo of two integers
    // Original code is modulo of HEALPix C++
    return (v1 >= 0) ? ((v1 < v2) ? v1 : (v1%v2)) : ((v1%v2) + v2);
  }
  
//______________________________________________________________________________
  Long_t Modulo(Long_t v1, Long_t v2)
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


#endif // T_HEAL_PIX
