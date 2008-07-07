// $Id: THealUtil.cxx,v 1.8 2008/07/07 07:00:14 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.

   This is a port of HEALPix C++ package to ROOT system.
   Original code is available at <http://healpix.jpl.nasa.gov> under GPL.
******************************************************************************/

#include <cmath>
#include <complex>
#include <vector>

#include <fitsio.h>

#include "TMath.h"

#include "THealAlm.h"
#include "THealFFT.h"
#include "THealPix.h"
#include "THealUtil.h"

namespace THealUtil {

Bool_t FitsReportError(Int_t status)
{
  if(status != 0){
      fits_report_error(stderr, status);
      return kFALSE;
  } // if
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t SaveToFits(const char* fname, const THealPix* hp)
{
  std::vector<const THealPix*> vec;
  vec.push_back((THealPix*)hp);

  return SaveToFits(fname, vec);
}

//______________________________________________________________________________
Bool_t SaveToFits(const char* fname, const std::vector<THealPix*>& hp)
{
  std::vector<const THealPix*> vec;
  for(UInt_t i = 0; i < vec.size(); i++){
    vec.push_back(hp[i]);
  } // i

  return SaveToFits(fname, vec);
}

//______________________________________________________________________________
Bool_t SaveToFits(const char* fname, const std::vector<const THealPix*>& hp)
{
  fitsfile* fptr;
  Int_t status = 0;
  
  fits_create_file(&fptr, fname, &status);
  if(!THealUtil::FitsReportError(status)){
    return kFALSE;
  } // if

  char** type = new char*[hp.size()];
  char** unit = new char*[hp.size()];
  char** form = new char*[hp.size()];
  for(UInt_t i = 0; i < hp.size(); i++){
    type[i] = new char[81];
    unit[i] = new char[81];
    form[i] = new char[81];
  } // i

  for(UInt_t i = 0; i < hp.size(); i++){
    strcpy(type[i], hp[i]->GetTitle());
    strcpy(unit[i], hp[i]->GetUnit().c_str());
    strcpy(form[i], Form("%d%s", hp[i]->GetNrows(), hp[i]->GetTypeString().c_str()));
  } // i

  fits_insert_btbl(fptr, 0, hp.size(), (char**)type, (char**)form, (char**)unit, "xtension", 0, &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(fptr, &status);
    return kFALSE;
  } // if

  fits_write_key(fptr, TSTRING, (char*)"PIXTYPE", (char*)"HEALPIX", (char*)"HEALPIX pixelisation", &status);
  fits_write_key(fptr, TSTRING, (char*)"ORDERING", (char*)hp[0]->GetSchemeString().c_str(), (char*)"Pixel ordering scheme, either RING or NESTED", &status);
  Long_t nside = hp[0]->GetNside();
  fits_write_key(fptr, TINT, (char*)"NSIDE", &nside, (char*)"Resolution parameter for HEALPIX", &status);
  Long_t first = 0;
  fits_write_key(fptr, TLONG, (char*)"FIRSTPIX", &first, (char*)"First pixel # (0 based)", &status);
  Long_t last = hp[0]->GetNpix() - 1;
  fits_write_key(fptr, TLONG, (char*)"LASTPIX", &last, (char*)"Last pixel # (0 based)", &status);
  fits_write_key(fptr, TSTRING, (char*)"INDXSCHM", (char*)"IMPLICIT", (char*)"Indexing: IMPLICIT or EXPLICIT", &status);
  Long_t grain = 0;
  fits_write_key(fptr, TLONG, (char*)"GRAIN", &grain, (char*)"Grain of pixel indexing", &status);
  if(!THealUtil::FitsReportError(status)){
    fits_close_file(fptr, &status);
    return kFALSE;
  } // if

  for(UInt_t i = 0; i < hp.size(); i++){
    if(hp[i]->GetType() == TDOUBLE){
      fits_write_col(fptr, hp[i]->GetType(), i + 1, 1, 1, hp[i]->GetNpix(),
		     ((THealPixD*)hp[i])->GetArray(), &status);
      Long_t nrows;
      fits_get_num_rows(fptr, &nrows, &status);
    } else if(hp[i]->GetType() == TFLOAT){
      fits_write_col(fptr, hp[i]->GetType(), i + 1, 1, 1, hp[i]->GetNpix(),
		     ((THealPixF*)hp[i])->GetArray(), &status);
      Long_t nrows;
      fits_get_num_rows(fptr, &nrows, &status);
    } else {
      // How should I do?
    } // if
    if(!THealUtil::FitsReportError(status)){
      fits_close_file(fptr, &status);
      return kFALSE;
    } // if
  } // i

  fits_close_file(fptr, &status);

  return kTRUE;
}

} // namespace THealUtil

namespace {
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
void read_work (const std::vector<std::complex<double> >& work, int nph, int mmax,
  bool shifted, const std::vector<std::complex<double> > &shiftarr,
  std::complex<double> *dataout)
  {
  int cnt2=0;
  for (int m=0; m<=mmax; ++m)
    {
    dataout[m] = work[cnt2];
    if (++cnt2==nph) cnt2=0;
    }
  if (shifted)
    for (int m=0; m<=mmax; ++m) dataout[m] *= shiftarr[m];
  }

//______________________________________________________________________________
template<typename T>
void fft_map2alm (int nph, int mmax, bool shifted,
  double weight, THealFFT::rfft &plan, const T *mapN, const T *mapS,
  std::complex<double> *phas_n, std::complex<double> *phas_s,
  const std::vector<std::complex<double> > &shiftarr, std::vector<std::complex<double> > &work)
  {
  for (int m=0; m<nph; ++m) work[m] = mapN[m]*weight;
  plan.forward_c(work);
  read_work (work, nph, mmax, shifted, shiftarr, phas_n);
  if (mapN!=mapS)
    {
    for (int m=0; m<nph; ++m) work[m] = mapS[m]*weight;
    plan.forward_c(work);
    read_work (work, nph, mmax, shifted, shiftarr, phas_s);
    }
  else
    for (int m=0; m<=mmax; ++m) phas_s[m]=0;
  }

//______________________________________________________________________________
void recalc_map2alm (int nph, int mmax, THealFFT::rfft &plan,
		     std::vector<std::complex<double> > &shiftarr)
{
  if (plan.size() == nph) return;
  plan.Set (nph);
  double f1 = TMath::Pi()/nph;
  for (int m=0; m<=mmax; ++m)
    {
    if (m<nph)
      shiftarr[m] = std::complex<double>(cos(m*f1),-sin(m*f1));
    else
      shiftarr[m]=-shiftarr[m-nph];
    }
  }

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
}

namespace THealUtil {

template<typename T> void Alm2Map(const THealAlm<T>& alm, THealPix& map)
{
  int lmax = alm.GetLmax(), mmax = alm.GetMmax(), nside = map.GetNside();

  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);

  std::complex<double> *b_north[chunksize], *b_south[chunksize];
  for(Int_t i = 0; i < chunksize; i++){
    b_north[i] = new std::complex<double>[mmax+1];
    b_south[i] = new std::complex<double>[mmax+1];
  } // i
  std::vector<double> cth(chunksize),sth(chunksize);

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=TMath::Min(llim+chunksize,2*nside);
    for (int ith=llim; ith<ulim; ++ith)
      {
      int nph, istart_north;
      bool shifted;
      map.GetRingInfo (ith+1,istart_north,nph, cth[ith-llim],sth[ith-llim],
                       shifted);
      }

#pragma omp parallel
{
    THealFFT::Ylmgen generator(lmax,mmax,1e-30);
    std::vector<double> Ylm;
    std::vector<std::complex<double> > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
      {
      for (int l=m; l<=lmax; ++l)
        alm_tmp[l]=alm(l,m);

      for (int ith=0; ith<ulim-llim; ++ith)
        {
        int l;
        generator.get_Ylm(cth[ith],sth[ith],m,Ylm,l);
        if (l<=lmax)
          {
	  std::complex<double> p1(0, 0), p2(0, 0);

          if ((l-m)&1) goto middle;
start:    p1.real() += alm_tmp[l].real()*Ylm[l]; p1.imag() += alm_tmp[l].imag()*Ylm[l];
          if (++l>lmax) goto end;
middle:   p2.real() += alm_tmp[l].real()*Ylm[l]; p2.imag() += alm_tmp[l].imag()*Ylm[l];
          if (++l<=lmax) goto start;
end:      b_north[ith][m] = p1+p2; b_south[ith][m] = p1-p2;
          }
        else
          {
          b_north[ith][m] = b_south[ith][m] = 0;
          }
        }
      }
} // end of parallel region

#pragma omp parallel
{
    std::vector<std::complex<double> > shiftarr(mmax+1), work(4*nside);
    THealFFT::rfft plan;
    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      int istart_north, istart_south, nph;
      double dum1, dum2;
      bool shifted;
      map.GetRingInfo (ith+1,istart_north,nph,dum1,dum2,shifted);
      istart_south = map.GetNpix()-istart_north-nph;
      recalc_alm2map (nph, mmax, plan, shiftarr);
      if(map.GetTypeString() == "D"){
      fft_alm2map (nph, mmax, shifted, plan,
        &(dynamic_cast<THealPixD&>(map))[istart_north],
        &(dynamic_cast<THealPixD&>(map))[istart_south],
        b_north[ith-llim], b_south[ith-llim], shiftarr, work);
      } else if(map.GetTypeString() == "F"){
      fft_alm2map (nph, mmax, shifted, plan,
        &(dynamic_cast<THealPixF&>(map))[istart_north],
        &(dynamic_cast<THealPixF&>(map))[istart_south],
        b_north[ith-llim], b_south[ith-llim], shiftarr, work);
      } // if
      }

} // end of parallel region
    }

  for(Int_t i = 0; i < chunksize; i++){
    delete b_north[i];
    delete b_south[i];
  } // i

  }

//______________________________________________________________________________
template<typename T>
void Map2Alm(const THealPix& map, THealAlm<T>& alm, Bool_t add_alm)
{
  std::vector<Double_t> weight(2*map.GetNside(), 1);
  Map2Alm(map, alm, weight, add_alm);
}

//______________________________________________________________________________
template<typename T>
void Map2Alm(const THealPix& map, THealAlm<T>& alm, const std::vector<Double_t>& weight, Bool_t add_alm)
{
  if(map.IsNested()){
    // to be modified
    return;
  } // if

  if((Int_t)weight.size() < 2*map.GetNside()){
    // to be modified
    return;
  } // if

  int lmax = alm.GetLmax(), mmax = alm.GetMmax(), nside = map.GetNside();

  int nchunks, chunksize;
  get_chunk_info(2*nside,nchunks,chunksize);

  std::complex<double> *phas_n[chunksize], *phas_s[chunksize];
  for(Int_t i = 0; i < chunksize; i++){
    phas_n[i] = new std::complex<double>[mmax+1];
    phas_s[i] = new std::complex<double>[mmax+1];
  } // i

  std::vector<double> cth(chunksize), sth(chunksize);
  double normfact = TMath::Pi()/(3*nside*nside);

  if (!add_alm) alm.SetToZero();

  for (int chunk=0; chunk<nchunks; ++chunk)
    {
    int llim=chunk*chunksize, ulim=TMath::Min(llim+chunksize,2*nside);
#pragma omp parallel
{
    std::vector<std::complex<double> > shiftarr(mmax+1), work(4*nside);
    THealFFT::rfft plan;

    int ith;
#pragma omp for schedule(dynamic,1)
    for (ith=llim; ith<ulim; ++ith)
      {
      int istart_north, istart_south, nph;
      bool shifted;
      map.GetRingInfo (ith+1,istart_north,nph,cth[ith-llim],sth[ith-llim],
                         shifted);
      istart_south = 12*nside*nside - istart_north - nph;
      recalc_map2alm (nph, mmax, plan, shiftarr);
      if(map.GetTypeString() == "D"){
	const THealPixD& tmp = dynamic_cast<const THealPixD&>(map);
	const Double_t* p1 = &(tmp.GetArray()[istart_north]);
	const Double_t* p2 = &(tmp.GetArray()[istart_south]);
      fft_map2alm (nph, mmax, shifted, weight[ith]*normfact, plan,
        p1, p2,	phas_n[ith-llim], phas_s[ith-llim], shiftarr, work);
      } else if(map.GetTypeString() == "F"){
	const THealPixF& tmp = dynamic_cast<const THealPixF&>(map);
	const Float_t* p1 = &(tmp.GetArray()[istart_north]);
	const Float_t* p2 = &(tmp.GetArray()[istart_south]);
      fft_map2alm (nph, mmax, shifted, weight[ith]*normfact, plan,
        p1, p2,	phas_n[ith-llim], phas_s[ith-llim], shiftarr, work);
      } // if
      }
} // end of parallel region

#pragma omp parallel
{
    THealFFT::Ylmgen generator(lmax,mmax,1e-30);
    std::vector<double> Ylm;
    std::vector<std::complex<double> > alm_tmp(lmax+1);
    int m;
#pragma omp for schedule(dynamic,1)
    for (m=0; m<=mmax; ++m)
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
} // end of parallel region
    }

  for(Int_t i = 0; i < chunksize; i++){
    delete phas_n[i];
    delete phas_s[i];
  } // i

  }

template void Alm2Map(const THealAlm<Float_t>& alm, THealPix& map);
template void Alm2Map(const THealAlm<Double_t>& alm, THealPix& map);
template void Map2Alm(const THealPix& map, THealAlm<Float_t>& alm, Bool_t add_alm);
template void Map2Alm(const THealPix& map, THealAlm<Double_t>& alm, Bool_t add_alm);
template void Map2Alm(const THealPix& map, THealAlm<Float_t>& alm, const std::vector<Double_t>& weight, Bool_t add_alm);
template void Map2Alm(const THealPix& map, THealAlm<Double_t>& alm, const std::vector<Double_t>& weight, Bool_t add_alm);
}
