// $Id: THealAlmMapTools.cxx,v 1.6 2008/07/05 23:24:34 oxon Exp $

#include <cmath>

#include "TMath.h"

#include "THealAlm.h"
#include "THealAlmMapTools.h"
#include "THealPix.h"

namespace {
void get_chunk_info (int nrings, int &nchunks, int &chunksize)
  {
  nchunks = nrings/TMath::Max(100,nrings/10) + 1;
  chunksize = (nrings+nchunks-1)/nchunks;
  }

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

namespace THealAlmMapTools {

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

template<typename T> void Map2Alm(const THealPix& map, THealAlm<T>& alm, Bool_t add_alm)
{
  std::vector<Double_t> weight(2*map.GetNside(), 1);
  Map2Alm(map, alm, weight, add_alm);
}

  template<typename T> void Map2Alm(const THealPix& map, THealAlm<T>& alm, const std::vector<Double_t>& weight, Bool_t add_alm)
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
