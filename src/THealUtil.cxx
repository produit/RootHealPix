// $Id: THealUtil.cxx,v 1.1 2008/06/24 08:16:44 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#include <cmath>
#include <vector>

#include <fitsio.h>

#include "TMath.h"

#include "THealPix.h"
#include "THealUtil.h"

namespace THealUtil {
void mk_xy2pix(int *x2pix, int *y2pix) {
  /* =======================================================================
   * subroutine mk_xy2pix
   * =======================================================================
   * sets the array giving the number of the pixel lying in (x,y)
   * x and y are in {1,128}
   * the pixel number is in {0,128**2-1}
   *
   * if  i-1 = sum_p=0  b_p * 2^p
   * then ix = sum_p=0  b_p * 4^p
   * iy = 2*ix
   * ix + iy in {0, 128**2 -1}
   * =======================================================================
   */
  int i, K,IP,I,J,ID;
  
  for(i = 0; i < 127; i++) x2pix[i] = 0;
  for( I=1;I<=128;I++ ) {
    J  = I-1;//            !pixel numbers
    K  = 0;//
    IP = 1;//
    truc : if( J==0 ) {
      x2pix[I-1] = K;
      y2pix[I-1] = 2*K;
    }
    else {
      ID = (int)fmod(J,2);
      J  = J/2;
      K  = IP*ID+K;
      IP = IP*4;
      goto truc;
    }
  }     
  
}


Int_t AngToPixNest(Long_t nside, Double_t theta, Double_t phi)
{
  // to be modified

  double z, za, z0, tt, tp, tmp;
  int    face_num,jp,jm;
  long   ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  double piover2 = 0.5*TMath::Pi(), pi = TMath::Pi(), twopi = 2.0*TMath::Pi();
  int    ns_max = 8192;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  if( theta<0. || theta>pi ) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  if( !setup_done ) {
    mk_xy2pix(x2pix,y2pix);
    setup_done = 1;
  }
  
  z  = cos(theta);
  za = fabs(z);
  z0 = 2./3.;
  if( phi>=twopi ) phi = phi - twopi;
  if( phi<0. )    phi = phi + twopi;
  tt = phi / piover2; /* in [0,4[ */
  
  if( za<=z0 ) { /* equatorial region */
    
    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */
    
    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;
    
    if( ifp==ifm ) face_num = (int)fmod(ifp,4) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = (int)fmod(ifp,4); /* (half-)faces 0 to 3 */
    else face_num = (int)fmod(ifm,4) + 8;           /* (half-)faces 8 to 11 */
    
    ix = (int)fmod(jm, ns_max);
    iy = ns_max - (int)fmod(jp, ns_max) - 1;
  }
  else { /* polar region, za > 2/3 */
    
    ntt = (int)floor(tt);
    if( ntt>=4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt( 3.*(1. - za) ); /* in ]0,1] */
    
    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floor( ns_max * tp          * tmp ); 

    /* that one goes away of the closest pole */
    jm = (int)floor( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);
    
    /* finds the face and pixel's (x,y) */
    if( z>=0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }
  
  ix_low = (int)fmod(ix,128);
  ix_hi  =     ix/128;
  iy_low = (int)fmod(iy,128);
  iy_hi  =     iy/128;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow(ns_max/nside,2));     /* in {0, nside**2 - 1} */
  return (long)( ipf + face_num*pow(nside,2)); /* in {0, 12*nside**2 - 1} */
}

//______________________________________________________________________________
Int_t AngToPixRing(Long_t nside, Double_t theta, Double_t phi)
{
  // to be modified

  int nl2, nl4, ncap, npix, jp, jm, ipix1;
  double  z, za, tt, tp, tmp;
  int ir, ip, kshift;
  
  double piover2 = 0.5*TMath::Pi();
  double PI=TMath::Pi();
  double twopi=2.0*TMath::Pi();
  double z0=2.0/3.0;
  long ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  
  if( theta<0. || theta>PI) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  
  z = cos(theta);
  za = fabs(z);
  if( phi >= twopi)  phi = phi - twopi;
  if (phi < 0.)     phi = phi + twopi;
  tt = phi / piover2;//  ! in [0,4)
  
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
  npix  = 12*nside*nside;
  
  if( za <= z0 ) {
    
    jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
    jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
    
    ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 0;
    if (fmod(ir,2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
    
    ip = (int)floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;// ! in {1,4n}
    if( ip>nl4 ) ip = ip - nl4;
    
    ipix1 = ncap + nl4*(ir-1) + ip ;
   }
  else {
    
    tp = tt - floor(tt);//      !MOD(tt,1.d0)
    tmp = sqrt( 3.*(1. - za) );
    
    jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
    jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
    
    ir = jp + jm + 1;//        ! ring number counted from the closest pole
    ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
    if( ip>4*ir ) ip = ip - 4*ir;
    
    ipix1 = 2*ir*(ir-1) + ip;
    if( z<=0. ) {
      ipix1 = npix - 2*ir*(ir+1) + ip;
    }
  }
  return ipix1 - 1;// ! in {0, npix-1}
}

//______________________________________________________________________________
Bool_t FitsReportError(Int_t status)
{
  if(status != 0){
      fits_report_error(stderr, status);
      return kFALSE;
  } // if
  
  return kTRUE;
}
  
//______________________________________________________________________________
Bool_t SaveToFits(const char* fname, const std::vector<THealPix*>& hp)
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
  fits_write_key(fptr, TSTRING, (char*)"ORDERING", (char*)hp[0]->GetOrderingTypeString().c_str(), (char*)"Pixel ordering scheme, either RING or NESTED", &status);
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
