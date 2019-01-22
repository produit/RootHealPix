// $Id: THealFFT.cxx,v 1.3 2008/07/07 07:19:05 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.

   This is a port of HEALPix C++ package to ROOT system.
   Original code is available at <http://healpix.jpl.nasa.gov> under GPL.
******************************************************************************/

#include <cmath>

#include "THealFFT.h"

// from fftpack.c of HEALPix C++
namespace {

void passf2(int ido, int l1, const double *cc, double *ch,
  const double *wa1)
  {
  int i, k, ah, ac;
  double ti2, tr2;

  if(ido<=2)
    {
    for(k=0; k<l1; k++)
      {
      ah=k*ido;
      ac=2*k*ido;
      ch[ah]=cc[ac]+cc[ac+ido];
      ch[ah+ido*l1]=cc[ac]-cc[ac+ido];
      ch[ah+1]=cc[ac+1]+cc[ac+ido+1];
      ch[ah+ido*l1+1]=cc[ac+1]-cc[ac+ido+1];
      }
    }
  else
    {
    for(k=0; k<l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ah=i+k*ido;
        ac=i+2*k*ido;
        ch[ah]=cc[ac]+cc[ac+ido];
        tr2=cc[ac]-cc[ac+ido];
        ch[ah+1]=cc[ac+1]+cc[ac+1+ido];
        ti2=cc[ac+1]-cc[ac+1+ido];
        ch[ah+l1*ido+1]=wa1[i]*ti2-wa1[i+1]*tr2;
        ch[ah+l1*ido]=wa1[i]*tr2+wa1[i+1]*ti2;
        }
      }
    }
  }

void passb2(int ido, int l1, const double *cc, double *ch,
  const double *wa1)
  {
  int i, k, ah, ac;
  double ti2, tr2;

  if(ido<=2)
    {
    for(k=0; k<l1; k++)
      {
      ah=k*ido;
      ac=2*k*ido;
      ch[ah]=cc[ac]+cc[ac+ido];
      ch[ah+ido*l1]=cc[ac]-cc[ac+ido];
      ch[ah+1]=cc[ac+1]+cc[ac+ido+1];
      ch[ah+ido*l1+1]=cc[ac+1]-cc[ac+ido+1];
      }
    }
  else
    {
    for(k=0; k<l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ah=i+k*ido;
        ac=i+2*k*ido;
        ch[ah]=cc[ac]+cc[ac+ido];
        tr2=cc[ac]-cc[ac+ido];
        ch[ah+1]=cc[ac+1]+cc[ac+1+ido];
        ti2=cc[ac+1]-cc[ac+1+ido];
        ch[ah+l1*ido+1]=wa1[i]*ti2+wa1[i+1]*tr2;
        ch[ah+l1*ido]=wa1[i]*tr2-wa1[i+1]*ti2;
        }
      }
    }
  }

void passf3(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2)
  {
  static const double taur=-0.5, taui=0.86602540378443864676;
  int i, k, ac, ah;
  double ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

  if(ido==2)
    {
    for(k=1; k<=l1; k++)
      {
      ac=(3*k-2)*ido;
      tr2=cc[ac]+cc[ac+ido];
      cr2=cc[ac-ido]+taur*tr2;
      ah=(k-1)*ido;
      ch[ah]=cc[ac-ido]+tr2;

      ti2=cc[ac+1]+cc[ac+ido+1];
      ci2=cc[ac-ido+1]+taur*ti2;
      ch[ah+1]=cc[ac-ido+1]+ti2;

      cr3=-taui*(cc[ac]-cc[ac+ido]);
      ci3=-taui*(cc[ac+1]-cc[ac+ido+1]);
      ch[ah+l1*ido]=cr2-ci3;
      ch[ah+2*l1*ido]=cr2+ci3;
      ch[ah+l1*ido+1]=ci2+cr3;
      ch[ah+2*l1*ido+1]=ci2-cr3;
      }
    }
  else
    {
    for(k=1; k<=l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ac=i+(3*k-2)*ido;
        tr2=cc[ac]+cc[ac+ido];
        cr2=cc[ac-ido]+taur*tr2;
        ah=i+(k-1)*ido;
        ch[ah]=cc[ac-ido]+tr2;
        ti2=cc[ac+1]+cc[ac+ido+1];
        ci2=cc[ac-ido+1]+taur*ti2;
        ch[ah+1]=cc[ac-ido+1]+ti2;
        cr3=-taui*(cc[ac]-cc[ac+ido]);
        ci3=-taui*(cc[ac+1]-cc[ac+ido+1]);
        dr2=cr2-ci3;
        dr3=cr2+ci3;
        di2=ci2+cr3;
        di3=ci2-cr3;
        ch[ah+l1*ido+1]=wa1[i]*di2-wa1[i+1]*dr2;
        ch[ah+l1*ido]=wa1[i]*dr2+wa1[i+1]*di2;
        ch[ah+2*l1*ido+1]=wa2[i]*di3-wa2[i+1]*dr3;
        ch[ah+2*l1*ido]=wa2[i]*dr3+wa2[i+1]*di3;
        }
      }
    }
  }

void passb3(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2)
  {
  static const double taur=-0.5, taui=0.86602540378443864676;
  int i, k, ac, ah;
  double ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

  if(ido==2)
    {
    for(k=1; k<=l1; k++)
      {
      ac=(3*k-2)*ido;
      tr2=cc[ac]+cc[ac+ido];
      cr2=cc[ac-ido]+taur*tr2;
      ah=(k-1)*ido;
      ch[ah]=cc[ac-ido]+tr2;

      ti2=cc[ac+1]+cc[ac+ido+1];
      ci2=cc[ac-ido+1]+taur*ti2;
      ch[ah+1]=cc[ac-ido+1]+ti2;

      cr3=taui*(cc[ac]-cc[ac+ido]);
      ci3=taui*(cc[ac+1]-cc[ac+ido+1]);
      ch[ah+l1*ido]=cr2-ci3;
      ch[ah+2*l1*ido]=cr2+ci3;
      ch[ah+l1*ido+1]=ci2+cr3;
      ch[ah+2*l1*ido+1]=ci2-cr3;
      }
    }
  else
    {
    for(k=1; k<=l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ac=i+(3*k-2)*ido;
        tr2=cc[ac]+cc[ac+ido];
        cr2=cc[ac-ido]+taur*tr2;
        ah=i+(k-1)*ido;
        ch[ah]=cc[ac-ido]+tr2;
        ti2=cc[ac+1]+cc[ac+ido+1];
        ci2=cc[ac-ido+1]+taur*ti2;
        ch[ah+1]=cc[ac-ido+1]+ti2;
        cr3=taui*(cc[ac]-cc[ac+ido]);
        ci3=taui*(cc[ac+1]-cc[ac+ido+1]);
        dr2=cr2-ci3;
        dr3=cr2+ci3;
        di2=ci2+cr3;
        di3=ci2-cr3;
        ch[ah+l1*ido+1]=wa1[i]*di2+wa1[i+1]*dr2;
        ch[ah+l1*ido]=wa1[i]*dr2-wa1[i+1]*di2;
        ch[ah+2*l1*ido+1]=wa2[i]*di3+wa2[i+1]*dr3;
        ch[ah+2*l1*ido]=wa2[i]*dr3-wa2[i+1]*di3;
        }
      }
    }
  }

void passf4(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3)
  {
  int i, k, ac, ah;
  double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
	
  if(ido==2)
    {
    for(k=0; k<l1; k++)
      {
      ac=4*k*ido+1;
      ti1=cc[ac]-cc[ac+2*ido];
      ti2=cc[ac]+cc[ac+2*ido];
      tr4=cc[ac+3*ido]-cc[ac+ido];
      ti3=cc[ac+ido]+cc[ac+3*ido];
      tr1=cc[ac-1]-cc[ac+2*ido-1];
      tr2=cc[ac-1]+cc[ac+2*ido-1];
      ti4=cc[ac+ido-1]-cc[ac+3*ido-1];
      tr3=cc[ac+ido-1]+cc[ac+3*ido-1];
      ah=k*ido;
      ch[ah]=tr2+tr3;
      ch[ah+2*l1*ido]=tr2-tr3;
      ch[ah+1]=ti2+ti3;
      ch[ah+2*l1*ido+1]=ti2-ti3;
      ch[ah+l1*ido]=tr1-tr4;
      ch[ah+3*l1*ido]=tr1+tr4;
      ch[ah+l1*ido+1]=ti1-ti4;
      ch[ah+3*l1*ido+1]=ti1+ti4;
      }
    }
  else
    {
    for(k=0; k<l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ac=i+1+4*k*ido;
        ti1=cc[ac]-cc[ac+2*ido];
        ti2=cc[ac]+cc[ac+2*ido];
        ti3=cc[ac+ido]+cc[ac+3*ido];
        tr4=cc[ac+3*ido]-cc[ac+ido];
        tr1=cc[ac-1]-cc[ac+2*ido-1];
        tr2=cc[ac-1]+cc[ac+2*ido-1];
        ti4=cc[ac+ido-1]-cc[ac+3*ido-1];
        tr3=cc[ac+ido-1]+cc[ac+3*ido-1];
        ah=i+k*ido;
        ch[ah]=tr2+tr3;
        cr3=tr2-tr3;
        ch[ah+1]=ti2+ti3;
        ci3=ti2-ti3;
        cr2=tr1-tr4;
        cr4=tr1+tr4;
        ci2=ti1-ti4;
        ci4=ti1+ti4;
        ch[ah+l1*ido]=wa1[i]*cr2+wa1[i+1]*ci2;
        ch[ah+l1*ido+1]=wa1[i]*ci2-wa1[i+1]*cr2;
        ch[ah+2*l1*ido]=wa2[i]*cr3+wa2[i+1]*ci3;
        ch[ah+2*l1*ido+1]=wa2[i]*ci3-wa2[i+1]*cr3;
        ch[ah+3*l1*ido]=wa3[i]*cr4+wa3[i+1]*ci4;
        ch[ah+3*l1*ido+1]=wa3[i]*ci4-wa3[i+1]*cr4;
        }
      }
    }
  }

void passb4(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3)
  {
  int i, k, ac, ah;
  double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
	
  if(ido==2)
    {
    for(k=0; k<l1; k++)
      {
      ac=4*k*ido+1;
      ti1=cc[ac]-cc[ac+2*ido];
      ti2=cc[ac]+cc[ac+2*ido];
      tr4=cc[ac+3*ido]-cc[ac+ido];
      ti3=cc[ac+ido]+cc[ac+3*ido];
      tr1=cc[ac-1]-cc[ac+2*ido-1];
      tr2=cc[ac-1]+cc[ac+2*ido-1];
      ti4=cc[ac+ido-1]-cc[ac+3*ido-1];
      tr3=cc[ac+ido-1]+cc[ac+3*ido-1];
      ah=k*ido;
      ch[ah]=tr2+tr3;
      ch[ah+2*l1*ido]=tr2-tr3;
      ch[ah+1]=ti2+ti3;
      ch[ah+2*l1*ido+1]=ti2-ti3;
      ch[ah+l1*ido]=tr1+tr4;
      ch[ah+3*l1*ido]=tr1-tr4;
      ch[ah+l1*ido+1]=ti1+ti4;
      ch[ah+3*l1*ido+1]=ti1-ti4;
      }
    }
  else
    {
    for(k=0; k<l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ac=i+1+4*k*ido;
        ti1=cc[ac]-cc[ac+2*ido];
        ti2=cc[ac]+cc[ac+2*ido];
        ti3=cc[ac+ido]+cc[ac+3*ido];
        tr4=cc[ac+3*ido]-cc[ac+ido];
        tr1=cc[ac-1]-cc[ac+2*ido-1];
        tr2=cc[ac-1]+cc[ac+2*ido-1];
        ti4=cc[ac+ido-1]-cc[ac+3*ido-1];
        tr3=cc[ac+ido-1]+cc[ac+3*ido-1];
        ah=i+k*ido;
        ch[ah]=tr2+tr3;
        cr3=tr2-tr3;
        ch[ah+1]=ti2+ti3;
        ci3=ti2-ti3;
        cr2=tr1+tr4;
        cr4=tr1-tr4;
        ci2=ti1+ti4;
        ci4=ti1-ti4;
        ch[ah+l1*ido]=wa1[i]*cr2-wa1[i+1]*ci2;
        ch[ah+l1*ido+1]=wa1[i]*ci2+wa1[i+1]*cr2;
        ch[ah+2*l1*ido]=wa2[i]*cr3-wa2[i+1]*ci3;
        ch[ah+2*l1*ido+1]=wa2[i]*ci3+wa2[i+1]*cr3;
        ch[ah+3*l1*ido]=wa3[i]*cr4-wa3[i+1]*ci4;
        ch[ah+3*l1*ido+1]=wa3[i]*ci4+wa3[i+1]*cr4;
        }
      }
    }
  }

void passf5(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3,
  const double *wa4)
  {
  static const double tr11= 0.3090169943749474241, ti11=0.95105651629515357212;
  static const double tr12=-0.8090169943749474241, ti12=0.58778525229247312917;
  int i, k, ac, ah;
  double ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4,
         ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

  if(ido==2)
    {
    for(k=1; k<=l1;++k)
      {
      ac=(5*k-4)*ido+1;
      ti5=cc[ac]-cc[ac+3*ido];
      ti2=cc[ac]+cc[ac+3*ido];
      ti4=cc[ac+ido]-cc[ac+2*ido];
      ti3=cc[ac+ido]+cc[ac+2*ido];
      tr5=cc[ac-1]-cc[ac+3*ido-1];
      tr2=cc[ac-1]+cc[ac+3*ido-1];
      tr4=cc[ac+ido-1]-cc[ac+2*ido-1];
      tr3=cc[ac+ido-1]+cc[ac+2*ido-1];
      ah=(k-1)*ido;
      ch[ah]=cc[ac-ido-1]+tr2+tr3;
      ch[ah+1]=cc[ac-ido]+ti2+ti3;
      cr2=cc[ac-ido-1]+tr11*tr2+tr12*tr3;
      ci2=cc[ac-ido]+tr11*ti2+tr12*ti3;
      cr3=cc[ac-ido-1]+tr12*tr2+tr11*tr3;
      ci3=cc[ac-ido]+tr12*ti2+tr11*ti3;
      cr5=-(ti11*tr5+ti12*tr4);
      ci5=-(ti11*ti5+ti12*ti4);
      cr4=-(ti12*tr5-ti11*tr4);
      ci4=-(ti12*ti5-ti11*ti4);
      ch[ah+l1*ido]=cr2-ci5;
      ch[ah+4*l1*ido]=cr2+ci5;
      ch[ah+l1*ido+1]=ci2+cr5;
      ch[ah+2*l1*ido+1]=ci3+cr4;
      ch[ah+2*l1*ido]=cr3-ci4;
      ch[ah+3*l1*ido]=cr3+ci4;
      ch[ah+3*l1*ido+1]=ci3-cr4;
      ch[ah+4*l1*ido+1]=ci2-cr5;
      }
    }
  else
    {
    for(k=1; k<=l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ac=i+1+(k*5-4)*ido;
        ti5=cc[ac]-cc[ac+3*ido];
        ti2=cc[ac]+cc[ac+3*ido];
        ti4=cc[ac+ido]-cc[ac+2*ido];
        ti3=cc[ac+ido]+cc[ac+2*ido];
        tr5=cc[ac-1]-cc[ac+3*ido-1];
        tr2=cc[ac-1]+cc[ac+3*ido-1];
        tr4=cc[ac+ido-1]-cc[ac+2*ido-1];
        tr3=cc[ac+ido-1]+cc[ac+2*ido-1];
        ah=i+(k-1)*ido;
        ch[ah]=cc[ac-ido-1]+tr2+tr3;
        ch[ah+1]=cc[ac-ido]+ti2+ti3;
        cr2=cc[ac-ido-1]+tr11*tr2+tr12*tr3;

        ci2=cc[ac-ido]+tr11*ti2+tr12*ti3;
        cr3=cc[ac-ido-1]+tr12*tr2+tr11*tr3;

        ci3=cc[ac-ido]+tr12*ti2+tr11*ti3;
        cr5=-(ti11*tr5+ti12*tr4);
        ci5=-(ti11*ti5+ti12*ti4);
        cr4=-(ti12*tr5-ti11*tr4);
        ci4=-(ti12*ti5-ti11*ti4);
        dr3=cr3-ci4;
        dr4=cr3+ci4;
        di3=ci3+cr4;
        di4=ci3-cr4;
        dr5=cr2+ci5;
        dr2=cr2-ci5;
        di5=ci2-cr5;
        di2=ci2+cr5;
        ch[ah+l1*ido]=wa1[i]*dr2+wa1[i+1]*di2;
        ch[ah+l1*ido+1]=wa1[i]*di2-wa1[i+1]*dr2;
        ch[ah+2*l1*ido]=wa2[i]*dr3+wa2[i+1]*di3;
        ch[ah+2*l1*ido+1]=wa2[i]*di3-wa2[i+1]*dr3;
        ch[ah+3*l1*ido]=wa3[i]*dr4+wa3[i+1]*di4;
        ch[ah+3*l1*ido+1]=wa3[i]*di4-wa3[i+1]*dr4;
        ch[ah+4*l1*ido]=wa4[i]*dr5+wa4[i+1]*di5;
        ch[ah+4*l1*ido+1]=wa4[i]*di5-wa4[i+1]*dr5;
        }
      }
    }
  }

void passb5(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3,
  const double *wa4)
  {
  static const double tr11= 0.3090169943749474241, ti11=0.95105651629515357212;
  static const double tr12=-0.8090169943749474241, ti12=0.58778525229247312917;
  int i, k, ac, ah;
  double ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4,
         ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

  if(ido==2)
    {
    for(k=1; k<=l1;++k)
      {
      ac=(5*k-4)*ido+1;
      ti5=cc[ac]-cc[ac+3*ido];
      ti2=cc[ac]+cc[ac+3*ido];
      ti4=cc[ac+ido]-cc[ac+2*ido];
      ti3=cc[ac+ido]+cc[ac+2*ido];
      tr5=cc[ac-1]-cc[ac+3*ido-1];
      tr2=cc[ac-1]+cc[ac+3*ido-1];
      tr4=cc[ac+ido-1]-cc[ac+2*ido-1];
      tr3=cc[ac+ido-1]+cc[ac+2*ido-1];
      ah=(k-1)*ido;
      ch[ah]=cc[ac-ido-1]+tr2+tr3;
      ch[ah+1]=cc[ac-ido]+ti2+ti3;
      cr2=cc[ac-ido-1]+tr11*tr2+tr12*tr3;
      ci2=cc[ac-ido]+tr11*ti2+tr12*ti3;
      cr3=cc[ac-ido-1]+tr12*tr2+tr11*tr3;
      ci3=cc[ac-ido]+tr12*ti2+tr11*ti3;
      cr5=ti11*tr5+ti12*tr4;
      ci5=ti11*ti5+ti12*ti4;
      cr4=ti12*tr5-ti11*tr4;
      ci4=ti12*ti5-ti11*ti4;
      ch[ah+l1*ido]=cr2-ci5;
      ch[ah+4*l1*ido]=cr2+ci5;
      ch[ah+l1*ido+1]=ci2+cr5;
      ch[ah+2*l1*ido+1]=ci3+cr4;
      ch[ah+2*l1*ido]=cr3-ci4;
      ch[ah+3*l1*ido]=cr3+ci4;
      ch[ah+3*l1*ido+1]=ci3-cr4;
      ch[ah+4*l1*ido+1]=ci2-cr5;
      }
    }
  else
    {
    for(k=1; k<=l1; k++)
      {
      for(i=0; i<ido-1; i+=2)
        {
        ac=i+1+(k*5-4)*ido;
        ti5=cc[ac]-cc[ac+3*ido];
        ti2=cc[ac]+cc[ac+3*ido];
        ti4=cc[ac+ido]-cc[ac+2*ido];
        ti3=cc[ac+ido]+cc[ac+2*ido];
        tr5=cc[ac-1]-cc[ac+3*ido-1];
        tr2=cc[ac-1]+cc[ac+3*ido-1];
        tr4=cc[ac+ido-1]-cc[ac+2*ido-1];
        tr3=cc[ac+ido-1]+cc[ac+2*ido-1];
        ah=i+(k-1)*ido;
        ch[ah]=cc[ac-ido-1]+tr2+tr3;
        ch[ah+1]=cc[ac-ido]+ti2+ti3;
        cr2=cc[ac-ido-1]+tr11*tr2+tr12*tr3;

        ci2=cc[ac-ido]+tr11*ti2+tr12*ti3;
        cr3=cc[ac-ido-1]+tr12*tr2+tr11*tr3;

        ci3=cc[ac-ido]+tr12*ti2+tr11*ti3;
        cr5=ti11*tr5+ti12*tr4;
        ci5=ti11*ti5+ti12*ti4;
        cr4=ti12*tr5-ti11*tr4;
        ci4=ti12*ti5-ti11*ti4;
        dr3=cr3-ci4;
        dr4=cr3+ci4;
        di3=ci3+cr4;
        di4=ci3-cr4;
        dr5=cr2+ci5;
        dr2=cr2-ci5;
        di5=ci2-cr5;
        di2=ci2+cr5;
        ch[ah+l1*ido]=wa1[i]*dr2-wa1[i+1]*di2;
        ch[ah+l1*ido+1]=wa1[i]*di2+wa1[i+1]*dr2;
        ch[ah+2*l1*ido]=wa2[i]*dr3-wa2[i+1]*di3;
        ch[ah+2*l1*ido+1]=wa2[i]*di3+wa2[i+1]*dr3;
        ch[ah+3*l1*ido]=wa3[i]*dr4-wa3[i+1]*di4;
        ch[ah+3*l1*ido+1]=wa3[i]*di4+wa3[i+1]*dr4;
        ch[ah+4*l1*ido]=wa4[i]*dr5-wa4[i+1]*di5;
        ch[ah+4*l1*ido+1]=wa4[i]*di5+wa4[i+1]*dr5;
        }
      }
    }
  }

void passfg(int *nac, int ido, int ip, int l1, int idl1,
  double *cc, double *ch, const double *wa)
  {
  int idij, idlj, ipph, i, j, k, l, jc, lc, ik, idj, idl, inc, idp, idx;
  double wai, war;

  ipph=(ip+1)/ 2;
  idp=ip*ido;
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    for(k=0; k<l1; k++)
      {
      for(i=0; i<ido; i++)
        {
        ch[i+(k+j*l1)*ido] = cc[i+(j+k*ip)*ido]+cc[i+(jc+k*ip)*ido];
        ch[i+(k+jc*l1)*ido]= cc[i+(j+k*ip)*ido]-cc[i+(jc+k*ip)*ido];
        }
      }
    }
  for(k=0; k<l1; k++)
    memcpy (ch+k*ido, cc+k*ip*ido, ido*sizeof(double));

  idl=2-ido;
  inc=0;
  for(l=1; l<ipph; l++)
    {
    lc=ip-l;
    idl+=ido;
    for(ik=0; ik<idl1; ik++)
      {
      cc[ik+l*idl1]=ch[ik]+wa[idl-2]*ch[ik+idl1];
      cc[ik+lc*idl1]=-wa[idl-1]*ch[ik+(ip-1)*idl1];
      }
    idlj=idl;
    inc+=ido;
    for(j=2; j<ipph; j++)
      {
      jc=ip-j;
      idlj+=inc;
      if(idlj>idp)
        idlj-=idp;
      war=wa[idlj-2];
      wai=wa[idlj-1];
      for(ik=0; ik<idl1; ik++)
        {
        cc[ik+l*idl1]+=war*ch[ik+j*idl1];
        cc[ik+lc*idl1]-=wai*ch[ik+jc*idl1];
        }
      }
    }
  for(j=1; j<ipph; j++)
    for(ik=0; ik<idl1; ik++)
      ch[ik]+=ch[ik+j*idl1];
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    for(ik=1; ik<idl1; ik+=2)
      {
      ch[ik-1+j *idl1]=cc[ik-1+j*idl1]-cc[ik  +jc*idl1];
      ch[ik-1+jc*idl1]=cc[ik-1+j*idl1]+cc[ik  +jc*idl1];
      ch[ik  +j *idl1]=cc[ik  +j*idl1]+cc[ik-1+jc*idl1];
      ch[ik  +jc*idl1]=cc[ik  +j*idl1]-cc[ik-1+jc*idl1];
      }
    }
  *nac=1;
  if(ido==2)
    return;
  *nac=0;
  for(ik=0; ik<idl1; ik++)
    cc[ik]=ch[ik];
  for(j=1; j<ip; j++)
    {
    for(k=0; k<l1; k++)
      {
      cc[(k+j*l1)*ido  ]=ch[(k+j*l1)*ido  ];
      cc[(k+j*l1)*ido+1]=ch[(k+j*l1)*ido+1];
      }
    }

  idj=2-ido;
  for(j=1; j<ip; j++)
    {
    idj+=ido;
    for(k=0; k<l1; k++)
      {
      idij=idj;
      for(i=3; i<ido; i+=2)
        {
        idij+=2;
        idx = (k+j*l1)*ido;
        cc[i-1+idx] = wa[idij-2]*ch[i-1+idx]+wa[idij-1]*ch[i  +idx];
        cc[i  +idx] = wa[idij-2]*ch[i  +idx]-wa[idij-1]*ch[i-1+idx];
        }
      }
    }
  }

void passbg(int *nac, int ido, int ip, int l1, int idl1,
  double *cc, double *ch, const double *wa)
  {
  int idij, idlj, ipph, i, j, k, l, jc, lc, ik, idj, idl, inc, idp, idx;
  double wai, war;

  ipph=(ip+1)/ 2;
  idp=ip*ido;
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    for(k=0; k<l1; k++)
      {
      for(i=0; i<ido; i++)
        {
        ch[i+(k+j*l1)*ido] = cc[i+(j+k*ip)*ido]+cc[i+(jc+k*ip)*ido];
        ch[i+(k+jc*l1)*ido]= cc[i+(j+k*ip)*ido]-cc[i+(jc+k*ip)*ido];
        }
      }
    }
  for(k=0; k<l1; k++)
    memcpy (ch+k*ido, cc+k*ip*ido, ido*sizeof(double));

  idl=2-ido;
  inc=0;
  for(l=1; l<ipph; l++)
    {
    lc=ip-l;
    idl+=ido;
    for(ik=0; ik<idl1; ik++)
      {
      cc[ik+l*idl1]=ch[ik]+wa[idl-2]*ch[ik+idl1];
      cc[ik+lc*idl1]=wa[idl-1]*ch[ik+(ip-1)*idl1];
      }
    idlj=idl;
    inc+=ido;
    for(j=2; j<ipph; j++)
      {
      jc=ip-j;
      idlj+=inc;
      if(idlj>idp)
        idlj-=idp;
      war=wa[idlj-2];
      wai=wa[idlj-1];
      for(ik=0; ik<idl1; ik++)
        {
        cc[ik+l*idl1]+=war*ch[ik+j*idl1];
        cc[ik+lc*idl1]+=wai*ch[ik+jc*idl1];
        }
      }
    }
  for(j=1; j<ipph; j++)
    for(ik=0; ik<idl1; ik++)
      ch[ik]+=ch[ik+j*idl1];
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    for(ik=1; ik<idl1; ik+=2)
      {
      ch[ik-1+j *idl1]=cc[ik-1+j*idl1]-cc[ik  +jc*idl1];
      ch[ik-1+jc*idl1]=cc[ik-1+j*idl1]+cc[ik  +jc*idl1];
      ch[ik  +j *idl1]=cc[ik  +j*idl1]+cc[ik-1+jc*idl1];
      ch[ik  +jc*idl1]=cc[ik  +j*idl1]-cc[ik-1+jc*idl1];
      }
    }
  *nac=1;
  if(ido==2)
    return;
  *nac=0;
  for(ik=0; ik<idl1; ik++)
    cc[ik]=ch[ik];
  for(j=1; j<ip; j++)
    {
    for(k=0; k<l1; k++)
      {
      cc[(k+j*l1)*ido  ]=ch[(k+j*l1)*ido  ];
      cc[(k+j*l1)*ido+1]=ch[(k+j*l1)*ido+1];
      }
    }

  idj=2-ido;
  for(j=1; j<ip; j++)
    {
    idj+=ido;
    for(k=0; k<l1; k++)
      {
      idij=idj;
      for(i=3; i<ido; i+=2)
        {
        idij+=2;
        idx = (k+j*l1)*ido;
        cc[i-1+idx] = wa[idij-2]*ch[i-1+idx]-wa[idij-1]*ch[i  +idx];
        cc[i  +idx] = wa[idij-2]*ch[i  +idx]+wa[idij-1]*ch[i-1+idx];
        }
      }
    }
  }


void radf2 (int ido, int l1, const double *cc, double *ch,
  const double *wa1)
  {
  int i, k, ic;
  double ti2, tr2;

  for(k=0; k<l1; k++)
    {
    ch[2*k*ido] = cc[k*ido]+cc[(k+l1)*ido];
    ch[(2*k+1)*ido+ido-1] = cc[k*ido]-cc[(k+l1)*ido];
    }
  if(ido<2)
    return;
  if(ido !=2)
    {
    for(k=0; k<l1; k++)
      {
      for(i=2; i<ido; i+=2)
        {
        ic=ido-i;
        tr2=wa1[i-2]*cc[i-1+(k+l1)*ido]+wa1[i-1]*cc[i+(k+l1)*ido];
        ti2=wa1[i-2]*cc[i+(k+l1)*ido]-wa1[i-1]*cc[i-1+(k+l1)*ido];
        ch[i+2*k*ido]=cc[i+k*ido]+ti2;
        ch[ic+(2*k+1)*ido]=ti2-cc[i+k*ido];
        ch[i-1+2*k*ido]=cc[i-1+k*ido]+tr2;
        ch[ic-1+(2*k+1)*ido]=cc[i-1+k*ido]-tr2;
        }
      }
    if(ido%2==1)
      return;
    }
  for(k=0; k<l1; k++)
    {
    ch[(2*k+1)*ido] = -cc[ido-1+(k+l1)*ido];
    ch[ido-1+2*k*ido] = cc[ido-1+k*ido];
    }
  }

void radb2(int ido, int l1, const double *cc, double *ch,
  const double *wa1)
  {
  int i, k, ic;
  double ti2, tr2;

  for(k=0; k<l1; k++)
    {
    ch[k*ido] = cc[2*k*ido]+cc[ido-1+(2*k+1)*ido];
    ch[(k+l1)*ido] = cc[2*k*ido]-cc[ido-1+(2*k+1)*ido];
    }
  if(ido<2)
    return;
  if(ido !=2)
    {
    for(k=0; k<l1;++k)
      {
      for(i=2; i<ido; i+=2)
        {
        ic=ido-i;
        ch[i-1+k*ido] = cc[i-1+2*k*ido]+cc[ic-1+(2*k+1)*ido];
        tr2 = cc[i-1+2*k*ido]-cc[ic-1+(2*k+1)*ido];
        ch[i+k*ido] = cc[i+2*k*ido]-cc[ic+(2*k+1)*ido];
        ti2 = cc[i+(2*k)*ido]+cc[ic+(2*k+1)*ido];
        ch[i-1+(k+l1)*ido] = wa1[i-2]*tr2-wa1[i-1]*ti2;
        ch[i+(k+l1)*ido] = wa1[i-2]*ti2+wa1[i-1]*tr2;
        }
      }
    if(ido%2==1)
      return;
    }
  for(k=0; k<l1; k++)
    {
    ch[ido-1+k*ido]=2*cc[ido-1+2*k*ido];
    ch[ido-1+(k+l1)*ido]=-2*cc[(2*k+1)*ido];
    }
  }

void radf3(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2)
  {
  static const double taur=-0.5, taui=0.86602540378443864676;
  int i, k, ic;
  double ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;

  for(k=0; k<l1; k++)
    {
    cr2=cc[(k+l1)*ido]+cc[(k+2*l1)*ido];
    ch[3*k*ido]=cc[k*ido]+cr2;
    ch[(3*k+2)*ido]=taui*(cc[(k+l1*2)*ido]-cc[(k+l1)*ido]);
    ch[ido-1+(3*k+1)*ido]=cc[k*ido]+taur*cr2;
    }
  if(ido==1)
    return;
  for(k=0; k<l1; k++)
    {
    for(i=2; i<ido; i+=2)
      {
      ic=ido-i;
      dr2=wa1[i-2]*cc[i-1+(k+l1)*ido]+
          wa1[i-1]*cc[i+(k+l1)*ido];
      di2=wa1[i-2]*cc[i+(k+l1)*ido]-wa1[i-1]*cc[i-1+(k+l1)*ido];
      dr3=wa2[i-2]*cc[i-1+(k+l1*2)*ido]+wa2[i-1]*cc[i+(k+l1*2)*ido];
      di3=wa2[i-2]*cc[i+(k+l1*2)*ido]-wa2[i-1]*cc[i-1+(k+l1*2)*ido];
      cr2=dr2+dr3;
      ci2=di2+di3;
      ch[i-1+3*k*ido]=cc[i-1+k*ido]+cr2;
      ch[i+3*k*ido]=cc[i+k*ido]+ci2;
      tr2=cc[i-1+k*ido]+taur*cr2;
      ti2=cc[i+k*ido]+taur*ci2;
      tr3=taui*(di2-di3);
      ti3=taui*(dr3-dr2);
      ch[i-1+(3*k+2)*ido]=tr2+tr3;
      ch[ic-1+(3*k+1)*ido]=tr2-tr3;
      ch[i+(3*k+2)*ido]=ti2+ti3;
      ch[ic+(3*k+1)*ido]=ti3-ti2;
      }
    }
  }

void radb3(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2)
  {
  static const double taur=-0.5, taui=0.86602540378443864676;
  int i, k, ic;
  double ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

  for(k=0; k<l1; k++)
    {
    tr2=2*cc[ido-1+(3*k+1)*ido];
    cr2=cc[3*k*ido]+taur*tr2;
    ch[k*ido]=cc[3*k*ido]+tr2;
    ci3=2*taui*cc[(3*k+2)*ido];
    ch[(k+l1)*ido]=cr2-ci3;
    ch[(k+2*l1)*ido]=cr2+ci3;
    }
  if(ido==1)
    return;
  for(k=0; k<l1; k++)
    {
    for(i=2; i<ido; i+=2)
      {
      ic=ido-i;
      tr2=cc[i-1+(3*k+2)*ido]+cc[ic-1+(3*k+1)*ido];
      cr2=cc[i-1+3*k*ido]+taur*tr2;
      ch[i-1+k*ido]=cc[i-1+3*k*ido]+tr2;
      ti2=cc[i+(3*k+2)*ido]-cc[ic+(3*k+1)*ido];
      ci2=cc[i+3*k*ido]+taur*ti2;
      ch[i+k*ido]=cc[i+3*k*ido]+ti2;
      cr3=taui*(cc[i-1+(3*k+2)*ido]-cc[ic-1+(3*k+1)*ido]);
      ci3=taui*(cc[i+(3*k+2)*ido]+cc[ic+(3*k+1)*ido]);
      dr2=cr2-ci3;
      dr3=cr2+ci3;
      di2=ci2+cr3;
      di3=ci2-cr3;
      ch[i-1+(k+l1)*ido]=wa1[i-2]*dr2-wa1[i-1]*di2;
      ch[i+(k+l1)*ido]=wa1[i-2]*di2+wa1[i-1]*dr2;
      ch[i-1+(k+2*l1)*ido]=wa2[i-2]*dr3-wa2[i-1]*di3;
      ch[i+(k+2*l1)*ido]=wa2[i-2]*di3+wa2[i-1]*dr3;
      }
    }
  }

void radf4(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3)
  {
  static const double hsqt2=0.70710678118654752440;
  int i, k, ic;
  double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

  for(k=0; k<l1; k++)
    {
    tr1=cc[(k+l1)*ido]+cc[(k+3*l1)*ido];
    tr2=cc[k*ido]+cc[(k+2*l1)*ido];
    ch[4*k*ido]=tr1+tr2;
    ch[ido-1+(4*k+3)*ido]=tr2-tr1;
    ch[ido-1+(4*k+1)*ido]=cc[k*ido]-cc[(k+2*l1)*ido];
    ch[(4*k+2)*ido]=cc[(k+3*l1)*ido]-cc[(k+l1)*ido];
    }
  if(ido<2)
    return;
  if(ido !=2)
    {
    for(k=0; k<l1; k++)
      {
      for(i=2; i<ido; i+=2)
        {
        ic=ido-i;
        cr2=wa1[i-2]*cc[i-1+(k+l1)*ido]+wa1[i-1]*cc[i+(k+l1)*ido];
        ci2=wa1[i-2]*cc[i+(k+l1)*ido]-wa1[i-1]*cc[i-1+(k+l1)*ido];
        cr3=wa2[i-2]*cc[i-1+(k+2*l1)*ido]+wa2[i-1]*cc[i+(k+2*l1)*ido];
        ci3=wa2[i-2]*cc[i+(k+2*l1)*ido]-wa2[i-1]*cc[i-1+(k+2*l1)*ido];
        cr4=wa3[i-2]*cc[i-1+(k+3*l1)*ido]+wa3[i-1]*cc[i+(k+3*l1)*ido];
        ci4=wa3[i-2]*cc[i+(k+3*l1)*ido]-wa3[i-1]*cc[i-1+(k+3*l1)*ido];
        tr1=cr2+cr4;
        tr4=cr4-cr2;
        ti1=ci2+ci4;
        ti4=ci2-ci4;
        ti2=cc[i+k*ido]+ci3;
        ti3=cc[i+k*ido]-ci3;
        tr2=cc[i-1+k*ido]+cr3;
        tr3=cc[i-1+k*ido]-cr3;
        ch[i-1+4*k*ido]=tr1+tr2;
        ch[ic-1+(4*k+3)*ido]=tr2-tr1;
        ch[i+4*k*ido]=ti1+ti2;
        ch[ic+(4*k+3)*ido]=ti1-ti2;
        ch[i-1+(4*k+2)*ido]=ti4+tr3;
        ch[ic-1+(4*k+1)*ido]=tr3-ti4;
        ch[i+(4*k+2)*ido]=tr4+ti3;
        ch[ic+(4*k+1)*ido]=tr4-ti3;
        }
      }
    if(ido%2==1)
      return;
    }
  for(k=0; k<l1; k++)
    {
    ti1=-hsqt2*(cc[ido-1+(k+l1)*ido]+cc[ido-1+(k+3*l1)*ido]);
    tr1=hsqt2*(cc[ido-1+(k+l1)*ido]-cc[ido-1+(k+3*l1)*ido]);
    ch[ido-1+4*k*ido]=tr1+cc[ido-1+k*ido];
    ch[ido-1+(4*k+2)*ido]=cc[ido-1+k*ido]-tr1;
    ch[(4*k+1)*ido]=ti1-cc[ido-1+(k+2*l1)*ido];
    ch[(4*k+3)*ido]=ti1+cc[ido-1+(k+2*l1)*ido];
    }
  }

void radb4(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3)
  {
  static const double sqrt2=1.41421356237309504880;
  int i, k, ic;
  double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

  for(k=0; k<l1; k++)
    {
    tr1=cc[4*k*ido]-cc[ido-1+(4*k+3)*ido];
    tr2=cc[4*k*ido]+cc[ido-1+(4*k+3)*ido];
    tr3=cc[ido-1+(4*k+1)*ido]+cc[ido-1+(4*k+1)*ido];
    tr4=cc[(4*k+2)*ido]+cc[(4*k+2)*ido];
    ch[k*ido]=tr2+tr3;
    ch[(k+l1)*ido]=tr1-tr4;
    ch[(k+2*l1)*ido]=tr2-tr3;
    ch[(k+3*l1)*ido]=tr1+tr4;
    }
  if(ido<2)
    return;
  if(ido !=2)
    {
    for(k=0; k<l1;++k)
      {
      for(i=2; i<ido; i+=2)
        {
        ic=ido-i;
        ti1=cc[i+4*k*ido]+cc[ic+(4*k+3)*ido];
        ti2=cc[i+4*k*ido]-cc[ic+(4*k+3)*ido];
        ti3=cc[i+(4*k+2)*ido]-cc[ic+(4*k+1)*ido];
        tr4=cc[i+(4*k+2)*ido]+cc[ic+(4*k+1)*ido];
        tr1=cc[i-1+4*k*ido]-cc[ic-1+(4*k+3)*ido];
        tr2=cc[i-1+4*k*ido]+cc[ic-1+(4*k+3)*ido];
        ti4=cc[i-1+(4*k+2)*ido]-cc[ic-1+(4*k+1)*ido];
        tr3=cc[i-1+(4*k+2)*ido]+cc[ic-1+(4*k+1)*ido];
        ch[i-1+k*ido]=tr2+tr3;
        cr3=tr2-tr3;
        ch[i+k*ido]=ti2+ti3;
        ci3=ti2-ti3;
        cr2=tr1-tr4;
        cr4=tr1+tr4;
        ci2=ti1+ti4;
        ci4=ti1-ti4;
        ch[i-1+(k+l1)*ido]=wa1[i-2]*cr2-wa1[i-1]*ci2;
        ch[i+(k+l1)*ido]=wa1[i-2]*ci2+wa1[i-1]*cr2;
        ch[i-1+(k+2*l1)*ido]=wa2[i-2]*cr3-wa2[i-1]*ci3;
        ch[i+(k+2*l1)*ido]=wa2[i-2]*ci3+wa2[i-1]*cr3;
        ch[i-1+(k+3*l1)*ido]=wa3[i-2]*cr4-wa3[i-1]*ci4;
        ch[i+(k+3*l1)*ido]=wa3[i-2]*ci4+wa3[i-1]*cr4;
        }
      }
    if(ido%2==1)
      return;
    }
  for(k=0; k<l1; k++)
    {
    ti1=cc[(4*k+1)*ido]+cc[(4*k+3)*ido];
    ti2=cc[(4*k+3)*ido]-cc[(4*k+1)*ido];
    tr1=cc[ido-1+4*k*ido]-cc[ido-1+(4*k+2)*ido];
    tr2=cc[ido-1+4*k*ido]+cc[ido-1+(4*k+2)*ido];
    ch[ido-1+k*ido]=tr2+tr2;
    ch[ido-1+(k+l1)*ido]=sqrt2*(tr1-ti1);
    ch[ido-1+(k+2*l1)*ido]=ti2+ti2;
    ch[ido-1+(k+3*l1)*ido]=-sqrt2*(tr1+ti1);
    }
  }

void radf5(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3, const double *wa4)
  {
  static const double tr11=0.3090169943749474241;
  static const double ti11=0.95105651629515357212;
  static const double tr12=-0.8090169943749474241;
  static const double ti12=0.58778525229247312917;
  int i, k, ic;
  double ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3,
         dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;

  for(k=0; k<l1; k++)
    {
    cr2=cc[(k+4*l1)*ido]+cc[(k+l1)*ido];
    ci5=cc[(k+4*l1)*ido]-cc[(k+l1)*ido];
    cr3=cc[(k+3*l1)*ido]+cc[(k+2*l1)*ido];
    ci4=cc[(k+3*l1)*ido]-cc[(k+2*l1)*ido];
    ch[5*k*ido]=cc[k*ido]+cr2+cr3;
    ch[ido-1+(5*k+1)*ido]=cc[k*ido]+tr11*cr2+tr12*cr3;
    ch[(5*k+2)*ido]=ti11*ci5+ti12*ci4;
    ch[ido-1+(5*k+3)*ido]=cc[k*ido]+tr12*cr2+tr11*cr3;
    ch[(5*k+4)*ido]=ti12*ci5-ti11*ci4;
    }
  if(ido==1)
    return;
  for(k=0; k<l1;++k)
    {
    for(i=2; i<ido; i+=2)
      {
      ic=ido-i;
      dr2=wa1[i-2]*cc[i-1+(k+l1)*ido]+wa1[i-1]*cc[i+(k+l1)*ido];
      di2=wa1[i-2]*cc[i+(k+l1)*ido]-wa1[i-1]*cc[i-1+(k+l1)*ido];
      dr3=wa2[i-2]*cc[i-1+(k+2*l1)*ido]+wa2[i-1]*cc[i+(k+2*l1)*ido];
      di3=wa2[i-2]*cc[i+(k+2*l1)*ido]-wa2[i-1]*cc[i-1+(k+2*l1)*ido];
      dr4=wa3[i-2]*cc[i-1+(k+3*l1)*ido]+wa3[i-1]*cc[i+(k+3*l1)*ido];
      di4=wa3[i-2]*cc[i+(k+3*l1)*ido]-wa3[i-1]*cc[i-1+(k+3*l1)*ido];
      dr5=wa4[i-2]*cc[i-1+(k+4*l1)*ido]+wa4[i-1]*cc[i+(k+4*l1)*ido];
      di5=wa4[i-2]*cc[i+(k+4*l1)*ido]-wa4[i-1]*cc[i-1+(k+4*l1)*ido];
      cr2=dr2+dr5;
      ci5=dr5-dr2;
      cr5=di2-di5;
      ci2=di2+di5;
      cr3=dr3+dr4;
      ci4=dr4-dr3;
      cr4=di3-di4;
      ci3=di3+di4;
      ch[i-1+5*k*ido]=cc[i-1+k*ido]+cr2+cr3;
      ch[i+5*k*ido]=cc[i+k*ido]+ci2+ci3;
      tr2=cc[i-1+k*ido]+tr11*cr2+tr12*cr3;
      ti2=cc[i+k*ido]+tr11*ci2+tr12*ci3;
      tr3=cc[i-1+k*ido]+tr12*cr2+tr11*cr3;
      ti3=cc[i+k*ido]+tr12*ci2+tr11*ci3;
      tr5=ti11*cr5+ti12*cr4;
      ti5=ti11*ci5+ti12*ci4;
      tr4=ti12*cr5-ti11*cr4;
      ti4=ti12*ci5-ti11*ci4;
      ch[i-1+(5*k+2)*ido]=tr2+tr5;
      ch[ic-1+(5*k+1)*ido]=tr2-tr5;
      ch[i+(5*k+2)*ido]=ti2+ti5;
      ch[ic+(5*k+1)*ido]=ti5-ti2;
      ch[i-1+(5*k+4)*ido]=tr3+tr4;
      ch[ic-1+(5*k+3)*ido]=tr3-tr4;
      ch[i+(5*k+4)*ido]=ti3+ti4;
      ch[ic+(5*k+3)*ido]=ti4-ti3;
      }
    }
  }

void radb5(int ido, int l1, const double *cc, double *ch,
  const double *wa1, const double *wa2, const double *wa3, const double *wa4)
  {
  static const double tr11=0.3090169943749474241;
  static const double ti11=0.95105651629515357212;
  static const double tr12=-0.8090169943749474241;
  static const double ti12=0.58778525229247312917;
  int i, k, ic;
  double ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4,
          ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

  for(k=0; k<l1; k++)
    {
    ti5=2*cc[(5*k+2)*ido];
    ti4=2*cc[(5*k+4)*ido];
    tr2=2*cc[ido-1+(5*k+1)*ido];
    tr3=2*cc[ido-1+(5*k+3)*ido];
    ch[k*ido]=cc[5*k*ido]+tr2+tr3;
    cr2=cc[5*k*ido]+tr11*tr2+tr12*tr3;
    cr3=cc[5*k*ido]+tr12*tr2+tr11*tr3;
    ci5=ti11*ti5+ti12*ti4;
    ci4=ti12*ti5-ti11*ti4;
    ch[(k+l1)*ido]=cr2-ci5;
    ch[(k+2*l1)*ido]=cr3-ci4;
    ch[(k+3*l1)*ido]=cr3+ci4;
    ch[(k+4*l1)*ido]=cr2+ci5;
    }
  if(ido==1)
    return;
  for(k=0; k<l1;++k)
    {
    for(i=2; i<ido; i+=2)
      {
      ic=ido-i;
      ti5=cc[i+(5*k+2)*ido]+cc[ic+(5*k+1)*ido];
      ti2=cc[i+(5*k+2)*ido]-cc[ic+(5*k+1)*ido];
      ti4=cc[i+(5*k+4)*ido]+cc[ic+(5*k+3)*ido];
      ti3=cc[i+(5*k+4)*ido]-cc[ic+(5*k+3)*ido];
      tr5=cc[i-1+(5*k+2)*ido]-cc[ic-1+(5*k+1)*ido];
      tr2=cc[i-1+(5*k+2)*ido]+cc[ic-1+(5*k+1)*ido];
      tr4=cc[i-1+(5*k+4)*ido]-cc[ic-1+(5*k+3)*ido];
      tr3=cc[i-1+(5*k+4)*ido]+cc[ic-1+(5*k+3)*ido];
      ch[i-1+k*ido]=cc[i-1+5*k*ido]+tr2+tr3;
      ch[i+k*ido]=cc[i+5*k*ido]+ti2+ti3;
      cr2=cc[i-1+5*k*ido]+tr11*tr2+tr12*tr3;

      ci2=cc[i+5*k*ido]+tr11*ti2+tr12*ti3;
      cr3=cc[i-1+5*k*ido]+tr12*tr2+tr11*tr3;

      ci3=cc[i+5*k*ido]+tr12*ti2+tr11*ti3;
      cr5=ti11*tr5+ti12*tr4;
      ci5=ti11*ti5+ti12*ti4;
      cr4=ti12*tr5-ti11*tr4;
      ci4=ti12*ti5-ti11*ti4;
      dr3=cr3-ci4;
      dr4=cr3+ci4;
      di3=ci3+cr4;
      di4=ci3-cr4;
      dr5=cr2+ci5;
      dr2=cr2-ci5;
      di5=ci2-cr5;
      di2=ci2+cr5;
      ch[i-1+(k+l1)*ido]=wa1[i-2]*dr2-wa1[i-1]*di2;
      ch[i+(k+l1)*ido]=wa1[i-2]*di2+wa1[i-1]*dr2;
      ch[i-1+(k+2*l1)*ido]=wa2[i-2]*dr3-wa2[i-1]*di3;
      ch[i+(k+2*l1)*ido]=wa2[i-2]*di3+wa2[i-1]*dr3;
      ch[i-1+(k+3*l1)*ido]=wa3[i-2]*dr4-wa3[i-1]*di4;
      ch[i+(k+3*l1)*ido]=wa3[i-2]*di4+wa3[i-1]*dr4;
      ch[i-1+(k+4*l1)*ido]=wa4[i-2]*dr5-wa4[i-1]*di5;
      ch[i+(k+4*l1)*ido]=wa4[i-2]*di5+wa4[i-1]*dr5;
      }
    }
  }

void radfg(int ido, int ip, int l1, int idl1,
  double *cc, double *ch, const double *wa)
  {
  static const double twopi=6.28318530717958647692;
  int idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is;
  double ai1, ai2, ar1, ar2, arg;
  double *csarr;
  int aidx;

  ipph=(ip+1)/ 2;
  if(ido !=1)
    {
    for(ik=0; ik<idl1; ik++)
      ch[ik]=cc[ik];
    for(j=1; j<ip; j++)
      for(k=0; k<l1; k++)
        ch[(k+j*l1)*ido]=cc[(k+j*l1)*ido];

    is=-ido;
    for(j=1; j<ip; j++)
      {
      is+=ido;
      for(k=0; k<l1; k++)
        {
        idij=is-1;
        for(i=2; i<ido; i+=2)
          {
          idij+=2;
          ch[i-1+(k+j*l1)*ido]=
            wa[idij-1]*cc[i-1+(k+j*l1)*ido]+wa[idij]*cc[i+(k+j*l1)*ido];
          ch[i+(k+j*l1)*ido]=
            wa[idij-1]*cc[i+(k+j*l1)*ido]-wa[idij]*cc[i-1+(k+j*l1)*ido];
          }
        }
      }

    for(j=1; j<ipph; j++)
      {
      jc=ip-j;
      for(k=0; k<l1; k++)
        {
        for(i=2; i<ido; i+=2)
          {
          cc[i-1+(k+j*l1)*ido]=ch[i-1+(k+j*l1)*ido]+ch[i-1+(k+jc*l1)*ido];
          cc[i-1+(k+jc*l1)*ido]=ch[i+(k+j*l1)*ido]-ch[i+(k+jc*l1)*ido];
          cc[i+(k+j*l1)*ido]=ch[i+(k+j*l1)*ido]+ch[i+(k+jc*l1)*ido];
          cc[i+(k+jc*l1)*ido]=ch[i-1+(k+jc*l1)*ido]-ch[i-1+(k+j*l1)*ido];
          }
        }
      }
    }
  else
    {                           /*now ido==1*/
    for(ik=0; ik<idl1; ik++)
      cc[ik]=ch[ik];
    }
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    for(k=0; k<l1; k++)
      {
      cc[(k+j*l1)*ido]=ch[(k+j*l1)*ido]+ch[(k+jc*l1)*ido];
      cc[(k+jc*l1)*ido]=ch[(k+jc*l1)*ido]-ch[(k+j*l1)*ido];
      }
    }

  csarr=(double *)malloc(2*ip*sizeof(double));
  arg=twopi / ip;
  csarr[0]=1.;
  csarr[1]=0.;
  csarr[2]=csarr[2*ip-2]=cos(arg);
  csarr[3]=sin(arg); csarr[2*ip-1]=-csarr[3];
  for (i=2; i<=ip/2; ++i)
    {
    csarr[2*i]=csarr[2*ip-2*i]=cos(i*arg);
    csarr[2*i+1]=sin(i*arg);
    csarr[2*ip-2*i+1]=-csarr[2*i+1];
    }
  for(l=1; l<ipph; l++)
    {
    lc=ip-l;
    ar1=csarr[2*l];
    ai1=csarr[2*l+1];
    for(ik=0; ik<idl1; ik++)
      {
      ch[ik+l*idl1]=cc[ik]+ar1*cc[ik+idl1];
      ch[ik+lc*idl1]=ai1*cc[ik+(ip-1)*idl1];
      }
    aidx=2*l;
    for(j=2; j<ipph; j++)
      {
      jc=ip-j;
      aidx+=2*l;
      if (aidx>=2*ip) aidx-=2*ip;
      ar2=csarr[aidx];
      ai2=csarr[aidx+1];
      for(ik=0; ik<idl1; ik++)
        {
        ch[ik+l*idl1]+=ar2*cc[ik+j*idl1];
        ch[ik+lc*idl1]+=ai2*cc[ik+jc*idl1];
        }
      }
    }
  free(csarr);

  for(j=1; j<ipph; j++)
    for(ik=0; ik<idl1; ik++)
      ch[ik]+=cc[ik+j*idl1];

  for(k=0; k<l1; k++)
    for(i=0; i<ido; i++)
      cc[i+k*ip*ido]=ch[i+k*ido];
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    j2=2*j;
    for(k=0; k<l1; k++)
      {
      cc[ido-1+(j2-1+k*ip)*ido] = ch[(k+j*l1)*ido];
      cc[(j2+k*ip)*ido] = ch[(k+jc*l1)*ido];
      }
    }
  if(ido==1)
    return;

  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    j2=2*j;
    for(k=0; k<l1; k++)
      {
      for(i=2; i<ido; i+=2)
        {
        ic=ido-i;
        cc[i-1+(j2+k*ip)*ido]=ch[i-1+(k+j*l1)*ido]+ch[i-1+(k+jc*l1)*ido];
        cc[ic-1+(j2-1+k*ip)*ido]=ch[i-1+(k+j*l1)*ido]-ch[i-1+(k+jc*l1)*ido];
        cc[i+(j2+k*ip)*ido]=ch[i+(k+j*l1)*ido]+ch[i+(k+jc*l1)*ido];
        cc[ic+(j2-1+k*ip)*ido]=ch[i+(k+jc*l1)*ido]-ch[i+(k+j*l1)*ido];
        }
      }
    }
  }

void radbg(int ido, int ip, int l1, int idl1,
  double *cc, double *ch, const double *wa)
  {
  static const double twopi=6.28318530717958647692;
  int     idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is;
  double ai1, ai2, ar1, ar2, arg;
  double *csarr;
  int aidx;

  ipph=(ip+1)/ 2;
  for(k=0; k<l1; k++)
    for(i=0; i<ido; i++)
      ch[i+k*ido]=cc[i+k*ip*ido];
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    j2=2*j;
    for(k=0; k<l1; k++)
      {
      ch[(k+j*l1)*ido]=cc[ido-1+(j2-1+k*ip)*ido]+cc[ido-1+(j2-1+k*ip)*ido];
      ch[(k+jc*l1)*ido]=cc[(j2+k*ip)*ido]+cc[(j2+k*ip)*ido];
      }
    }

  if(ido !=1)
    {
    for(j=1; j<ipph; j++)
      {
      jc=ip-j;
      for(k=0; k<l1; k++)
        {
        for(i=2; i<ido; i+=2)
          {
          ic=ido-i;
          ch[i-1+(k+j*l1)*ido] =
            cc[i-1+(2*j+k*ip)*ido]+cc[ic-1+(2*j-1+k*ip)*ido];
          ch[i-1+(k+jc*l1)*ido] =
            cc[i-1+(2*j+k*ip)*ido]-cc[ic-1+(2*j-1+k*ip)*ido];
          ch[i+(k+j*l1)*ido]=cc[i+(2*j+k*ip)*ido]-cc[ic+(2*j-1+k*ip)*ido];
          ch[i+(k+jc*l1)*ido]=cc[i+(2*j+k*ip)*ido]+cc[ic+(2*j-1+k*ip)*ido];
          }
        }
      }
    }

  csarr=(double *)malloc(2*ip*sizeof(double));
  arg=twopi / ip;
  csarr[0]=1.;
  csarr[1]=0.;
  csarr[2]=csarr[2*ip-2]=cos(arg);
  csarr[3]=sin(arg); csarr[2*ip-1]=-csarr[3];
  for (i=2; i<=ip/2; ++i)
    {
    csarr[2*i]=csarr[2*ip-2*i]=cos(i*arg);
    csarr[2*i+1]=sin(i*arg);
    csarr[2*ip-2*i+1]=-csarr[2*i+1];
    }
  for(l=1; l<ipph; l++)
    {
    lc=ip-l;
    ar1=csarr[2*l];
    ai1=csarr[2*l+1];
    for(ik=0; ik<idl1; ik++)
      {
      cc[ik+l*idl1]=ch[ik]+ar1*ch[ik+idl1];
      cc[ik+lc*idl1]=ai1*ch[ik+(ip-1)*idl1];
      }
    aidx=2*l;
    for(j=2; j<ipph; j++)
      {
      jc=ip-j;
      aidx+=2*l;
      if (aidx>=2*ip) aidx-=2*ip;
      ar2=csarr[aidx];
      ai2=csarr[aidx+1];
      for(ik=0; ik<idl1; ik++)
        {
        cc[ik+l*idl1]+=ar2*ch[ik+j*idl1];
        cc[ik+lc*idl1]+=ai2*ch[ik+jc*idl1];
        }
      }
    }
  free(csarr);

  for(j=1; j<ipph; j++)
    for(ik=0; ik<idl1; ik++)
      ch[ik]+=ch[ik+j*idl1];

  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    for(k=0; k<l1; k++)
      {
      ch[(k+j*l1)*ido]=cc[(k+j*l1)*ido]-cc[(k+jc*l1)*ido];
      ch[(k+jc*l1)*ido]=cc[(k+j*l1)*ido]+cc[(k+jc*l1)*ido];
      }
    }

  if(ido==1)
    return;
  for(j=1; j<ipph; j++)
    {
    jc=ip-j;
    for(k=0; k<l1; k++)
      {
      for(i=2; i<ido; i+=2)
        {
        ch[i-1+(k+j*l1)*ido]=cc[i-1+(k+j*l1)*ido]-cc[i+(k+jc*l1)*ido];
        ch[i-1+(k+jc*l1)*ido]=cc[i-1+(k+j*l1)*ido]+cc[i+(k+jc*l1)*ido];
        ch[i+(k+j*l1)*ido]=cc[i+(k+j*l1)*ido]+cc[i-1+(k+jc*l1)*ido];
        ch[i+(k+jc*l1)*ido]=cc[i+(k+j*l1)*ido]-cc[i-1+(k+jc*l1)*ido];
        }
      }
    }
  for(ik=0; ik<idl1; ik++)
    cc[ik]=ch[ik];
  for(j=1; j<ip; j++)
    for(k=0; k<l1; k++)
      cc[(k+j*l1)*ido]=ch[(k+j*l1)*ido];

  is=-ido;
  for(j=1; j<ip; j++)
    {
    is+=ido;
    for(k=0; k<l1; k++)
      {
      idij=is-1;
      for(i=2; i<ido; i+=2)
        {
        idij+=2;
        cc[i-1+(k+j*l1)*ido]=
          wa[idij-1]*ch[i-1+(k+j*l1)*ido]-wa[idij]*ch[i+(k+j*l1)*ido];
        cc[i+(k+j*l1)*ido]=
          wa[idij-1]*ch[i+(k+j*l1)*ido]+wa[idij]*ch[i-1+(k+j*l1)*ido];
        }
      }
    }
  }


/*----------------------------------------------------------------------
   cfftf1, cfftb1, cfftf, cfftb, cffti1, cffti. Complex FFTs.
  ----------------------------------------------------------------------*/

void cfftf1(int n, double c[], double ch[], const double wa[],
  const int ifac[])
  {
  int idot, k1, l1, l2, na, nf, ip, iw, nac, ido, idl1;
  double *p1, *p2;

  nf=ifac[1];
  na=0;
  l1=1;
  iw=0;
  for(k1=2; k1<=nf+1; k1++)
    {
    ip=ifac[k1];
    l2=ip*l1;
    ido=n / l2;
    idot=ido+ido;
    idl1=idot*l1;
    p1 = (na==0) ? c : ch;
    p2 = (na==0) ? ch : c;
    if(ip==4)
      passf4(idot, l1, p1, p2, wa+iw, wa+iw+idot, wa+iw+2*idot);
    else if(ip==2)
      passf2(idot, l1, p1, p2, wa+iw);
    else if(ip==3)
      passf3(idot, l1, p1, p2, wa+iw, wa+iw+idot);
    else if(ip==5)
      passf5(idot, l1, p1, p2, wa+iw, wa+iw+idot, wa+iw+2*idot, wa+iw+3*idot);
    else
      {
      passfg(&nac, idot, ip, l1, idl1, p1, p2, &wa[iw]);
      if(nac==0)
        na=1-na;
      }
    na=1-na;
    l1=l2;
    iw+=(ip-1)*idot;
    }
  if(na!=0)
    memcpy (c,ch,2*n*sizeof(double));
  }

void cfftb1(int n, double c[], double ch[], const double wa[],
  const int ifac[])
  {
  int idot, k1, l1, l2, na, nf, ip, iw, nac, ido, idl1;
  double *p1, *p2;

  nf=ifac[1];
  na=0;
  l1=1;
  iw=0;
  for(k1=2; k1<=nf+1; k1++)
    {
    ip=ifac[k1];
    l2=ip*l1;
    ido=n / l2;
    idot=ido+ido;
    idl1=idot*l1;
    p1 = (na==0) ? c : ch;
    p2 = (na==0) ? ch : c;
    if(ip==4)
      passb4(idot, l1, p1, p2, wa+iw, wa+iw+idot, wa+iw+2*idot);
    else if(ip==2)
      passb2(idot, l1, p1, p2, wa+iw);
    else if(ip==3)
      passb3(idot, l1, p1, p2, wa+iw, wa+iw+idot);
    else if(ip==5)
      passb5(idot, l1, p1, p2, wa+iw, wa+iw+idot, wa+iw+2*idot, wa+iw+3*idot);
    else
      {
      passbg(&nac, idot, ip, l1, idl1, p1, p2, &wa[iw]);
      if(nac==0)
        na=1-na;
      }
    na=1-na;
    l1=l2;
    iw+=(ip-1)*idot;
    }
  if(na!=0)
    memcpy (c,ch,2*n*sizeof(double));
  }

void cfftf(int n, double c[], double wsave[])
  {
  if(n!=1)
    cfftf1(n, c, wsave, wsave+2*n,(int*)(wsave+4*n));
  }

void cfftb(int n, double c[], double wsave[])
  {
  if(n!=1)
    cfftb1(n, c, wsave, wsave+2*n,(int*)(wsave+4*n));
  }

void cffti1(int n, double wa[], int ifac[])
  {
  static const int ntryh[4]= {3, 4, 2, 5};
  static const double twopi=6.28318530717958647692;
  double argh, argld, arg, fi;
  int idot, ntry=0, i, j, i1, k1, l1, l2, ib;
  int ld, ii, nf, ip, nl, nq, nr, ido, ipm;

  nl=n;
  nf=0;
  j=0;
startloop:
  j++;
  if(j<=4)
    ntry=ntryh[j-1];
  else
    ntry+=2;
  do
    {
    nq=nl / ntry;
    nr=nl-ntry*nq;
    if(nr !=0)
      goto startloop;
    nf++;
    ifac[nf+1]=ntry;
    nl=nq;
    if(ntry==2 && nf !=1)
      {
      for(i=2; i<=nf; i++)
        {
        ib=nf-i+2;
        ifac[ib+1]=ifac[ib];
        }
      ifac[2]=2;
      }
    }
  while(nl !=1);
  ifac[0]=n;
  ifac[1]=nf;
  argh=twopi /(double)n;
  i=1;
  l1=1;
  for(k1=1; k1<=nf; k1++)
    {
    ip=ifac[k1+1];
    ld=0;
    l2=l1*ip;
    ido=n / l2;
    idot=ido+ido+2;
    ipm=ip-1;
    for(j=1; j<=ipm; j++)
      {
      i1=i;
      wa[i-1]=1;
      wa[i]=0;
      ld+=l1;
      fi=0;
      argld=ld*argh;
      for(ii=4; ii<=idot; ii+=2)
        {
        i+=2;
        fi+=1;
        arg=fi*argld;
        wa[i-1]=cos(arg);
        wa[i]=sin(arg);
        }
      if(ip>5)
        {
        wa[i1-1]=wa[i-1];
        wa[i1]=wa[i];
        }
      }
    l1=l2;
    }
  }

void cffti(int n, double wsave[])
  {
  if (n!=1)
    cffti1(n, wsave+2*n,(int*)(wsave+4*n));
  }


/*----------------------------------------------------------------------
   rfftf1, rfftb1, rfftf, rfftb, rffti1, rffti. Real FFTs.
  ----------------------------------------------------------------------*/

void rfftf1(int n, double c[], double ch[], const double wa[],
  const int ifac[])
  {
  int i, k1, l1, l2, na, kh, nf, ip, iw, ido, idl1;
  double *p1, *p2;

  nf=ifac[1];
  na=1;
  l2=n;
  iw=n-1;
  for(k1=1; k1<=nf;++k1)
    {
    kh=nf-k1;
    ip=ifac[kh+2];
    l1=l2 / ip;
    ido=n / l2;
    idl1=ido*l1;
    iw-=(ip-1)*ido;
    na=1-na;
    p1 = (na==0) ? c : ch;
    p2 = (na==0) ? ch : c;
    if(ip==4)
      radf4(ido, l1, p1, p2, wa+iw, wa+iw+ido, wa+iw+2*ido);
    else if(ip==2)
      radf2(ido, l1, p1, p2, wa+iw);
    else if(ip==3)
      radf3(ido, l1, p1, p2, wa+iw, wa+iw+ido);
    else if(ip==5)
      radf5(ido, l1, p1, p2, wa+iw, wa+iw+ido, wa+iw+2*ido, wa+iw+3*ido);
    else
      {
      if(ido==1)
        na=1-na;
      if(na==0)
        radfg(ido, ip, l1, idl1, c, ch, wa+iw);
      else
        radfg(ido, ip, l1, idl1, ch, c, wa+iw);
      na=1-na;
      }
    l2=l1;
    }
  if(na==1)
    return;
  for(i=0; i<n; i++)
    c[i]=ch[i];
}

void rfftb1(int n, double c[], double ch[], const double wa[],
  const int ifac[])
  {
  int k1, l1, l2, na, nf, ip, iw, ido, idl1;
  double *p1, *p2;

  nf=ifac[1];
  na=0;
  l1=1;
  iw=0;
  for(k1=1; k1<=nf; k1++)
    {
    ip=ifac[k1+1];
    l2=ip*l1;
    ido=n / l2;
    idl1=ido*l1;
    p1 = (na==0) ? c : ch;
    p2 = (na==0) ? ch : c;
    if(ip==4)
      radb4(ido, l1, p1, p2, wa+iw, wa+iw+ido, wa+iw+2*ido);
    else if(ip==2)
      radb2(ido, l1, p1, p2, wa+iw);
    else if(ip==3)
      radb3(ido, l1, p1, p2, wa+iw, wa+iw+ido);
    else if(ip==5)
      radb5(ido, l1, p1, p2, wa+iw, wa+iw+ido, wa+iw+2*ido, wa+iw+3*ido);
    else
      {
      radbg(ido, ip, l1, idl1, p1, p2, wa+iw);
      if(ido!=1)
        na=1-na;
      }
    na=1-na;
    l1=l2;
    iw+=(ip-1)*ido;
    }
  if(na!=0)
    memcpy (c,ch,n*sizeof(double));
  }

void rfftf(int n, double r[], double wsave[])
  {
  if(n!=1)
    rfftf1(n, r, wsave, wsave+n,(int*)(wsave+2*n));
  }

void rfftb(int n, double r[], double wsave[])
  {
  if(n!=1)
    rfftb1(n, r, wsave, wsave+n,(int*)(wsave+2*n));
  }

void rffti1(int n, double wa[], int ifac[])
  {
  static const int ntryh[4]={4, 2, 3, 5};
  static const double twopi=6.28318530717958647692;
  double argh, argld, arg, fi;
  int ntry=0, i, j, k1, l1, l2, ib, ld, ii, nf, ip, nl, is, nq, nr;
  int ido, ipm, nfm1;

  nl=n;
  nf=0;
  j=0;
startloop:
 ++j;
  if(j<=4)
    ntry=ntryh[j-1];
  else
    ntry+=2;
  do
    {
    nq=nl / ntry;
    nr=nl-ntry*nq;
    if(nr !=0)
      goto startloop;
    ++nf;
    ifac[nf+1]=ntry;
    nl=nq;
    if(ntry==2 && nf !=1)
      {
      for(i=2; i<=nf; i++)
        {
        ib=nf-i+2;
        ifac[ib+1]=ifac[ib];
        }
      ifac[2]=2;
      }
    }
  while(nl !=1);
  ifac[0]=n;
  ifac[1]=nf;
  argh=twopi /(double)(n);
  is=0;
  nfm1=nf-1;
  l1=1;
  if(nfm1==0)
    return;
  for(k1=1; k1<=nfm1; k1++)
    {
    ip=ifac[k1+1];
    ld=0;
    l2=l1*ip;
    ido=n / l2;
    ipm=ip-1;
    for(j=1; j<=ipm;++j)
      {
      ld+=l1;
      i=is;
      argld=(double)ld*argh;

      fi=0;
      for(ii=3; ii<=ido; ii+=2)
        {
        i+=2;
        fi+=1;
        arg=fi*argld;
        wa[i-2]=cos(arg);
        wa[i-1]=sin(arg);
        }
      is+=ido;
      }
    l1=l2;
    }
  }

void rffti(int n, double wsave[])
  {
  if (n!=1)
    rffti1(n, wsave+n,(int*)(wsave+2*n));
  }
}
// fftpack.c end

// from bluestein.c of HEALPix C++
namespace {
/* returns the largest prime factor of n */
/*
  int largest_prime_factor (int n)
  {
  int maxdiv=0,x,limit,tmp;
  while (((tmp=(n>>1))<<1)==n)
    { maxdiv=2; n=tmp; }

  limit=(int)sqrt(n+0.01);
  for (x=3; x<=limit; x+=2)
  while ((tmp=(n/x))*x==n)
    {
    maxdiv=x;
    n=tmp;
    limit=(int)sqrt(n+0.01);
    }
  if (n>maxdiv) maxdiv=n;

  return maxdiv;
  }
*/
  
/* returns the sum of all prime factors of n */
int prime_factor_sum (int n)
  {
  int result=0,x,limit,tmp;
  while (((tmp=(n>>1))<<1)==n)
    { result+=2; n=tmp; }

  limit=(int)sqrt(n+0.01);
  for (x=3; x<=limit; x+=2)
  while ((tmp=(n/x))*x==n)
    {
    result+=x;
    n=tmp;
    limit=(int)sqrt(n+0.01);
    }
  if (n>1) result+=n;

  return result;
  }

/* returns the smallest power of 2 which is >= n */
/*
int nextpow2 (int n)
  {
  int cnt=1;
  while ((n>>cnt)>0)
    ++cnt;
  return 1<<cnt;
  }
*/
  
/* returns the smallest composite of 2, 3 and 5 which is >= n */
int good_size(int n)
  {
  int maxfactors=1, i, j, k, f2=1, f3, f5, bestfac, guessfac;
  while ((n>>maxfactors)>0)
    ++maxfactors;
  bestfac=1<<maxfactors;

  for (i=0; i<maxfactors; ++i)
    {
    f3=1;
    for (j=0; j<maxfactors-i; ++j)
      {
      if (f2*f3>bestfac) break;
      f5=1;
      for (k=0; k<maxfactors-i-j; ++k)
        {
        guessfac = f2*f3*f5;
        if (guessfac>=bestfac) break;
        if ((guessfac>=n) && (guessfac<bestfac))
          bestfac=guessfac;
        f5*=5;
        }
      f3*=3;
      }
    f2*=2;
    }
  return bestfac;
  }

void bluestein_i (int n, double **tstorage)
  {
  static const double pi=3.14159265358979323846;
  int n2=good_size(n*2-1);
  int m, coeff;
  double angle, xn2;
  double *bk, *bkf, *work;
  double pibyn=pi/n;
  *tstorage = (double *)malloc (sizeof(double)*(1+2*n+8*n2+15));
  ((int *)(*tstorage))[0]=n2;
  bk  = *tstorage+1;
  bkf = *tstorage+1+2*n;
  work= *tstorage+1+2*(n+n2);

/* initialize b_k */
  bk[0] = 1;
  bk[1] = 0;

  coeff=0;
  for (m=1; m<n; ++m)
    {
    coeff+=2*m-1;
    if (coeff>=2*n) coeff-=2*n;
    angle = pibyn*coeff;
    bk[2*m] = cos(angle);
    bk[2*m+1] = sin(angle);
    }

/* initialize the zero-padded, Fourier transformed b_k. Add normalisation.  */
  xn2 = 1./n2;
  bkf[0] = bk[0]*xn2;
  bkf[1] = bk[1]*xn2;
  for (m=2; m<2*n; m+=2)
    {
    bkf[m]   = bkf[2*n2-m]   = bk[m]   *xn2;
    bkf[m+1] = bkf[2*n2-m+1] = bk[m+1] *xn2;
    }
  for (m=2*n;m<=(2*n2-2*n+1);++m)
    bkf[m]=0.;
  cffti (n2,work);
  cfftf (n2,bkf,work);
  }

void bluestein (int n, double *data, double *tstorage, int isign)
  {
  int n2=*((int *)tstorage);
  int m;
  double *bk, *bkf, *akf, *work;
  bk  = tstorage+1;
  bkf = tstorage+1+2*n;
  work= tstorage+1+2*(n+n2);
  akf = tstorage+1+2*n+6*n2+15;  

/* initialize a_k and FFT it */
  if (isign>0)
    for (m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   - data[m+1]*bk[m+1];
      akf[m+1] = data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  else
    for (m=0; m<2*n; m+=2)
      {
      akf[m]   = data[m]*bk[m]   + data[m+1]*bk[m+1];
      akf[m+1] =-data[m]*bk[m+1] + data[m+1]*bk[m];
      }
  for (m=2*n; m<2*n2; ++m)
    akf[m]=0;

  cfftf (n2,akf,work);

/* do the convolution */
  if (isign>0)
    for (m=0; m<2*n2; m+=2)
      {
      double im = -akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  =  akf[m]*bkf[m]   + akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }
  else
    for (m=0; m<2*n2; m+=2)
      {
      double im = akf[m]*bkf[m+1] + akf[m+1]*bkf[m];
      akf[m  ]  = akf[m]*bkf[m]   - akf[m+1]*bkf[m+1];
      akf[m+1]  = im;
      }


/* inverse FFT */
  cfftb (n2,akf,work);

/* multiply by b_k* */
  if (isign>0)
    for (m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] - bk[m+1]*akf[m+1];
      data[m+1] = bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  else
    for (m=0; m<2*n; m+=2)
      {
      data[m]   = bk[m]  *akf[m] + bk[m+1]*akf[m+1];
      data[m+1] =-bk[m+1]*akf[m] + bk[m]  *akf[m+1];
      }
  }
}
// bluestein.c end

// ls_fft.c of HEALPix C++
namespace {
  using namespace THealFFT;

complex_plan make_complex_plan (int length)
  {
  complex_plan plan = (complex_plan) malloc(sizeof(complex_plan_i));
  int pfsum = prime_factor_sum(length);
  double comp1 = length*pfsum;
  double comp2 = 2*3*length*log(3.*length);
  plan->length=length;
  plan->bluestein = (comp2<comp1);
  if (plan->bluestein)
    bluestein_i (length,&(plan->work));
  else
    {
    plan->work=(double *)malloc((4*length+15)*sizeof(double));
    cffti(length, plan->work);
    }
  return plan;
  }

void kill_complex_plan (complex_plan plan)
  {
  free(plan->work);
  free(plan);
  }
}

void THealFFT::complex_plan_forward (complex_plan plan, double *data)
  {
  if (plan->bluestein)
    bluestein (plan->length, data, plan->work, -1);
  else
    cfftf (plan->length, data, plan->work);
  }

void THealFFT::complex_plan_backward (complex_plan plan, double *data)
  {
  if (plan->bluestein)
    bluestein (plan->length, data, plan->work, 1);
  else
    cfftb (plan->length, data, plan->work);
  }

namespace {

  using namespace THealFFT;

real_plan make_real_plan (int length)
  {
  real_plan plan = (real_plan) malloc(sizeof(real_plan_i));
  int pfsum = prime_factor_sum(length);
  double comp1 = .5*length*pfsum;
  double comp2 = 2*3*length*log(3.*length);
  plan->length=length;
  plan->bluestein = (comp2<comp1);
  if (plan->bluestein)
    bluestein_i (length,&(plan->work));
  else
    {
    plan->work=(double *)malloc((2*length+15)*sizeof(double));
    rffti(length, plan->work);
    }
  return plan;
  }

void kill_real_plan (real_plan plan)
  {
  free(plan->work);
  free(plan);
  }
}

void THealFFT::real_plan_forward_fftpack (real_plan plan, double *data)
  {
  if (plan->bluestein)
    {
    int m;
    int n=plan->length;
    double *tmp = (double *)malloc(2*n*sizeof(double));
    for (m=0; m<n; ++m)
      {
      tmp[2*m] = data[m];
      tmp[2*m+1] = 0.;
      }
    bluestein(n,tmp,plan->work,-1);
    data[0] = tmp[0];
    memcpy (data+1, tmp+2, (n-1)*sizeof(double));
    free (tmp);
    }
  else
    rfftf (plan->length, data, plan->work);
  }

/*
namespace {
void fftpack2halfcomplex (double *data, int n)
  {
  int m;
  double *tmp = (double *)malloc(n*sizeof(double));
  tmp[0]=data[0];
  for (m=1; m<(n+1)/2; ++m)
    {
    tmp[m]=data[2*m-1];
    tmp[n-m]=data[2*m];
    }
  if (!(n&1))
    tmp[n/2]=data[n-1];
  memcpy (data,tmp,n*sizeof(double));
  free(tmp);
  }

void halfcomplex2fftpack (double *data, int n)
  {
  int m;
  double *tmp = (double *)malloc(n*sizeof(double));
  tmp[0]=data[0];
  for (m=1; m<(n+1)/2; ++m)
    {
    tmp[2*m-1]=data[m];
    tmp[2*m]=data[n-m];
    }
  if (!(n&1))
    tmp[n-1]=data[n/2];
  memcpy (data,tmp,n*sizeof(double));
  free(tmp);
  }

void real_plan_forward_fftw (real_plan plan, double *data)
  {
  real_plan_forward_fftpack (plan, data);
  fftpack2halfcomplex (data,plan->length);
  }
}
*/

void THealFFT::real_plan_backward_fftpack (real_plan plan, double *data)
  {
  if (plan->bluestein)
    {
    int m;
    int n=plan->length;
    double *tmp = (double *)malloc(2*n*sizeof(double));
    tmp[0]=data[0];
    tmp[1]=0.;
    memcpy (tmp+2,data+1, (n-1)*sizeof(double));
    if ((n&1)==0) tmp[n+1]=0.;
    for (m=2; m<n; m+=2)
      {
      tmp[2*n-m]=tmp[m];
      tmp[2*n-m+1]=-tmp[m+1];
      }
    bluestein (n, tmp, plan->work, 1);
    for (m=0; m<n; ++m)
      data[m] = tmp[2*m];
    free (tmp);
    }
  else
    rfftb (plan->length, data, plan->work);
  }

/*
namespace {
  void real_plan_backward_fftw (real_plan plan, double *data)
  {
  halfcomplex2fftpack (data,plan->length);
  real_plan_backward_fftpack (plan, data);
  }
}
*/

void THealFFT::real_plan_forward_c (real_plan plan, double *data)
  {
  int m;
  int n=plan->length;

  if (plan->bluestein)
    {
    for (m=1; m<2*n; m+=2)
      data[m]=0;
    bluestein (plan->length, data, plan->work, -1);
    data[1]=0;
    for (m=2; m<n; m+=2)
      {
      double avg;
      avg = 0.5*(data[2*n-m]+data[m]);
      data[2*n-m] = data[m] = avg;
      avg = 0.5*(data[2*n-m+1]-data[m+1]);
      data[2*n-m+1] = avg;
      data[m+1] = -avg;
      }
    if ((n&1)==0) data[n+1] = 0.;
    }
  else
    {
    for (m=0; m<n; ++m) data[m+1] = data[2*m];
    rfftf (n, data+1, plan->work);
    data[0] = data[1];
    data[1] = 0;
    for (m=2; m<n; m+=2)
      {
      data[2*n-m]   =  data[m];
      data[2*n-m+1] = -data[m+1];
      }
    if ((n&1)==0) data[n+1] = 0.;
    }
  }

void THealFFT::real_plan_backward_c (real_plan plan, double *data)
  {
  int m;
  int n=plan->length;

  if (plan->bluestein)
    {
    data[1]=0;
    for (m=2; m<n; m+=2)
      {
      double avg;
      avg = 0.5*(data[2*n-m]+data[m]);
      data[2*n-m] = data[m] = avg;
      avg = 0.5*(data[2*n-m+1]-data[m+1]);
      data[2*n-m+1] = avg;
      data[m+1] = -avg;
      }
    if ((n&1)==0) data[n+1] = 0.;
    bluestein (plan->length, data, plan->work, 1);
    for (m=1; m<2*n; m+=2)
      data[m]=0;
    }
  else
    {
    data[1] = data[0];
    rfftb (n, data+1, plan->work);
    for (m=n-1; m>=0; --m)
      {
      data[2*m]   = data[m+1];
      data[2*m+1] = 0.;
      }
    }
  }
// ls_fft.c end

namespace THealFFT {
  cfft::cfft(int size_)
    : n(size_), plan(make_complex_plan(size_)) {}

  cfft::~cfft ()
  { if (plan!=0) kill_complex_plan (plan); }

  void cfft::Set (int size_)
  {
    if (plan!=0) kill_complex_plan (plan);
    n=size_;
    plan=make_complex_plan(size_);
  }

  rfft::rfft (int size_)
    : n(size_), plan(make_real_plan(size_)) {}

  rfft::~rfft ()
  { if (plan!=0) kill_real_plan (plan); }

  void rfft::Set (int size_)
  {
    if (plan!=0) kill_real_plan (plan);
    n=size_;
    plan=make_real_plan(size_);
  }

} // THealFFT
