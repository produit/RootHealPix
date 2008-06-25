// $Id: CalcTest.C,v 1.1 2008/06/25 07:59:48 oxon Exp $
// Author: Akira Okumura 2008/06/25

#include <iostream>

#include "TSystem.h"

#ifndef __CINT__
#include "THealPix.h"
#include "THealUtil.h"
#endif // __CINT__

void CalcTest()
{
  gSystem->Load("libRootHealPix");

  THealPixF hpf1("hpf1", "", 1);
  THealPixF hpf2("hpf2", "", 1);
  THealPixD hpd1("hpd1", "", 1);
  THealPixD hpd2("hpd2", "", 1);

  const Int_t kNface = 12;
  for(Int_t i = 0; i < kNface; i++){
    hpf1.SetBinContent(i, i);
    hpf2.SetBinContent(i, kNface - i);
    hpd1.SetBinContent(i, i);
    hpd2.SetBinContent(i, kNface - i);
  } // i

  THealPixF hpf_bad1("hpf_bad1", "", 2);
  THealPixF hpf_bad2("hpf_bad2", "", 1, kTRUE);
  THealPixD hpd_bad1("hpd_bad1", "", 2);
  THealPixD hpd_bad2("hpd_bad2", "", 1, kTRUE);

  THealPixF hpf_res1 = 2*hpf1;
  THealPixF hpf_res2 = hpf1*2;
  THealPixF hpf_res3 = hpf1 + hpf2;
  THealPixF hpf_res4 = hpf1 - hpf2;
  THealPixF hpf_res5 = hpf1 * hpf2;
  THealPixF hpf_res6 = hpf1 / hpf2;
  
  for(Int_t i = 0; i < kNface; i++){
    std::cout << hpf_res1.GetBinContent(i) << "\t"
	      << hpf_res2.GetBinContent(i) << "\t"
	      << hpf_res3.GetBinContent(i) << "\t"
	      << hpf_res4.GetBinContent(i) << "\t"
	      << hpf_res5.GetBinContent(i) << "\t"
	      << hpf_res6.GetBinContent(i) << "\n";
  } // i

  THealPixF tmpf;
  tmpf = hpf1 + hpf_bad1;
  tmpf = hpf1 + hpf_bad2;
  tmpf = hpf1 - hpf_bad1;
  tmpf = hpf1 - hpf_bad2;
  tmpf = hpf1 * hpf_bad1;
  tmpf = hpf1 * hpf_bad2;
  tmpf = hpf1 / hpf_bad1;
  tmpf = hpf1 / hpf_bad2;

  THealPixD hpd_res1 = 2*hpd1;
  THealPixD hpd_res2 = hpd1*2;
  THealPixD hpd_res3 = hpd1 + hpd2;
  THealPixD hpd_res4 = hpd1 - hpd2;
  THealPixD hpd_res5 = hpd1 * hpd2;
  THealPixD hpd_res6 = hpd1 / hpd2;
  
  for(Int_t i = 0; i < kNface; i++){
    std::cout << hpd_res1.GetBinContent(i) << "\t"
	      << hpd_res2.GetBinContent(i) << "\t"
	      << hpd_res3.GetBinContent(i) << "\t"
	      << hpd_res4.GetBinContent(i) << "\t"
	      << hpd_res5.GetBinContent(i) << "\t"
	      << hpd_res6.GetBinContent(i) << "\n";
  } // i

  THealPixD tmpd;
  tmpd = hpd1 + hpd_bad1;
  tmpd = hpd1 + hpd_bad2;
  tmpd = hpd1 - hpd_bad1;
  tmpd = hpd1 - hpd_bad2;
  tmpd = hpd1 * hpd_bad1;
  tmpd = hpd1 * hpd_bad2;
  tmpd = hpd1 / hpd_bad1;
  tmpd = hpd1 / hpd_bad2;
}
