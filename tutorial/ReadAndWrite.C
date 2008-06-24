// $Id: ReadAndWrite.C,v 1.1 2008/06/24 08:51:15 oxon Exp $
// Author: Akira Okumura 2008/06/24

#include <vector>
#include <iostream>

#include "TSystem.h"

#ifndef __CINT__
#include "THealPix.h"
#include "THealUtil.h"
#endif // __CINT__

void ReadAndWrite()
{
  gSystem->Load("libRootHealPix");

  THealPixD* hpd1 = new THealPixD("hpd1", "data1", 5);
  THealPixD* hpd2 = new THealPixD("hpd2", "data2", 4);
  hpd1->SetUnit("MeV"); // Set the unit of entries
  hpd2->SetUnit("GeV");
  hpd1->SetDegree(); // change the unit of Fill(theta, phi)
  hpd2->SetDegree();

  for(Int_t i = 0; i<1800; i++){
    for(Int_t j = 0; j<3600; j++){
      hpd1->Fill(i/10., j/10., 1);
      hpd2->Fill(i/10., j/10., 5);
    } // j
  } // i

  std::vector<THealPix*> vec;
  vec.push_back(hpd1);
  vec.push_back(hpd2);

  THealUtil::SaveToFits("!test.fits", vec); // hpd1/2 are saved in one file

  delete hpd1;
  hpd1 = 0;
  delete hpd2;
  hpd2 = 0;

  // Read only "data2" from FITS file
  THealPixD* hpd = THealPixD::ReadFits("test.fits", "data2");

  std::cout << "Order: " << hpd->GetOrder() << std::endl;
  std::cout << "Nside: " << hpd->GetNside() << std::endl;
  std::cout << "Npix : " << hpd->GetNpix()  << std::endl;
  std::cout << "Unit : " << hpd->GetUnit()  << std::endl;
  std::cout << "Title: " << hpd->GetTitle() << std::endl;

}
