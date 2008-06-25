// $Id: Rebin.C,v 1.1 2008/06/25 04:31:30 oxon Exp $
// Author: Akira Okumura 2008/06/24

#include <vector>
#include <iostream>

#include "TSystem.h"
#include "TH1D.h"
#include "TCanvas.h"

#ifndef __CINT__
#include "THealPix.h"
#include "THealUtil.h"
#endif // __CINT__

void Rebin()
{
  gSystem->Load("libRootHealPix");

  THealPixD* hpd1 = new THealPixD("hpd1", "unknown", 5); // RING
  hpd1->SetDegree(); // change the unit of Fill(theta, phi)
 
  THealPixD* hpd2 = new THealPixD("hpd2", "", 5, kTRUE); // NESTED
  hpd2->SetDegree(); // change the unit of Fill(theta, phi)

  for(Int_t i = 0; i<1800; i++){
    for(Int_t j = 0; j<3600; j++){
      hpd1->Fill(i/10., j/10., 1);
      hpd2->Fill(i/10., j/10., 1);
    } // j
  } // i

  TH1D* h1[3] = {new TH1D("h11", "", 50, 0, 1000), new TH1D("h12", "", 50, 0, 250),
		 new TH1D("h13", "", 100, 0, 0)};
  TH1D* h2[3] = {new TH1D("h21", "", 50, 0, 1000), new TH1D("h22", "", 50, 0, 250),
		 new TH1D("h23", "", 100, 0, 0)};
  TH1D* h3[3] = {new TH1D("h31", "", 50, 0, 4000), new TH1D("h32", "", 50, 0, 1000),
		 new TH1D("h33", "", 100, 0, 0)};
  TH1D* h4[3] = {new TH1D("h41", "", 50, 0, 4000), new TH1D("h42", "", 50, 0, 1000),
		 new TH1D("h43", "", 100, 0, 0)};

  THealPixD* hpd3 = (THealPixD*)hpd1->Rebin(6, "hpd3"); // one step higher level
  THealPixD* hpd4 = (THealPixD*)hpd2->Rebin(6, "hpd4"); // one step higher level
  THealPixD* hpd5 = (THealPixD*)hpd1->Rebin(4, "hpd5"); // one step lower level
  THealPixD* hpd6 = (THealPixD*)hpd2->Rebin(4, "hpd6"); // one step lower level

  for(Int_t i = 0; i < hpd1->GetNpix(); i++){
    // Check RING 5->6 rebinning
    Int_t inest = hpd1->Ring2Nest(i);
    Double_t content = hpd1->GetBinContent(i);
    h1[0]->Fill(content);
    for(Int_t j = 0; j < 4; j++){
      Double_t rebinned = hpd3->GetBinContent(hpd3->Nest2Ring(inest*4 + j));
      h1[1]->Fill(rebinned);
      content -= rebinned;
    } // j
    h1[2]->Fill(content);

    // Check NESTED 5->6 rebinning
    content = hpd2->GetBinContent(i);
    h2[0]->Fill(content);
    for(Int_t j = 0; j < 4; j++){
      Double_t rebinned = hpd4->GetBinContent(i*4 + j);
      h2[1]->Fill(rebinned);
      content -= rebinned;
    } // j
    h2[2]->Fill(content);
  } // i

  for(Int_t i = 0; i < hpd5->GetNpix(); i++){
    // Check RING 5->4 rebinning
    Int_t inest = hpd5->Ring2Nest(i);
    Double_t content = hpd5->GetBinContent(i);
    h3[0]->Fill(content);
    for(Int_t j = 0; j < 4; j++){
      Double_t rebinned = hpd1->GetBinContent(hpd1->Nest2Ring(inest*4 + j));
      h3[1]->Fill(rebinned);
      content -= rebinned;
    } // j
    h3[2]->Fill(content);

    // Check NESTED 5->4 rebinning
    content = hpd6->GetBinContent(i);
    h4[0]->Fill(content);
    for(Int_t j = 0; j < 4; j++){
      Double_t rebinned = hpd2->GetBinContent(i*4 + j);
      h4[1]->Fill(rebinned);
      content -= rebinned;
    } // j
    h4[2]->Fill(content);
  } // i

  TCanvas* can = new TCanvas("can", "", 600, 800);
  can->Divide(3, 4);
  for(Int_t i = 0; i < 3; i++){
    can->cd(i + 1);
    h1[i]->Draw();
    can->cd(i + 4);
    h2[i]->Draw();
    can->cd(i + 7);
    h3[i]->Draw();
    can->cd(i + 10);
    h4[i]->Draw();
  } // i
}
