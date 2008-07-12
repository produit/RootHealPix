// $Id: THealPixCube.h,v 1.2 2008/07/12 22:07:07 oxon Exp $
// Author: Akira Okumura 2008/07/11

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_PIX_CUBE
#define T_HEAL_PIX_CUBE

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THealPixCube                                                         //
//                                                                      //
// Cubic HEALPix Image class                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "THealPix.h"

class THealPixCube : public TNamed {
private:
  THealPixCube& operator=(const THealPixCube&); // Not implemented
  THealPixCube(const THealPixCube& cube); // Not implemented

protected:
  TAxis         fWaxis;   //W axis descriptor
  Double_t      fMaximum; //Maximum value for plotting
  Double_t      fMinimum; //Minimum value for plotting
  Int_t         fN;       //Number of HEALPix
  THealPix**    fHeals;   //[fN]Pointer to array of THealPix

protected:
  THealPixCube();
  THealPixCube(const char* name, const char* title, Int_t nbins,
	       Double_t wlow, Double_t wup);

public:
  virtual ~THealPixCube();
  
  virtual Int_t    GetN() const {return fN;}
          TAxis*   GetWaxis() const {return &((THealPixCube*)this)->fWaxis;}
  
  THealPix* operator[](Int_t n) const;

  ClassDef(THealPixCube, 1);
};

//______________________________________________________________________________
class THealPixCubeF : public THealPixCube {
private:
protected:
  THealPixCubeF();

public:
  THealPixCubeF(const char* name, const char* title, Int_t order, Int_t nbins,
		Double_t wlow, Double_t wup, Bool_t nested = kFALSE);
  virtual ~THealPixCubeF();

  THealPixCubeF* ReadFits(const char* fname, const char* colname);

  ClassDef(THealPixCubeF, 1);
};

//______________________________________________________________________________
class THealPixCubeD : public THealPixCube {
private:
protected:
  THealPixCubeD();

public:
  THealPixCubeD(const char* name, const char* title, Int_t order, Int_t nbins,
		Double_t wlow, Double_t wup, Bool_t nested = kFALSE);
  virtual ~THealPixCubeD();

  THealPixCubeD* ReadFits(const char* fname, const char* colname);

  ClassDef(THealPixCubeD, 1);
};

#endif // T_HEAL_PIX_CUBE
