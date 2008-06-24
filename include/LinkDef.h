// $Id: LinkDef.h,v 1.1 2008/06/24 08:16:43 oxon Exp $
// Author: Akira Okumura 2008/06/20

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class THealPix+;
#pragma link C++ class std::vector<THealPix*>+;
#pragma link C++ class THealPixD+;

#pragma link C++ namespace THealUtil;
#pragma link C++ function THealUtil::SaveToFits(const char*, const std::vector<THealPix*>&);

#endif
