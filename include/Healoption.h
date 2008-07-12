// $Id: Healoption.h,v 1.1 2008/07/12 02:28:49 oxon Exp $
// Author: Akira Okumura 2008/07/10

/*****************************************************************************
   Copyright (C) 2008-, Akira Okumura
   All rights reserved.
******************************************************************************/

#ifndef T_HEAL_OPTION
#define T_HEAL_OPTION

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Healoption                                                           //
//                                                                      //
// HEALPix option structure.                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

typedef struct Healoption_t {
   // chopt may be the concatenation of the following options:

   int Axis;        // "A"  Axis are not drawn around the graph.
   int Off;         // "][" With H option, the first and last vertical lines are not drawn.
   int Keep;        // "K"  The status of the HEALPix is kept in memory
   int Same;        // "S"  Histogram is plotted in the current PAD.
   int Update;      // "U"  Update histogram previously plotted with option K
   int Color;       // "COL"   Draw 2D plot with Colored boxes.
   int Contour;     // "CONT"  Draw 2D plot as a Contour plot.
   int Func;        // "FUNC"  Draw only the function (for example in case of fit).
   int Heal;        // "HEAL"  Draw only the histogram.
   int Lego;        // "LEGO"  Draw as a Lego plot(LEGO,Lego=1, LEGO1,Lego1=11, LEGO2,Lego=12).
   int Scat;        // "SCAT"  Draw 2D plot a Scatter plot.
   int Surf;        // "SURF"  Draw as a Surface (SURF,Surf=1, SURF1,Surf=11, SURF2,Surf=12)
   int Tri;         // "TRI"   Draw 2D plot with Delaunay triangles.
   int System;      // type of coordinate system(1=car,2=pol,3=cyl,4=sph,5=psr)
   int Zscale;      // "Z"   to display the Z scale (color palette)
   int FrontBox;    //  = 0 to suppress the front box
   int BackBox;     //  = 0 to suppress the back box
   int List;        //  = 1 to generate the TObjArray "contours"
   int HighRes;     //  = 1 to select high resolution
   int Proj;        //  = 1 to get an Aitoff projection
                    //  = 2 to get a Lambert projection
   int AxisPos;     //  Axis position
   int Zero;        // if selected with any LEGO option the empty are not drawn.

   // the following structure members are set to 1 if the corresponding option
   // in the current style is selected.
   int Logz;        // log scale in Z. Also set by histogram option

} Healoption_t;

#endif // HEAL_OPTION
