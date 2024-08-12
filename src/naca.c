/* This file is part of LLT, a library for Lifting Line Theory
 * calculations
 *
 * Copyright (C) 2021 Michael Carley
 *
 * LLT is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  LLT is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LLT.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <amc.h>

#include <llt.h>

/*
 * Based on Ladson, Charles L., Brooks, Cuyler W., Jr, Hill, Acquilla
 * S., and Sproles, Darrell W., Computer Program to Obtain Ordinates
 * for NACA Airfoils, NASA Technical Memorandum 4741, 1996.
 */

/**
 * @{
 * @ingroup shapes
 *
 */

/** 
 * Evaluate NACA four digit aerofoil shape
 * 
 * Evaluate NACA four digit aerofoil shape using formulae of Ladson,
 * Charles L., Brooks, Cuyler W., Jr., Hill, Acquilla S., Sproles,
 * Darrell W., `Computer Program to Obtain Ordinates for NACA
 * Airfoils', NASA-TM-4741, 1996
 * 
 * https://ntrs.nasa.gov/citations/19970008124
 * 
 * Resulting section is NACApmxx, e.g. NACA2312 would have
 * \f$t=0.12\f$, \f$p=0.02\f$, \f$m=0.3\f$. 
 * 
 * @param t thickness to chord ratio;
 * @param p maximum camber as fraction of chord;
 * @param m location of maximum camber as fraction of chord;
 * @param x ordinate \f$-1\leq x\leq1\f$.
 * 
 * @return abcissa \f$y(|x|)\f$, lower surface for \f$x<0\f$, upper
 * surface otherwise.
 */

gdouble llt_naca_four(gdouble t, gdouble p, gdouble m, gdouble x)
  
{
  gdouble y, sgn = 1 ;
  gdouble a[] = {0.2969, -0.1260, -0.3516, 0.2843, -0.1015} ;

  g_assert(x >= -1.0 && x <= 1.0) ;
  
  if ( x < 0 ) { x = -x ; sgn = -1 ; }

  y = a[0]*sqrt(x) ;
  y += x*(a[1] + x*(a[2] + x*(a[3] + x*a[4]))) ;
  
  y *= sgn*t/0.2 ;

  if ( p == 0.0 || m == 0.0 ) return y ;
  /*camber distribution*/
  if ( x <= m ) {
    y += p/m/m*(x*(2*m - x)) ;
  } else {
    y += p/(1.0-m)/(1.0-m)*((1.0-2*m) + x*(2*m - x)) ;
  }
  
  return y ;
}

/**
 * @}
 */
