/* This file is part of LLT, a library for Lifting Line Theory
 * calculations
 *
 * Copyright (C) 2024 Michael Carley
 *
 * LLT is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. LLT is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LLT.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <amc.h>

#include "llt.h"
#include "llt-private.h"

gint llt_lifting_line_wing_rectangular(llt_lifting_line_t *L,
				       llt_section_t *s,
				       gdouble y0, gdouble y1,
				       gdouble ch0, gdouble ch1,
				       gdouble tw0, gdouble tw1,
				       gint n)

/*
 * leading edge at negative x (assuming flow towards positive x,
 * trailing edge)
 */

{
  gdouble tw, ch, *xs0, *xs1, *xc ;
  gint i ;
  
  if ( n > llt_lifting_line_point_number_max(L) ) {
    g_error("%s: not enough points in line (%d) for %d points",
	    __FUNCTION__, llt_lifting_line_point_number_max(L), n) ;
  }

  /*set segment end and wake-shedding points*/
  for ( i = 0 ; i <= n ; i ++ ) {
    ch = ch0 + (ch1 - ch0)*i/n ;
    xs0 = llt_lifting_line_bound_point(L, i) ;
    xs0[0] = 0.0 ;
    xs0[1] = y0 + (y1 - y0)*i/n ;
    xs0[2] = 0.0 ;
    xs0 = llt_lifting_line_wake_point(L, i) ;
    xs0[0] = 0.0 + 0.75*ch ;
    xs0[1] = y0 + (y1 - y0)*i/n ;
    xs0[2] = 0.0 ;
  }

  /*control points*/
  for ( i = 0 ; i < n ; i ++ ) {
    xs0 = llt_lifting_line_bound_point(L, i  ) ;
    xs1 = llt_lifting_line_bound_point(L, i+1) ;
    xc = llt_lifting_line_point(L, i) ;
    xc[0] = 0.5*(xs0[0] + xs1[0]) ;
    xc[1] = 0.5*(xs0[1] + xs1[1]) ;
    xc[2] = 0.5*(xs0[2] + xs1[2]) ;

    ch = ch0 + (ch1 - ch0)*(i+0.5)/n ;
    llt_lifting_line_chord(L, i) = ch ;

    tw = tw0 + (tw1 - tw0)*(i+0.5)/n ;

    xc = llt_lifting_line_normal(L, i) ;
    xc[0] = sin(tw) ;
    xc[1] = 0.0 ;
    xc[2] = cos(tw) ;
    xc = llt_lifting_line_unit_chord(L, i) ;
    xc[0] = cos(tw) ;
    xc[1] = 0.0 ;
    xc[2] = -sin(tw) ;    
    llt_lifting_line_point_section(L, i) = s ;
  }

  llt_lifting_line_point_number(L) = n ;
  
  return 0 ;
}

gint llt_lifting_line_wing_elliptical(llt_lifting_line_t *L,
				      llt_section_t *s,
				      gdouble a, gdouble b,
				      gint n)

/*
 * leading edge at negative x (assuming flow towards positive x,
 * trailing edge)
 */

{
  gdouble tw, ch, *xs0, *xs1, *xc, y ;
  gint i ;
  
  if ( n > llt_lifting_line_point_number_max(L) ) {
    g_error("%s: not enough points in line (%d) for %d points",
	    __FUNCTION__, llt_lifting_line_point_number_max(L), n) ;
  }

  /*set segment end and wake-shedding points*/
  for ( i = 0 ; i <= n ; i ++ ) {
    y = -0.5*b + b*i/n ;
    ch = a*sqrt(1.0-4.0*y*y/b/b) ;
    xs0 = llt_lifting_line_bound_point(L, i) ;
    xs0[0] = 0.0 ;
    xs0[1] = y ;
    xs0[2] = 0.0 ;
    xs0 = llt_lifting_line_wake_point(L, i) ;
    xs0[0] = 0.0 + 0.75*ch ;
    xs0[1] = y ;
    xs0[2] = 0.0 ;
  }

  /*control points*/
  for ( i = 0 ; i < n ; i ++ ) {
    xs0 = llt_lifting_line_bound_point(L, i  ) ;
    xs1 = llt_lifting_line_bound_point(L, i+1) ;
    xc = llt_lifting_line_point(L, i) ;
    y = 0.5*(xs0[1] + xs1[1]) ;
    xc[0] = 0.5*(xs0[0] + xs1[0]) ;
    xc[1] = y ;
    xc[2] = 0.5*(xs0[2] + xs1[2]) ;

    ch = a*sqrt(1.0-4.0*y*y/b/b) ;
    
    llt_lifting_line_chord(L, i) = ch ;
    tw = 0 ;
    
    /* tw = tw0 + (tw1 - tw0)*(i+0.5)/n ; */

    xc = llt_lifting_line_normal(L, i) ;
    xc[0] = sin(tw) ;
    xc[1] = 0.0 ;
    xc[2] = cos(tw) ;
    xc = llt_lifting_line_unit_chord(L, i) ;
    xc[0] = cos(tw) ;
    xc[1] = 0.0 ;
    xc[2] = -sin(tw) ;
    llt_lifting_line_point_section(L, i) = s ;
  }

  llt_lifting_line_point_number(L) = n ;
  
  return 0 ;
}

gint llt_lifting_line_wing_helical(llt_lifting_line_t *L,
				   llt_section_t *s,
				   gdouble r0, gdouble r1,
				   gdouble ch,
				   gdouble U, gdouble Om,
				   gint n)

/*
 * leading edge at negative x (assuming flow towards positive x,
 * trailing edge)
 */

{
  gdouble tw, *xs0, *xs1, *xc, r ;
  gint i ;
  
  if ( n > llt_lifting_line_point_number_max(L) ) {
    g_error("%s: not enough points in line (%d) for %d points",
	    __FUNCTION__, llt_lifting_line_point_number_max(L), n) ;
  }

  /*set segment end and wake-shedding points*/
  for ( i = 0 ; i <= n ; i ++ ) {
    xs0 = llt_lifting_line_bound_point(L, i) ;
    r = r0 + (r1 - r0)*i/n ;
    xs0[0] = 0.0 ;
    xs0[1] = r ;
    xs0[2] = 0.0 ;
    xs0 = llt_lifting_line_wake_point(L, i) ;

    tw = -atan2(U,Om*r) ;
    
    xs0[0] = 0.0 + 0.75*ch*cos(tw) ;
    xs0[1] = r ; 
    xs0[2] = -0.75*ch*sin(tw) ;
  }

  /*control points*/
  for ( i = 0 ; i < n ; i ++ ) {
    xs0 = llt_lifting_line_bound_point(L, i  ) ;
    xs1 = llt_lifting_line_bound_point(L, i+1) ;
    xc = llt_lifting_line_point(L, i) ;
    r = 0.5*(xs0[1] + xs1[1]) ;
    xc[0] = 0.5*(xs0[0] + xs1[0]) ;
    xc[1] = r ;
    xc[2] = 0.5*(xs0[2] + xs1[2]) ;

    llt_lifting_line_chord(L, i) = ch ;

    tw = -atan2(U,Om*r) ;
    
    xc = llt_lifting_line_normal(L, i) ;
    xc[0] = sin(tw) ;
    xc[1] = 0.0 ;
    xc[2] = cos(tw) ;
    xc = llt_lifting_line_unit_chord(L, i) ;
    xc[0] = cos(tw) ;
    xc[1] = 0.0 ;
    xc[2] = -sin(tw) ;    
    llt_lifting_line_point_section(L, i) = s ;
  }

  llt_lifting_line_point_number(L) = n ;
  
  return 0 ;
}
