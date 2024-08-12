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

llt_section_t *llt_section_alloc(gint nRe, gint npts)

{
  llt_section_t *s ;

  s = (llt_section_t *)g_malloc0(sizeof(llt_section_t)) ;

  s->iRe = (gint *)g_malloc0((nRe+1)*sizeof(gint)) ;
  s->nRemax = nRe ;

  s->data = (gdouble *)g_malloc0(npts*LLT_SECTION_DATA_SIZE*sizeof(gdouble)) ;
  s->nptsmax = npts ;
  
  return s ;
}

gint llt_section_write(FILE *f, llt_section_t *s)

{
  gint i, j ;

  if ( llt_section_point_number(s) !=
       llt_section_Re_index(s,llt_section_Re_number(s))) {
    g_error("%s: final index (%d) does not match number of data points (%d)",
	    __FUNCTION__, llt_section_Re_index(s,llt_section_Re_number(s)),
	    llt_section_point_number(s)) ;
  }
  
  fprintf(f, "%s\n", llt_section_name(s)) ;

  fprintf(f, "%d:", llt_section_Re_number(s)) ;
  for ( i = 0 ; i <= llt_section_Re_number(s) ; i ++ ) {
    fprintf(f, " %d", llt_section_Re_index(s,i)) ;
  }
  fprintf(f, "\n") ;

  for ( i = 0 ; i < llt_section_point_number(s) ; i ++ ) {
    fprintf(f, "%e", s->data[i*LLT_SECTION_DATA_SIZE+0]) ;
    for ( j = 1 ; j < LLT_SECTION_DATA_SIZE ; j ++ ) {
      fprintf(f, " %e", s->data[i*LLT_SECTION_DATA_SIZE+j]) ;
    }
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

gint llt_section_read(FILE *f, llt_section_t *s)

{
  gchar buf[256] ;
  gint n, i ;
  
  fscanf(f, "%[^\n]s", buf) ;
  llt_section_name(s) = g_strdup(buf) ;

  fscanf(f, "%d", &n) ;

  if ( n <= 0 ) {
    g_error("%s: invalid number (%d) of Reynolds numbers", __FUNCTION__, n) ;
  }
  if ( n > llt_section_Re_number_max(s) ) {
    g_error("%s: not enough space (%d) for %d Reynolds numbers",
	    __FUNCTION__, llt_section_Re_number_max(s), n) ;
  }
  llt_section_Re_number(s) = n ;
  
  fscanf(f, "%*c") ;
  for ( i = 0 ; i <= llt_section_Re_number(s) ; i ++ ) {
    fscanf(f, "%d", &(llt_section_Re_index(s,i))) ;
  }

  if ( llt_section_point_number_max(s) <
       llt_section_Re_index(s,llt_section_Re_number(s))) {
    g_error("%s: not enough space (%d) for %d data points",
	    __FUNCTION__, llt_section_point_number_max(s),
	    llt_section_Re_index(s,llt_section_Re_number(s))) ;
  }

  llt_section_point_number(s) =
    llt_section_Re_index(s,llt_section_Re_number(s)) ;
  for ( i = 0 ;
	i < llt_section_point_number(s)*LLT_SECTION_DATA_SIZE ;
	i ++ ) {
    fscanf(f, "%lg", &(s->data[i])) ;
  }
  
  return 0 ;
}

llt_section_t *llt_section_read_alloc(FILE *f, gint m)

{

  gchar buf[256] ;
  gint nRe, i, iRe[128] ;
  llt_section_t *s ;

  if ( m < 0 ) {
    g_error("%s: margin (%d) must be non-negative", __FUNCTION__, m) ;
  }
  
  fscanf(f, "%[^\n]s", buf) ;
  
  fscanf(f, "%d", &nRe) ;

  if ( nRe <= 0 ) {
    g_error("%s: invalid number (%d) of Reynolds numbers", __FUNCTION__, nRe) ;
  }
  if ( nRe > 127 ) {
    g_error("%s: not enough buffer space for %d Reynolds numbers",
	    __FUNCTION__, nRe) ;
  }
  fscanf(f, "%*c") ;
  for ( i = 0 ; i <= nRe ; i ++ ) {
    fscanf(f, "%d", &(iRe[i])) ;
  }

  s = llt_section_alloc(nRe, iRe[nRe]+m) ;
  llt_section_name(s) = g_strdup(buf) ;
  
  llt_section_Re_number(s) = nRe ;  
  for ( i = 0 ; i <= nRe ; i ++ ) {
    llt_section_Re_index(s,i) = iRe[i] ;
  }

  llt_section_point_number(s) =
    llt_section_Re_index(s,llt_section_Re_number(s)) ;
  for ( i = 0 ;
	i < llt_section_point_number(s)*LLT_SECTION_DATA_SIZE ;
	i ++ ) {
    fscanf(f, "%lg", &(s->data[i])) ;
  }

  llt_section_interp_Re(s)    = LLT_INTERPOLATE_LINEAR ;
  llt_section_interp_alpha(s) = LLT_INTERPOLATE_LINEAR ;
  
  return s ;
}

static void interp_alpha(gdouble *al0, gint astr,
			 gdouble *f, gint fstr,
			 gint na,
			 gdouble al,
			 gdouble *fi, gint nf,
			 llt_interpolate_t interp,
			 llt_extrapolate_t elow,
			 llt_extrapolate_t ehigh)

{
  gdouble t ;
  gint i, j ;
  
  g_assert(interp == LLT_INTERPOLATE_LINEAR) ;

  if ( al < al0[0] ) {
    /*hold C_L constant below minimum defined \alpha*/
    for ( j = 0 ; j < nf ; j ++ ) {
      fi[j] = f[0*fstr+j] ;
    }
    return ;
  }
  if ( al > al0[(na-1)*astr] ) {
    /*hold C_L constant above maximum defined \alpha*/
    for ( j = 0 ; j < nf ; j ++ ) {
      fi[j] = f[(na-1)*fstr+j] ;
    }
    return ;
  }
  
  for ( i = 1 ; (i < na) && (al > al0[i*astr]) ; i ++ ) ;

  g_assert(i < na) ;
  i -- ;
  g_assert(al >= al0[i*astr]) ;
  g_assert(al < al0[(i+1)*astr]) ;
  
  t = (al - al0[i*astr])/(al0[(i+1)*astr]-al0[i*astr]) ;

  g_assert((t >= 0) && (t <= 1)) ;
  
  for ( j = 0 ; j < nf ; j ++ ) {
    fi[j] = f[i*fstr+j]*(1-t) + f[(i+1)*fstr+j]*t ;
  }
  
  return ;
}

static void index_locate(gdouble *data, gint dstr, gint *idx, gint ni,
			 gdouble x,
			 llt_extrapolate_t elow,
			 llt_extrapolate_t ehigh,
			 gint *i0, gint *i1,
			 gdouble *t)

{
  gint i ;
  gdouble x0, x1 ;

  if ( ni == 1 ) {
    *i0 = *i1 = 0 ; *t = 0 ;
    return ;
  }
  
  if ( x < data[0] ) {
    if ( elow == LLT_EXTRAPOLATE_CONSTANT ) {
      *i0 = *i1 = 0 ; *t = 0 ;
      return ;
    } else {
      *i0 = 0 ; *i1 = 1 ;
      x0 = data[idx[*i0]*dstr] ; x1 = data[idx[*i1]*dstr] ;
      *t = (x - x0)/(x1 - x0) ;
      return ;
    }
  }
  if ( x > data[idx[ni-1]] ) {
    if ( ehigh == LLT_EXTRAPOLATE_CONSTANT ) {
      *i0 = *i1 = ni-1 ; *t = 0 ;
      return ;
    } else {
      *i0 = ni-2 ; *i1 = ni-1 ;
      x0 = data[idx[*i0]*dstr] ; x1 = data[idx[*i1]*dstr] ;
      *t = (x - x0)/(x1 - x0) ;
      return ;
    }
  }

  
  for ( (*i1) = 1 ; ((*i1) < ni) && (x > data[idx[(*i1)]*dstr]) ;
	(*i1) ++ ) ;
  /* { */
  /*   fprintf(stderr, "%d: %e %e %e\n", */
  /* 	    (*i0), data[idx[(*i0)]*dstr], x, data[idx[(*i1)]*dstr]) ; */
  /* } */
  *i0 = *i1 - 1 ;
  for ( i = 0 ; (i < ni-1) && (x > data[idx[i]*dstr]) ;
	i ++ ) ;
  
  x0 = data[idx[*i0]*dstr] ; x1 = data[idx[*i1]*dstr] ;
  *t = (x - x0)/(x1 - x0) ;
  g_assert(x0 <= x) ;
  g_assert(x1 > x) ;
  
  return ;
}

gint llt_section_interp(llt_section_t *s, gdouble Re, gdouble al,
			gdouble *cl, gdouble *cd)

{
  gint ia, na, i0, i1 ;
  gdouble data[8]={0}, t ;

  index_locate(s->data, 4, s->iRe, s->nRe, Re,
	       llt_section_extrap_Re_low(s),
	       llt_section_extrap_Re_high(s),
	       &i0, &i1, &t) ;
  
  ia = llt_section_Re_index(s,i0  ) ;
  na = llt_section_Re_index(s,i0+1) - llt_section_Re_index(s,i0) ;
  interp_alpha(&(s->data[ia*LLT_SECTION_DATA_SIZE+1]), LLT_SECTION_DATA_SIZE,
	       &(s->data[ia*LLT_SECTION_DATA_SIZE+2]), LLT_SECTION_DATA_SIZE,
	       na,
	       al,
	       data, LLT_SECTION_DATA_SIZE-2,
	       llt_section_interp_alpha(s),
	       llt_section_extrap_alpha_low(s),
	       llt_section_extrap_alpha_high(s)) ;

  *cl = data[0] ; *cd = data[1] ;

  if ( t == 0 ) return 0 ;  

  ia = llt_section_Re_index(s,i1  ) ;
  na = llt_section_Re_index(s,i1+1) - llt_section_Re_index(s,i1) ;
  interp_alpha(&(s->data[ia*LLT_SECTION_DATA_SIZE+1]), LLT_SECTION_DATA_SIZE,
	       &(s->data[ia*LLT_SECTION_DATA_SIZE+2]), LLT_SECTION_DATA_SIZE,
	       na,
	       al,
	       data, LLT_SECTION_DATA_SIZE-2,
	       llt_section_interp_alpha(s),
	       llt_section_extrap_alpha_low(s),
	       llt_section_extrap_alpha_high(s)) ;

  *cl = *cl + (data[0] - *cl)*t;
  *cd = *cd + (data[1] - *cd)*t;
  
  return 0 ;
}

gint llt_section_linear(llt_section_t *s, gdouble a)

{
  if ( llt_section_point_number_max(s) < 2 ) {
    g_error("%s: not enough space (%d) for two-point interpolation",
	    __FUNCTION__, llt_section_point_number_max(s)) ;
  }
  llt_section_point_number(s) = 2 ;
  llt_section_Re_number(s) = 1 ;
  llt_section_Re_index(s,0) = 0 ;
  llt_section_Re_index(s,1) = 2 ;  

  s->data[0] = 1 ; s->data[1] = -180 ; s->data[2] = -a*M_PI ; s->data[3] = 0 ;
  s->data[4] = 1 ; s->data[5] =  180 ; s->data[6] =  a*M_PI ; s->data[7] = 0 ;

  llt_section_interp_Re(s)    = LLT_INTERPOLATE_LINEAR ;
  llt_section_interp_alpha(s) = LLT_INTERPOLATE_LINEAR ;

  llt_section_extrap_Re_low(s)     = LLT_EXTRAPOLATE_CONSTANT ;
  llt_section_extrap_Re_high(s)    = LLT_EXTRAPOLATE_CONSTANT ;
  llt_section_extrap_alpha_low(s)  = LLT_EXTRAPOLATE_EXTEND ;
  llt_section_extrap_alpha_high(s) = LLT_EXTRAPOLATE_EXTEND ;
  
  return 0 ;
}

gint llt_section_entry_append(llt_section_t *s,
			      gdouble Re, gdouble al,
			      gdouble cl, gdouble cd)

{
  gint n ;
  
  if ( llt_section_point_number(s) >= llt_section_point_number_max(s) ) {
    g_error("%s: available space (%d entries) exceeded (%d)\n",
	    __FUNCTION__,
	    llt_section_point_number_max(s),
	    llt_section_point_number(s)) ;
  }  

  n = llt_section_point_number(s) ;
  llt_section_Re(s,n)    = Re ; 
  llt_section_alpha(s,n) = al ; 
  llt_section_cl(s,n)    = cl ; 
  llt_section_cd(s,n)    = cd ; 
  llt_section_point_number(s) ++ ;
  
  return 0 ;
}

static gint compare_section(const void *a, const void *b)

{
  const gdouble *A, *B ;

  A = (gdouble *)a ; B = (gdouble *)b ; 

  /*sort on Reynolds number first*/
  if ( A[LLT_SECTION_DATA_RE] < B[LLT_SECTION_DATA_RE] ) return -1 ;
  if ( A[LLT_SECTION_DATA_RE] > B[LLT_SECTION_DATA_RE] ) return  1 ;

  /*same Reynolds number, sort on \alpha*/
  if ( A[LLT_SECTION_DATA_ALPHA] < B[LLT_SECTION_DATA_ALPHA] ) return -1 ;
  if ( A[LLT_SECTION_DATA_ALPHA] > B[LLT_SECTION_DATA_ALPHA] ) return  1 ;
  
  return 0 ;
}

gint llt_section_sort(llt_section_t *s)

{
  qsort(s->data, llt_section_point_number(s),
	LLT_SECTION_DATA_SIZE*sizeof(gdouble), compare_section) ;
  
  return 0 ;
}

gint llt_section_index(llt_section_t *s)

{
  gint i ;

  llt_section_Re_index(s,0) = 0 ;
  llt_section_Re_number(s)  = 1 ;

  for ( i = 1 ; i < llt_section_point_number(s) ; i ++ ) {
    if ( llt_section_Re(s,i) < llt_section_Re(s,i-1) ) {
      g_error("%s: section data must be sorted before indexing (entry %d)",
	      __FUNCTION__, i) ;
    }
    if ( llt_section_Re(s,i) > llt_section_Re(s,i-1) ) {
      llt_section_Re_index(s,llt_section_Re_number(s)) = i ;
      llt_section_Re_number(s) ++ ;
    }
  }

  llt_section_Re_index(s,llt_section_Re_number(s)) =
    llt_section_point_number(s) ;
  
  return 0 ;
}
