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

llt_assembly_t *llt_assembly_alloc(void)

{
  llt_assembly_t *a ;

  a = (llt_assembly_t *)g_malloc0(sizeof(llt_assembly_t)) ;

  llt_assembly_lifting_line_number(a) = 0 ;
  
  return a ;
}

gint llt_assembly_lifting_line_add(llt_assembly_t *A, llt_lifting_line_t *L,
				   amc_transform_chain_t *C)

{
  gint n ;
  
  if ( (n = llt_assembly_lifting_line_number(A)) >=
       LLT_ASSEMBLY_LINE_NUMBER_MAX ) {
    g_error("%s: maximum number of lines (%d) exceeded",
	    __FUNCTION__, LLT_ASSEMBLY_LINE_NUMBER_MAX) ;
  }

  llt_assembly_lifting_line(A,n) = L ;
  llt_assembly_transform(A,n) = amc_transform_alloc(3, 1) ;
  llt_assembly_transform_chain(A,n) = C ;

  llt_assembly_lifting_line_number(A) ++ ;
  
  return 0 ;
}

gint llt_assembly_initialise(llt_assembly_t *A, gdouble t)
  
{
  gint i, np, nw ;
  llt_lifting_line_t *L, *L0 ;
  
  /*count the lifting lines and initialise the indices*/
  A->indices[0] = 0 ; nw = 0 ;
  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    L = llt_assembly_lifting_line(A,i) ;
    A->indices[i+1] = A->indices[i] + llt_lifting_line_point_number(L) ;
    llt_assembly_lifting_line_current(A,i) =
      llt_lifting_line_alloc(llt_lifting_line_point_number(L)) ;
    llt_assembly_wake(A,i) =
      llt_lattice_alloc(llt_lifting_line_point_number(L)+1,
			llt_assembly_time_steps(A)) ;
    /*total number of wake points to allocate wake velocity buffer*/
    nw += llt_lattice_point_number_max(llt_assembly_wake(A,i)) ;
  }

  np = A->indices[llt_assembly_lifting_line_number(A)] ;

  /*allocate lifting line and solver data and connect circulations to
   lifting lines*/
  A->G = (gdouble *)g_malloc0((np+         /*current \Gamma*/
			       np+         /*updated \Gamma*/
			       3*np+       /*time derivative of normals*/
			       3*np+       /*externally imposed velocity*/
			       2*np+       /*in-plane velocity*/
			       3*np+       /*induced velocity*/
			       np*3*np+    /*induced velocity matrix*/
			       np+         /*local \alpha*/
			       np+         /*c_l*/
			       np+         /*c_d*/
			       np)         /*Re*/
			      *sizeof(gdouble)) ;
  A->Gu = &(A->G[np]) ;
  A->dn = &(A->Gu[np]) ;
  A->ue = &(A->dn[3*np]) ;
  A->vs = &(A->ue[3*np]) ;
  A->ui = &(A->vs[2*np]) ;
  A->U  = &(A->ui[3*np]) ;
  A->al = &(A->U[3*np*np]) ;
  A->cl = &(A->al[np]) ;
  A->cd = &(A->cl[np]) ;
  A->Re = &(A->cd[np]) ;

  A->uw = (gdouble *)g_malloc0(3*nw*sizeof(gdouble)) ;
  
  /*set current lifting lines to position at time t, and generate the first
   row of the wake lattices*/
  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    L0 = llt_assembly_lifting_line(A,i) ;
    L  = llt_assembly_lifting_line_current(A,i) ;
    L->G = &(A->G[A->indices[i]]) ;
    if ( llt_assembly_transform_chain(A,i) != NULL ) {
      amc_transform_chain_matrices_evaluate(llt_assembly_transform_chain(A,i),
					    t,
					    llt_assembly_transform(A,i), 1) ;
    }
    llt_lifting_line_move(L0, llt_assembly_transform(A,i),
			  L, NULL, 0, NULL, 0) ;
    llt_lattice_initialise(llt_assembly_wake(A,i), L, t, 0) ;
  }
  
  return 0 ;
}

gint llt_assembly_lifting_line_force(llt_assembly_t *A, gint i, gdouble *F)

{
  gdouble *cl, *cd, len, *c, *n, ch, u2, *vs ;
  gint *indices, j, idx ;
  llt_lifting_line_t *L ;
  
  indices = A->indices ;
  vs = A->vs ; cl = A->cl ; cd = A->cd ;

  F[0] = F[1] = F[2] = 0.0 ;

  L = llt_assembly_lifting_line_current(A,i) ;
  for ( j = 0 ; j < llt_lifting_line_point_number(L) ; j ++ ) {
    len = llt_lifting_line_segment_length(L, j) ;
    idx = indices[i] + j ;
    c = llt_lifting_line_unit_chord(L,j) ;
    n = llt_lifting_line_normal(L,j) ;
    u2 = vs[2*idx+0]*vs[2*idx+0] + vs[2*idx+1]*vs[2*idx+1] ;
    ch = llt_lifting_line_chord(L,j) ;
    F[0] += cl[idx]*0.5*u2*ch*len*n[0] ;
    F[1] += cl[idx]*0.5*u2*ch*len*n[1] ;
    F[2] += cl[idx]*0.5*u2*ch*len*n[2] ;
    F[0] += cd[idx]*0.5*u2*ch*len*c[0] ;
    F[1] += cd[idx]*0.5*u2*ch*len*c[1] ;
    F[2] += cd[idx]*0.5*u2*ch*len*c[2] ;
  }
  
  return 0 ;
}

gint llt_assembly_lifting_line_force_moment(llt_assembly_t *A, gint i,
					    gdouble *F,
					    gdouble *x, gdouble *M)

{
  gdouble *cl, *cd, len, *c, *n, ch, u2, *vs, df[3], dM[3], *xc, r[3] ;
  gint *indices, j, idx ;
  llt_lifting_line_t *L ;
  
  indices = A->indices ;
  vs = A->vs ; cl = A->cl ; cd = A->cd ;

  F[0] = F[1] = F[2] = 0.0 ;
  M[0] = M[1] = M[2] = 0.0 ;

  L = llt_assembly_lifting_line_current(A,i) ;
  for ( j = 0 ; j < llt_lifting_line_point_number(L) ; j ++ ) {
    len = llt_lifting_line_segment_length(L, j) ;
    xc = llt_lifting_line_point(L, j) ;
    llt_vector_diff(r,xc,x) ;
    idx = indices[i] + j ;
    c = llt_lifting_line_unit_chord(L,j) ;
    n = llt_lifting_line_normal(L,j) ;
    u2 = vs[2*idx+0]*vs[2*idx+0] + vs[2*idx+1]*vs[2*idx+1] ;
    ch = llt_lifting_line_chord(L,j) ;
    df[0] = cl[idx]*n[0] + cd[idx]*c[0] ; 
    df[1] = cl[idx]*n[1] + cd[idx]*c[1] ; 
    df[2] = cl[idx]*n[2] + cd[idx]*c[2] ;
    df[0] *= 0.5*u2*ch*len ;
    df[1] *= 0.5*u2*ch*len ;
    df[2] *= 0.5*u2*ch*len ;
    F[0] += df[0] ; F[1] += df[1] ; F[2] += df[2] ;
    llt_vector_cross(dM,r,df) ;
    M[0] += dM[0] ; M[1] += dM[1] ; M[2] += dM[2] ; 
  }
  
  return 0 ;
}

gint llt_assembly_force(llt_assembly_t *A, gdouble *F)

{
  gint i ;
  gdouble dF[3] ;

  F[0] = F[1] = F[2] = 0.0 ;
  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    llt_assembly_lifting_line_force(A, i, dF) ;
    F[0] += dF[0] ; F[1] += dF[1] ; F[2] += dF[2] ; 
  }
  
  return 0 ;
}

gint llt_assembly_force_moment(llt_assembly_t *A, gdouble *F,
			       gdouble *x, gdouble *M)

{
  gint i ;
  gdouble dF[3], dM[3] ;

  F[0] = F[1] = F[2] = 0.0 ;
  M[0] = M[1] = M[2] = 0.0 ;
  
  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    llt_assembly_lifting_line_force_moment(A, i, dF, x, dM) ;
    F[0] += dF[0] ; F[1] += dF[1] ; F[2] += dF[2] ; 
    M[0] += dM[0] ; M[1] += dM[1] ; M[2] += dM[2] ; 
  }
  
  return 0 ;
}
