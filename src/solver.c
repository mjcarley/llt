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
#ifdef HAVE_STRING_H
#include <string.h>
#endif /*HAVE_STRING_H*/

#include <glib.h>

#include <blaswrap.h>

#include <amc.h>

#include "llt.h"
#include "llt-private.h"

#include "amc.h"

static gint llt_solver_step_relaxation(llt_assembly_t *A, llt_solver_t *s,
				       gdouble dt)

{
  gint iter, *indices, idx, i, j, nrows, ncols, i1 = 1, np ;
  gdouble err, *U, *ui, *G, *Gu, *ue, *cl, *cd, *Re, al, bt, relax, nu, *vs ;
  gdouble *dn ;
  llt_lifting_line_t *L ;
  
  np = llt_assembly_point_number(A) ;
  indices = A->indices ;
  U = A->U ; ui = A->ui ; G = A->G ; Gu = A->Gu ;
  ue = A->ue ; vs = A->vs ; dn = A->dn ;
  cl = A->cl ; cd = A->cd ; Re = A->Re ;
  nu = llt_solver_viscosity(s) ;
  /*relaxation solver*/
  relax = llt_solver_relaxation_factor(s) ;

  memcpy(Gu, G, np*sizeof(gdouble)) ;
  
  for ( (iter = 0), (err = G_MAXDOUBLE) ;
	(iter < llt_solver_iteration_number_max(s)) &&
	  (err > llt_solver_tolerance(s)) ;
	iter ++ ) {
    err = 0 ;
    /*evaluate self-induced velocity*/
    al = 1.0 ; bt = 0.0 ; nrows = 3*np ; ncols = np ;
    blaswrap_dgemv(FALSE, nrows, ncols, al, U, ncols, Gu, i1, bt, ui, i1) ;
    /*evaluate section quantities*/
    for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
      L = llt_assembly_lifting_line_current(A,i) ;
      for ( j = 0 ; j < llt_lifting_line_point_number(L) ; j ++ ) {
	gdouble ul[3], u2, *c, *n, ch, F, Fl, Fd, ds[3], tmp[3] ;
	gdouble t1, t2, t3, t4, g, slen ;
	llt_lifting_line_segment_unit(ds, L, j, &slen) ;
	/*local velocity*/
	idx = indices[i] + j ;
	ul[0] = ui[3*idx+0] + ue[3*idx+0] ;
	ul[1] = ui[3*idx+1] + ue[3*idx+1] ;
	ul[2] = ui[3*idx+2] + ue[3*idx+2] ;
	c = llt_lifting_line_unit_chord(L,j) ;
	n = llt_lifting_line_normal(L,j) ;
	vs[2*idx+0] = llt_vector_scalar(ul,c) ;
	vs[2*idx+1] = llt_vector_scalar(ul,n) ;
	ch = llt_lifting_line_chord(L,j) ;
	al = atan2(vs[2*idx+1], vs[2*idx+0]) ;
	u2 = vs[2*idx+0]*vs[2*idx+0] + vs[2*idx+1]*vs[2*idx+1] ;
	Re[idx] = ch*sqrt(u2)/nu ;
	if ( isnan(al) ) {
	  g_error("%s: NaN incidence", __FUNCTION__) ;
	}
	llt_section_interp(llt_lifting_line_point_section(L,j),
			   Re[idx], al*180.0/M_PI, &(cl[idx]), &(cd[idx])) ;
	Fl = 0.5*u2*cl[idx]*ch ;
	Fd = 0.5*u2*cd[idx]*ch ;
	F = sqrt(Fl*Fl + Fd*Fd) ;
	llt_vector_cross(tmp,ul,ds) ;
	t1 = llt_vector_scalar(tmp,n) ;
	t2 = llt_vector_scalar(tmp,c) ;
	/*Sugar-Gabor unsteady terms*/
	t3 = llt_vector_scalar(&(dn[3*idx]), n) ;
	t4 = llt_vector_scalar(&(dn[3*idx]), c) ;
	g = (F + ch*G[idx]/dt)/(sqrt(t1*t1 + t2*t2) + ch/dt +
				ch*sqrt(t3*t3 + t4*t4)) ;
	/*(quasi-) steady assumption*/
	/* g = F/sqrt(t1*t1 + t2*t2) ; */
	
	Gu[idx] += relax*(g - Gu[idx]) ;
	if ( fabs(g - Gu[idx]) > 100 ) {
	  g_error("%s: large change in G", __FUNCTION__) ;
	}
	err = MAX(err, fabs((g - Gu[idx]))) ;
	if ( isnan(Gu[idx]) ) {
	  g_error("%s: NaN at entry %d of line %d (%d)",
		  __FUNCTION__, j, i, idx) ;
	}
      }
    }
    fflush(stdout) ;
  }

  memcpy(G, Gu, np*sizeof(gdouble)) ;  
  
  llt_solver_iteration_number(s) = iter ;
  llt_solver_error(s) = err ;
  
  return 0 ;
}

gint llt_solver_step(llt_assembly_t *A, llt_solver_t *solver, gdouble t,
		     gdouble dt, gboolean solve)

{
  gint i, j, np, *indices ;
  llt_lifting_line_t *L, *T, *L0 ;
  llt_lattice_t *W ;
  gdouble *U, *ue, *dn, d ;
  
  if ( solve ) {
    ue = A->ue ;
    dn = A->dn ;
  } else {
    ue = dn = NULL ;
  }
  /*set new lifting line positions and add trailing edge to wake*/
  indices = A->indices ;
  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    L0 = llt_assembly_lifting_line(A,i) ;
    L  = llt_assembly_lifting_line_current(A,i) ;
    if ( llt_assembly_transform_chain(A,i) != NULL ) {
      amc_transform_chain_matrices_evaluate(llt_assembly_transform_chain(A,i),
					    t,
					    llt_assembly_transform(A,i), 1) ;
    }
    if ( ue != NULL ) 
      llt_lifting_line_move(L0, llt_assembly_transform(A,i), L,
			    &(ue[3*indices[i]]), 3, &(dn[3*indices[i]]), 3) ;
    else 
      llt_lifting_line_move(L0, llt_assembly_transform(A,i), L,
			    NULL, 3, NULL, 0) ;
    llt_lattice_wake_shed(llt_assembly_wake(A,i), L, t) ;
  }

  if ( !solve ) return 0 ;

  /*
   * solver: 
   * per time step:
   *   -generate interaction matrix for bound vortices using all LLs
   *   -evaluate externally imposed velocity (motion plus wake)
   *   -iteratively solve for new bound circulation
   */

  U = A->U ;
  np = llt_assembly_point_number(A) ;
  d = llt_assembly_regularisation(A) ;
  
  memset(U, 0, np*3*np*sizeof(gdouble)) ;
  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    L = llt_assembly_lifting_line_current(A,i) ;
    for ( j = 0 ; j < llt_assembly_lifting_line_number(A) ; j ++ ) {
      T = llt_assembly_lifting_line_current(A,j) ;
      llt_lifting_line_velocity_matrix(L, d,
				       llt_lifting_line_point(T,0),
				       LLT_LIFTING_LINE_POINT_SIZE,
				       llt_lifting_line_point_number(T),
				       &(U[3*indices[j]*np+indices[i]]), np) ;
    }
  }

  for ( i = 0 ; i < 3*np*np ; i ++ ) {
    if ( isnan(U[i]) ) {
      g_error("%s: NaN at entry %d of velocity matrix", __FUNCTION__, i) ;
    }
  }
  
  /*reverse direction of velocity due to motion*/
  for ( i = 0 ; i < 3*np ; i ++ ) ue[i] *= -1 ;
  
  /*add wake contributions to externally-imposed velocity*/
  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    L = llt_assembly_lifting_line_current(A,i) ;
    for ( j = 0 ; j < llt_assembly_lifting_line_number(A) ; j ++ ) {
      if ( (W = llt_assembly_wake(A,j)) != NULL ) {
	llt_lattice_velocity_lifting_line(W, L, d, &(ue[3*indices[i]]), 3) ;
      }
    }
  }  

  switch ( llt_solver_method(solver) ) {
  default:
    g_error("%s: unrecognised solver method %d",
	    __FUNCTION__, llt_solver_method(solver)) ;
    break ;
  case LLT_SOLVER_RELAXATION:
    return llt_solver_step_relaxation(A, solver, dt) ;
    break ;
  }    
  
  return 0 ;
}

gint llt_solver_wakes_update(llt_assembly_t *A, gdouble dt)

{
  gint i, j, k, idx, *indices, off, ustr ;
  gdouble *u, al, bt, d, *x, *G, *y1, *y2, *y3, *y4, len ;
  llt_lattice_t *W ;
  
  u = A->uw ; ustr = 3 ;
  d = llt_assembly_regularisation(A) ;
  G = A->G ; indices = A->indices ;
  for ( i = off = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    W = llt_assembly_wake(A,i) ;
    al = 1.0 ; bt = 0.0 ;
    for ( j = 0 ; j < llt_assembly_lifting_line_number(A) ; j ++ ) {
      llt_lifting_line_velocity_lattice(llt_assembly_lifting_line_current(A,j),
					&(G[indices[j]]), 1,
					W, d,
					al, bt, &(u[off*ustr]), ustr) ;
      bt = 1.0 ;
      llt_lattice_velocity_lattice(llt_assembly_wake(A,j), W, d,
				   al, bt, &(u[off*ustr]), ustr) ;
    }
    off += llt_lattice_point_number(W)*llt_lattice_time_steps(W) ;
  }

  for ( i = off = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    W = llt_assembly_wake(A,i) ;
    /*time_steps-1 so that the new wake points remain attached to the
      trailing edge
     */
    for ( j = 0 ;
	  j < llt_lattice_point_number(W)*(llt_lattice_time_steps(W)-1) ;
	  j ++ ) {
      x = &(W->x[LLT_LATTICE_POINT_SIZE*j]) ;
      u[(off+j)*ustr+0] = u[(off+j)*ustr+0]*dt + x[0] ; 
      u[(off+j)*ustr+1] = u[(off+j)*ustr+1]*dt + x[1] ; 
      u[(off+j)*ustr+2] = u[(off+j)*ustr+2]*dt + x[2] ; 
    }
    off += llt_lattice_point_number(W)*llt_lattice_time_steps(W) ;
  }

  /*u now contains updated wake point positions: use them to compute
   vortex stretching effect on circulations*/
  for ( i = off = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    W = llt_assembly_wake(A,i) ;
    for ( j = 1 ; j < llt_lattice_time_steps(W)-1 ; j ++ ) {
      for ( k = 0 ; k < llt_lattice_point_number(W) - 1 ; k ++ ) {
	idx = llt_lattice_index(W,j+0,k+0) ; y1 = &(u[(off+idx)*ustr]) ;
	idx = llt_lattice_index(W,j+0,k+1) ; y2 = &(u[(off+idx)*ustr]) ;
	idx = llt_lattice_index(W,j-1,k+1) ; y3 = &(u[(off+idx)*ustr]) ;
	idx = llt_lattice_index(W,j-1,k+0) ; y4 = &(u[(off+idx)*ustr]) ;
	len = llt_ring_length(y1,y2,y3,y4) ;
	len = llt_lattice_ring_length(W,j,k)/len ;
	llt_lattice_gamma(W,j,k) *= len ;
      }
    }
    /* off += llt_lattice_point_number(W)*llt_lattice_time_steps(W) ; */
  }
  
  for ( i = off = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    W = llt_assembly_wake(A,i) ;
    for ( j = 0 ; j < llt_lattice_time_steps(W)-1 ; j ++ ) {
      for ( k = 0 ; k < llt_lattice_point_number(W) ; k ++ ) {
	idx = llt_lattice_index(W,j,k) ;
	x = llt_lattice_point(W,j,k) ;
	x[0] = u[(off+idx)*ustr+0] ;
	x[1] = u[(off+idx)*ustr+1] ;
	x[2] = u[(off+idx)*ustr+2] ;
      }
    }
    off += llt_lattice_point_number(W)*llt_lattice_time_steps(W) ;
  }
  
  return 0 ;
}
