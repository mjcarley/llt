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


/**
 * @addtogroup lattice
 * @{
 */

/** 
 * Allocate a lattice to store wake data
 * 
 * @param np (maximum) number of points on lattice edge: this should be equal
 * to the number of wake (trailing edge) points on the lifting line shedding
 * the wake;
 * @param nt (maximum) number of time steps (shed rows) in lattice.
 * 
 * @return newly allocated ::llt_lattice_t

 */

llt_lattice_t *llt_lattice_alloc(gint np, gint nt)
  
{
  llt_lattice_t *l ;

  if ( np <= 0 ) {
    g_error("%s: invalid number (%d) of wake points", __FUNCTION__, np) ;
  }

  if ( nt <= 0 ) {
    g_error("%s: invalid number (%d) of time points", __FUNCTION__, nt) ;
  }
  
  l = (llt_lattice_t *)g_malloc0(sizeof(llt_lattice_t)) ;

  /*nodes and circulations on lattice*/
  l->x = (gdouble *)g_malloc0(LLT_LATTICE_POINT_SIZE*np*nt*sizeof(gdouble)) ;
  llt_lattice_point_number_max(l) = np*nt ;
  llt_lattice_point_number(l) = 0 ;
  llt_lattice_time_steps(l) = 0 ;
  
  return l ;
}

/** 
 * Evaluate velocity induced by lattice at field point
 * 
 * @param l vortex lattice;
 * @param d regularization parameter for velocity evaluation;
 * @param x field point for evaluation of velocity;
 * @param al scaling factor for lattice velocity;
 * @param bt scaling factor for initial velocity;
 * @param u on entry contains \f$\mathbf{u}_{0}\f$, on exit contains
 * \f$\beta \mathbf{u}_{0} + \alpha\mathbf{u}_{l}\f$ with
 * \f$\mathbf{u}_{l}\f$ velocity generated by lattice \a l.
 * 
 * @return 0 on success.
 */

gint llt_lattice_velocity(llt_lattice_t *l, gdouble d, gdouble *x,
			  gdouble al, gdouble bt, gdouble *u)

{
  gdouble du[3], *y1, *y2, G ;
  gint i, j ;
  
  u[0] *= bt ; u[1] *= bt ; u[2] *= bt ;

  for ( i = 1 ; i < llt_lattice_time_steps(l) ; i ++ ) {
    for ( j = 0 ; j < llt_lattice_point_number(l)-1 ; j ++ ) {
      G = llt_lattice_gamma(l,i,j) ;
      y1 = llt_lattice_point(l,i+0,j+0) ;
      y2 = llt_lattice_point(l,i+0,j+1) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      u[0] += al*G*du[0] ; u[1] += al*G*du[1] ; u[2] += al*G*du[2] ;
      y1 = llt_lattice_point(l,i+0,j+1) ;
      y2 = llt_lattice_point(l,i-1,j+1) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      u[0] += al*G*du[0] ; u[1] += al*G*du[1] ; u[2] += al*G*du[2] ;
      y1 = llt_lattice_point(l,i-1,j+1) ;
      y2 = llt_lattice_point(l,i-1,j+0) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      u[0] += al*G*du[0] ; u[1] += al*G*du[1] ; u[2] += al*G*du[2] ;
      y1 = llt_lattice_point(l,i-1,j+0) ;
      y2 = llt_lattice_point(l,i+0,j+0) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      u[0] += al*G*du[0] ; u[1] += al*G*du[1] ; u[2] += al*G*du[2] ;
    }
  }

  return 0 ;
}

/** 
 * Evaluate lattice-induced velocity at each collocation point of a lifting
 * line
 * 
 * @param W lattice (wake);
 * @param L lifting line;
 * @param d regularization parameter for velocity evaluation;
 * @param u pointer to array of velocities;
 * @param ustr stride between successive velocities.
 * 
 * @return 0 on success.
 */

gint llt_lattice_velocity_lifting_line(llt_lattice_t *W,
				       llt_lifting_line_t *L,
				       gdouble d,
				       gdouble *u, gint ustr)

{
  gint i ;
  gdouble *x ;

  if ( ustr < 3 ) {
    g_error("%s: velocity stride (%d) must be at least three",
	    __FUNCTION__, ustr) ;
  }
  
  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_point(L,i) ;
    llt_lattice_velocity(W, d, x, 1.0, 1.0, &(u[i*ustr])) ;
  }
  
  return 0 ;
}

/** 
 * Write vortex lattice data to output as nodes of lattice and
 * circulations
 * 
 * @param f output file stream;
 * @param W lattice to write.
 * 
 * @return 0 on success.
 */

gint llt_lattice_write(FILE *f, llt_lattice_t *W)

{
  gint i, j ;
  gdouble *x ;

  for ( i = 0 ; i < llt_lattice_time_steps(W) ; i ++ ) {
    for ( j = 0 ; j < llt_lattice_point_number(W) ; j ++ ) {
      x = llt_lattice_point(W,i,j) ;
      fprintf(f, "%e %e %e %e ",
	      x[0], x[1], x[2], llt_lattice_gamma(W,i,j)) ;
    }
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

/** 
 * Initialise a vortex lattice by connection to a lifting line
 * 
 * @param W vortex lattice;
 * @param L lifting line from which vortex lattice in \a W is being shed;
 * @param t initial time;
 * @param G0 initial circulation.
 * 
 * @return 0 on success.
 */

gint llt_lattice_initialise(llt_lattice_t *W, llt_lifting_line_t *L,
			    gdouble t, double G0)

{
  gint i ;
  gdouble *xl, *xw ;
  
  llt_lattice_time_steps(W) = 0 ;
  llt_lattice_point_number(W) = llt_lifting_line_point_number(L) + 1 ;
  for ( i = 0 ; i < llt_lattice_point_number(W) ; i ++ ) {
    xw = llt_lifting_line_wake_point(L,i) ;
    xl = llt_lattice_point(W,0,i) ;
    xl[0] = xw[0] ; xl[1] = xw[1] ; xl[2] = xw[2] ;
    llt_lattice_gamma(W,0,i) = G0 ;
  }
  llt_lattice_time_shed(W,0) = t ;
  llt_lattice_time_steps(W) = 1 ;
  
  return 0 ;
}

/** 
 * Update a vortex lattice by adding shed vorticity from a lifting line
 * 
 * @param W wake lattice;
 * @param L lifting line generating wake in \a W;
 * @param t shedding time (used to calculate viscous effects).
 * 
 * @return 0 on success.
 */

gint llt_lattice_wake_shed(llt_lattice_t *W, llt_lifting_line_t *L, gdouble t)

{
  gint i, nt, np ;
  gdouble *xw, *xl ;
  
  nt = llt_lattice_time_steps(W) ;
  np = llt_lattice_point_number(W) ;

  if ( (nt+1)*np > llt_lattice_point_number_max(W) ) {
    g_error("%s: not enough space for %d (%dx%d) points",
	    __FUNCTION__, nt*np, nt, np) ;
  }

  if ( np != llt_lifting_line_point_number(L)+1 ) {
    g_error("%s: wake point number (%d) does not match "
	    "lifting line shedding point number (%d)",
	    __FUNCTION__, np, llt_lifting_line_point_number(L)+1) ;
  }
  
  for ( i = 0 ; i < np ; i ++ ) {
    xw = llt_lattice_point(W,nt,i) ;
    xl = llt_lifting_line_wake_point(L,i) ;
    xw[0] = xl[0] ; xw[1] = xl[1] ; xw[2] = xl[2] ; 
  }

  llt_lattice_time_shed(W,nt) = t ;
  llt_lattice_time_steps(W) ++ ;

  /*if no circulation is defined on the lifting line, return*/
  if ( L->G == NULL ) return 0 ;

  for ( i = 0 ; i < np-1 ; i ++ ) {
    llt_lattice_gamma(W,nt,i) = L->G[i] ;
  }
    
  return 0 ;
}

/** 
 * Evaluate velocity "induced" by one vortex lattice on another lattice
 * 
 * @param S source lattice;
 * @param T target lattice;
 * @param d regularisation parameter;
 * @param al \f$\alpha\f$;
 * @param bt \f$\beta\f$;
 * @param u on entry contains \f$\mathbf{u}_{0}\f$, on exit contains
 * \f$\beta \mathbf{u}_{0} + \alpha\mathbf{u}_{S}\f$ with 
 * \f$\mathbf{u}_{S}\f$ velocity generated by lattice \a S;
 * @param ustr stride in \a u.
 * 
 * @return 0 on success.
 */

gint llt_lattice_velocity_lattice(llt_lattice_t *S, llt_lattice_t *T,
				  gdouble d, gdouble al, gdouble bt,
				  gdouble *u, gint ustr)

{
  gint i, j, idx ;
  gdouble *x ;
  
  for ( i = 0 ; i < llt_lattice_time_steps(T) ; i ++ ) {
    for ( j = 0 ; j < llt_lattice_point_number(T) ; j ++ ) {
      idx = llt_lattice_index(T, i, j) ;
      u[idx*ustr+0] *= bt ; u[idx*ustr+1] *= bt ; u[idx*ustr+2] *= bt ; 
      x = llt_lattice_point(T, i, j) ;
      llt_lattice_velocity(S, d, x, al, 1.0, &(u[idx*ustr])) ;
    }
  }

  return 0 ;
}

/** 
 * Convenience function for length of vortex ring in lattice
 * 
 * @param l an ::llt_lattice_t;
 * @param i time step when vortex ring was shed;
 * @param j index of point from which vortex ring was shed.
 * 
 * @return length of vortex ring with corner at point \f$(i,j)\f$.
 */

gdouble llt_lattice_ring_length(llt_lattice_t *l, gint i, gint j)

{
  gdouble len, *y1, *y2, *y3, *y4 ;
  
  if ( i <= 0 || i >= llt_lattice_time_steps(l) ) {
    g_error("%s: time index (%d) out of range (1,%d)",
	    __FUNCTION__, i, llt_lattice_time_steps(l)-1) ;
  }

  if ( j < 0 || j >= llt_lattice_point_number(l)-1 ) {
    g_error("%s: point index (%d) out of range (0,%d)",
	    __FUNCTION__, j, llt_lattice_point_number(l)-1) ;
  }

  y1 = llt_lattice_point(l,i+0,j+0) ;
  y2 = llt_lattice_point(l,i+0,j+1) ;
  y3 = llt_lattice_point(l,i-1,j+1) ;
  y4 = llt_lattice_point(l,i-1,j+0) ;

  len = llt_ring_length(y1,y2,y3,y4) ;
  
  return len ;
}

gint llt_lattice_truncate(llt_lattice_t *l, gint nt)

{
  gint n ;
  
  if ( nt >= llt_lattice_time_steps(l) ) {
    /*wipe the whole lattice and return*/
    llt_lattice_time_steps(l) = 0 ;

    return 0 ;
  }

  /*move the lattice data to remove nt time steps*/
  n = llt_lattice_time_steps(l) - nt ;
  memmove(l->x,
	  &(l->x[nt*LLT_LATTICE_POINT_SIZE*llt_lattice_point_number(l)]),
	  n*LLT_LATTICE_POINT_SIZE*llt_lattice_point_number(l)*
	  sizeof(gdouble)) ;
  
  llt_lattice_time_steps(l) -= nt ;
  
  return 0 ;
}

/**
 * @}
 */
