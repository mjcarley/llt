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
 * Allocate a new lifting line
 * 
 * @param n maximum number of control points (or segments) on line.
 * 
 * @return newly allocated ::llt_lifting_line_t
 */

llt_lifting_line_t *llt_lifting_line_alloc(gint n)
  
{
  llt_lifting_line_t *l ;

  l = (llt_lifting_line_t *)g_malloc0(sizeof(llt_lifting_line_t)) ;

  /*control points where velocity etc is evaluated*/
  l->xc = (gdouble *)g_malloc0(n*LLT_LIFTING_LINE_POINT_SIZE*sizeof(gdouble)) ;
  /*segment endpoints and wake shedding points*/
  l->xs = (gdouble *)g_malloc0((n+1)*LLT_LIFTING_LINE_BOUND_SIZE
			       *sizeof(gdouble)) ;
  l->n = l->nmax = n ;
  l->s = (llt_section_t **)g_malloc0(n*sizeof(llt_section_t *)) ;
  
  return l ;
}

/** 
 * Deep copy one lifting line into another \f$M:=L\f$
 * 
 * @param M on exit contains copy of all data from \a L;
 * @param L lifting line to be copied.
 * 
 * @return 0 on success
 */

gint llt_lifting_line_copy(llt_lifting_line_t *M, llt_lifting_line_t *L)

{
  gint i ;
  gdouble *xl, *xm ;
  
  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    xl = llt_lifting_line_point(L,i) ;
    xm = llt_lifting_line_point(M,i) ;
    xm[0] = xl[0] ; xm[1] = xl[1] ; xm[2] = xl[2] ;

    xl = llt_lifting_line_normal(L,i) ;
    xm = llt_lifting_line_normal(M,i) ;
    xm[0] = xl[0] ; xm[1] = xl[1] ; xm[2] = xl[2] ;

    xl = llt_lifting_line_unit_chord(L,i) ;
    xm = llt_lifting_line_unit_chord(M,i) ;
    xm[0] = xl[0] ; xm[1] = xl[1] ; xm[2] = xl[2] ;

    xl = llt_lifting_line_bound_point(L,i) ;
    xm = llt_lifting_line_bound_point(M,i) ;
    xm[0] = xl[0] ; xm[1] = xl[1] ; xm[2] = xl[2] ;

    llt_lifting_line_point_section(M,i) =
      llt_lifting_line_point_section(L,i) ;
  }

  i = llt_lifting_line_point_number(L) ;
  xl = llt_lifting_line_bound_point(L,i) ;
  xm = llt_lifting_line_bound_point(M,i) ;
  xm[0] = xl[0] ; xm[1] = xl[1] ; xm[2] = xl[2] ;
  
  return 0 ;
}  

/** 
 * Apply motion transformation to lifting line, evaluating velocities if
 * required
 * 
 * @param L basic lifting line in reference position; 
 * @param T transform to be applied to \a L;
 * @param M on exit contains lifting line in transformed position;
 * @param u if not NULL, on exit contains velocity of control points of \a M; 
 * @param ustr stride between successive velocity entries in \a u;
 * @param n if not NULL, on exit contains derivatives of normals of \a M; 
 * @param nstr stride between successive entries in \a n.
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_move(llt_lifting_line_t *L,
			   amc_transform_t *T,
			   llt_lifting_line_t *M,
			   gdouble *u, gint ustr,
			   gdouble *n, gint nstr)

{
  gint i ;

  if ( llt_lifting_line_point_number_max(M) >
       llt_lifting_line_point_number(L) ) {
    g_error("%s: not enough elements for %d points", __FUNCTION__,
	    llt_lifting_line_point_number(L)) ;
  }

  if ( T == NULL ) {
    llt_lifting_line_copy(M, L) ;
    if ( u == NULL ) return 0 ;
    for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
      u[i*ustr+0] = u[i*ustr+1] = u[i*ustr+2] = 0 ;
    }
    
    return 0 ;
  }
  
  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    amc_transform_matrix_apply(T, 0,
			       llt_lifting_line_point(L,i),
			       llt_lifting_line_point(M,i)) ;
    amc_transform_matrix_apply(T, 0,
			       llt_lifting_line_bound_point(L,i),
			       llt_lifting_line_bound_point(M,i)) ;    
    amc_transform_matrix_apply(T, 0,
			       llt_lifting_line_wake_point(L,i),
			       llt_lifting_line_wake_point(M,i)) ;
    amc_transform_matrix_apply_vector(T, 0,
				      /* llt_lifting_line_point(L,i), */
				      llt_lifting_line_normal(L,i),
				      llt_lifting_line_normal(M,i)) ;
    amc_transform_matrix_apply_vector(T, 0,
				      /* llt_lifting_line_point(L,i), */
				      llt_lifting_line_unit_chord(L,i),
				      llt_lifting_line_unit_chord(M,i)) ;
    llt_lifting_line_chord(M,i) = llt_lifting_line_chord(L,i) ;
    llt_lifting_line_point_section(M,i) =
      llt_lifting_line_point_section(L,i) ;
  }
  i = llt_lifting_line_point_number(L) ;
  
  amc_transform_matrix_apply(T, 0,
			     llt_lifting_line_bound_point(L,i),
			     llt_lifting_line_bound_point(M,i)) ;    
  amc_transform_matrix_apply(T, 0,
			     llt_lifting_line_wake_point(L,i),
			     llt_lifting_line_wake_point(M,i)) ;    

  llt_lifting_line_point_number(M) = llt_lifting_line_point_number(L) ;

  if ( u == NULL ) return 0 ;

  if ( ustr < 3 ) {
    g_error("%s: ustr (%d) must be at least 3", __FUNCTION__, ustr) ;
  }
  
  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    amc_transform_matrix_apply(T, 1,
			       llt_lifting_line_point(L,i),
			       &(u[i*ustr])) ;
  }

  if ( n == NULL ) return 0 ;
  
  if ( nstr < 3 ) {
    g_error("%s: nstr (%d) must be at least 3", __FUNCTION__, nstr) ;
  }

  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    amc_transform_matrix_apply_vector(T, 1,
				      /* llt_lifting_line_point(L,i), */
				      llt_lifting_line_normal(L,i),
				      &(n[i*nstr])) ;
  }

  return 0 ;
}

/** 
 * Write control points of a lifting line to file.
 * 
 * @param f output file stream;
 * @param L lifting line to be written to file.
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_points_write(FILE *f, llt_lifting_line_t *L)

{
  gint i ;
  gdouble *x ;

  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_point(L,i) ;
    fprintf(f, "%e %e %e ", x[0], x[1], x[2]) ;
  }
  fprintf(f, "\n") ;
  
  return 0 ;
}

/** 
 * Write bound and wake segment endpoints of lifting line to file
 * 
 * @param f output file stream;
 * @param L lifting line to write.
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_geometry_write(FILE *f, llt_lifting_line_t *L)

{
  gint i ;
  gdouble *x ;

  for ( i = 0 ; i <= llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_bound_point(L,i) ;
    fprintf(f, "%e %e %e ", x[0], x[1], x[2]) ;
  }
  fprintf(f, "\n") ;
  for ( i = 0 ; i <= llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_wake_point(L,i) ;
    fprintf(f, "%e %e %e ", x[0], x[1], x[2]) ;
  }
  fprintf(f, "\n") ;
  
  return 0 ;
}

/** 
 * Evaluate velocity induced by vortex segment of unit circulation,
 * using equation 9 of van Garrel,
 * https://publications.ecn.nl/WIN/2003/ECN-C--03-079
 * 
 * @param x field point where velocity is to be evaluated;
 * @param y1 first end point of segment;
 * @param y2 second end point of segment;
 * @param d regularisation parameter;
 * @param u on exit contains velocity at \a x.
 * 
 * @return 0 on success.
 */

gint llt_segment_velocity(gdouble *x, gdouble *y1, gdouble *y2, gdouble d,
			  gdouble *u)

{
  gdouble r1[3], r2[3], dr[3], R1, R2, sc, len2 ;

  llt_vector_diff(dr,y2,y1) ;
  len2 = llt_vector_scalar(dr,dr) ;

  /* if ( len2 < 1e-9 ) { u[0] = u[1] = u[2] = 0.0 ; return 0 ; } */
  
  llt_vector_diff(r1,x,y1) ;
  llt_vector_diff(r2,x,y2) ;

  R1 = llt_vector_length(r1) ;
  R2 = llt_vector_length(r2) ;

  llt_vector_cross(u,r1,r2) ;

  sc = R1*R2*(R1*R2 + llt_vector_scalar(r1,r2)) + d*d*len2 ;
  sc = (R1 + R2)/sc/4.0/M_PI ;

  u[0] *= sc ; u[1] *= sc ; u[2] *= sc ; 
  
  return 0 ;
}

/** 
 * Find unit vector aligned with bound segment of a lifting line
 * 
 * @param ds on exit contains unit vector aligned with segment \a i of \a L;
 * @param L lifting line containing bound segment;
 * @param i index of segment whose unit vector is to be evaluated;
 * @param len if not NULL, contains length of segment on exit.
 *
 * @return 0 on success.
 */

gint llt_lifting_line_segment_unit(gdouble *ds, llt_lifting_line_t *L, gint i,
				   gdouble *len)

{
  gdouble s ;
  
  if ( i >= llt_lifting_line_point_number(L) ) {
    g_error("%s: point index (%d) out of range on %d point lifting line",
	    __FUNCTION__, i, llt_lifting_line_point_number(L)) ;
  }

  llt_vector_diff(ds,
		  llt_lifting_line_bound_point(L,i+1),
		  llt_lifting_line_bound_point(L,i+0)) ;
  s = llt_vector_length(ds) ;

  ds[0] /= s ; ds[1] /= s ; ds[2] /= s ;

  if ( len != NULL ) *len = s ;
  
  return 0 ;
}

/** 
 * Evaluate matrix for velocity induced by lifting line at arbitrary
 * points (in practice control points on the same, or another, lifting
 * line). After generation, velocities at field points can be
 * evaluated as \f$\mathbf{u}=[A]\Gamma\f$, where \f$\Gamma\f$ is the
 * vector of circulations on segments of \a L.
 * 
 * @param L lifting line inducing velocity;
 * @param d regularisation parameter;
 * @param xcp array of control points for velocity evaluation;
 * @param xstr stride between successive control points;
 * @param ncp number of control points;
 * @param A allocated array for velocity matrix;
 * @param lda leading dimension of \a A (used when generating matrix for
 * multiple lifting lines).
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_velocity_matrix(llt_lifting_line_t *L, gdouble d,
				      gdouble *xcp, gint xstr, gint ncp,
				      gdouble *A, gint lda)

{
  gint i, j, off0, off1, off2 ;
  gdouble *x, *y1, *y2, du[3] ;
  
  for ( i = 0 ; i < ncp ; i ++ ) {
    x = &(xcp[i*xstr]) ;
    /*offsets into appropriate rows of A*/
    off0 = 3*i*lda + 0*lda ;
    off1 = 3*i*lda + 1*lda ;
    off2 = 3*i*lda + 2*lda ;
    for ( j = 0 ; j < llt_lifting_line_point_number(L) ; j ++ ) {
      y1 = llt_lifting_line_bound_point(L,j+0) ;
      y2 = llt_lifting_line_bound_point(L,j+1) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      A[off0+j] += du[0] ; A[off1+j] += du[1] ; A[off2+j] += du[2] ; 
      y1 = llt_lifting_line_bound_point(L,j+1) ;
      y2 = llt_lifting_line_wake_point(L,j+1) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      A[off0+j] += du[0] ; A[off1+j] += du[1] ; A[off2+j] += du[2] ; 
      y1 = llt_lifting_line_wake_point(L,j+1) ;
      y2 = llt_lifting_line_wake_point(L,j+0) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      A[off0+j] += du[0] ; A[off1+j] += du[1] ; A[off2+j] += du[2] ; 
      y1 = llt_lifting_line_wake_point(L,j+0) ;
      y2 = llt_lifting_line_bound_point(L,j+0) ;
      llt_segment_velocity(x, y1, y2, d, du) ;
      A[off0+j] += du[0] ; A[off1+j] += du[1] ; A[off2+j] += du[2] ; 
    }
  }
  
  return 0 ;
}

/** 
 * Write data for a lifting line to file
 * 
 * @param f output file stream;
 * @param L lifting line to write.
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_write(FILE *f, llt_lifting_line_t *L)

{
  gint i, j ;
  gdouble *x ;
  
  fprintf(f, "LINE %d\n", llt_lifting_line_point_number(L)) ;

  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_point(L,i) ;
    fprintf(f, "%e", x[0]) ;
    for ( j = 1 ; j < LLT_LIFTING_LINE_POINT_SIZE ; j ++ ) {
      fprintf(f, " %e", x[j]) ;
    }
    fprintf(f, "\n") ;
  }

  for ( i = 0 ; i <= llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_bound_point(L,i) ;
    fprintf(f, "%e", x[0]) ;
    for ( j = 1 ; j < LLT_LIFTING_LINE_BOUND_SIZE ; j ++ ) {
      fprintf(f, " %e", x[j]) ;
    }
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

/** 
 * Read a lifting line written to file with ::llt_lifting_line_write
 * 
 * @param f input file stream;
 * @param L on exit contains lifting line read from \a f; a check is 
 * performed to ensure that \a L has enough space allocated to contain
 * the input data.
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_read(FILE *f, llt_lifting_line_t *L)

{
  gint n, i, j ;
  gdouble *x ;

  fscanf(f, "%*[^ ]s") ;
  fscanf(f, "%d", &n) ;

  if ( llt_lifting_line_point_number_max(L) < n ) {
    g_error("%s: not enough space (%d points) for %d elements",
	    __FUNCTION__, llt_lifting_line_point_number_max(L), n) ;
  }

  llt_lifting_line_point_number(L) = n ;

  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_point(L,i) ;
    for ( j = 0 ; j < LLT_LIFTING_LINE_POINT_SIZE ; j ++ ) {
      fscanf(f, "%lg", &(x[j])) ;
    }
  }

  for ( i = 0 ; i <= llt_lifting_line_point_number(L) ; i ++ ) {
    x = llt_lifting_line_bound_point(L,i) ;
    for ( j = 0 ; j < LLT_LIFTING_LINE_BOUND_SIZE ; j ++ ) {
      fscanf(f, "%lg", &(x[j])) ;
    }
  }
  
  return 0 ;
}

/** 
 * Find the length of a lifting line segment
 * 
 * @param L an ::llt_lifting_line_t;
 * @param i index of segment of \a L.
 * 
 * @return length of segment \a i of \a L.
 */

gdouble llt_lifting_line_segment_length(llt_lifting_line_t *L, gint i)

{
  gdouble s, ds[3] ;

  if ( i >= llt_lifting_line_point_number(L) ) {
    g_error("%s: point index (%d) out of range on %d point lifting line",
	    __FUNCTION__, i, llt_lifting_line_point_number(L)) ;
  }

  llt_vector_diff(ds,
		  llt_lifting_line_bound_point(L,i+1),
		  llt_lifting_line_bound_point(L,i+0)) ;
  s = llt_vector_length(ds) ;

  return s ;
}

/** 
 * Evaluate velocity "induced" by a lifting line
 * 
 * @param L a lifting line;
 * @param G circulation at each segment of \a L;
 * @param gstr stride between elements of \a G;
 * @param d regularisation parameter for velocity evaluation;
 * @param x evaluation point;
 * @param al \f$\alpha\f$;
 * @param bt \f$\beta\f$;
 * @param u on entry contains \f$\mathbf{u}_{0}\f$, on exit contains
 * \f$\beta \mathbf{u}_{0} + \alpha\mathbf{u}_{l}\f$ with
 * \f$\mathbf{u}_{l}\f$ velocity generated by lifting line \a L.
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_velocity(llt_lifting_line_t *L, gdouble *G, gint gstr,
			       gdouble d, gdouble *x, gdouble al, gdouble bt,
			       gdouble *u)
{
  gint i ;
  gdouble du[3], *y1, *y2 ;
  
  u[0] *= bt ; u[1] *= bt ; u[2] *= bt ;

  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    gdouble ut[3] = {0} ;
    y1 = llt_lifting_line_bound_point(L,i+0) ;
    y2 = llt_lifting_line_bound_point(L,i+1) ;
    llt_segment_velocity(x, y1, y2, d, du) ;
    ut[0] += du[0] ; ut[1] += du[1] ; ut[2] += du[2] ; 

    y1 = llt_lifting_line_bound_point(L,i+1) ;
    y2 = llt_lifting_line_wake_point(L,i+1) ;
    llt_segment_velocity(x, y1, y2, d, du) ;
    ut[0] += du[0] ; ut[1] += du[1] ; ut[2] += du[2] ; 

    y1 = llt_lifting_line_wake_point(L,i+1) ;
    y2 = llt_lifting_line_wake_point(L,i+0) ;
    llt_segment_velocity(x, y1, y2, d, du) ;
    ut[0] += du[0] ; ut[1] += du[1] ; ut[2] += du[2] ; 

    y1 = llt_lifting_line_wake_point(L,i+0) ;
    y2 = llt_lifting_line_bound_point(L,i+0) ;
    llt_segment_velocity(x, y1, y2, d, du) ;
    ut[0] += du[0] ; ut[1] += du[1] ; ut[2] += du[2] ; 

    u[0] += al*G[i*gstr]*ut[0] ;
    u[1] += al*G[i*gstr]*ut[1] ;
    u[2] += al*G[i*gstr]*ut[2] ; 
  }
  
  return 0 ;
}

/** 
 * Evaluate velocity "induced" by a lifting line on a lattice
 * 
 * @param L lifting line;
 * @param G circulation at each segment of \a L;
 * @param gstr stride between elements of \a G;
 * @param T target lattice;
 * @param d regularisation parameter;
 * @param al \f$\alpha\f$;
 * @param bt \f$\beta\f$;
 * @param u on entry contains \f$\mathbf{u}_{0}\f$, on exit contains
 * \f$\beta \mathbf{u}_{0} + \alpha\mathbf{u}_{S}\f$ with \f$\mathbf{u}_{S}\f$
 * velocity generated by lattice \a S;
 * @param ustr stride in \a u.
 * 
 * @return 0 on success.
 */

gint llt_lifting_line_velocity_lattice(llt_lifting_line_t *L,
				       gdouble *G, gint gstr,
				       llt_lattice_t *T,
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
      llt_lifting_line_velocity(L, G, gstr, d, x, al, 1.0, &(u[idx*ustr])) ;
    }
  }

  return 0 ;
}
