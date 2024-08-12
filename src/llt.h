/* This file is part of LLT, a library for Lifting Line Theory
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

#ifndef __LLT_H_INCLUDED__
#define __LLT_H_INCLUDED__

#include <glib.h>


#define LLT_SECTION_DATA_SIZE  4
#define LLT_SECTION_DATA_RE    0
#define LLT_SECTION_DATA_ALPHA 1
#define LLT_SECTION_DATA_CL    2
#define LLT_SECTION_DATA_CD    3

typedef enum {
  LLT_INTERPOLATE_CONSTANT = 0,
  LLT_INTERPOLATE_LINEAR = 1
} llt_interpolate_t ;

typedef enum {
  LLT_EXTRAPOLATE_CONSTANT = 0,
  LLT_EXTRAPOLATE_EXTEND = 1
} llt_extrapolate_t ;

typedef struct _llt_section_t llt_section_t ;
struct _llt_section_t {
  gchar *name ;
  gint *iRe, nRe, nRemax, npts, nptsmax ;
  gdouble *data ;
  llt_interpolate_t interp_Re, interp_al ;
  llt_extrapolate_t extrap_Re_low, extrap_Re_high,
    extrap_al_low, extrap_al_high ;
} ;  

#define llt_section_name(_s)              ((_s)->name)
#define llt_section_Re_index(_s,_i)       ((_s)->iRe[(_i)])
#define llt_section_Re_number(_s)         ((_s)->nRe)
#define llt_section_Re_number_max(_s)     ((_s)->nRemax)
#define llt_section_point_number(_s)      ((_s)->npts)
#define llt_section_point_number_max(_s)  ((_s)->nptsmax)
#define llt_section_interp_Re(_s)         ((_s)->interp_Re)
#define llt_section_interp_alpha(_s)      ((_s)->interp_al)
#define llt_section_extrap_Re_low(_s)     ((_s)->extrap_Re_low)
#define llt_section_extrap_alpha_low(_s)  ((_s)->extrap_al_low)
#define llt_section_extrap_Re_high(_s)    ((_s)->extrap_Re_high)
#define llt_section_extrap_alpha_high(_s) ((_s)->extrap_al_high)
#define llt_section_Re(_s,_i)					\
  ((_s)->data[(_i)*LLT_SECTION_DATA_SIZE+LLT_SECTION_DATA_RE])
#define llt_section_alpha(_s,_i)					\
  ((_s)->data[(_i)*LLT_SECTION_DATA_SIZE+LLT_SECTION_DATA_ALPHA])
#define llt_section_cl(_s,_i)					\
  ((_s)->data[(_i)*LLT_SECTION_DATA_SIZE+LLT_SECTION_DATA_CL])
#define llt_section_cd(_s,_i)					\
  ((_s)->data[(_i)*LLT_SECTION_DATA_SIZE+LLT_SECTION_DATA_CD])

#define LLT_LIFTING_LINE_POINT_SIZE   10
#define LLT_LIFTING_LINE_X             0
#define LLT_LIFTING_LINE_NORMAL        3
#define LLT_LIFTING_LINE_UNIT          6
#define LLT_LIFTING_LINE_CHORD         9
#define LLT_LIFTING_LINE_BOUND_SIZE    6
#define LLT_LIFTING_LINE_BOUND_X       0
#define LLT_LIFTING_LINE_WAKE_X        3

typedef struct _llt_lifting_line_t llt_lifting_line_t ;
struct  _llt_lifting_line_t {
  gint n, nmax ;
  gdouble *xc, *xs, *G ;
  llt_section_t **s ;
} ;

#define llt_lifting_line_point_number(_L)  ((_L)->n)
#define llt_lifting_line_point_number_max(_L)  ((_L)->nmax)
#define llt_lifting_line_point_section(_L,_i) ((_L)->s[(_i)])
#define llt_lifting_line_point(_L,_i)				\
  (&((_L)->xc[(_i)*LLT_LIFTING_LINE_POINT_SIZE+LLT_LIFTING_LINE_X]))
#define llt_lifting_line_normal(_L,_i)					\
  (&((_L)->xc[(_i)*LLT_LIFTING_LINE_POINT_SIZE+LLT_LIFTING_LINE_NORMAL]))
#define llt_lifting_line_unit_chord(_L,_i)				\
  (&((_L)->xc[(_i)*LLT_LIFTING_LINE_POINT_SIZE+LLT_LIFTING_LINE_UNIT]))
#define llt_lifting_line_chord(_L,_i)				\
  ((_L)->xc[(_i)*LLT_LIFTING_LINE_POINT_SIZE+LLT_LIFTING_LINE_CHORD])
#define llt_lifting_line_bound_point(_L,_i)		\
  (&((_L)->xs[(_i)*LLT_LIFTING_LINE_BOUND_SIZE+LLT_LIFTING_LINE_BOUND_X]))
#define llt_lifting_line_wake_point(_L,_i)		\
  (&((_L)->xs[(_i)*LLT_LIFTING_LINE_BOUND_SIZE+LLT_LIFTING_LINE_WAKE_X]))

#define LLT_LATTICE_POINT_SIZE 4

typedef struct _llt_lattice_t llt_lattice_t ;
struct _llt_lattice_t {
  gint npmax, nt, np ;
  gdouble *x ;
} ;

#define llt_lattice_point_number_max(_l) ((_l)->npmax)
#define llt_lattice_point_number(_l)     ((_l)->np)
#define llt_lattice_time_steps(_l)       ((_l)->nt)
#define llt_lattice_index(_l,_i,_j)	 ((_i)*(_l)->np+(_j))
#define llt_lattice_point(_l,_i,_j)				\
  (&((_l)->x[LLT_LATTICE_POINT_SIZE*(llt_lattice_index((_l),(_i),(_j)))+0]))
#define llt_lattice_gamma(_l,_i,_j)				\
  ((_l)->x[LLT_LATTICE_POINT_SIZE*(llt_lattice_index((_l),(_i),(_j)))+3])
  /* ((_l)->x[LLT_LATTICE_POINT_SIZE*((_i)*(_l)->np+(_j))+3]) */
#define llt_lattice_time_shed(_l,_i)		\
  llt_lattice_gamma((_l),(_i),(llt_lattice_point_number((_l))-1))

#define LLT_ASSEMBLY_LINE_NUMBER_MAX 32

typedef struct _llt_assembly_t llt_assembly_t ;
struct _llt_assembly_t {
  gint nl, nt, indices[LLT_ASSEMBLY_LINE_NUMBER_MAX+1] ;
  llt_lifting_line_t *current[LLT_ASSEMBLY_LINE_NUMBER_MAX],
    *lines[LLT_ASSEMBLY_LINE_NUMBER_MAX] ;
  llt_lattice_t *wakes[LLT_ASSEMBLY_LINE_NUMBER_MAX] ;
  amc_transform_chain_t *chains[LLT_ASSEMBLY_LINE_NUMBER_MAX] ;
  amc_transform_t *transforms[LLT_ASSEMBLY_LINE_NUMBER_MAX] ;
  gdouble reg, *G, *Gu, *dn, *ue, *vs, *ui, *U, *al, *cl, *cd, *Re ;
  gdouble *uw ;
} ;

#define llt_assembly_time_steps(_a) ((_a)->nt)
#define llt_assembly_lifting_line_number(_a) ((_a)->nl)
#define llt_assembly_regularisation(_a)      ((_a)->reg)
#define llt_assembly_lifting_line(_a,_i)     ((_a)->lines[(_i)])
#define llt_assembly_lifting_line_current(_a,_i) ((_a)->current[(_i)])
#define llt_assembly_wake(_a,_i)             ((_a)->wakes[(_i)])
#define llt_assembly_transform_chain(_a,_i)  ((_a)->chains[(_i)])
#define llt_assembly_transform(_a,_i)        ((_a)->transforms[(_i)])
#define llt_assembly_point_number(_a)		\
  ((_a)->indices[llt_assembly_lifting_line_number((_a))])

typedef enum {
  LLT_SOLVER_UNDEFINED = 0,
  LLT_SOLVER_RELAXATION = 1
} llt_solver_method_t ;

typedef struct _llt_solver_t llt_solver_t ;

struct _llt_solver_t {
  llt_solver_method_t method ;
  gdouble s, tol, err, nu ;
  gint iter, iter_max ;
} ;

#define llt_solver_method(_s)               ((_s)->method)
#define llt_solver_relaxation_factor(_s)    ((_s)->s)
#define llt_solver_tolerance(_s)            ((_s)->tol)
#define llt_solver_error(_s)                ((_s)->err)
#define llt_solver_iteration_number(_s)     ((_s)->iter)
#define llt_solver_iteration_number_max(_s) ((_s)->iter_max)
#define llt_solver_viscosity(_s)            ((_s)->nu)

typedef enum {
  LLT_SECTION_SOLVER_UNDEFINED = 0,
  LLT_SECTION_SOLVER_XFOIL     = 1
} llt_section_solver_t ;

typedef struct _llt_script_t llt_script_t ;

struct _llt_script_t {
  llt_section_solver_t type ;
  gchar *name, *exec, *path, *dcty ;
  gint osize, ia, icl, icd ;
} ;

#define llt_script_name(_s) ((_s)->name)
#define llt_script_exec(_s) ((_s)->exec)
#define llt_script_path(_s) ((_s)->path)
#define llt_script_dcty(_s) ((_s)->dcty)
#define llt_script_type(_s) ((_s)->type)
#define llt_script_output_size(_s)        ((_s)->osize)
#define llt_script_output_index_alpha(_s) ((_s)->ia)
#define llt_script_output_index_cl(_s)    ((_s)->icl)
#define llt_script_output_index_cd(_s)    ((_s)->icd)

llt_lifting_line_t *llt_lifting_line_alloc(gint n) ;
llt_lattice_t *llt_lattice_alloc(gint np, gint nt) ;
gint llt_lifting_line_move(llt_lifting_line_t *L,
			   amc_transform_t *T,
			   llt_lifting_line_t *M,
			   gdouble *u, gint ustr,
			   gdouble *n, gint nstr) ;
gint llt_lifting_line_wing_rectangular(llt_lifting_line_t *L,
				       llt_section_t *s,
				       gdouble y0, gdouble y1,
				       gdouble ch0, gdouble ch1,
				       gdouble tw0, gdouble tw1,
				       gint n) ;
gint llt_lifting_line_wing_elliptical(llt_lifting_line_t *L,
				      llt_section_t *s,
				      gdouble a, gdouble b, gint n) ;
gint llt_lifting_line_wing_helical(llt_lifting_line_t *L,
				   llt_section_t *s,
				   gdouble r0, gdouble r1,
				   gdouble ch,
				   gdouble U, gdouble Om,
				   gint n) ;

gint llt_lifting_line_velocity_matrix(llt_lifting_line_t *L, gdouble d,
				      gdouble *xcp, gint xstr, gint nx,
				      gdouble *A, gint lda) ;
gint llt_lifting_line_velocity(llt_lifting_line_t *L, gdouble *G, gint gstr,
			       gdouble d, gdouble *x, gdouble al, gdouble bt,
			       gdouble *u) ;
gint llt_lifting_line_velocity_lattice(llt_lifting_line_t *L,
				       gdouble *G, gint gstr,
				       llt_lattice_t *T,
				       gdouble d, gdouble al, gdouble bt,
				       gdouble *u, gint ustr) ;

gint llt_lifting_line_points_write(FILE *f, llt_lifting_line_t *L) ;
gint llt_lifting_line_geometry_write(FILE *f, llt_lifting_line_t *L) ;
gint llt_lattice_write(FILE *f, llt_lattice_t *W) ;
gint llt_segment_velocity(gdouble *x, gdouble *y1, gdouble *y2, gdouble d,
			  gdouble *u) ;
gint llt_assembly_lifting_line_add(llt_assembly_t *A, llt_lifting_line_t *L,
				   amc_transform_chain_t *C) ;
gint llt_lattice_velocity(llt_lattice_t *l, gdouble d, gdouble *x,
			  gdouble al, gdouble bt, gdouble *u) ;
gint llt_lattice_velocity_lattice(llt_lattice_t *T, llt_lattice_t *S,
				  gdouble d, gdouble al, gdouble bt,
				  gdouble *u, gint ustr) ;
gdouble llt_lattice_ring_length(llt_lattice_t *l, gint i, gint j) ;
gint llt_lattice_truncate(llt_lattice_t *l, gint nt) ;

gint llt_lifting_line_copy(llt_lifting_line_t *M, llt_lifting_line_t *L) ;
gint llt_lifting_line_read(FILE *f, llt_lifting_line_t *L) ;
gint llt_lifting_line_write(FILE *f, llt_lifting_line_t *L) ;

llt_assembly_t *llt_assembly_alloc(void) ;
gint llt_assembly_initialise(llt_assembly_t *A, gdouble t) ;
gint llt_assembly_force(llt_assembly_t *A, gdouble *F) ;
gint llt_assembly_force_moment(llt_assembly_t *A, gdouble *F,
			       gdouble *x, gdouble *M) ;
gint llt_assembly_lifting_line_force(llt_assembly_t *A, gint i, gdouble *F) ;
gint llt_assembly_lifting_line_force_moment(llt_assembly_t *A, gint i,
					    gdouble *F,
					    gdouble *x, gdouble *M) ;

gint llt_lattice_velocity_lifting_line(llt_lattice_t *W,
				       llt_lifting_line_t *L,
				       gdouble d,
				       gdouble *u, gint ustr) ;
gint llt_lattice_initialise(llt_lattice_t *W, llt_lifting_line_t *L,
			    gdouble t, double G0) ;
gint llt_lattice_wake_shed(llt_lattice_t *W, llt_lifting_line_t *L, gdouble t) ;
gint llt_lifting_line_segment_unit(gdouble *ds, llt_lifting_line_t *L, gint i,
				   gdouble *len) ;
gdouble llt_lifting_line_segment_length(llt_lifting_line_t *L, gint i) ;

gint llt_solver_step(llt_assembly_t *A, llt_solver_t *s, gdouble t,
		     gdouble dt, gboolean solve) ;
gint llt_solver_wakes_update(llt_assembly_t *A, gdouble dt) ;

llt_section_t *llt_section_alloc(gint nRe, gint npts) ;
gint llt_section_write(FILE *f, llt_section_t *s) ;
gint llt_section_read(FILE *f, llt_section_t *s) ;
llt_section_t *llt_section_read_alloc(FILE *f, gint m) ;
gint llt_section_interp(llt_section_t *s, gdouble Re, gdouble al,
			gdouble *cl, gdouble *cd) ;
gint llt_section_linear(llt_section_t *s, gdouble a) ;
gint llt_section_entry_append(llt_section_t *s,
			      gdouble Re, gdouble al,
			      gdouble cl, gdouble cd) ;
gint llt_section_sort(llt_section_t *s) ;
gint llt_section_index(llt_section_t *s) ;

llt_script_t *llt_script_alloc(gchar *name) ;
gchar *llt_script_locate(gchar *script, gchar *evar) ;
gint llt_script_init(llt_script_t *s, llt_section_solver_t type,
		     gchar *dcty, gboolean dcty_force) ;
gint llt_script_section_make(llt_script_t *script, gchar *aerofoil,
			     gchar *polar,
			     gdouble Re, gdouble amin, gdouble amax,
			     gdouble da) ;
gint llt_script_section_add(llt_script_t *script, gchar *aerofoil,
			    gdouble Re, gdouble amin, gdouble amax,
			    gdouble da, llt_section_t *s) ;

gdouble llt_naca_four(gdouble t, gdouble p, gdouble m, gdouble x) ;

#endif /*__LLT_H_INCLUDED__*/
