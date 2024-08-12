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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <amc.h>

#include "llt.h"
#include "llt-private.h"

gchar *progname ;

static FILE *file_open(gchar *file, gchar *mode,
		       gchar *file_default, FILE *fdefault)

{
  FILE *f ;

  if ( file_default != NULL ) {
    if ( strcmp(file, file_default) == 0 ) return fdefault ;
  }
  
  f = fopen(file, mode) ;

  if ( f == NULL ) {
    fprintf(stderr, "%s: cannot open file \"%s\"", progname, file) ;
    exit(1) ;
  }
  
  return f ;
}

static void file_close(FILE *f)

{
  if ( f == stdin || f == stdout || f == stderr ) return ;

  fclose(f) ;
  
  return ;
}

static gint parse_test(gchar *str)

{
  gchar *tests[] = {"wake", "solve", "section", "script", NULL} ;
  gint i ;

  for ( i = 0 ; tests[i] != NULL ; i ++ ) {
    if ( strcmp(str, tests[i]) == 0 ) return i ;
  }
  
  return -1 ;
}

static void wake_test(llt_assembly_t *W, gdouble t0, gdouble t1, gint nt)

{
  gint i ;
  FILE *f ;
  gchar fname[256] ;
  llt_solver_t solver ;
  
  llt_assembly_time_steps(W) = nt+1 ;
  llt_assembly_initialise(W, t0) ;

  for ( i = 1 ; i < nt ; i ++ ) {
    llt_solver_step(W, &solver, t0 + (t1 - t0)*i/(nt-1), 0, FALSE) ;
  }

  for ( i = 0 ; i < llt_assembly_lifting_line_number(W) ; i ++ ) {
    sprintf(fname, "line-%d.dat", i) ;
    f = fopen(fname, "w") ;
    llt_lifting_line_geometry_write(f, llt_assembly_lifting_line_current(W,i)) ;
    fclose(f) ;
    if ( llt_assembly_wake(W,i) != NULL ) {
      sprintf(fname, "wake-%d.dat", i) ;
      f = fopen(fname, "w") ;
      llt_lattice_write(f, llt_assembly_wake(W,i)) ;
      fclose(f) ;
    }
  }
  
  return ;
}

static void solver_test(llt_assembly_t *W, gdouble t0, gdouble t1, gint nt)

{
  gint i, j, idx ;
  FILE *f ;
  gchar fname[256] ;
  gdouble t, *x ;
  llt_lifting_line_t *L ;
  llt_solver_t solver ;
  
  llt_assembly_time_steps(W) = nt+1 ;
  llt_assembly_initialise(W, t0) ;
  llt_assembly_regularisation(W) = 0.0125 ;

  llt_solver_method(&solver) = LLT_SOLVER_RELAXATION ;
  llt_solver_relaxation_factor(&solver) = 0.0625 ;
  llt_solver_tolerance(&solver) = 1e-3 ;
  llt_solver_iteration_number_max(&solver) = 40 ;
  llt_solver_viscosity(&solver) = 1.5e-5 ;
  
  for ( i = 1 ; i < nt ; i ++ ) {
    t = t0 + (t1 - t0)*i/(nt-1) ;
    /* fprintf(stderr, "time step %d: t = %lg\n", i, t) ; */
    llt_solver_step(W, &solver, t, (t1 - t0)/(nt-1), TRUE) ;
  }

  for ( i = 0 ; i < llt_assembly_lifting_line_number(W) ; i ++ ) {
    sprintf(fname, "line-%d.dat", i) ;
    f = fopen(fname, "w") ;
    llt_lifting_line_geometry_write(f, llt_assembly_lifting_line_current(W,i)) ;
    fclose(f) ;
    if ( llt_assembly_wake(W,i) != NULL ) {
      sprintf(fname, "wake-%d.dat", i) ;
      f = fopen(fname, "w") ;
      llt_lattice_write(f, llt_assembly_wake(W,i)) ;
      fclose(f) ;
    }
  }

  for ( i = 0 ; i < llt_assembly_lifting_line_number(W) ; i ++ ) {
    L = llt_assembly_lifting_line(W, i) ;
    for ( j = 0 ; j < llt_lifting_line_point_number(L) ; j ++ ) {
      idx = W->indices[i] + j ;
      x = llt_lifting_line_point(L,j) ;
      fprintf(stdout, "%e %e %e %e %e\n",
	      x[0], x[1], x[2], W->cl[idx], W->Re[idx]) ;
    }
  }
  
  return ;
}

static void assembly_wing(llt_assembly_t *A, gint np, gdouble a, gdouble b,
			  gdouble *U)

{
  llt_lifting_line_t *L ;
  llt_section_t *s ;
  gdouble tw0, tw1 ;
  amc_transform_t *T ;
  amc_transform_chain_t *C ;
  
  L = llt_lifting_line_alloc(np) ;
  tw0 = 3.0*M_PI/180 ; tw1 = 3.0*M_PI/180 ;
  /* tw0 = tw1 = 0.0 ; */

  fprintf(stderr, "rectangular wing geometry\n") ;
  fprintf(stderr, "=========================\n") ;
  fprintf(stderr, "span:  b = %lg\n", b) ;
  fprintf(stderr, "chord: a = %lg\n", a) ;

  s = llt_section_alloc(1, 2) ;
  llt_section_linear(s, 2.0*M_PI) ;
  
  llt_lifting_line_wing_rectangular(L, s, -0.5*b, 0.5*b, a, a, tw0, tw1, np) ;

  /*motion definition*/
  C = amc_transform_chain_alloc(2) ;

  T = amc_transform_alloc(3, 1) ; T->order = 1 ;
  amc_transform_variable_add(T, "U", U) ;
  amc_transform_translation(T, 0, "-U*t", 0, NULL, 0, "-0.1*U*t", 1) ;
  amc_transform_expressions_compile(T) ;

  amc_transform_chain_transform_add(C, T) ;

  llt_assembly_lifting_line_add(A, L, C) ;

  return ;
}

static void assembly_propeller(llt_assembly_t *A, gint np,
			       gdouble a, gdouble b,
			       gdouble *U, gdouble *Om)

{
  llt_lifting_line_t *L ;
  llt_section_t *s ;
  gdouble tw0, tw1 ; 
  amc_transform_t *T, *R, *R2 ;
  amc_transform_chain_t *C1, *C2 ;
  
  L = llt_lifting_line_alloc(np) ;
  tw0 = 13.0*M_PI/180 ; tw1 = 13.0*M_PI/180 ;
  
  s = llt_section_alloc(1, 2) ;
  llt_section_linear(s, 2.0*M_PI) ;

  fprintf(stderr, "%s: making lifting lines for propeller\n", __FUNCTION__) ;
  llt_lifting_line_wing_rectangular(L, s, 0.1*b, b, a, a, tw0, tw1, np) ;
  llt_lifting_line_wing_rectangular(L, s, 0.1*b, b, a, a, tw0, tw1, np) ;

  /*motion definition*/
  C1 = amc_transform_chain_alloc(8) ;
  C2 = amc_transform_chain_alloc(8) ;

  /*forward flight parallel to z axis*/
  T = amc_transform_alloc(3, 1) ; T->order = 1 ;
  amc_transform_variable_add(T, "U", U) ;
  amc_transform_translation(T, 0, NULL, 0, NULL, 0, "U*t", 1) ;
  amc_transform_expressions_compile(T) ;

  /*rotation about z axis at velocity Om*/
  R = amc_transform_alloc(3, 1) ; R->order = 1 ;
  amc_transform_variable_add(R, "Omega", Om) ;
  amc_transform_rotation_z(R, 0, "Omega*t", 1) ;
  amc_transform_expressions_compile(R) ;

  /*rotation about z axis by 180 degrees for second blade*/
  R2 = amc_transform_alloc(3, 1) ; R2->order = 1 ;
  amc_transform_rotation_z(R2, 0, "pi", 1) ;
  amc_transform_expressions_compile(R2) ;
  
  amc_transform_chain_transform_add(C1, R) ;
  amc_transform_chain_transform_add(C1, T) ;

  amc_transform_chain_transform_add(C2, R2) ;
  amc_transform_chain_transform_add(C2, R) ;
  amc_transform_chain_transform_add(C2, T) ;

  fprintf(stderr, "%s: assembling lifting line(s)\n", __FUNCTION__) ;
  llt_assembly_lifting_line_add(A, L, C1) ;
  llt_assembly_lifting_line_add(A, L, C2) ;

  return ;
}

static void assembly_elliptical(llt_assembly_t *A, gint np,
				gdouble a, gdouble b,
				gdouble *U)

{
  llt_lifting_line_t *L ;
  llt_section_t *s ;
  amc_transform_t *T ;
  amc_transform_chain_t *C ;
  
  L = llt_lifting_line_alloc(np) ;
  
  s = llt_section_alloc(1, 2) ;
  llt_section_linear(s, 2.0*M_PI) ;

  fprintf(stderr, "%s: making lifting line for wing\n", __FUNCTION__) ;
  llt_lifting_line_wing_elliptical(L, s, a, b, np) ;

  /*motion definition*/
  C = amc_transform_chain_alloc(2) ;

  T = amc_transform_alloc(3, 1) ; T->order = 1 ;
  amc_transform_variable_add(T, "U", U) ;
  amc_transform_translation(T, 0, "-U*t", 0, NULL, 0, "-0.1*U*t", 1) ;
  amc_transform_expressions_compile(T) ;

  amc_transform_chain_transform_add(C, T) ;

  fprintf(stderr, "%s: assembling lifting line(s)\n", __FUNCTION__) ;
  llt_assembly_lifting_line_add(A, L, C) ;

  return ;
}

static void assembly_van_garrel(llt_assembly_t *A, gint np,
				gdouble a, gdouble b,
				gdouble *U, gdouble *Om)

{
  llt_lifting_line_t *L1 ;
  llt_section_t *s ;
  gdouble lm ; 
  amc_transform_t *T, *R ;
  amc_transform_chain_t *C1 ;

  /*
   * test case from van Garrel 2003, helical rotor blade aligned with
   * flow: lift should be zero at all stations at design speed
   */

  /*velocity ratio*/
  lm = 10.0 ;
  /* *Om = lm*(10)/b ; */
  
  L1 = llt_lifting_line_alloc(np) ;
  
  s = llt_section_alloc(1, 2) ;
  llt_section_linear(s, 2.0*M_PI) ;

  fprintf(stderr, "%s: making lifting lines for propeller\n", __FUNCTION__) ;
  llt_lifting_line_wing_helical(L1, s, 0.2*b, b, a, 10, lm*10/b, np) ;

  /*motion definition*/
  C1 = amc_transform_chain_alloc(8) ;

  /*forward flight parallel to z axis*/
  T = amc_transform_alloc(3, 1) ; T->order = 1 ;
  amc_transform_variable_add(T, "U", U) ;
  amc_transform_translation(T, 0, NULL, 0, NULL, 0, "-U*t", 1) ;
  amc_transform_expressions_compile(T) ;

  /*rotation about z axis at velocity Om*/
  R = amc_transform_alloc(3, 1) ; R->order = 1 ;
  amc_transform_variable_add(R, "Omega", Om) ;
  amc_transform_rotation_z(R, 0, "Omega*t", 1) ;
  amc_transform_expressions_compile(R) ;

  amc_transform_chain_transform_add(C1, R) ;
  amc_transform_chain_transform_add(C1, T) ;

  fprintf(stderr, "%s: assembling lifting line(s)\n", __FUNCTION__) ;
  llt_assembly_lifting_line_add(A, L1, C1) ;

  /* *Om *= 60/2.0/M_PI ; */
  /* *Om = lm*(*U)/b*60/2.0/M_PI ; */
  
  return ;
}

static void geometry_set(llt_assembly_t *A, gchar *g, gint np,
			 gdouble a, gdouble b,
			 gdouble *U, gdouble *Om)

{
  if ( strcmp(g, "rectangular") == 0 ) {
    assembly_wing(A, np, a, b, U) ;

    return ;
  }

  if ( strcmp(g, "prop") == 0 ) {
    assembly_propeller(A, np, a, b, U, Om) ;

    return ;
  }

  if ( strcmp(g, "elliptical") == 0 ) {
    assembly_elliptical(A, np, a, b, U) ;

    return ;
  }

  if ( strcmp(g, "van-garrel") == 0 ) {
    assembly_van_garrel(A, np, a, b, U, Om) ;

    return ;
  }
  
  g_error("unrecognised geometry \"%s\"", g) ;
  
  return ;
}

static void section_test(void)

{
  llt_section_t *s, *t ;
  gdouble a, cl, cd, al, cref ;
  FILE *f ;
  
  fprintf(stderr, "aerofoil section test\n") ;
  fprintf(stderr, "=====================\n") ;

  s = llt_section_alloc(1, 2) ;

  a = 2.0*M_PI ;

  llt_section_linear(s, a) ;

  al = 0.1 ;

  cref = al*a ;

  llt_section_interp(s, 1e5, al*180/M_PI, &cl, &cd) ;

  fprintf(stderr, "%lg %lg (%lg)\n", cl, cref, fabs(cl-cref)) ;

  f = file_open("section-test.dat", "w", NULL, NULL) ;
  llt_section_write(f, s) ;
  file_close(f) ;
  
  f = file_open("section-test.dat", "r", NULL, NULL) ;
  t = llt_section_read_alloc(f, 0) ;
  file_close(f) ;
  f = file_open("section-check.dat", "w", NULL, NULL) ;
  llt_section_write(f, t) ;
  file_close(f) ;
  
  return ;
}

static void script_location_test(void)

{
  gchar *script = "llt-xfoil" ;
  gchar *evar = "LLTXFOIL" ;
  gchar *path ;

  fprintf(stderr, "script location test\n") ;
  fprintf(stderr, "====================\n") ;

  fprintf(stderr, "searching for \"%s\"\n", script) ;

  path = llt_script_locate(script, evar) ;
  fprintf(stderr, "path to \"%s\" is \"%s\"\n", script, path) ;
  
  return ;
}

/* static void kernel_test(void) */

/* { */
/*   gdouble x[3], y[3], K[3], K1[3], K2[3], dK[9], dKc[9], err, ee, d ; */
/*   gint i ; */
  
/*   fprintf(stderr, "kernel gradient check\n") ; */
/*   fprintf(stderr, "=====================\n") ; */

/*   ee = 0.1 ; d = 1e-3 ; */
  
/*   x[0] = 0.3 ; x[1] = -0.9 ; x[2] = 1.7 ; */
/*   y[0] = 1.8 ; y[1] = 1.5 ; y[2] = -0.3 ; */

/*   llt_vorticity_kernel_WL(x, y, ee, K, dK) ; */

/*   x[0] += 0.5*d ; llt_vorticity_kernel_WL(x, y, ee, K1, NULL) ; */
/*   x[0] -=     d ; llt_vorticity_kernel_WL(x, y, ee, K2, NULL) ; */
/*   x[0] += 0.5*d ; */

/*   llt_vector_diff(&(dKc[0]), K1, K2) ; */
/*   dKc[0] /= d ; dKc[1] /= d ; dKc[2] /= d ;  */

/*   x[1] += 0.5*d ; llt_vorticity_kernel_WL(x, y, ee, K1, NULL) ; */
/*   x[1] -=     d ; llt_vorticity_kernel_WL(x, y, ee, K2, NULL) ; */
/*   x[1] += 0.5*d ; */

/*   llt_vector_diff(&(dKc[3]), K1, K2) ; */
/*   dKc[3] /= d ; dKc[4] /= d ; dKc[5] /= d ;  */

/*   x[2] += 0.5*d ; llt_vorticity_kernel_WL(x, y, ee, K1, NULL) ; */
/*   x[2] -=     d ; llt_vorticity_kernel_WL(x, y, ee, K2, NULL) ; */
/*   x[2] += 0.5*d ; */

/*   llt_vector_diff(&(dKc[6]), K1, K2) ; */
/*   dKc[6] /= d ; dKc[7] /= d ; dKc[8] /= d ;  */

/*   err = 0 ; */
/*   for ( i = 0 ; i < 9 ; i ++ ) { */
/*     fprintf(stderr, "%lg %lg (%lg)\n", dK[i], dKc[i], fabs(dK[i]-dKc[i])) ; */
/*     err = MAX(err, fabs(dK[i]-dKc[i])) ; */
/*   } */

/*   fprintf(stderr, "maximum error = %lg\n", err) ; */
  
/*   return ; */
/* } */

/* static void vortex_test(void) */

/* { */
/*   gdouble *u, th, ph, r, a, dV, V, x[3], w[3], w0, *p, s, e, dt, U, G, rho, g ; */
/*   gdouble xw, yw, zw, rw, tw ; */
/*   llt_vorticity_t *R ; */
/*   gint np, nph, nth, ns, i, j, k, n, nx, ny, nz ; */

/*   G = 1.5 ; r = 1.7 ; a = 0.3 ; */
/*   U = G/4.0/M_PI/r*(log(8.0*r/a) - 0.25) ; */
/*   w0 = G/M_PI/a/a ; */
  
/*   fprintf(stderr, "vortex particle check\n") ; */
/*   fprintf(stderr, "=====================\n") ; */
/*   fprintf(stderr, "circulation = %lg\n", G) ; */
/*   fprintf(stderr, "radius      = %lg\n", r) ; */
/*   fprintf(stderr, "core radius = %lg\n", a) ; */
/*   fprintf(stderr, "vorticity   = %lg\n", w0) ; */
/*   fprintf(stderr, "velocity    = %lg\n", U) ; */
  
/*   nth = 1024 ; nph = 8 ; ns = 4 ;  */
/*   e = 0.025 ; dt = 0.025 ; */
/*   np = nth*(nph*ns + 1) ; */

/*   np = 65536 ; */
/*   u = (gdouble *)g_malloc0(np*12*sizeof(gdouble)) ; */

/*   R = llt_vorticity_alloc(np) ; */

/*   llt_vortex_particle_number(R) = 0 ; */

/*   nx = 64 ; ny = 64 ; nz = 16 ; */

/* #if 1   */
/*   dV = 2*(r+a)/(nx-1)*2*(r+a)/(ny-1)*2*a/(nz-1) ; */
/*   g = 0.0 ; V = 0.0 ; */
/*   for ( i = 0 ; i <= nx ; i ++ ) { */
/*     xw = -(r+a) + 2*(r+a)*i/nx ; */
/*     for ( j = 0 ; j <= ny ; j ++ ) { */
/*       yw = -(r+a) + 2*(r+a)*j/ny ; */
/*       rw = sqrt(xw*xw + yw*yw) ; */
/*       tw = atan2(yw, xw) ; */
/*       for ( k = 0 ; k <= nz ; k ++ ) { */
/* 	zw = -a + 2*a*k/nz ; */
/* 	/\*distance from vortex centreline*\/ */
/* 	s = (rw - r)*(rw - r) + zw*zw ; */
/* 	if ( s < 1.1*a*a ) { */
/* 	  s = sqrt(s) ; */
/* 	  x[0] = xw ; x[1] = yw ; x[2] = zw ; */
/* 	  w[0] = -w0*sin(tw) ; w[1] = w0*cos(tw) ; w[2] = 0.0 ; */
/* 	  /\* w[0] *= exp(-s*s/a/a) ; *\/ */
/* 	  /\* w[1] *= exp(-s*s/a/a) ; *\/ */
/* 	  /\* w[2] *= exp(-s*s/a/a) ; *\/ */
/* 	  llt_vortex_particle_add(R, x, w, dV) ; */
/* 	  /\* g += dV ; *\/ */
/* 	  V += dV ; */
/* 	} */
/*       } */
/*     } */
/*   } */
/* #endif */
  
/* #if 0   */
/*   s = 2.0*M_PI/nth*r ; */
/*   /\* V = 4.0/3.0*M_PI*s*s*s ; *\/ */
/*   V = M_PI*a*a*s ; */
/*   /\* V = M_PI*a*a*2*M_PI*r/nth ; *\/ */
/*   dV = M_PI*a*a*2*M_PI*r/nth/ns/nph ; */
/*   g = 0.0 ; V = 0.0 ; */
/*   for ( i = 0 ; i < nth ; i ++ ) { */
/*     th = 2.0*M_PI*i/nth ; */
/*     s = 0 ; ph = 0 ; */
/*     rho = r + s*cos(ph) ; */
/*     x[0] = rho*cos(th) ; x[1] = rho*sin(th) ; x[2] = s*sin(ph) ; */
/*     w[0] = -w0*sin(th) ; w[1] = w0*cos(th) ; w[2] = 0 ; */
/*     llt_vortex_particle_add(R, x, w, dV) ; */
/*     g += w0*dV ; V += dV ; */
    
/*     for ( j = 0 ; j < nph ; j ++ ) { */
/*       ph = j*2.0*M_PI/nph ; */
/*       for ( k = 1 ; k < ns ; k ++ ) { */
/* 	s = a*k/(ns-1) ; */
/* 	rho = r + s*cos(ph) ; */
/* 	x[0] = rho*cos(th) ; x[1] = rho*sin(th) ; x[2] = s*sin(ph) ; */
/* 	w[0] = -w0*sin(th) ; w[1] = w0*cos(th) ; w[2] = 0 ; */
/* 	llt_vortex_particle_add(R, x, w, dV) ; */
/* 	g += w0*dV ; V += dV ; */
/*       } */
/*     } */
    
/*     /\* x[0] = r*cos(th) ; x[1] = r*sin(th) ; x[2] = 0 ; *\/ */
/*     /\* w[0] = -w0*sin(th) ; w[1] = w0*cos(th) ; w[2] = 0 ; *\/ */
/*     /\* llt_vortex_particle_add(R, x, w, V) ; *\/ */
/*   } */
/* #endif */

/* #if 0 */
/*   dV = r*2.0*M_PI/nth*a*a*M_PI ; */
/*   for ( i = 0 ; i < nth ; i ++ ) { */
/*     th = 2.0*M_PI*i/nth ; */
/*     rho = r ; */
/*     x[0] = rho*cos(th) ; x[1] = rho*sin(th) ; x[2] = s*sin(ph) ; */
/*     w[0] = -w0*sin(th) ; w[1] = w0*cos(th) ; w[2] = 0 ; */
/*     llt_vortex_particle_add(R, x, w, dV) ; */
/*   } */
/* #endif */
  
/*   fprintf(stderr, "circulation = %lg\n", g) ; */
/*   fprintf(stderr, "volume      = %lg\n", V) ; */

/*   np = llt_vortex_particle_number(R) ; */
/*   fprintf(stderr, "particles   = %d\n", np) ; */
/*   for ( i = 0 ; i < np ; i ++ ) { */
/*     p = llt_vortex_particle(R, i) ; */
/*     fprintf(stdout, "%e %e %e %e %e %e\n", */
/* 	    p[0], p[1], p[2], p[3], p[4], p[5]) ; */
/*   } */
  
/*   for ( n = 0 ; n < 160 ; n ++ ) { */
/*     memset(u, 0, 12*np*sizeof(gdouble)) ; */
/*     for ( i = 0 ; i < np ; i ++ ) { */
/*       llt_vorticity_velocity_gradient(R, e, LLT_KERNEL_WINCKELMANS_LEONARD, */
/* 				      llt_vortex_particle(R, i), */
/* 				      &(u[12*i+0]), &(u[12*i+3])) ; */
/*     } */
    
/*     for ( i = 0 ; i < np ; i ++ ) { */
/*       p = llt_vortex_particle(R, i) ; */
/*       p[0] += u[12*i+0]*dt ; p[1] += u[12*i+1]*dt ; p[2] += u[12*i+2]*dt ;  */
/*     } */
/*   } */
  
/*   for ( i = 0 ; i < np ; i ++ ) { */
/*     p = llt_vortex_particle(R, i) ; */
/*     fprintf(stdout, "%e %e %e %e %e %e\n", */
/* 	    p[0], p[1], p[2], p[3], p[4], p[5]) ; */
/*   } */
  
/*   return ; */
/* } */

gint main(gint argc, gchar **argv)

{
  gint np, test, nt ;
  llt_assembly_t *A ;
  gdouble U, Om, t0, t1, a, b ;
  gchar ch ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  np = 20 ;

  t0 = 0.0 ; t1 = 2.0 ; nt = 128 ;
  
  test = 0 ;

  U = 1.0 ; Om = 100 ;
  a = 1 ; b = 5 ;
  A = llt_assembly_alloc() ;

  while ( (ch = getopt(argc, argv, "a:b:g:n:O:p:s:T:t:U:")) != EOF ) {
    switch (ch ) {
    default: g_assert_not_reached() ; break ;
    case 'a': a = atof(optarg) ; break ;
    case 'b': b = atof(optarg) ; break ;
    case 'g': geometry_set(A, optarg, np, a, b, &U, &Om) ; break ;
    case 'n': nt = atoi(optarg) ; break ;
    case 'O': Om = atof(optarg) ; break ;
    case 'p': np = atoi(optarg) ; break ;
    case 's': if ( (test = parse_test(optarg)) == -1 ) {
	fprintf(stderr, "unrecognised test \"%s\"\n", optarg) ;
	return 1 ;
      }
      break ;
    case 't': t0 = atof(optarg) ; break ;
    case 'T': t1 = atof(optarg) ; break ;
    case 'U': U  = atof(optarg) ; break ;
    }
  }

  Om *= 2.0*M_PI/60.0 ;
  
  if ( test == 0 ) {
    wake_test(A, t0, t1, nt) ;

    return 0 ;
  }
  
  if ( test == 1 ) {
    solver_test(A, t0, t1, nt) ;

    return 0 ;
  }

  if ( test == 2 ) {
    section_test() ;

    return 0 ;
  }

  if ( test == 3 ) {
    script_location_test() ;

    return 0 ;
  }
  
  return 0 ;
}
