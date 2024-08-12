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

static void rotor_assemble(llt_assembly_t *A, llt_lifting_line_t *L,
			   llt_section_t *s,
			   gint nblade, gdouble *Om, gdouble *U)

{
  amc_transform_t *T, *R, *B[32] ;
  amc_transform_chain_t *C[32] ;
  gint i ;
  gchar rot[32] ;
  
  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    L->s[i] = s ;
  }    
  
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
  
  for ( i = 0 ; i < nblade ; i ++ ) {
    B[i] = amc_transform_alloc(3, 1) ; B[i]->order = 1 ;
    sprintf(rot, "%d*2*pi/%d", i, nblade) ;
    amc_transform_rotation_z(B[i], 0, rot, 1) ;
    amc_transform_expressions_compile(B[i]) ;
    C[i] = amc_transform_chain_alloc(8) ;
    amc_transform_chain_transform_add(C[i], B[i]) ;
    amc_transform_chain_transform_add(C[i], R) ;
    amc_transform_chain_transform_add(C[i], T) ;
    llt_assembly_lifting_line_add(A, L, C[i]) ;
  }

  return ;
}

static void wing_assemble(llt_assembly_t *A, llt_lifting_line_t *L,
			  llt_section_t *s,
			  gdouble *U, gdouble *Om, gdouble *h,
			  gdouble *al0, gdouble *al) 

{
  amc_transform_t *T, *P ;
  amc_transform_chain_t *C ;
  gint i ;
  
  for ( i = 0 ; i < llt_lifting_line_point_number(L) ; i ++ ) {
    L->s[i] = s ;
  }    
  
  /*forward flight parallel to x axis with vertical oscillation*/
  T = amc_transform_alloc(3, 1) ; T->order = 1 ;
  amc_transform_variable_add(T, "U", U) ;
  amc_transform_variable_add(T, "Omega", Om) ;
  amc_transform_variable_add(T, "h", h) ;
  if ( (*Om) != 0 ) {
    amc_transform_translation(T, 0, "-U*t", 0, NULL, 0, "h*sin(Omega*t)", 1) ;
  } else {
    amc_transform_translation(T, 0, "-U*t", 0, NULL, 0, NULL, 1) ;
  }
  amc_transform_expressions_compile(T) ;

  /*rotation about y axis (pitch) at velocity Om*/
  P = amc_transform_alloc(3, 1) ; P->order = 1 ;
  amc_transform_variable_add(P, "Omega", Om) ;
  amc_transform_variable_add(P, "alpha", al0) ;
  amc_transform_variable_add(P, "dalpha", al) ;
  amc_transform_rotation_y(P, 0, "alpha + dalpha*sin(Omega*t)", 1) ;
  amc_transform_expressions_compile(P) ;

  C = amc_transform_chain_alloc(8) ;
  amc_transform_chain_transform_add(C, P) ;
  amc_transform_chain_transform_add(C, T) ;
  llt_assembly_lifting_line_add(A, L, C) ;

  return ;
}

static void print_help(FILE *f, gdouble al0, gdouble al, gint nblade,
		       gdouble h, gint nt, gdouble Om,
		       gdouble t0, gdouble t1, gdouble U)

{
  fprintf(stderr,
	  "Usage:\n\n"
	  "  %s [options]\n\n"
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -A # wing incidence (%lg)\n"
	  "  -a # wing pitch oscillation amplitude (%lg)\n"
	  "  -B # number of blades in rotor problem (%d)\n"
	  "  -c # maximum number of time steps in wake before cutting off\n"
	  "  -d # filename stub for output data\n"
	  "  -H # wing plunge amplitude (%lg)\n"
	  "  -L # lifting line file name\n"
	  "  -n # number of time steps (%d)\n"
	  "  -O # angular velocity in rad/s (%lg)\n"
	  "  -R # angular velocity in rpm (%lg)\n"
	  "  -r select rotor problem\n"
	  "  -s # aerofoil section file\n"
	  "  -t # start time for simulation (%lg)\n"
	  "  -T # end time for simulation (%lg)\n"
	  "  -U # linear velocity (%lg)\n"
	  "  -W update wakes at each time step\n"
	  "  -w select wing problem\n",
	  progname, al0, al, nblade, h, nt, Om, Om*60/2.0/M_PI, t0, t1, U) ;
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gint nt, nblade, i, j, idx, cutwake ;
  llt_assembly_t *A ;
  llt_lifting_line_t *L ;
  llt_section_t *s ;
  llt_solver_t solver ;
  gdouble U, Om, t0, t1, t, *x, F[3], M[3], x0[]={0,0,0}, h, al, al0, dt ;
  gboolean update_wakes, wing, rotor, truncate_wakes ;
  gchar ch, *lfile, *sfile, fname[128], *dstub, dfile[512] ;
  FILE *f ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  
  t0 = 0.0 ; t1 = 2.0 ; nt = 128 ;
  update_wakes = FALSE ;
  wing = FALSE ; rotor = FALSE ;
  truncate_wakes = FALSE ; cutwake = 0 ;
  
  U = 1.0 ; Om = 0 ; nblade = 1 ;
  h = 0.0 ; al = 0.0 ; al0 = 0 ;
  A = llt_assembly_alloc() ;
  s = NULL ;
  L = llt_lifting_line_alloc(512) ;

  sfile = lfile = NULL ;
  dstub = NULL ;
  
  while ( (ch = getopt(argc, argv, "hA:a:c:d:B:H:L:n:O:R:rs:t:T:U:Ww")) !=
	  EOF ) {
    switch (ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help(stderr, al0, al, nblade, h, nt, Om, t0, t1, U) ;
      return 0 ;
      break ;
    case 'A': al0 = atof(optarg) ; break ;
    case 'a': al = atof(optarg) ; break ;
    case 'B': nblade = atoi(optarg) ; break ;
    case 'c':
      cutwake = atoi(optarg) ;
      truncate_wakes = TRUE ;
      break ;
    case 'd': dstub = g_strdup(optarg) ; break ;
    case 'H': h = atof(optarg) ; break ;
    case 'L': lfile = g_strdup(optarg) ; break ;
    case 'n': nt = atoi(optarg) ; break ;
    case 'O': Om = atof(optarg) ; break ;
    case 'R': Om = atof(optarg)*2.0*M_PI/60.0 ; break ;
    case 'r': rotor = TRUE ; break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    case 't': t0 = atof(optarg) ; break ;
    case 'T': t1 = atof(optarg) ; break ;
    case 'U': U  = atof(optarg) ; break ;
    case 'W': update_wakes = TRUE ; break ;
    case 'w': wing = TRUE ; break ;
    }
  }

  if ( lfile == NULL ) {
    fprintf(stderr, "%s: lifting line file name not specified\n",
	    progname) ;
    return 1 ;
  } else {
    fprintf(stderr, "%s: reading lifting line data from \"%s\"\n",
	    progname, lfile) ;
    f = file_open(lfile, "r", "-", stdin) ;
    llt_lifting_line_read(f, L) ;
    file_close(f) ;
  }

  if ( sfile == NULL ) {
    fprintf(stderr,
	    "%s: no section data file specified, using thin aerofoil\n",
	    progname) ;
    s = llt_section_alloc(1, 2) ;

    llt_section_linear(s, 2.0*M_PI) ;
  } else {
    fprintf(stderr, "%s: reading section data from \"%s\"\n",
	    progname, sfile) ;
    f = file_open(sfile, "r", "-", stdin) ;
    s = llt_section_read_alloc(f, 0) ;
    file_close(f) ;
  }

  if ( !rotor && !wing ) {
    fprintf(stderr, "%s: select wing or rotor problem\n", progname) ;
    return 1 ;
  }
  
  if ( rotor ) {
    if ( wing ) {
      fprintf(stderr, "%s: select wing or rotor problem, not both\n",
	      progname) ;
      return 1 ;
    }
    /*rotor problem*/
    if ( Om == 0 ) {
      fprintf(stderr, "%s: rotor problem requires finite rotation speed\n",
	      progname) ;
      return 1 ;
    }

    rotor_assemble(A, L, s, nblade, &Om, &U) ;
  }

  if ( wing ) {
    if ( rotor ) {
      fprintf(stderr, "%s: select wing or rotor problem, not both\n",
	      progname) ;
      return 1 ;
    }
    /*wing problem*/

    wing_assemble(A, L, s, &U, &Om, &h, &al0, &al) ;
  }
  
  llt_assembly_time_steps(A) = nt+1 ;
  llt_assembly_initialise(A, t0) ;
  llt_assembly_regularisation(A) = 0.025 ;

  llt_solver_method(&solver) = LLT_SOLVER_RELAXATION ;
  llt_solver_relaxation_factor(&solver) = 0.03125 ;
  llt_solver_tolerance(&solver) = 1e-3 ;
  llt_solver_iteration_number_max(&solver) = 200 ;
  llt_solver_viscosity(&solver) = 1.5e-5 ;

  for ( i = 1 ; i < nt ; i ++ ) {
    t = t0 + (t1 - t0)*i/(nt-1) ;
    dt = (t1 - t0)/(nt-1) ;
    fprintf(stderr, "%s: time step %d: t = %lg\n", progname, i, t) ;

    llt_solver_step(A, &solver, t, dt, TRUE) ;
    fprintf(stderr, "     %d iterations; error = %lg; \n",
	    llt_solver_iteration_number(&solver),
	    llt_solver_error(&solver)) ;
    if ( update_wakes ) llt_solver_wakes_update(A, (t1 - t0)/(nt-1)) ;
    llt_assembly_force_moment(A, F, x0, M) ;
    /* fprintf(stdout, "%e %e %e %e\n", t, F[0], F[1], F[2]) ; */
    /* fflush(stdout) ; */

    if ( dstub != NULL ) {
      sprintf(dfile, dstub, i) ;
      f = file_open(dfile, "w", "-", stdout) ;
      L = llt_assembly_lifting_line(A, 0) ;
      for ( j = 0 ; j < llt_lifting_line_point_number(L) ; j ++ ) {
	idx = A->indices[0] + j ;
	x = llt_lifting_line_point(L,j) ;
	fprintf(f, "%e %e %e %e %e %e %e\n",
		x[0], x[1], x[2], A->cl[idx], A->Re[idx],
		A->vs[2*idx+0], A->vs[2*idx+1]) ;
      }
      file_close(f) ;
    }
    if ( truncate_wakes ) {
      if ( i >= cutwake ) {
	for ( j = 0 ; j < llt_assembly_lifting_line_number(A) ; j ++ ) {
	  llt_lattice_truncate(llt_assembly_wake(A,j), 1) ;
	}
      }
    }    
    
  }

  for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) {
    sprintf(fname, "line-%d.dat", i) ;
    f = fopen(fname, "w") ;
    llt_lifting_line_geometry_write(f, llt_assembly_lifting_line_current(A,i)) ;
    fclose(f) ;
    if ( llt_assembly_wake(A,i) != NULL ) {
      sprintf(fname, "wake-%d.dat", i) ;
      f = fopen(fname, "w") ;
      llt_lattice_write(f, llt_assembly_wake(A,i)) ;
      fclose(f) ;
    }
  }

  /* for ( i = 0 ; i < llt_assembly_lifting_line_number(A) ; i ++ ) { */
  /*   L = llt_assembly_lifting_line(A, i) ; */
  /*   for ( j = 0 ; j < llt_lifting_line_point_number(L) ; j ++ ) { */
  /*     idx = A->indices[i] + j ; */
  /*     x = llt_lifting_line_point(L,j) ; */
  /*     fprintf(stdout, "%e %e %e %e %e\n", */
  /* 	      x[0], x[1], x[2], A->cl[idx], A->Re[idx]) ; */
  /*   } */
  /* } */

  /* llt_assembly_lifting_line_force_moment(A, 0, F, x0, M) ; */
  llt_assembly_force_moment(A, F, x0, M) ;
  fprintf(stderr, "force: %e %e %e\n", F[0], F[1], F[2]) ;
  fprintf(stderr, "moment: %e %e %e\n", M[0], M[1], M[2]) ;
  
  return 0 ;
}
