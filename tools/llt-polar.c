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

/* static gint compare_double(const void *a, const void *b) */

/* { */
/*   if ( *((gdouble *)a) < *((gdouble *)b) ) return -1 ; */
/*   if ( *((gdouble *)a) > *((gdouble *)b) ) return  1 ; */

/*   return 0 ; */
/* } */

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

static void print_help(FILE *f, gdouble alpha_min, gdouble alpha_max,
		       gdouble delta_alpha)

{
  fprintf(stderr,
	  "Usage:\n\n"
	  "  %s [options]\n\n"
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -a # minimum alpha (%lg)\n"
	  "  -A # maximum alpha (%lg)\n"	  
	  "  -d # delta alpha (%lg)\n",
	  progname, alpha_min, alpha_max, delta_alpha) ;
  
  return ;
}

static gint parse_reynolds(gchar *arg, gdouble *Re, gint *nRe)

{
  gchar **tok, *nptr ;

  tok = g_strsplit(arg, ",", 0) ;

  for ( *nRe = 0 ; tok[(*nRe)] != NULL ; (*nRe) ++ ) {
    Re[(*nRe)] = strtod(tok[(*nRe)], &nptr) ;
    if ( Re[(*nRe)] == 0 ) {
      /*check if conversion was performed*/
      if ( nptr == tok[(*nRe)] ) {
	fprintf(stderr, "%s: cannot parse \"%s\"\n", progname, tok[(*nRe)]) ;
	return 1 ;
      }
    }
    if ( Re[(*nRe)] <= 0 ) {
      fprintf(stderr, "%s: Reynolds numbers must be strictly positive (%s)\n",
	      progname, tok[(*nRe)]) ;
      return 1 ;
    }
  }
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gchar ch, *sfile, *gfile, *name ;
  gdouble alpha_min, alpha_max, delta_alpha, Re[64] ;
  gint nRe, i, na ;
  llt_section_t *s ;
  llt_script_t *script ;
  FILE *f ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  alpha_min = -10 ; alpha_max = 10 ; delta_alpha = 0.25 ;
  sfile = NULL ; gfile = NULL ; name = NULL ;
  
  while ( (ch = getopt(argc, argv, "ha:A:d:g:n:R:s:")) != EOF ) {
    switch (ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help(stderr, alpha_min, alpha_max, delta_alpha) ;
      return 0 ;
      break ;
    case 'a': alpha_min = atof(optarg) ; break ;
    case 'A': alpha_max = atof(optarg) ; break ;
    case 'd': delta_alpha = atof(optarg) ; break ;
    case 'g': gfile = g_strdup(optarg) ; break ;
    case 'n': name = g_strdup(optarg) ; break ;
    case 'R':
      if ( parse_reynolds(optarg, Re, &nRe) != 0 ) return 1 ;
      break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    }
  }

  if ( sfile == NULL ) sfile = g_strdup("-") ;
  if ( name == NULL ) name = g_strdup("LLT") ;

  if ( gfile == NULL ) {
    fprintf(stderr, "%s: no geometry file specified\n", progname) ;
    return 1 ;
  }    
  
  /* qsort(Re, nRe, sizeof(gdouble), compare_double) ; */
  
  fprintf(stderr, "%s: %d Reynolds numbers\n", progname, nRe) ;
  for ( i = 0 ; i < nRe ; i ++ ) {
    fprintf(stderr, "  %lg\n", Re[i]) ;
  }    

  na = (gint)ceil((alpha_max - alpha_min)/delta_alpha) + 1 ;
  fprintf(stderr, "%s: %d angles of incidence %lg--%lg (%lg)\n",
	  progname, na, alpha_min, alpha_max, delta_alpha) ;

  s = llt_section_alloc(nRe, na*nRe) ;
  llt_section_name(s) = name ;
  script = llt_script_alloc("aerofoils") ;

  llt_script_init(script, LLT_SECTION_SOLVER_XFOIL, "./Xfoil", TRUE) ;

  fprintf(stderr, "%s: section solver \"%s\"\n",
	  progname, llt_script_name(script)) ;
  fprintf(stderr, "  running \"%s\" in directory \"%s\"\n",
	  llt_script_path(script), llt_script_dcty(script)) ;

  for ( i = 0 ; i < nRe ; i ++ ) {
    fprintf(stderr, "%s: generating section data for Re=%lg\n",
	    progname, Re[i]) ;
    llt_script_section_add(script, gfile, Re[i],
			   alpha_min, alpha_max, delta_alpha, s) ;
  }

  llt_section_sort(s) ;
  llt_section_index(s) ;

  f = file_open(sfile, "w", "-", stdout) ;
  llt_section_write(f, s) ;
  file_close(f) ;
  
  return 0 ;
}
