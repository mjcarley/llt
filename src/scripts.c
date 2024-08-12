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
#include <stdlib.h>
#include <errno.h>

#include <glib.h>

#include <amc.h>

#include "llt.h"
#include "llt-private.h"

llt_script_t *llt_script_alloc(gchar *name)

{
  llt_script_t *s ;

  s = (llt_script_t *)g_malloc0(sizeof(llt_script_t)) ;

  llt_script_name(s) = g_strdup(name) ;
  
  return s ;
}

gchar *llt_script_locate(gchar *script, gchar *evar)

{
  gchar *p ;
  const gchar *e ;

  if ( evar != NULL ) {
    e = g_getenv(evar) ;
    if ( e != NULL ) {
      p = g_find_program_in_path(e) ;
      if ( p != NULL ) return p ;
    }
  }

  if ( script == NULL ) {
    g_error("%s: script and path cannot both be NULL", __FUNCTION__) ;
  }
  p = g_find_program_in_path(script) ;
  if ( p != NULL ) return p ;
  
  return NULL ;
}

gint llt_script_init(llt_script_t *s, llt_section_solver_t type,
		     gchar *dcty, gboolean dcty_force)


{
  gint mode ;

  mode = 0755 ;
  
  if ( type == LLT_SECTION_SOLVER_UNDEFINED ) {
    g_error("%s: cannot initialise script of undefined type", __FUNCTION__) ;
  }

  switch ( type ) {
  default:
    g_error("%s: unknown script type %d", __FUNCTION__, llt_script_type(s)) ;
    break ;
  case LLT_SECTION_SOLVER_XFOIL:
    /*find path to executable script*/
    if ( llt_script_exec(s) == NULL )
      llt_script_exec(s) = g_strdup("llt-xfoil") ;
    llt_script_path(s) = llt_script_locate(llt_script_exec(s), "LLTXFOIL") ;
    llt_script_output_size(s) = 7 ;
    llt_script_output_index_alpha(s) = 0 ;
    llt_script_output_index_cl(s)    = 1 ;
    llt_script_output_index_cd(s)    = 2 ;
    break ;
  }

  llt_script_type(s) = type ;

  if ( llt_script_path(s) == NULL ) {
    g_error("%s: cannot locate executable script \"%s\" for \"%s\"",
	    __FUNCTION__, llt_script_exec(s), llt_script_name(s)) ;
  }

  /*expand directory path and check it*/
  llt_script_dcty(s) = realpath(dcty, NULL) ;

  if ( llt_script_dcty(s) == NULL ) {
    /*check errno for non-existent directory and make it if necessary*/
    if ( errno == ENOENT ) {
      if ( dcty_force ) {
	g_mkdir_with_parents(dcty, mode) ;
	llt_script_dcty(s) = realpath(dcty, NULL) ;
      }
    }
  }
  
  return 0 ;
}  

/** 
 * Call a script to generate aerofoil section data
 * 
 * @param s an ::llt_script_t containing the script information;
 * @param aerofoil aerofoil geometry file;
 * @param polar on output file containing aerofoil aerodynamic data;
 * @param Re Reynolds number;
 * @param amin minimum \f$\alpha\f$ for polar data;
 * @param amax maximum \f$\alpha\f$ for polar data;
 * @param da increment in incidence \f$\Delta\alpha\f$.
 * 
 * @return 0 on success.
 */

gint llt_script_section_make(llt_script_t *s, gchar *aerofoil,
			     gchar *polar,
			     gdouble Re, gdouble amin, gdouble amax,
			     gdouble da)

{
  gchar cmd[2048] ;
  gint status ;
  
  if ( llt_script_type(s) == LLT_SECTION_SOLVER_UNDEFINED ) {
    g_error("%s: cannot run uninitialised script", __FUNCTION__) ;
  }

  sprintf(cmd, "cp %s %s/llt-aerofoil.dat", aerofoil, llt_script_dcty(s)) ;
  status = system(cmd) ;
  if ( status != 0 ) {
    g_error("%s: system error on \"%s\"", __FUNCTION__, cmd) ;
  }

  switch ( llt_script_type(s) ) {
  default:
    g_error("%s: unknown script type %d",
	    __FUNCTION__, llt_script_type(s)) ;
    break ;
  case LLT_SECTION_SOLVER_XFOIL:
    sprintf(cmd,
	    "%s --directory=%s --aerofoil=llt-aerofoil.dat --Reynolds=%lg "
	    "--polar-file=llt-polar.dat --alpha-min=%lg --alpha-max=%lg "
	    "--delta-alpha=%lg",
	    llt_script_path(s), llt_script_dcty(s), Re, amin, amax, da) ;
    break ;
  }

  status = system(cmd) ;
  if ( status != 0 ) {
    g_error("%s: system error on \"%s\"", __FUNCTION__, cmd) ;
  }
  
  sprintf(cmd, "cp %s/llt-polar.dat %s", llt_script_dcty(s), polar) ;
  status = system(cmd) ;
  if ( status != 0 ) {
    g_error("%s: system error on \"%s\"", __FUNCTION__, cmd) ;
  }

  return 0 ;
}

gint llt_script_section_add(llt_script_t *script, gchar *aerofoil,
			    gdouble Re, gdouble amin, gdouble amax,
			    gdouble da, llt_section_t *s)

{
  gchar *polar = "llt-polar-tmp.dat" ;
  FILE *f ;
  gint i, j, nlines ;
  gdouble data[32], cl, cd, al ;
  
  llt_script_section_make(script, aerofoil, polar, Re, amin, amax, da) ;

  f = fopen(polar, "r") ;
  if ( f == NULL ) {
    g_error("%s: cannot open temporary polar file \"%s\"",
	    __FUNCTION__, polar) ;
  }

  fscanf(f, "%d", &nlines) ;

  for ( i = 0 ; i < nlines ; i ++ ) {
    for ( j = 0 ; j < llt_script_output_size(script) ; j ++ ) {
      fscanf(f, "%lg", &(data[j])) ;
    }
    al = cl = cd = 0.0 ;
    if ( llt_script_output_index_alpha(script) >= 0 )
      al = data[llt_script_output_index_alpha(script)] ;
    if ( llt_script_output_index_cl(script) >= 0 )
      cl = data[llt_script_output_index_cl(script)] ;
    if ( llt_script_output_index_cd(script) >= 0 )
      cd = data[llt_script_output_index_cd(script)] ;
    llt_section_entry_append(s, Re, al, cl, cd) ;
  }
  
  fclose(f) ;
  
  return 0 ;
}
