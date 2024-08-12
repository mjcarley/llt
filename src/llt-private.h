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

#ifndef __LLT_PRIVATE_H_INCLUDED__
#define __LLT_PRIVATE_H_INCLUDED__

#include <glib.h>

/*r = x - y*/
#define llt_vector_diff(_r,_x,_y)		\
  do {						\
    (_r)[0] = (_x)[0] - (_y)[0] ;		\
    (_r)[1] = (_x)[1] - (_y)[1] ;		\
    (_r)[2] = (_x)[2] - (_y)[2] ;		\
  } while (0)

#define llt_vector_scalar(_x,_y)		\
  ((_x)[0]*(_y)[0]+(_x)[1]*(_y)[1]+(_x)[2]*(_y)[2])

#define llt_vector_cross(_u,_x,_y)					\
  do {									\
    (_u)[0] = (_x)[1]*(_y)[2] - (_x)[2]*(_y)[1] ;			\
    (_u)[1] = (_x)[2]*(_y)[0] - (_x)[0]*(_y)[2] ;			\
    (_u)[2] = (_x)[0]*(_y)[1] - (_x)[1]*(_y)[0] ;			\
  } while (0)

#define llt_vector_length(_x) sqrt(llt_vector_scalar(_x,_x))

#define llt_vector_distance2(_x,_y)		\
  (((_x)[0] - (_y)[0])*((_x)[0] - (_y)[0])	\
   + ((_x)[1] - (_y)[1])*((_x)[1] - (_y)[1])	\
   + ((_x)[2] - (_y)[2])*((_x)[2] - (_y)[2]))

#define llt_ring_length(_y1,_y2,_y3,_y4)	\
  (sqrt(llt_vector_distance2((_y1),(_y2))) +	\
   sqrt(llt_vector_distance2((_y2),(_y3))) +	\
   sqrt(llt_vector_distance2((_y3),(_y4))) +	\
   sqrt(llt_vector_distance2((_y4),(_y1))))

#endif /*__LLT_PRIVATE_H_INCLUDED__*/

