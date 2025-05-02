/* GNU Ocrad - Optical Character Recognition program
   Copyright (C) 2003-2024 Antonio Diaz Diaz.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "common.h"
#include "user_filter.h"

int verbosity = 0;

bool Ocrad::similar( const int a, const int b,
                     const int percent_dif, const int abs_dif )
  {
  int difference = std::abs( a - b );
  if( percent_dif > 0 && difference <= abs_dif ) return true;
  int max_abs = std::max( std::abs( a ), std::abs( b ) );
  return ( difference * 100 <= max_abs * percent_dif );
  }


Control::~Control()
  {
  for( unsigned f = filters.size(); f > 0; --f )
    if( filters[f-1].user_filterp )
      delete filters[f-1].user_filterp;
  }
