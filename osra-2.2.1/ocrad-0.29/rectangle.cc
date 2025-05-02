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
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "common.h"
#include "rectangle.h"


namespace {

void error( const char * const msg ) { throw Ocrad::Internal( msg ); }

} // end namespace


Rectangle::Rectangle( const int l, const int t, const int r, const int b )
  {
  if( r < l || b < t )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "l = %d, t = %d, r = %d, b = %d\n", l, t, r, b );
    error( "bad parameter building a Rectangle." );
    }
  left_ = l; top_ = t; right_ = r; bottom_ = b;
  }

Rectangle::Rectangle( const int h, const int w )
  {
  if( h <= 0 || w <= 0 )
    {
    if( verbosity >= 0 ) std::fprintf( stderr, "h = %d, w = %d\n", h, w );
    error( "bad parameter building a Rectangle." );
    }
  left_ = 0; top_ = 0; right_ = w - 1; bottom_ = h - 1;
  }


void Rectangle::left( const int l )
  {
  if( l <= right_ ) left_ = l;
  else error( "left: bad parameter resizing a Rectangle." );
  }

void Rectangle::top( const int t )
  {
  if( t <= bottom_ ) top_ = t;
  else error( "top: bad parameter resizing a Rectangle." );
  }

void Rectangle::right( const int r )
  {
  if( r >= left_ ) right_ = r;
  else error( "right: bad parameter resizing a Rectangle." );
  }

void Rectangle::bottom( const int b )
  {
  if( b >= top_ ) bottom_ = b;
  else error( "bottom: bad parameter resizing a Rectangle." );
  }

void Rectangle::height( const int h )
  {
  if( h > 0 ) bottom_ = top_ + h - 1;
  else error( "height: bad parameter resizing a Rectangle." );
  }

void Rectangle::width( const int w )
  {
  if( w > 0 ) right_ = left_ + w - 1;
  else error( "width: bad parameter resizing a Rectangle." );
  }


int Rectangle::v_overlap_percent( const Rectangle & re ) const
  {
  int ov = std::min( bottom_, re.bottom_ ) - std::max( top_, re.top_ ) + 1;
  if( ov > 0 )
    ov = std::max( 1, ( ov * 100 ) / std::min( height(), re.height() ) );
  else ov = 0;
  return ov;
  }


bool Rectangle::is_hcentred_in( const Rectangle & re ) const
  {
  if( this->h_includes( re.hcenter() ) ) return true;
  int w = std::min( re.height(), re.width() ) / 2;
  if( width() < w )
    {
    int d = ( w + 1 ) / 2;
    if( hcenter() - d <= re.hcenter() && hcenter() + d >= re.hcenter() )
      return true;
    }
  return false;
  }


bool Rectangle::is_vcentred_in( const Rectangle & re ) const
  {
  if( this->v_includes( re.vcenter() ) ) return true;
  int h = std::min( re.height(), re.width() ) / 2;
  if( height() < h )
    {
    int d = ( h + 1 ) / 2;
    if( vcenter() - d <= re.vcenter() && vcenter() + d >= re.vcenter() )
      return true;
    }
  return false;
  }


bool Rectangle::precedes( const Rectangle & re ) const
  {
  if( right_ < re.left_ ) return true;
  if( this->h_overlaps( re ) &&
      ( ( top_ < re.top_ ) || ( top_ == re.top_ && left_ < re.left_ ) ) )
    return true;
  return false;
  }


bool Rectangle::v_precedes( const Rectangle & re ) const
  {
  if( bottom_ < re.vcenter() || vcenter() < re.top_ ) return true;
  if( this->includes_vcenter( re ) && re.includes_vcenter( *this ) )
    return this->h_precedes( re );
  return false;
  }


int Rectangle::hypoti( const int c1, const int c2 )
  {
  long long temp = c1; temp *= temp;
  long long target = c2; target *= target; target += temp;
  int lower = std::max( std::abs( c1 ), std::abs( c2 ) );
  int upper = std::abs( c1 ) + std::abs( c2 );
  while( upper - lower > 1 )
    {
    int m = ( lower + upper ) / 2;
    temp = m; temp *= temp;
    if( temp < target ) lower = m; else upper = m;
    }
  temp = lower; temp *= temp; target *= 2; target -= temp;
  temp = upper; temp *= temp;
  if( target < temp ) return lower;
  else return upper;
  }
