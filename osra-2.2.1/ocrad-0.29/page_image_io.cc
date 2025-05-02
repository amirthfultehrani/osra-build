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
#include <cctype>
#include <climits>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>

#include "common.h"
#include "page_image.h"


namespace {

uint8_t pnm_getrawbyte( FILE * const f )
  {
  int ch = std::fgetc( f );

  if( ch == EOF )
    throw Page_image::Error( "End-of-file reading pnm file." );

  return ch;
  }


uint8_t pnm_getc( FILE * const f )
  {
  uint8_t ch = pnm_getrawbyte( f );
  if( ch == '#' )			// comment
    { do ch = pnm_getrawbyte( f ); while( ch != '\n' ); }
  return ch;
  }


int pnm_getint( FILE * const f )
  {
  uint8_t ch;
  int i = 0;

  do ch = pnm_getc( f ); while( std::isspace( ch ) );
  if( !std::isdigit( ch ) )
    throw Page_image::Error( "Junk in pnm file where an integer should be." );
  do {
    if( ( INT_MAX - (ch - '0') ) / 10 < i )
      throw Page_image::Error( "Number too big in pnm file." );
    i = (i * 10) + (ch - '0');
    ch = pnm_getc( f );
    }
  while( std::isdigit( ch ) );
  return i;
  }


uint8_t pbm_getbit( FILE * const f )
  {
  uint8_t ch;

  do ch = pnm_getc( f ); while( std::isspace( ch ) );

  if( ch == '0' ) return 0;
  if( ch == '1' ) return 1;
  throw Page_image::Error( "Junk in pbm file where bits should be." );
  }

} // end namespace


void Page_image::test_size( const int rows, const int cols )
  {
  if( rows < 3 || cols < 3 )
    throw Error( "Image too small. Minimum size is 3x3." );
  if( INT_MAX / rows < cols )
    throw Error( "Image too big. 'int' will overflow." );
  }


void Page_image::read_p1( FILE * const f, const int cols, const bool invert )
  {
  maxval_ = 1; threshold_ = 0;
  const int rows = height();
  if( !invert )
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; ++col )
        data[row].push_back( 1 - pbm_getbit( f ) );
  else
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; ++col )
        data[row].push_back( pbm_getbit( f ) );
  }


void Page_image::read_p4( FILE * const f, const int cols, const bool invert )
  {
  maxval_ = 1; threshold_ = 0;
  const int rows = height();
  if( !invert )
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; )
        {
        uint8_t byte = pnm_getrawbyte( f );
        for( uint8_t mask = 0x80; mask > 0 && col < cols; mask >>= 1, ++col )
          data[row].push_back( ( byte & mask ) ? 0 : 1 );
        }
  else
    for( int row = 0; row < rows; ++row )
      for( int col = 0; col < cols; )
        {
        uint8_t byte = pnm_getrawbyte( f );
        for( uint8_t mask = 0x80; mask > 0 && col < cols; mask >>= 1, ++col )
          data[row].push_back( ( byte & mask ) ? 1 : 0 );
        }
  }


void Page_image::read_p2( FILE * const f, const int cols, const bool invert )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Error( "Zero maxval in pgm file." );
  maxval_ = std::min( maxval, 255 );
  threshold_ = maxval_ / 2;
  const int rows = height();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      int val = pnm_getint( f );
      if( val > maxval ) throw Error( "Pixel value > maxval in pgm file." );
      if( invert ) val = maxval - val;
      if( maxval > 255 ) { val *= 255; val /= maxval; }
      data[row].push_back( val );
      }
  }


void Page_image::read_p5( FILE * const f, const int cols, const bool invert )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Error( "Zero maxval in pgm file." );
  if( maxval > 255 ) throw Error( "maxval > 255 in pgm \"P5\" file." );
  maxval_ = maxval;
  threshold_ = maxval_ / 2;
  const int rows = height();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      uint8_t val = pnm_getrawbyte( f );
      if( val > maxval_ ) throw Error( "Pixel value > maxval in pgm file." );
      if( invert ) val = maxval_ - val;
      data[row].push_back( val );
      }
  }


void Page_image::read_p3( FILE * const f, const int cols, const bool invert )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Error( "Zero maxval in ppm file." );
  maxval_ = std::min( maxval, 255 );
  threshold_ = maxval_ / 2;
  const int rows = height();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      const int r = pnm_getint( f );			// Red value
      const int g = pnm_getint( f );			// Green value
      const int b = pnm_getint( f );			// Blue value
      if( r > maxval || g > maxval || b > maxval )
        throw Error( "Pixel value > maxval in ppm file." );
      int val;
      if( !invert ) val = std::min( r, std::min( g, b ) );
      else val = maxval - std::max( r, std::max( g, b ) );
      if( maxval > 255 ) { val *= 255; val /= maxval; }
      data[row].push_back( val );
      }
  }


void Page_image::read_p6( FILE * const f, const int cols, const bool invert )
  {
  const int maxval = pnm_getint( f );
  if( maxval == 0 ) throw Error( "Zero maxval in ppm file." );
  if( maxval > 255 ) throw Error( "maxval > 255 in ppm \"P6\" file." );
  maxval_ = maxval;
  threshold_ = maxval_ / 2;
  const int rows = height();

  for( int row = 0; row < rows; ++row )
    for( int col = 0; col < cols; ++col )
      {
      const uint8_t r = pnm_getrawbyte( f );	// Red value
      const uint8_t g = pnm_getrawbyte( f );	// Green value
      const uint8_t b = pnm_getrawbyte( f );	// Blue value
      if( r > maxval_ || g > maxval_ || b > maxval_ )
        throw Error( "Pixel value > maxval in ppm file." );
      uint8_t val;
      if( !invert ) val = std::min( r, std::min( g, b ) );	// darkest color
      else val = maxval_ - std::max( r, std::max( g, b ) );	// lightest color
      data[row].push_back( val );
      }
  }


/* Create a Page_image from a png, pbm, pgm, or ppm file.
   "P1" (pbm), "P4" (pbm RAWBITS), "P2" (pgm), "P5" (pgm RAWBITS),
   "P3" (ppm), "P6" (ppm RAWBITS) file formats are recognized.
*/
Page_image::Page_image( FILE * const f, const bool invert )
  {
  unsigned char filetype = 0;
  const unsigned char first_byte = pnm_getrawbyte( f );

  if( first_byte == 'P' )
    {
    unsigned char ch = pnm_getrawbyte( f );
    if( ch >= '1' && ch <= '6' ) filetype = ch;
    }
  else if( read_check_png_sig8( f, first_byte ) )
    { read_png( f, 8, invert ); return; }
  if( filetype == 0 )
    throw Error( "Bad magic number - not a png, pbm, pgm, or ppm file." );

  const int cols = pnm_getint( f );
  if( cols == 0 ) throw Error( "Zero width in pnm file." );
  const int rows = pnm_getint( f );
  if( rows == 0 ) throw Error( "Zero height in pnm file." );
  test_size( rows, cols );
  data.resize( rows );
  for( int row = 0; row < rows; ++row ) data[row].reserve( cols );

  switch( filetype )
    {
    case '1': read_p1( f, cols, invert ); break;
    case '2': read_p2( f, cols, invert ); break;
    case '3': read_p3( f, cols, invert ); break;
    case '4': read_p4( f, cols, invert ); break;
    case '5': read_p5( f, cols, invert ); break;
    case '6': read_p6( f, cols, invert ); break;
    }

  if( verbosity >= 1 )
    {
    std::fprintf( stderr, "file type is P%c\n", filetype );
    std::fprintf( stderr, "file size is %dw x %dh\n", width(), height() );
    }
  }


bool Page_image::save( FILE * const f, const char filetype ) const
  {
  if( filetype == '7' ) { write_png( f, 1 ); return true; }
  if( filetype == '8' ) { write_png( f, 8 ); return true; }
  if( filetype < '1' || filetype > '6' ) return false;
  std::fprintf( f, "P%c\n%d %d\n", filetype, width(), height() );
  if( filetype != '1' && filetype != '4' )
    std::fprintf( f, "%d\n", maxval_ );

  if( filetype == '1' )					// pbm
    for( int row = 0; row < height(); ++row )
      {
      for( int col = 0; col < width(); ++col )
        std::putc( get_bit( row, col ) ? '1' : '0', f );
      std::putc( '\n', f );
      }
  else if( filetype == '4' )				// pbm RAWBITS
    for( int row = 0; row < height(); ++row )
      {
      uint8_t byte = 0, mask = 0x80;
      for( int col = 0; col < width(); ++col )
        {
        if( get_bit( row, col ) ) byte |= mask;
        mask >>= 1;
        if( mask == 0 ) { std::putc( byte, f ); byte = 0; mask = 0x80; }
        }
      if( mask != 0x80 ) std::putc( byte, f ); // incomplete byte at end of row
      }
  else if( filetype == '2' )				// pgm
    for( int row = 0; row < height(); ++row )
      {
      for( int col = 0; col < width() - 1; ++col )
        std::fprintf( f, "%d ", data[row][col] );
      std::fprintf( f, "%d\n", data[row][width()-1] );
      }
  else if( filetype == '5' )				// pgm RAWBITS
    for( int row = 0; row < height(); ++row )
      for( int col = 0; col < width(); ++col )
        std::fprintf( f, "%c", data[row][col] );
  else if( filetype == '3' )				// ppm
    for( int row = 0; row < height(); ++row )
      {
      for( int col = 0; col < width() - 1; ++col )
        {
        const uint8_t d = data[row][col];
        std::fprintf( f, "%d %d %d ", d, d, d );
        }
      const uint8_t d = data[row][width()-1];
      std::fprintf( f, "%d %d %d\n", d, d, d );
      }
  else if( filetype == '6' )				// ppm RAWBITS
    for( int row = 0; row < height(); ++row )
      for( int col = 0; col < width(); ++col )
        {
        const uint8_t d = data[row][col];
        std::fprintf( f, "%c%c%c", d, d, d );
        }
  return true;
  }
