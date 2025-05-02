/* GNU Ocrad - Optical Character Recognition program
   Copyright (C) 2022-2024 Antonio Diaz Diaz.

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

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <string>
#include <vector>
#include <stdint.h>	// uint8_t
#include <png.h>

#include "common.h"
#include "page_image.h"

#if PNG_LIBPNG_VER < 10209
#error "Wrong libpng version. At least 1.2.9 is required."
#endif


/* Check whether a file is a PNG file by using png_sig_cmp().
 * png_sig_cmp() returns zero if the image is a PNG, and nonzero otherwise.
 *
 * If this call is successful, and you are going to keep the file open,
 * you should call png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK); once
 * you have created the png_ptr, so that libpng knows your application
 * has read that many bytes from the start of the file.
 */
bool read_check_png_sig8( FILE * const f, const int first_byte )
  {
  enum { sig_size = 8 };
  png_byte buf[sig_size];
  const unsigned i = first_byte >= 0 && first_byte <= 0xFF;
  if( i ) buf[0] = first_byte;

  // read the signature bytes
  if( std::fread( buf + i, 1, sig_size - i, f ) != sig_size - i ) return false;
  return !png_sig_cmp( buf, 0, sig_size );
  }

/* Create a Page_image with elements in range [0, maxval] from a PNG file.
   'sig_read' indicates the number of magic bytes already read (maybe 0). */
void Page_image::read_png( FILE * const f, const int sig_read, const bool invert )
  {
  /* Create and initialize the png_struct with the desired error handler
   * functions.  If you want to use the default stderr and longjump method,
   * you can supply NULL for the last three parameters.  We also supply the
   * compiler header file version, so that we know if the application
   * was compiled with a compatible version of the library.  REQUIRED.
   */
  png_structp png_ptr =
    png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
  if( !png_ptr ) throw std::bad_alloc();

  /* Allocate/initialize the memory for image information.  REQUIRED. */
  png_infop info_ptr = png_create_info_struct( png_ptr );
  if( !info_ptr )
    {
    png_destroy_read_struct( &png_ptr, NULL, NULL );	// avoid memory leak
    throw std::bad_alloc();
    }

  if( setjmp( png_jmpbuf( png_ptr ) ) )		// Set error handling
    {
    /* Free all of the memory associated with the png_ptr and info_ptr. */
    png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
    throw Error( "Error reading PNG image." );
    }

  /* Set up the input control if you are using standard C streams. */
  png_init_io( png_ptr, f );

  /* If we have already read some of the signature */
  png_set_sig_bytes( png_ptr, sig_read );

  /* Read the entire image (including pixels) into the info structure with a
     call to png_read_png, equivalent to png_read_info(), followed by the
     set of transformations indicated by the transform mask, then
     png_read_image(), and finally png_read_end():
     PNG_TRANSFORM_PACKING does not expand 1-bit samples to bytes.
     (It seems to do nothing).
  */
  png_read_png( png_ptr, info_ptr,
                PNG_TRANSFORM_STRIP_ALPHA | PNG_TRANSFORM_EXPAND, NULL );

  /* At this point you have read the entire image.
     Now let's extract the data. */

  const int rows = png_get_image_height( png_ptr, info_ptr );
  const int cols = png_get_image_width( png_ptr, info_ptr );
  test_size( rows, cols );
  const unsigned bit_depth  = png_get_bit_depth( png_ptr, info_ptr );
  const unsigned maxval = ( 1 << bit_depth ) - 1;
  const unsigned color_type = png_get_color_type( png_ptr, info_ptr );
  const unsigned channels   = png_get_channels( png_ptr, info_ptr );
  /* bit_depth  - holds the bit depth of one of the image channels. (valid
                  values are 1, 2, 4, 8, 16 and depend also on color_type.
     color_type - describes which color/alpha channels are present.
                  PNG_COLOR_TYPE_GRAY (bit depths 1, 2, 4, 8, 16)
                  PNG_COLOR_TYPE_GRAY_ALPHA (bit depths 8, 16)
                  PNG_COLOR_TYPE_PALETTE (bit depths 1, 2, 4, 8)
                  PNG_COLOR_TYPE_RGB (bit_depths 8, 16)
                  PNG_COLOR_TYPE_RGB_ALPHA (bit_depths 8, 16)
                  PNG_COLOR_MASK_PALETTE
                  PNG_COLOR_MASK_COLOR
                  PNG_COLOR_MASK_ALPHA
     channels   - number of channels of info for the color type
                  (valid values are 1 (GRAY, PALETTE), 2 (GRAY_ALPHA),
                  3 (RGB), 4 (RGB_ALPHA or RGB + filler byte)) */

  if( ( color_type != PNG_COLOR_TYPE_GRAY && color_type != PNG_COLOR_TYPE_RGB ) ||
      ( channels != 1 && channels != 3 ) || bit_depth > 8 )
    throw Error( "Unsupported type of PNG image." );

  const png_bytepp row_pointers = png_get_rows( png_ptr, info_ptr );

  data.resize( rows );
  for( int row = 0; row < rows; ++row ) data[row].reserve( cols );
  maxval_ = maxval;
  threshold_ = maxval_ / 2;
  if( channels == 1 )
    for( int row = 0; row < rows; ++row )
      {
      const png_byte * ptr = row_pointers[row];
      for( int col = 0; col < cols; ++col )
        {
        const uint8_t val = *ptr++;
        data[row].push_back( invert ? maxval - val : val );
        }
      }
  else if( channels == 3 )
    for( int row = 0; row < rows; ++row )
      {
      const png_byte * ptr = row_pointers[row];
      for( int col = 0; col < cols; ++col )
        {
        const uint8_t r = *ptr++;
        const uint8_t g = *ptr++;
        const uint8_t b = *ptr++;
        uint8_t val;
        if( !invert ) val = std::min( r, std::min( g, b ) );
        else val = maxval_ - std::max( r, std::max( g, b ) );
        data[row].push_back( val );
        }
      }

  /* Clean up after the read, and free any memory allocated.  REQUIRED. */
  png_destroy_read_struct( &png_ptr, &info_ptr, NULL );

  if( verbosity >= 1 )
    {
    std::fprintf( stderr, "file type is PNG %s\n",
                  ( channels == 1 ) ? "greyscale" : "color" );
    std::fprintf( stderr, "file size is %dw x %dh\n", width(), height() );
    }
  }


// write greyscale images with bit depths of 1 or 8 only
//
void Page_image::write_png( FILE * const f, const unsigned bit_depth ) const
  {
  if( bit_depth != 1 && bit_depth != 8 )
    throw Error( "Invalid bit depth writing PNG image." );
  /* row_bytes is the width x number of channels x ceil(bit_depth / 8) */
  const unsigned row_bytes = width();		// width() * 1 * 1
  png_byte * const png_pixels = (png_byte *)std::malloc( height() * row_bytes );
  if( !png_pixels ) throw std::bad_alloc();
  png_byte ** const row_pointers =
    (png_byte **)std::malloc( height() * sizeof row_pointers[0] );
  if( !row_pointers ) { std::free( png_pixels ); throw std::bad_alloc(); }

  // fill png_pixels[] and set row_pointers
  int idx = 0;					// index in png_pixels[]
  if( bit_depth == 1 )
    for( int row = 0; row < height(); ++row )
      for( int col = 0; col < width(); ++col )
        png_pixels[idx++] = !get_bit( row, col );
  else if( maxval_ == 1 )			// expand to bit_depth == 8
    for( int row = 0; row < height(); ++row )
      for( int col = 0; col < width(); ++col )
        png_pixels[idx++] = data[row][col] ? 255 : 0;
  else						// bit_depth == 8
    for( int row = 0; row < height(); ++row )
      for( int col = 0; col < width(); ++col )
        png_pixels[idx++] = data[row][col];
  for( int i = 0; i < height(); ++i )
    row_pointers[i] = png_pixels + ( i * row_bytes );

  /* Create and initialize the png_struct with the desired error handler
   * functions.  If you want to use the default stderr and longjump method,
   * you can supply NULL for the last three parameters.  We also check that
   * the library version is compatible with the one used at compile time,
   * in case we are using dynamically linked libraries.  REQUIRED.
   */
  png_structp png_ptr =
    png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
  if( !png_ptr ) { std::free( row_pointers ); std::free( png_pixels );
                   throw std::bad_alloc(); }

  /* Allocate/initialize the image information data.  REQUIRED. */
  png_infop info_ptr = png_create_info_struct( png_ptr );
  if( !info_ptr )
    {
    png_destroy_write_struct( &png_ptr,  NULL );
    std::free( row_pointers ); std::free( png_pixels ); throw std::bad_alloc();
    }

  /* Set up error handling. */
  if( setjmp( png_jmpbuf( png_ptr ) ) )
    {
    /* If we get here, we had a problem writing the file. */
    png_destroy_write_struct( &png_ptr, &info_ptr );
    std::free( row_pointers ); std::free( png_pixels ); throw std::bad_alloc();
    }

  /* Set up the output control if you are using standard C streams. */
  png_init_io( png_ptr, f );

  png_set_IHDR( png_ptr, info_ptr, width(), height(), bit_depth,
                PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE );

  png_set_rows( png_ptr, info_ptr, row_pointers );

  // write the PNG image
  png_write_png( png_ptr, info_ptr, PNG_TRANSFORM_PACKING, NULL );

  /* Clean up after the write, and free any allocated memory. */
  png_destroy_write_struct( &png_ptr, &info_ptr );
  std::free( row_pointers );
  std::free( png_pixels );
  }


void show_png_info( FILE * const f, const char * const input_filename,
                    const int sig_read )
  {
  if( verbosity >= 0 ) std::printf( "%s\n", input_filename );
  png_structp png_ptr =
    png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL, NULL, NULL );
  if( !png_ptr ) throw std::bad_alloc();

  png_infop info_ptr = png_create_info_struct( png_ptr );
  if( !info_ptr )
    {
    png_destroy_read_struct( &png_ptr, NULL, NULL );	// avoid memory leak
    throw std::bad_alloc();
    }

  if( setjmp( png_jmpbuf( png_ptr ) ) )		// Set error handling
    {
    /* Free all of the memory associated with the png_ptr and info_ptr. */
    png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
    throw Page_image::Error( "Error reading PNG image." );
    }

  /* Set up the input control if you are using standard C streams. */
  png_init_io( png_ptr, f );

  /* If we have already read some of the signature */
  png_set_sig_bytes( png_ptr, sig_read );

  png_read_info( png_ptr, info_ptr );		// read image info

  /* Now let's print the data. */

  const unsigned height     = png_get_image_height( png_ptr, info_ptr );
  const unsigned width      = png_get_image_width( png_ptr, info_ptr );
  const long size           = height * width;
  const unsigned bit_depth  = png_get_bit_depth( png_ptr, info_ptr );
  const unsigned color_type = png_get_color_type( png_ptr, info_ptr );
  const unsigned channels   = png_get_channels( png_ptr, info_ptr );
  const unsigned interlace_type = png_get_interlace_type( png_ptr, info_ptr );
  const char * ct;
  if( color_type == PNG_COLOR_TYPE_GRAY_ALPHA ) ct = "greyscale with alpha channel";
  else if( color_type == PNG_COLOR_TYPE_GRAY )      ct = "greyscale";
  else if( color_type == PNG_COLOR_TYPE_PALETTE )   ct = "colormap";
  else if( color_type == PNG_COLOR_TYPE_RGB )       ct = "RGB";
  else if( color_type == PNG_COLOR_TYPE_RGB_ALPHA ) ct = "RGB with alpha channel";
  else if( color_type == PNG_COLOR_MASK_PALETTE )   ct = "mask colormap";
  else if( color_type == PNG_COLOR_MASK_COLOR )     ct = "mask color";
  else if( color_type == PNG_COLOR_MASK_ALPHA )     ct = "mask alpha";
  else ct = "unknown color_type";

  if( verbosity >= 0 )
    std::printf( "  PNG image %4u x %4u (%5.2f megapixels), "
                 "%2u-bit %s, %u channel(s), %sinterlaced\n",
                 width, height, size / 1000000.0, bit_depth, ct, channels,
                 ( interlace_type == PNG_INTERLACE_NONE ) ? "non-" : "" );

  png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
  }
