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

struct OCRAD_Pixmap;
class Mask;
class Rational;
class Rectangle;
class Track;
class Transformation;

class Page_image
  {
public:
  struct Error
    {
    const char * const msg;
    explicit Error( const char * const s ) : msg( s ) {}
    };

private:
  std::vector< std::vector< uint8_t > > data;	// 256 level greymap
  uint8_t maxval_, threshold_;			// x > threshold == white

  static void test_size( const int rows, const int cols );
  void read_p1( FILE * const f, const int cols, const bool invert );
  void read_p4( FILE * const f, const int cols, const bool invert );
  void read_p2( FILE * const f, const int cols, const bool invert );
  void read_p5( FILE * const f, const int cols, const bool invert );
  void read_p3( FILE * const f, const int cols, const bool invert );
  void read_p6( FILE * const f, const int cols, const bool invert );
  void read_png( FILE * const f, const int sig_read, const bool invert );
  void write_png( FILE * const f, const unsigned bit_depth ) const;

public:
  // Create a Page_image from a png, pbm, pgm, or ppm file
  Page_image( FILE * const f, const bool invert );

  // Create a Page_image from a OCRAD_Pixmap
  Page_image( const OCRAD_Pixmap & image, const bool invert );

  // Create a reduced Page_image
  Page_image( const Page_image & source, const int scale );

  bool get_bit( const int row, const int col ) const
    { return data[row][col] <= threshold_; }
  bool get_bit( const int row, const int col, const uint8_t th ) const
    { return data[row][col] <= th; }
  void set_bit( const int row, const int col, const bool bit )
    { data[row][col] = ( bit ? 0 : maxval_ ); }

  int height() const { return data.size(); }
  int width() const { return data.empty() ? 0 : data[0].size(); }
  uint8_t maxval() const { return maxval_; }
  uint8_t threshold() const { return threshold_; }
  void threshold( const Rational & th );	// 0 <= th <= 1, else auto
  void threshold( const int th );		// 0 <= th <= 255, else auto

  bool cut( const Rational ltwh[4] );
  void draw_mask( const Mask & m );
  void draw_rectangle( const Rectangle & re );
  void draw_track( const Track & tr );
  bool save( FILE * const f, const char filetype ) const;
  bool change_scale( int n );
  void transform( const Transformation & t );
  };


// defined in png_io.cc
bool read_check_png_sig8( FILE * const f, const int first_byte = -1 );
void show_png_info( FILE * const f, const char * const input_filename,
                    const int sig_read );
