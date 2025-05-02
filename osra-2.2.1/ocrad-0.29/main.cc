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
/*
   Exit status: 0 for a normal exit, 1 for environmental problems
   (file not found, invalid command-line options, I/O errors, etc), 2 to
   indicate a corrupt or invalid input file, 3 for an internal consistency
   error (e.g., bug) which caused ocrad to panic.
*/

#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <string>
#include <vector>
#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#if defined __MSVCRT__ || defined __OS2__ || defined __DJGPP__
#include <fcntl.h>
#include <io.h>
#endif

#include "arg_parser.h"
#include "common.h"
#include "rational.h"
#include "user_filter.h"
#include "page_image.h"
#include "textpage.h"


namespace {

const char * const program_name = "ocrad";
const char * const program_year = "2024";
const char * invocation_name = program_name;		// default value


void show_help()
  {
  std::printf( "GNU Ocrad is an OCR (Optical Character Recognition) program based on a\n"
               "feature extraction method. It reads images in png or pnm formats and\n"
               "produces text in byte (8-bit) or UTF-8 formats. The formats pbm (bitmap),\n"
               "pgm (greyscale), and ppm (color) are collectively known as pnm.\n"
               "\nOcrad includes a layout analyser able to separate the columns or blocks\n"
               "of text normally found on printed pages.\n"
               "\nFor best results the characters should be at least 20 pixels high. If\n"
               "they are smaller, try the option --scale. Scanning the image at 300 dpi\n"
               "usually produces a character size good enough for ocrad.\n"
               "Merged, very bold, or very light (broken) characters are normally not\n"
               "recognized correctly. Try to avoid them.\n"
               "\nUsage: %s [options] [files]\n", invocation_name );
  std::printf( "\nOptions:\n"
               "  -h, --help                display this help and exit\n"
               "  -V, --version             output version information and exit\n"
               "  -a, --append              append text to output file\n"
               "  -c, --charset=<name>      try '--charset=help' for a list of names\n"
               "  -e, --filter=<name>       try '--filter=help' for a list of names\n"
               "  -E, --user-filter=<file>  user-defined filter, see manual for format\n"
               "  -f, --force               force overwrite of output file\n"
               "  -F, --format=<fmt>        output format (byte, utf8)\n"
               "  -i, --invert              invert image levels (white on black)\n"
               "  -l, --layout              perform layout analysis\n"
               "  -o, --output=<file>       place the output into <file>\n"
               "  -q, --quiet               suppress all messages\n"
               "  -s, --scale=[-]<n>        scale input image by [1/]<n>\n"
               "  -t, --transform=<name>    try '--transform=help' for a list of names\n"
               "  -T, --threshold=<n>       threshold for binarization (0.0-1.0)\n"
               "  -u, --cut=<l,t,w,h>       cut the input image by the rectangle given\n"
               "  -v, --verbose             be verbose\n"
               "  -x, --export=<file>       export results in ORF format to <file>\n" );
  if( verbosity >= 1 )
    {
    std::printf( "  -1..8                     pnm or png output file type (debug)\n"
                 "  -C, --copy                copy input to output (debug)\n"
                 "  -D, --debug=<level>       (0-100) output intermediate data (debug)\n"
                 "  -I, --info                show png image information\n" );
    }
  std::printf( "\nIf no files are specified, or if a file is '-', ocrad reads the image\n"
               "from standard input.\n"
               "If the option -o is not specified, ocrad sends text to standard output.\n"
               "\nExit status: 0 for a normal exit, 1 for environmental problems\n"
               "(file not found, invalid command-line options, I/O errors, etc), 2 to\n"
               "indicate a corrupt or invalid input file, 3 for an internal consistency\n"
               "error (e.g., bug) which caused ocrad to panic.\n"
               "\nReport bugs to bug-ocrad@gnu.org\n"
               "Ocrad home page: http://www.gnu.org/software/ocrad/ocrad.html\n"
               "General help using GNU software: http://www.gnu.org/gethelp\n" );
  }


void show_version()
  {
  std::printf( "GNU %s %s\n", program_name, PROGVERSION );
  std::printf( "Copyright (C) %s Antonio Diaz Diaz.\n", program_year );
  std::printf( "License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n"
               "This is free software: you are free to change and redistribute it.\n"
               "There is NO WARRANTY, to the extent permitted by law.\n" );
  }


struct Input_control
  {
  Transformation transformation;
  int scale;
  Rational threshold, ltwh[4];
  bool copy, cut, invert, layout;

  Input_control()
    : scale( 1 ), threshold( -1 ),
      copy( false ), cut( false ), invert( false ), layout( false ) {}

  void parse_cut_rectangle( const char * const arg, const char * const pn );
  void parse_scale( const char * const arg, const char * const pn );
  void parse_threshold( const char * const arg, const char * const pn );
  };


void show_error( const char * const msg, const int errcode = 0,
                 const bool help = false )
  {
  if( verbosity < 0 ) return;
  if( msg && msg[0] )
    std::fprintf( stderr, "%s: %s%s%s\n", program_name, msg,
                  ( errcode > 0 ) ? ": " : "",
                  ( errcode > 0 ) ? std::strerror( errcode ) : "" );
  if( help )
    std::fprintf( stderr, "Try '%s --help' for more information.\n",
                  invocation_name );
  }


void show_file_error( const char * const filename, const char * const msg,
                      const int errcode = 0 )
  {
  if( verbosity >= 0 )
    std::fprintf( stderr, "%s: %s: %s%s%s\n", program_name, filename, msg,
                  ( errcode > 0 ) ? ": " : "",
                  ( errcode > 0 ) ? std::strerror( errcode ) : "" );
  }


void internal_error( const char * const msg )
  {
  if( verbosity >= 0 )
    std::fprintf( stderr, "%s: internal error: %s\n", program_name, msg );
  std::exit( 3 );
  }


void show_option_error( const char * const arg, const char * const msg,
                        const char * const option_name )
  {
  if( verbosity >= 0 )
    std::fprintf( stderr, "%s: '%s': %s option '%s'.\n",
                  program_name, arg, msg, option_name );
  }


void Input_control::parse_cut_rectangle( const char * const arg, const char * const pn )
  {
  int c = ltwh[0].parse( arg );				// left
  if( c && arg[c] == ',' && ltwh[0] >= -1 )
    {
    int i = c + 1;
    c = ltwh[1].parse( &arg[i] );			// top
    if( c && arg[i+c] == ',' && ltwh[1] >= -1 )
      {
      i += c + 1; c = ltwh[2].parse( &arg[i] );		// width
      if( c && arg[i+c] == ',' && ltwh[2] > 0 )
        {
        i += c + 1; c = ltwh[3].parse( &arg[i] );	// height
        if( c && ltwh[3] > 0 ) { cut = true; return; }
        }
      }
    }
  show_option_error( arg, "Invalid cut rectangle in", pn );
  std::exit( 1 );
  }


void Input_control::parse_scale( const char * const arg, const char * const pn )
  {
  errno = 0;
  const long n = std::strtol( arg, 0, 0 );
  if( errno || n == 0 || n > INT_MAX || n < -INT_MAX )
    { show_option_error( arg, "Invalid scale factor in", pn ); std::exit( 1 ); }
  scale = n;
  }


void Input_control::parse_threshold( const char * const arg, const char * const pn )
  {
  Rational th;
  if( th.parse( arg ) && th >= 0 && th <= 1 ) { threshold = th; return; }
  show_option_error( arg, "Threshold out of limits [0.0, 1.0] in", pn );
  std::exit( 1 );
  }


bool make_dirs( const char * const name, const int len )
  {
  int i = len;
  while( i > 0 && name[i-1] != '/' ) --i;	// remove last component
  while( i > 0 && name[i-1] == '/' ) --i;	// remove slash(es)
  const int dirsize = i;	// size of dirname without trailing slash(es)

  for( i = 0; i < dirsize; )	// if dirsize == 0, dirname is '/' or empty
    {
    while( i < dirsize && name[i] == '/' ) ++i;
    const int first = i;
    while( i < dirsize && name[i] != '/' ) ++i;
    if( first < i )
      {
      const std::string partial( name, 0, i );
      const mode_t mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
      struct stat st;
      if( stat( partial.c_str(), &st ) == 0 )
        { if( !S_ISDIR( st.st_mode ) ) { errno = ENOTDIR; return false; } }
      else if( mkdir( partial.c_str(), mode ) != 0 && errno != EEXIST )
        return false;		// if EEXIST, another process created the dir
      }
    }
  return true;
  }


FILE * open_outstream( const char * const outfile_name, const bool append,
                       const bool force )
  {
  const int len = strlen( outfile_name );
  if( len > 0 && outfile_name[len-1] == '/' ) errno = EISDIR;
  else {
    if( !make_dirs( outfile_name, len ) )
      { show_file_error( outfile_name, "Error creating intermediate directory",
                         errno ); return 0; }
    const char * const mode = append ? "a" : force ? "w" : "wx";
    FILE * outfile = std::fopen( outfile_name, mode );
    if( outfile ) return outfile;
    if( errno == EEXIST )
      { show_file_error( outfile_name, "Output file already exists, skipping." );
        return 0; }
    }
  show_file_error( outfile_name, "Can't create output file", errno );
  return 0;
  }


const char * my_basename( const char * filename )
  {
  const char * c = filename;
  while( *c ) { if( *c == '/' ) { filename = c + 1; } ++c; }
  return filename;
  }


int process_file( FILE * const infile, const char * const infile_name,
                  const Input_control & input_control,
                  const Control & control )
  {
  Page_image page_image( infile, input_control.invert );

  if( input_control.cut )
    {
    if( !page_image.cut( input_control.ltwh ) )
      { show_file_error( infile_name, "File totally cut away." ); return 1; }
    if( verbosity >= 1 )
      std::fprintf( stderr, "%s: File cut to %dw x %dh\n", infile_name,
                    page_image.width(), page_image.height() );
    }

  page_image.transform( input_control.transformation );
  if( !page_image.change_scale( input_control.scale ) )
    { show_file_error( infile_name, "Scale factor too small for file." );
      return 1; }
  page_image.threshold( input_control.threshold );
  if( verbosity >= 1 )
    {
    const Rational th( page_image.threshold(), page_image.maxval() );
    std::fprintf( stderr, "maxval = %d, threshold = %d (%s)\n",
                  page_image.maxval(), page_image.threshold(),
                  th.to_decimal( 1, -3 ).c_str() );
    }

  if( input_control.copy )
    {
    if( control.outfile ) page_image.save( control.outfile, control.filetype );
    return 0;
    }

  Textpage textpage( page_image, my_basename( infile_name ), control,
                     input_control.layout );
  if( control.debug_level == 0 )
    {
    if( control.outfile ) textpage.print( control );
    if( control.exportfile ) textpage.xprint( control );
    }
  if( verbosity >= 1 ) std::fputc( '\n', stderr );
  return 0;
  }

int process_full_file( FILE * const infile, const char * const infile_name,
                       const Input_control & input_control,
                       const Control & control )
  {
  if( verbosity >= 1 )
    std::fprintf( stderr, "Processing file '%s'\n", infile_name );
  int retval = 0;
  while( true )
    {
    const int tmp = process_file( infile, infile_name, input_control, control );
    if( retval < tmp ) retval = tmp;
    if( control.outfile ) std::fflush( control.outfile );
    if( control.exportfile ) std::fflush( control.exportfile );
    if( infile != stdin ) break;
    if( tmp <= 1 )		// detect EOF
      {
      int ch;
      do ch = std::fgetc( infile ); while( ch == 0 || std::isspace( ch ) );
      std::ungetc( ch, infile );
      }
    if( tmp > 1 || std::feof( infile ) || std::ferror( infile ) ) break;
    }
  return retval;
  }

} // end namespace


void Charset::enable( const char * const arg, const char * const pn )
  {
  const Value charset_value[charsets] = { ascii, iso_8859_9, iso_8859_15 };
  const char * const charset_name[charsets] =
    { "ascii", "iso-8859-9", "iso-8859-15" };

  for( int i = 0; i < charsets; ++i )
    if( std::strcmp( arg, charset_name[i] ) == 0 )
      { charset_ |= charset_value[i]; return; }
  if( verbosity < 0 ) std::exit( 1 );
  if( arg && std::strcmp( arg, "help" ) != 0 )
    show_option_error( arg, "Unknown charset name in", pn );
  std::fputs( "  Valid charset names:", stderr );
  for( int i = 0; i < charsets; ++i )
    std::fprintf( stderr, "  %s", charset_name[i] );
  std::fputc( '\n', stderr );
  std::exit( 1 );
  }


void Transformation::set( const char * const arg, const char * const pn )
  {
  struct T_entry { const char * const name; Transformation::Type type; };
  const T_entry T_table[] = {
    { "none",      Transformation::none },
    { "rotate90",  Transformation::rotate90 },
    { "rotate180", Transformation::rotate180 },
    { "rotate270", Transformation::rotate270 },
    { "mirror_lr", Transformation::mirror_lr },
    { "mirror_tb", Transformation::mirror_tb },
    { "mirror_d1", Transformation::mirror_d1 },
    { "mirror_d2", Transformation::mirror_d2 },
    { 0, Transformation::none } };

  for( int i = 0; T_table[i].name != 0; ++i )
    if( std::strcmp( arg, T_table[i].name ) == 0 )
      { type_ = T_table[i].type; return; }
  if( verbosity < 0 ) std::exit( 1 );
  if( arg && std::strcmp( arg, "help" ) != 0 )
    show_option_error( arg, "Unknown bitmap trasformation name in", pn );
  std::fputs( "  Valid transformation names:", stderr );
  for( int i = 0; T_table[i].name != 0; ++i )
    std::fprintf( stderr, "  %s", T_table[i].name );
  std::fputs( "\n  Rotations are made counter-clockwise.\n", stderr );
  std::exit( 1 );
  }


void Control::add_filter( const char * const arg, const char * const pn )
  {
  struct F_entry { const char * const name; Filter::Type type; };
  const F_entry F_table[] = {
    { "letters",        Filter::letters },
    { "letters_only",   Filter::letters_only },
    { "numbers",        Filter::numbers },
    { "numbers_only",   Filter::numbers_only },
    { "same_height",    Filter::same_height },
    { "text_block",     Filter::text_block },
    { "upper_num",      Filter::upper_num },
    { "upper_num_mark", Filter::upper_num_mark },
    { "upper_num_only", Filter::upper_num_only },
    { 0, Filter::user } };

  for( int i = 0; F_table[i].name != 0; ++i )
    if( std::strcmp( arg, F_table[i].name ) == 0 )
      { filters.push_back( Filter( F_table[i].type ) ); return; }
  if( verbosity < 0 ) std::exit( 1 );
  if( arg && std::strcmp( arg, "help" ) != 0 )
    show_option_error( arg, "Unknown filter name in", pn );
  std::fputs( "  Valid filter names:", stderr );
  for( int i = 0; F_table[i].name != 0; ++i )
    std::fprintf( stderr, "  %s", F_table[i].name );
  std::fputc( '\n', stderr );
  std::exit( 1 );
  }


void Control::add_user_filter( const char * const file_name )
  {
  const User_filter * const user_filterp = new User_filter( file_name );
  const int retval = user_filterp->retval();
  if( retval == 0 ) { filters.push_back( Filter( user_filterp ) ); return; }
  if( user_filterp->error().empty() )
    show_file_error( file_name, "Can't open filter file", errno );
  else
    show_file_error( file_name, user_filterp->error().c_str() );
  delete user_filterp;
  std::exit( retval );
  }


void Control::set_format( const char * const arg, const char * const pn )
  {
  if( std::strcmp( arg, "byte" ) == 0 ) { utf8 = false; return; }
  if( std::strcmp( arg, "utf8" ) == 0 ) { utf8 = true; return; }
  show_option_error( arg, "Unknown format name in", pn );
  std::exit( 1 );
  }


/* 'infile' contains the scanned image to be converted to text.
   'outfile' is the destination for the text version of the scanned image
   (or for an image file if debugging).
   'exportfile' is the Ocr Results File.
*/
int main( const int argc, const char * const argv[] )
  {
  Input_control input_control;
  Control control;
  const char * outfile_name = 0, * exportfile_name = 0;
  bool append = false;
  bool force = false;
  bool info_mode = false;
  if( argc > 0 ) invocation_name = argv[0];

  const Arg_parser::Option options[] =
    {
    { '1', 0,             Arg_parser::no  },		// P1 to P6
    { '2', 0,             Arg_parser::no  },
    { '3', 0,             Arg_parser::no  },
    { '4', 0,             Arg_parser::no  },
    { '5', 0,             Arg_parser::no  },
    { '6', 0,             Arg_parser::no  },
    { '7', 0,             Arg_parser::no  },		// png mono
    { '8', 0,             Arg_parser::no  },		// png greyscale
    { 'a', "append",      Arg_parser::no  },
    { 'c', "charset",     Arg_parser::yes },
    { 'C', "copy",        Arg_parser::no  },
    { 'D', "debug",       Arg_parser::yes },
    { 'e', "filter",      Arg_parser::yes },
    { 'E', "user-filter", Arg_parser::yes },
    { 'f', "force",       Arg_parser::no  },
    { 'F', "format",      Arg_parser::yes },
    { 'h', "help",        Arg_parser::no  },
    { 'i', "invert",      Arg_parser::no  },
    { 'I', "info",        Arg_parser::no  },
    { 'l', "layout",      Arg_parser::no  },
    { 'o', "output",      Arg_parser::yes },
    { 'q', "quiet",       Arg_parser::no  },
    { 's', "scale",       Arg_parser::yes },
    { 't', "transform",   Arg_parser::yes },
    { 'T', "threshold",   Arg_parser::yes },
    { 'u', "cut",         Arg_parser::yes },
    { 'v', "verbose",     Arg_parser::no  },
    { 'V', "version",     Arg_parser::no  },
    { 'x', "export",      Arg_parser::yes },
    {  0 , 0,             Arg_parser::no  } };

  const Arg_parser parser( argc, argv, options );
  if( parser.error().size() )				// bad option
    { show_error( parser.error().c_str(), 0, true ); return 1; }

  int retval = 0;
  int argind = 0;
  for( ; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    if( !code ) break;					// no more options
    const char * const pn = parser.parsed_name( argind ).c_str();
    const char * const arg = parser.argument( argind ).c_str();
    switch( code )
      {
      case '1': case '2': case '3': case '4': case '5':
      case '6': case '7': case '8': control.filetype = code; break;
      case 'a': append = true; break;
      case 'c': control.charset.enable( arg, pn ); break;
      case 'C': input_control.copy = true; break;
      case 'D': control.debug_level = std::strtol( arg, 0, 0 ); break;
      case 'e': control.add_filter( arg, pn ); break;
      case 'E': control.add_user_filter( arg ); break;
      case 'f': force = true; break;
      case 'F': control.set_format( arg, pn ); break;
      case 'h': show_help(); return 0;
      case 'i': input_control.invert = true; break;
      case 'I': info_mode = true; break;
      case 'l': input_control.layout = true; break;
      case 'o': outfile_name = arg; break;
      case 'q': verbosity = -1; break;
      case 's': input_control.parse_scale( arg, pn ); break;
      case 't': input_control.transformation.set( arg, pn ); break;
      case 'T': input_control.parse_threshold( arg, pn ); break;
      case 'u': input_control.parse_cut_rectangle( arg, pn ); break;
      case 'v': if( verbosity < 4 ) ++verbosity; break;
      case 'V': show_version(); return 0;
      case 'x': exportfile_name = arg; break;
      default: internal_error( "uncaught option." );
      }
    } // end process options

#if defined __MSVCRT__ || defined __OS2__ || defined __DJGPP__
  setmode( STDIN_FILENO, O_BINARY );
  setmode( STDOUT_FILENO, O_BINARY );
#endif

  if( !info_mode && outfile_name && std::strcmp( outfile_name, "-" ) != 0 )
    {
    control.outfile = open_outstream( outfile_name, append, force );
    if( !control.outfile ) return 1;
    }

  if( !info_mode && exportfile_name && control.debug_level == 0 &&
      !input_control.copy )
    {
    if( std::strcmp( exportfile_name, "-" ) == 0 )
      { control.exportfile = stdout; if( !outfile_name ) control.outfile = 0; }
    else
      {
      control.exportfile = open_outstream( exportfile_name, false, true );
      if( !control.exportfile ) return 1;
      }
    std::fprintf( control.exportfile,
                  "# Ocr Results File. Created by GNU Ocrad version %s\n",
                  PROGVERSION );
    }

  // process any remaining command-line arguments (input files)
  bool stdin_used = false;
  for( bool first = true; first || argind < parser.arguments(); first = false )
    {
    const char * infile_name = ( argind < parser.arguments() ) ?
                               parser.argument( argind++ ).c_str() : "-";
    FILE * infile;
    if( std::strcmp( infile_name, "-" ) == 0 )
      {
      if( stdin_used ) continue; else stdin_used = true;
      infile = stdin;
      }
    else
      {
      infile = std::fopen( infile_name, "rb" );
      if( !infile )
        {
        show_file_error( infile_name, "Can't open input file", errno );
        if( retval < 1 ) retval = 1;
        continue;
        }
      }
    int tmp;
    try {
      if( info_mode )
        {
        if( !read_check_png_sig8( infile ) )
          { tmp = 2; show_file_error( infile_name, "Not a PNG file." ); }
        else { tmp = 0; show_png_info( infile, infile_name, 8 ); }
        }
      else
        tmp = process_full_file( infile, infile_name, input_control, control );
      }
    catch( std::bad_alloc & )
      { show_file_error( infile_name, "Not enough memory." ); tmp = 1; }
    catch( Page_image::Error & e )
      { show_file_error( infile_name, e.msg ); tmp = 2; }
    catch( Ocrad::Internal & e )
      { show_file_error( infile_name, e.msg ); std::exit( 3 ); }
    if( std::fclose( infile ) != 0 )
      { show_file_error( infile_name, "Error closing input file", errno );
        if( tmp < 1 ) tmp = 1; }
    if( retval < tmp ) retval = tmp;
    }
  if( control.outfile ) std::fclose( control.outfile );
  if( control.exportfile ) std::fclose( control.exportfile );
  return retval;
  }
