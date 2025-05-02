/* Ocrcheck - A test program for the library ocradlib
   Copyright (C) 2009-2024 Antonio Diaz Diaz.

   This program is free software: you have unlimited permission to
   copy, distribute, and modify it.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>
#include <unistd.h>

#include "arg_parser.h"
#include "ocradlib.h"


namespace {

const char * const program_name = "ocradcheck";
const char * invocation_name = program_name;		// default value


void show_help()
  {
  std::printf( "Ocradcheck is a test program for the library ocradlib. It reads the image\n"
               "files specified, feeds them to the OCR engine, and sends the resulting text\n"
               "to stdout.\n"
               "\nUsage: %s [options] [files]\n", invocation_name );
  std::printf( "\nOptions:\n"
               "  -h, --help                display this help and exit\n"
               "  -V, --version             output version information and exit\n"
               "  -i, --invert              invert image levels (white on black)\n"
               "  -l, --layout              perform layout analysis\n"
               "  -u, --utf8                output text in UTF-8 format [default byte]\n"
               "\nIf no files are specified, or if a file is '-', ocradcheck reads the image\n"
               "from standard input.\n" );
  }


void show_version()
  {
  std::printf( "%s %s\n", program_name, PROGVERSION );
  std::printf( "Using ocradlib %s\n", OCRAD_version() );
  }


void show_error( const char * const msg, const int errcode = 0,
                 const bool help = false )
  {
  if( msg && msg[0] )
    std::fprintf( stderr, "%s: %s%s%s\n", program_name, msg,
                  ( errcode > 0 ) ? ": " : "",
                  ( errcode > 0 ) ? std::strerror( errcode ) : "" );
  if( help )
    std::fprintf( stderr, "Try '%s --help' for more information.\n",
                  invocation_name );
  }

} // end namespace


int main( const int argc, const char * const argv[] )
  {
  bool invert = false;
  bool layout = false;
  bool utf8 = false;
  if( argc > 0 ) invocation_name = argv[0];

  const Arg_parser::Option options[] =
    {
    { 'h', "help",        Arg_parser::no  },
    { 'i', "invert",      Arg_parser::no  },
    { 'l', "layout",      Arg_parser::no  },
    { 'u', "utf8",        Arg_parser::no  },
    { 'V', "version",     Arg_parser::no  },
    {  0 , 0,             Arg_parser::no  } };

  const Arg_parser parser( argc, argv, options );
  if( parser.error().size() )				// bad option
    { show_error( parser.error().c_str(), 0, true ); return 1; }

  int argind = 0;
  for( ; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    if( !code ) break;					// no more options
    switch( code )
      {
      case 'h': show_help(); return 0;
      case 'i': invert = true; break;
      case 'l': layout = true; break;
      case 'u': utf8 = true; break;
      case 'V': show_version(); return 0;
      default: std::fprintf( stderr, "%s: internal error: uncaught option.\n",
                             program_name ); return 3;
      }
    } // end process options

  if( OCRAD_version()[0] != OCRAD_version_string[0] )
    { std::fputs( "wrong library version.\n", stderr ); return 3; }

  if( std::strcmp( PROGVERSION, OCRAD_version_string ) != 0 )
    { std::fputs( "wrong library version_string.\n", stderr ); return 3; }

  // process any remaining command-line arguments (input files)
  bool stdin_used = false;
  for( bool first = true; first || argind < parser.arguments(); first = false )
    {
    const char * infile_name;
    if( argind < parser.arguments() )
      infile_name = parser.argument( argind++ ).c_str();
    else
      { infile_name = "-"; if( stdin_used ) continue; else stdin_used = true; }

    OCRAD_Descriptor * const ocrdes = OCRAD_open();
    if( !ocrdes || OCRAD_get_errno( ocrdes ) != OCRAD_ok )
      {
      OCRAD_close( ocrdes );
      std::fputs( "Not enough memory.\n", stderr );
      return 1;
      }

    if( OCRAD_set_image_from_file( ocrdes, infile_name, invert ) < 0 )
      {
      const OCRAD_Errno ocr_errno = OCRAD_get_errno( ocrdes );
      OCRAD_close( ocrdes );
      if( ocr_errno == OCRAD_mem_error )
        std::fputs( "Not enough memory.\n", stderr );
      else
        std::fprintf( stderr, "%s: Can't open file for reading.\n", infile_name );
      return 1;
      }
//    std::fprintf( stderr, "ocradcheck: testing file '%s'\n", infile_name );

    if( ( utf8 && OCRAD_set_utf8_format( ocrdes, true ) < 0 ) ||
        OCRAD_set_threshold( ocrdes, -1 ) < 0 ||	// auto threshold
        OCRAD_recognize( ocrdes, layout ) < 0 )
      {
      const OCRAD_Errno ocr_errno = OCRAD_get_errno( ocrdes );
      OCRAD_close( ocrdes );
      if( ocr_errno == OCRAD_mem_error )
        { std::fputs( "Not enough memory.\n", stderr ); return 1; }
      std::fprintf( stderr, "%s: internal error: invalid argument.\n",
                    program_name ); return 3;
      }

    const int blocks = OCRAD_result_blocks( ocrdes );
    int chars_total_by_block = 0;
    int chars_total_by_line = 0;
    int chars_total_by_count = 0;
    for( int b = 0; b < blocks; ++b )
      {
      const int lines = OCRAD_result_lines( ocrdes, b );
      chars_total_by_block += OCRAD_result_chars_block( ocrdes, b );
      for( int l = 0; l < lines; ++l )
        {
        const char * const s = OCRAD_result_line( ocrdes, b, l );
        chars_total_by_line += OCRAD_result_chars_line( ocrdes, b, l );
        if( s && s[0] )
          {
          std::fputs( s, stdout );
          const int len = std::strlen( s ) - 1;
          if( !utf8 )
            chars_total_by_count += len;
          else
            for( int i = 0; i < len; ++i )
              if( (uint8_t)s[i] < 128 || (uint8_t)s[i] >= 0xC0 )
                ++chars_total_by_count;
          }
        }
      std::fputc( '\n', stdout );
      }
    const int chars_total = OCRAD_result_chars_total( ocrdes );
    if( chars_total_by_block != chars_total ||
        chars_total_by_line != chars_total ||
        chars_total_by_count != chars_total )
      {
      std::fprintf( stderr, "library_error: character counts differ.\n"
                    "%d  %d  %d  %d\n", chars_total, chars_total_by_block,
                    chars_total_by_line, chars_total_by_count );
      return 1;
      }
    OCRAD_close( ocrdes );
    }
  return 0;
  }
