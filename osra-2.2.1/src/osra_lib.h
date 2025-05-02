/******************************************************************************
 OSRA: Optical Structure Recognition Application

 Created by Igor Filippov, 2007-2013 (igor.v.filippov@gmail.com)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

// Header: osra_lib.h
//
// Defines types and functions of OSRA library.
//

#include <string> // std::string
#include <ostream> // std:ostream
#include <map>

//
// Section: Functions
//

// Function: osra_process_image()
//
// Parameters:
//      image_data - the binary image
//
// Returns:
//      0, if processing was completed successfully
int osra_process_image(
#ifdef OSRA_LIB
  const char *image_data,
  int image_length,
  std::ostream &structure_output_stream,
#else
  const std::string input_file,
  const std::string output_file,
#endif
  int rotate,
  bool invert,
  int input_resolution,
  int pdf_input_resolution,
  double threshold,
  int do_unpaper,
  bool jaggy,
  bool adaptive,
  bool keep_option,
  const std::string output_format,
  const std::string embedded_format,
  bool show_confidence,
  bool show_resolution_guess,
  bool show_page,
  bool show_coordinates,
  bool show_avg_bond_length,
  bool show_learning,
  const std::map<std::string, std::string> spelling,
  const std::map<std::string, std::string> superatom,
  bool debug,
  bool verbose,
  const std::string output_image_file_prefix,
  const std::string resize,
  const std::string preview,
  const std::string recognized_chars,
  bool v3000
);

int load_superatom_spelling_maps(
				 std::map<std::string, std::string> &spelling,
				 std::map<std::string, std::string> &superatom,
				 const std::string &osra_dir,
				 const std::string &spelling_file,
				 const std::string &superatom_file,
				 bool verbose);
enum class eProcessing
{
  Normal,
  NoThinning,
  Scaling
};
