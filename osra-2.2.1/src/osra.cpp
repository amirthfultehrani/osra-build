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
#include <string.h> // strncpy()
#include <libgen.h> // dirname()
#include <pthread.h>
#include <signal.h>
#include <unistd.h>

#include <filesystem>
#include <chrono>
#include <future>
#include <csignal>
#include <shared_mutex>

#include <openbabel/oberror.h>
#include <tclap/CmdLine.h>

#include "ctpl_stl.h"
#include "osra_lib.h"
#include "osra.h"
#include "config.h" // PACKAGE_VERSION

#ifndef SIGUSR1
#define SIGUSR1 16
#endif

// This will hide the output "Warning: non-positive median line gap" from GOCR. Remove after this is fixed:
#ifdef _WIN32
#include <windows.h>
#define NULL_DEVICE "NUL:"
#else
#define NULL_DEVICE "/dev/null"
#endif

#ifndef SCHED_IDLE
#define SCHED_IDLE SCHED_OTHER
#endif

void signal_handler(int /*signum*/)
{
  std::this_thread::sleep_for(std::chrono::hours(10 * 365 * 24));
}

int main(int argc,
         char **argv
        )
{
  TCLAP::CmdLine cmd("OSRA: Optical Structure Recognition Application, created by Igor Filippov, 2013", ' ',
                     PACKAGE_VERSION);

  //
  // Image pre-processing options
  //
  TCLAP::ValueArg<double> rotate_option("R", "rotate", "Rotate image clockwise by specified number of degrees", false, 0,
                                        "0..360");
  cmd.add(rotate_option);

  TCLAP::SwitchArg invert_option("n", "negate", "Invert color (white on black)", false);
  cmd.add(invert_option);

  TCLAP::ValueArg<int> resolution_option("r", "resolution", "Resolution in dots per inch", false, 0, "default: auto");
  cmd.add(resolution_option);

  TCLAP::ValueArg<int> pdf_resolution_option("", "pdf", "Resolution in dots per inch for PDF rendering", false, 300, "default: 300");
  cmd.add(pdf_resolution_option);

  TCLAP::ValueArg<double> threshold_option("t", "threshold", "Gray level threshold", false, 0, "0.2..0.8");
  cmd.add(threshold_option);

  TCLAP::ValueArg<int> do_unpaper_option("u", "unpaper", "Pre-process image with unpaper algorithm, rounds", false, 0,
                                         "default: 0 rounds");
  cmd.add(do_unpaper_option);

  TCLAP::SwitchArg jaggy_option("j", "jaggy", "Additional thinning/scaling down of low quality documents", false);
  cmd.add(jaggy_option);

  TCLAP::SwitchArg adaptive_option("i", "adaptive", "Adaptive thresholding pre-processing, useful for low light/low contrast images", false);
  cmd.add(adaptive_option);

  TCLAP::SwitchArg keep_option("k", "keep", "Keep image unsegmented, do not separate molecules from text", false);
  cmd.add(keep_option);

  TCLAP::ValueArg<int> timeout_option("", "timeout", "Timeout in seconds per file processing", false, 0, "default: no timeout");
  cmd.add(timeout_option);

  //
  // Output format options
  //
  TCLAP::ValueArg<std::string> output_format_option("f", "format", "Output format", false, "can", "can/smi/sdf");
  cmd.add(output_format_option);

  TCLAP::SwitchArg v3000_option("", "v3000", "Use V3000 format for MDL MOL and SDF output", false);
  cmd.add(v3000_option);

  TCLAP::ValueArg<std::string> embedded_format_option("", "embedded-format", "Embedded format", false, "", "inchi/smi/can");
  cmd.add(embedded_format_option);

  TCLAP::SwitchArg show_confidence_option("p", "print", "Print out confidence estimate", false);
  cmd.add(show_confidence_option);

  TCLAP::SwitchArg show_resolution_guess_option("g", "guess", "Print out resolution guess", false);
  cmd.add(show_resolution_guess_option);

  TCLAP::SwitchArg show_page_option("e", "page", "Show page number for PDF/PS/TIFF documents (only for SDF/SMI/CAN output format)", false);
  cmd.add(show_page_option);

  TCLAP::SwitchArg show_coordinates_option("c", "coordinates", "Show surrounding box coordinates (only for SDF/SMI/CAN output format)", false);
  cmd.add(show_coordinates_option);

  TCLAP::SwitchArg show_avg_bond_length_option("b", "bond", "Show average bond length in pixels (only for SDF/SMI/CAN output format)", false);
  cmd.add(show_avg_bond_length_option);

  //
  // Dictionaries options
  //
  TCLAP::ValueArg<std::string> spelling_file_option("l", "spelling", "Spelling correction dictionary", false, "", "configfile");
  cmd.add(spelling_file_option);

  TCLAP::ValueArg<std::string> superatom_file_option("a", "superatom", "Superatom label map to SMILES", false, "", "configfile");
  cmd.add(superatom_file_option);

  TCLAP::ValueArg<std::string> recognized_chars_option("", "ocr", "OCR character filter", false, "", "oOcCnNHFsSBuUgMeEXYZRPp23456789AmThDGQ");
  cmd.add(recognized_chars_option);

  //
  // Debugging options
  //
  TCLAP::SwitchArg debug_option("d", "debug", "Print out debug information on spelling corrections", false);
  cmd.add(debug_option);

  TCLAP::SwitchArg verbose_option("v", "verbose", "Be verbose and print the program flow", false);
  cmd.add(verbose_option);

  TCLAP::ValueArg<std::string> output_image_file_prefix_option("o", "output", "Write recognized structures to image files with given prefix", false, "", "filename prefix");
  cmd.add(output_image_file_prefix_option);

  TCLAP::ValueArg<std::string> resize_option("s", "size", "Resize image on output", false, "", "dimensions, 300x400");
  cmd.add(resize_option);

  TCLAP::ValueArg<std::string> preview_option("", "preview", "Preview Image", false, "", "filename");
  cmd.add(preview_option);
  //
  // Input-output options
  //
  TCLAP::ValueArg<std::string> output_file_option("w", "write", "Write recognized structures to a file or folder", false, "", "filename");
  cmd.add(output_file_option);

  TCLAP::SwitchArg show_learning_option("", "learn", "Print out all structure guesses with confidence parameters", false);
  cmd.add(show_learning_option);
  
  TCLAP::UnlabeledMultiArg<std::string> input_file_option("in", "input file(s) or a single folder", true, "", "filename(s)");
  cmd.add(input_file_option);
  
  cmd.parse(argc, argv);

  freopen(NULL_DEVICE, "w", stderr);
  OpenBabel::obErrorLog.StopLogging();
  
  // Calculating the current dir:
  char progname[1024];
  strncpy(progname, cmd.getProgramName().c_str(), sizeof(progname) - 1);
  progname[sizeof(progname) - 1] = '\0';
  std::string osra_dir = dirname(progname);

  int result = 0;
  std::vector<std::string> input_files = input_file_option.getValue();
  std::string output_file = output_file_option.getValue();
  std::filesystem::path out_path(output_file); 

  if (input_files.size() != 1 && output_file != "" && !std::filesystem::exists(out_path))
    {
      std::filesystem::create_directories(out_path);
    }

  if (input_files.size() == 1 && std::filesystem::is_directory(input_files[0]))
    {
      std::filesystem::path in_path(input_files[0]);
      input_files.clear();
      for (auto const& dir_entry : std::filesystem::directory_iterator(in_path))
	if (dir_entry.is_regular_file())
	  {
	    input_files.emplace_back(dir_entry.path().string());
	  }
    }
  
  std::vector<std::string> output_files(input_files.size(), output_file);
  if (std::filesystem::is_directory(out_path))
    { 
      output_files.clear();     
      std::string ext = output_format_option.getValue();
      std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
      for (const auto &input_file : input_files)
	{
	  std::filesystem::path p(out_path);
	  p /= std::filesystem::path(input_file).filename();
	  p.replace_extension(ext);
	  output_files.push_back(p.string());
	}
    }
  else if (input_files.size() != 1 && output_file != "")
    {
      std::cout << "Multiple input files given, but a single regular file as an output" << std::endl;
      return ERROR_ILLEGAL_ARGUMENT_COMBINATION;
    }

  std::string output_format =  output_format_option.getValue();					     
  std::string embedded_format = embedded_format_option.getValue();
      
  std::transform(output_format.begin(), output_format.end(), output_format.begin(), ::tolower);
  std::transform(embedded_format.begin(), embedded_format.end(), embedded_format.begin(), ::tolower);

  std::map<std::string, std::string> spelling, superatom;
  int err = load_superatom_spelling_maps(spelling, superatom, osra_dir, spelling_file_option.getValue(), superatom_file_option.getValue(), verbose_option.getValue());
  if (err != 0) return err;

  int timeout = timeout_option.getValue();

  auto rotate = rotate_option.getValue();
  auto invert = invert_option.getValue();
  auto resolution = resolution_option.getValue();
  auto pdf_resolution = pdf_resolution_option.getValue();
  auto threshold = threshold_option.getValue();
  auto do_unpaper = do_unpaper_option.getValue();
  auto jaggy = jaggy_option.getValue();
  auto adaptive = adaptive_option.getValue();
  auto keep = keep_option.getValue();
  auto show_confidence = show_confidence_option.getValue();
  auto show_resolution_guess = show_resolution_guess_option.getValue();
  auto show_page = show_page_option.getValue();
  auto show_coordinates = show_coordinates_option.getValue();
  auto show_avg_bond_length = show_avg_bond_length_option.getValue();
  auto show_learning = show_learning_option.getValue();
  auto debug = debug_option.getValue();
  auto verbose = verbose_option.getValue();
  auto output_image_file_prefix = output_image_file_prefix_option.getValue();
  auto resize = resize_option.getValue();
  auto preview = preview_option.getValue();
  auto recognized_chars = recognized_chars_option.getValue();
  auto v3000 = v3000_option.getValue();

  if (timeout == 0)
    {
      for (size_t file_counter = 0; file_counter < input_files.size(); ++file_counter)
	{
	  const auto &input_file = input_files[file_counter];
	  const auto &output_file = output_files[file_counter];
	  int single_result = osra_process_image(
						 input_file,
						 output_file,
						 rotate,
						 invert,
						 resolution,
						 pdf_resolution,
						 threshold,
						 do_unpaper,
						 jaggy,
						 adaptive,
						 keep,
						 output_format,
						 embedded_format,
						 show_confidence,
						 show_resolution_guess,
						 show_page,
						 show_coordinates,
						 show_avg_bond_length,
						 show_learning,
						 spelling,
						 superatom,
						 debug,
						 verbose,
						 output_image_file_prefix,
						 resize,
						 preview,
						 recognized_chars,
						 v3000);
	  result = std::max(result, single_result);    
	}
    }
  else
    {
#ifndef _WIN32
      struct sigaction action;
      action.sa_handler = signal_handler;
      sigemptyset(&action.sa_mask);
      action.sa_flags = 0;
      sigaction(SIGUSR1, &action, NULL);
#endif
      std::vector< std::future<int> > futures;
      ctpl::thread_pool p(1);   
      struct sched_param param;
      param.sched_priority = 0;  
      size_t start = 0;
      while (start < input_files.size())
	{
	  for (size_t file_counter = start; file_counter < input_files.size(); ++file_counter)
	    {
	      const auto &input_file = input_files[file_counter];
	      const auto &output_file = output_files[file_counter];
	      futures.emplace_back( p.push([input_file,
					    output_file,
					    rotate,
					    invert,
					    resolution,
					    pdf_resolution,
					    threshold,
					    do_unpaper,
					    jaggy,
					    adaptive,
					    keep,
					    output_format,
					    embedded_format,
					    show_confidence,
					    show_resolution_guess,
					    show_page,
					    show_coordinates,
					    show_avg_bond_length,
					    show_learning,
					    spelling,
					    superatom,
					    debug,
					    verbose,
					    output_image_file_prefix,
					    resize,
					    preview,
					    recognized_chars,
					    v3000]
					   (int) {
					     return  osra_process_image(
					   input_file,
					   output_file,
					   rotate,
					   invert,
					   resolution,
					   pdf_resolution,
					   threshold,
					   do_unpaper,
					   jaggy,
					   adaptive,
					   keep,
					   output_format,
					   embedded_format,
					   show_confidence,
					   show_resolution_guess,
					   show_page,
					   show_coordinates,
					   show_avg_bond_length,
					   show_learning,
					   spelling,
					   superatom,
					   debug,
					   verbose,
					   output_image_file_prefix,
					   resize,
					   preview,
					   recognized_chars,
					   v3000);
	      }));
	    }
	  for (auto &f : futures)
	    {
	      ++start;
	      if (f.wait_for(std::chrono::seconds(timeout)) == std::future_status::ready)
		{
		  result = std::max(result, f.get());
		}
	      else
		{
		  auto nh = p.get_thread(0).native_handle();
#ifdef _WIN32
		  HANDLE win_nh =  pthread_gethandle(nh);
		  if (win_nh)
		    SuspendThread(win_nh);
#else
		  if (typeid(nh) == typeid(pthread_t))
		    pthread_setschedparam(nh, SCHED_IDLE, &param);
#endif
		  p.resize(0);
		  p.clear_queue();
		  futures.clear();
		  if (typeid(nh) == typeid(pthread_t))
		    pthread_kill(nh, SIGUSR1);

		  std::cout << "Timeout reached: " << input_files[start - 1] << std::endl;	  
		  p.resize(1);
		  break;
		}
	    }
	}
      std::raise(SIGTERM);
    }

  return result;
}
