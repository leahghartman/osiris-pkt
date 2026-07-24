cmake_minimum_required(VERSION 3.19)

#
# TODO:
#
#   There is a way I found to 'maske' much faster... Go into the CmakeFiles directory at the file 'Makefile2'
#     If you take out att the lines that end \requires and save that file as Makefile_fdst
#     
#     I think you have to run 'make' the slow wasy once. Then 'make -f Makefile_fdst target'
#       Where target is 'CMakeFiles\osiris-2D.e.dir\all' (put in correct executable name/dimension as needed)
#
#     TODO: make a quick python script/utility to do this for ytou automatically.. its much faster.
#           or use Ninja.. its great but wierd and takes a little getting used to.
#


#------------------------------------------------------------------------------
# Establish the root of the Osiris repository we are using.
#   Note: The following assumes this CMakeLists.txt file in the directory:
#         OSIRIS_ROOT/tools/misc/cmake
#------------------------------------------------------------------------------
get_filename_component( OSIRIS_ROOT_DIR "${CMAKE_SOURCE_DIR}/../../.." ABSOLUTE )
message("The Osiris package Root Directory is at:")
message("    ${OSIRIS_ROOT_DIR}")
message("")
SET_PROPERTY(GLOBAL PROPERTY OSIRIS_ROOT_DIR "${OSIRIS_ROOT_DIR}")

#------------------------------------------------------------------------------
# Configure CMake itself
#------------------------------------------------------------------------------

# Tell Cmake the directory where any addtional 'modules' can be found. We
#   need these additional modules in order to find and configure Libraries 
#   that CMake does not have builit-in support for (yet)
# Currently use external modules to find/configure: 
#   - the library FFTW
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${OSIRIS_ROOT_DIR}/tools/misc/cmake/modules/")

#------------------------------------------------------------------------------
# some setup boilerplate:
#------------------------------------------------------------------------------

# 1: create a variable the tell us if this is the initial CMAKE setup 
#    (i.e. the first time through this script) in this directory. 
#    This is used mainly as a mark so that we wont overwrite user specified values
if(NOT DEFINED IS_INITAL_SETUP)
    set( IS_INITAL_SETUP "1"  CACHE INTERNAL "")
endif(NOT DEFINED IS_INITAL_SETUP)

# 2: For sanity, impose a Build Type if not set
if(NOT CMAKE_BUILD_TYPE)
    Message("Note: The build type is unset, so setting the build type to 'Release'")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Most common build options: Release and Debug")
endif(NOT CMAKE_BUILD_TYPE)

#------------------------------------------------------------------------------
# start setting up our project...
#   TODO: HDF5 'Find' package is a MESS and UNUSABLE. They have even admitted it. It has changed several times. Nail down the most modern one once and  for all.
#------------------------------------------------------------------------------
project (osiris)
enable_language (Fortran)
include(CmakeUtitlities.cmake)
include(OsirisBuildUtils.cmake)
FIND_PACKAGE(MPI REQUIRED)

if(UNIX)
    # On unix, get the HDF5 library. Hints my be given by setting the HDF5_ROOT environment variable 
    FIND_PACKAGE(HDF5 COMPONENTS C Fortran)
    set( HDF5_LIBRARIES "${HDF5_Fortran_LIBRARY_hdf5_fortran}" 
                        "${HDF5_Fortran_LIBRARY_hdf5}" 
                        "${HDF5_Fortran_LIBRARY_z}" 
                        "${HDF5_hdf5_LIBRARY_RELEASE}"
                        "${HDF5_hdf5_fortran_LIBRARY_RELEASE}")

    # Posix compilers uses "-" for compiler switches.q
    set(COMP_SWITCH "-")
endif(UNIX)

# Do Windows Visual Studio options, if needed
WINDOWS_VS_CONFIG()


#------------------------------------------------------------------------------
# Handle the Fortran <--> C naming convention housekeeping. This is set by the user
#   on a system-by-sytem basis in normal GNU Make, but it is a nice feature of CMake
#   this is figures can be set automatically. I still leave it as an option in the
#   User Interface GUI if this gets things wrong for some reason.
#   (Note: ${FortranCInterface_GLOBAL_CASE} may be of some use one day. It has value of either "UPPER" or "LOWER")
include(FortranCInterface)
string(LENGTH "${FortranCInterface_GLOBAL__PREFIX}" NUM_UNDERSCORES)
set(Osiris_NUMBER_UNDERSCORES_FORTRAN   "${NUM_UNDERSCORES}" CACHE STRING    "The number of unserscores in FORTRAN external symbols (used in FORTRAN/C calls). Should be set automaticaly to correct value." )


# ---------------------------------------------------
# Input parameters
#  Note: They are named this way because naything before the '_' till determine which node the option appears in
#         in the Cmake-GUI or CCMAKE user interfaces.
# ---------------------------------------------------
set(OsirisOptions_WRAPACIOUS             0                  CACHE BOOL      "Build Osiris with Wrapacious")
set(Osiris_P_X_DIM                       2                  CACHE STRING    "Number of spatial dimensions (1,2, or 3)")
set(Osiris_SIMD                         ""                  CACHE STRING    "Build Osiris with support Intel SIMD Support. Valid values are SSE, AVX, BGP, BGQ.  Set to empty to disable." )
set(Osiris_USE_DOUBLE                    1                  CACHE BOOL      "Compile with Double Percision Numbers (if not set, single percision is used)")
set(Osiris_NUMBER_UNDERSCORES_FORTRAN   "1"                 CACHE STRING    "The number of unserscores in FORTRAN external symbols (used in FORTRAN/C calls). Should be set automaticaly to correct value." )
set(OsirisOptions_FFTW		               0                  CACHE BOOL      "Build Osiris with support for FFTW." )
set(OsirisOptions_PARALLEL_IO            0                  CACHE BOOL      "Build Osiris with support for parallel I/O." )
set(OsirisOptions_ALLOC_EXT              0                  CACHE BOOL      "Build Osiris with support for new Fortran2000-like allocation extensions.")
set(OsirisOptions_NO_FILE_FLUSH          0                  CACHE BOOL      "Bulid Osiris with such that no file flushes occur.")
set(OsirisOptions_USE_LOG                0                  CACHE BOOL      "Build Osiris with log file output.")
set(OsirisOptions_Scion_RST_IO           0                  CACHE BOOL      "IF this is check, restarts will use the Scion Checkpointsystem. Otherwise, normal binary restarts will be used")
set(OsirisOptions_System_Timer          "__POSIX_TIMER__"   CACHE STRING    "Select the type of timer. If in doubt, use the default _POSIX_TIMER_ (note: _POSIX_TIMER_ works on Windows too)")
set(OsirisOptions_NO_NML_INTFILE         0                  CACHE BOOL      "Check this if your Fortran compiler does NOT support internal namelists (for older compilers.. for exmaple ifort 11 on windows needs this.")
set(OsirisOptions_DISABLE_PGC            0                  CACHE BOOL      "Compile the code WITHOUT PGC support.")
set(OsirisOptions_DISABLE_TILES          0                  CACHE BOOL      "Compile the code WITHOUT Tiles (dynamic load balance) support.")
set(OsirisOptions_DISABLE_OVERDENSE      0                  CACHE BOOL      "Compile the code WITHOUT Overdense support.")
set(OsirisOptions_DISABLE_QED            0                  CACHE BOOL      "Compile the code WITHOUT QED support.")
set(OsirisOptions_DISABLE_SHEAR          0                  CACHE BOOL      "Compile the code WITHOUT SHEAR support.")
set(OsirisOptions_DISABLE_CYLMODES       0                  CACHE BOOL      "Compile the code WITHOUT Cylrindrical Modes support.")
set(OsirisOptions_DISABLE_RAD            1                  CACHE BOOL      "Compile the code WITHOUT Radiation support.")
set(OsirisOptions_DISABLE_FEI            0                  CACHE BOOL      "Compile the code WITHOUT Fei solver support.")
set(OsirisOptions_DISABLE_QEDCYL         0                  CACHE BOOL      "Compile the code WITHOUT QEDCYL support.")
set(OsirisOptions_DISABLE_GR             0                  CACHE BOOL      "Compile the code WITHOUT GR module.")
set(OsirisOptions_ENABLE_CUDA            0                  CACHE BOOL      "Build Osiris with CUDA support.")
set(OsirisOptions_ENABLE_Coroutines      0                  CACHE BOOL      "Build Osiris with Support for Co-Routines.")


# these lines make some of the fields into listboxes in the CMAKE GUI (to constrain choices in UI to allowed values only)
set_property(CACHE Osiris_SIMD                              PROPERTY STRINGS " " SSE AVX BGP BGQ)
set_property(CACHE OsirisOptions_System_Timer               PROPERTY STRINGS _POSIX_TIMER_ __MPI_TIMER__ __AIX_TIMER__ __MACH_TIMER__  __HRTIME_TIMER__)
set_property(CACHE Osiris_NUMBER_UNDERSCORES_FORTRAN        PROPERTY STRINGS "0" "1" "2" )


# ---------------------------------------------------
# Process the input variables into compiler flags....
# ---------------------------------------------------
# dimension
#       add_flag_restrict_to_list_or_error( 'input var', 'var_name', 'target list(s) that will be appended',
#                                           'allowed_values', 'error message if invalid value' )
#   A name/value pair 'var_name'='an_allowed_value' will be added to the targest list(s).
#   If the 'var_name' is empty then '-Dan_allowed_value' will be added to the targest list(s).
#
add_flag_restrict_to_list_or_error(     
                                        "Osiris_P_X_DIM" 
                                        "P_X_DIM"
                                        "OSIRIS_FPP_FLAGS;CMAKE_C_FLAGS;CMAKE_CXX_FLAGS"
                                        "1;2;3" 
                                        "Fatal Error: P_X_DIM (the number of spatial dimension) MUST BE SET to 1,2 or 3!")

# number of underscore in Fortran external symbols
#       add_flag_restrict_to_list_or_error( 'input var', 'var_name', 'target list(s) that will be appended', 
#                                           'allowed_values','values written for allowed values', 'error message if invalid value' )
#   A name/value pair 'var_name'='a_value_written_for_allowed_value' will be added to the of the targest list(s).
#   If the 'var_name' is empty then '-Da_value_written_for_allowed_value' will be added to the of the targest list(s).
#           Example:    add_flag_restrict_to_list_or_error( 'num_dims_string_input',"dims",CMAKE_C_FLAGS,"one:two:three",'1:2:3', "You have an invalid number of dimensions.")
#           Now say 'num_dims_string_input'='two'. Then CMAKE_C_FLAGS would have 'num_dims=2; added to the end
#
add_flag_restrict_to_list__map_values__or_error( 
                                        "Osiris_NUMBER_UNDERSCORES_FORTRAN" 
                                        ""
                                        "OSIRIS_FPP_FLAGS;CMAKE_C_FLAGS;CMAKE_CXX_FLAGS"  
                                        "0;1;2" "FORTRANNOUNDERSCORE;FORTRANSINGLEUNDERSCORE;FORTRANDOUBLEUNDERSCORE"
                                        "Fatal Error: You must set the number of underscores to 0, 1, or 2." )

#   SIMD configuration
#       append_compile_flag( 'compile flag', 'target list(s) where flag will be appended'    )
#
STRING(STRIP "${Osiris_SIMD}" Osiris_SIMD)
if (NOT "${Osiris_SIMD}" STREQUAL "")
    append_compile_flag( "SIMD"                     "OSIRIS_FPP_FLAGS;CMAKE_C_FLAGS;CMAKE_CXX_FLAGS")
    append_compile_flag( "SIMD_${Osiris_SIMD}"      "OSIRIS_FPP_FLAGS;CMAKE_C_FLAGS;CMAKE_CXX_FLAGS")
endif (NOT "${Osiris_SIMD}" STREQUAL "")

#   Various true/false flags
#       bool_yes_no__flag('input var',  'target list(s) where values appended', 'value written if input true', 'value written if input false')
#
#                           'input var'               'target list(s)'    'value if true'       'value if false'
#
bool_yes_no__flag("OsirisOptions_PARALLEL_IO"       "OSIRIS_FPP_FLAGS"  "Osiris_PARALLEL_IO"        ""                  )
bool_yes_no__flag("OsirisOptions_ALLOC_EXT"         "OSIRIS_FPP_FLAGS"  "__HAS_ALLOC_EXT__"         ""                  )
bool_yes_no__flag("OsirisOptions_System_Timer"      "OSIRIS_FPP_FLAGS"  "${Osiris_System_Timer}"    ""                  )
bool_yes_no__flag("OsirisOptions_NO_FILE_FLUSH"     "OSIRIS_FPP_FLAGS"  "__NO_FILE_FLUSH__"         ""                  )
bool_yes_no__flag("OsirisOptions_USE_LOG"           "OSIRIS_FPP_FLAGS"  "__USE_LOG__"               ""                  )
bool_yes_no__flag("OsirisOptions_Scion_RST_IO"      "OSIRIS_FPP_FLAGS"  "__RST_SION__"              "__RST_BINARY__"    )
bool_yes_no__flag("OsirisOptions_NO_NML_INTFILE"    "OSIRIS_FPP_FLAGS"  "__NO_NML_INTFILE__"        ""                  )
bool_yes_no__flag("OsirisOptions_FFTW"              "OSIRIS_FPP_FLAGS;CMAKE_C_FLAGS;CMAKE_CXX_FLAGS"	
                                                                        "FFTW_ENABLED"              ""                  )
bool_yes_no__flag("Osiris_USE_DOUBLE"               "OSIRIS_FPP_FLAGS;CMAKE_C_FLAGS;CMAKE_CXX_FLAGS"
                                                                        "PRECISION_DOUBLE"          "PRECISION_SINGLE"  )
bool_yes_no__flag("OsirisOptions_DISABLE_PGC"       "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_PGC"        )
bool_yes_no__flag("OsirisOptions_DISABLE_QED"       "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_QED"        )
bool_yes_no__flag("OsirisOptions_DISABLE_SHEAR"     "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_SHEAR"      )
bool_yes_no__flag("OsirisOptions_DISABLE_CYLMODES"  "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_CYLMODES"   )
bool_yes_no__flag("OsirisOptions_DISABLE_RAD"       "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_RAD"        )
bool_yes_no__flag("OsirisOptions_DISABLE_FEI"       "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_FEI"        )
bool_yes_no__flag("OsirisOptions_DISABLE_FEI"       "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_FEI"        )
bool_yes_no__flag("OsirisOptions_DISABLE_QEDCYL"    "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_QEDCYL"     )
bool_yes_no__flag("OsirisOptions_DISABLE_TILES"     "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_TILES"      )
bool_yes_no__flag("OsirisOptions_DISABLE_OVERDENSE" "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_OVERDENSE"  )
bool_yes_no__flag("OsirisOptions_DISABLE_GR"        "OSIRIS_FPP_FLAGS"  ""                          "ENABLE_GR")
bool_yes_no__flag("OsirisOptions_ENABLE_Coroutines" "OSIRIS_FPP_FLAGS"  "ENABLE_ENABLE_COROUTINE"   ""                  )
bool_yes_no__flag("OsirisOptions_ENABLE_CUDA"       "OSIRIS_FPP_FLAGS"  "ENABLE_CUDA"               ""                  )



# add compile flags that are set independent of user input.
append_compile_flag( "DOS_REV"                      "3"                 "OSIRIS_FPP_FLAGS")
append_compile_flag( "HDF5"                         "YES"               "OSIRIS_FPP_FLAGS")


#-----------------------------------------------------------------
# Source Files
# For better or Worse, I decided to user the reverse of the ususal way fo doing things.
#   Here, all files under ./source are considered as part of the source code UNLESS
#   specifically excludeded.
#-----------------------------------------------------------------

SET(OBJS_C  "${OSIRIS_ROOT_DIR}/source/system/os-sys-multi-c.c" 
            "${OSIRIS_ROOT_DIR}/source/simd/os-spec-push-c.c" 
            "${OSIRIS_ROOT_DIR}/source/simd/os-spec-current-c.c"
            "${OSIRIS_ROOT_DIR}/source/zdf/zdf-parallel.c"
            "${OSIRIS_ROOT_DIR}/source/zdf/zdf.c"
            )

# I take the route ozf adding all files in the source directlory to the build automatically then
#   culling out the unwanted. We can always switch back to explictly listing the source files.
file( GLOB_RECURSE  
      OBJSALL 
      RELATIVE  "${OSIRIS_ROOT_DIR}/source" 
                "${OSIRIS_ROOT_DIR}/source/*.f03" 
                "${OSIRIS_ROOT_DIR}/source/*.f90" )

#---------- Add lines like the following to exclude files from being built -------------------------  
remove_matching_from_list(OBJSALL "meta/")                    # <--- removes any file under the meta/ directory from the build
remove_matching_from_list(OBJSALL "os-spline.f90")
remove_matching_from_list(OBJSALL "os-spec-diag-def.f03")     # This is #included inside another file not compiled directly
remove_matching_from_list(OBJSALL "os-psource-def.f03")       # This is #included inside another file not compiled directly
remove_matching_from_list(OBJSALL "ioperf.f90")               # This is it's own standalone executable (TODO: add a target for it in this CmakeFiles.txt)
remove_matching_from_list(OBJSALL "ioperf.f03")
remove_matching_from_list(OBJSALL "os-restart-io-binary.f90")
remove_matching_from_list(OBJSALL "os-restart-io-sion.f90")
remove_matching_from_list(OBJSALL "os-main.f03")              # We seperate this out and explicitly add it into the list of files.
remove_matching_from_list(OBJSALL "cuda/")                    # 
remove_matching_from_list(OBJSALL "os-emf-solver-define.f03")
remove_matching_from_list(OBJSALL "os-emf-bnd-define.f03")
remove_matching_from_list(OBJSALL "os-emf-marder.f90")
remove_matching_from_list(OBJSALL "build/")                  # the GNU make version puts temp build files in source/build.. ignore these if present.


if( OsirisOptions_DISABLE_PGC )
  remove_matching_from_list(OBJSALL "pgc/")
endif()
if( OsirisOptions_DISABLE_QED )
  remove_matching_from_list(OBJSALL "qed/")
endif()
if( OsirisOptions_DISABLE_SHEAR )
  remove_matching_from_list(OBJSALL "shear/")
endif()
if( OsirisOptions_DISABLE_CYLMODES )
  remove_matching_from_list(OBJSALL "cyl_modes/")
endif()
if( OsirisOptions_DISABLE_RAD )
  remove_matching_from_list(OBJSALL "radio/")
endif()
if( OsirisOptions_DISABLE_FEI )
  remove_matching_from_list(OBJSALL "fei/")
endif()
if( OsirisOptions_DISABLE_TILES )
  remove_matching_from_list(OBJSALL "tiles/")
endif()
if( OsirisOptions_DISABLE_OVERDENSE )
  remove_matching_from_list(OBJSALL "overdense/")
endif()
if( OsirisOptions_DISABLE_QEDCYL )
  remove_matching_from_list(OBJSALL "qed-cyl/")
endif()
if( OsirisOptions_DISABLE_GR )
  remove_matching_from_list(OBJSALL "gr/")
endif()

# Specify the file with the program's 'main'
SET(OBJSMAIN "os-main.f03")


# ---------------------------------------------------
# ---------------------------------------------------
# Deal with optional 3rd party source code/library packages
# ---------------------------------------------------
# ---------------------------------------------------

# -------- FFTW ------
if( OsirisOptions_FFTW)
  message("-------- Using the FFTW Library -----------")

  # Tell cmake to configure FFTW
  #  This package looks around the system to find FFTW.. if
  #  you wish to 'guide' where to look, set the 'FFTW_ROOT'
  #  environment variable to the path of the FFTW installation you want to use!
  #
  #  Also usefull that I may add: If the 'FFTW_USE_STATIC_LIBS' is True then
  #  the static libraries are used. Which is probably the more usefull case for us.
  #
  # If sucessful, the 'find_package' sets:
  #     FFTW_FOUND
  #     FFTW_INCLUDES
  #     FFTW_LIBRARIES
  find_package(FFTW QUIET)

  if(NOT FFTW_FOUND)
    message( "*****************************************")
    message( "              Error!                     ")
    message( "*****************************************")
    message( "You have requested that Osiris be built with support for FFTW.")
    message( "But the FFTW Library could not be found! ")
    message( "                                         ")
    message( "Please set your 'FFTW_ROOT' environment variable to the path of ")
    message( "a valid FFTW installation and re-reun Cmake.")
    message( "                                         ")
    message( FATAL_ERROR "" )
  endif(NOT FFTW_FOUND)

  #set(FFTW_DIR "C:/lab_code/third_party/fftw")
  #set(FFTW_LIBRARIES "${FFTW_DIR}/libfftw3-3.lib" "${FFTW_DIR_LIB}/libfftw3f-3.lib")
  #set(FFTW_INCLUDES "C:/lab_code/third_party/fftw")
endif( OsirisOptions_FFTW)


# -------- CUDA ------
SET(CUDA_TOOLS_INCLUDE_DIRS "")
if( OsirisOptions_ENABLE_CUDA )
  # Tell cmake to configure CUDA
  FIND_PACKAGE(CUDA)

  # Add in any extra source code    
  SET(OBJS_CUDA_FORTRAN cuda/os-simulation-cuda-define.f03 cuda/os-simulation-cuda.f03 cuda/os-cuda-serialize-state.f03)
  SET(OBJS_CUDA_C source/cuda/os-cuda-input-deck.cpp 
                  source/cuda/os-cuda-state.cpp
                  source/cuda/simple-cpp-profiler.cpp
                  source/cuda/os-cuda-device-runtime-interaction.cpp
                  source/cuda/os-cuda-main.cpp )

  SET(CUDA_TOOLS_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/opencontainers_1_7_6/include")

endif(OsirisOptions_ENABLE_CUDA)

# -------- Opencontainers ------
#  (Used for the CUDA input deck)
if( OsirisOptions_ENABLE_CUDA )
    message("Untaring needed source code libraries....")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar "xf" "${CMAKE_SOURCE_DIR}/tools/misc/cuda/opencontainers.tar")
    message("Finished Untaring needed source code libraries....")
endif( OsirisOptions_ENABLE_CUDA )

# -------- Boost (part of Boost) ------
#  (Used for Coroutine/Continuations)
if( OsirisOptions_ENABLE_Coroutines )
  # unzip the boost .zip file (I zipped boost so that it doesn't pollute the repor with 500 extra files)
  message("Untaring needed Boost source code libraries....")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar "xf" "${OSIRIS_ROOT_DIR}/tools/misc/boost/boost.tar")

  # set some compiler defintions to fake the ones Boost natually does.
  add_definitions(-DBOOST_CONTEXT_EXPORT="")
  ADD_DEFINITIONS("-DBOOST_ALL_NO_LIB")

  # get all the source code files from Boost
  file( GLOB_RECURSE OBJSBOOST RELATIVE "${OSIRIS_ROOT_DIR}" "${CMAKE_BINARY_DIR}/boost/libs/context/*.cpp")

  if(MSVC)
    remove_matching_from_list(OBJSBOOST "context/src/posix")  
    add_definitions(-DBOOST_USE_WINFIB)
  else(MSVC)
    add_definitions(-DBOOST_USE_UCONTEXT)
    remove_matching_from_list(OBJSBOOST "context/src/windows")
    #set_source_files_properties(foo.s PROPERTIES COMPILE_FLAGS "-x assembler-with-cpp")
  endif(MSVC)

  include_directories( "${CMAKE_BINARY_DIR}/boost" )
endif(OsirisOptions_ENABLE_Coroutines)


# ---------------------------------------------------
# ---------------------------------------------------
# Deal with any specific flags  needed for specifc compilers
# ---------------------------------------------------
# ---------------------------------------------------
if(UNIX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

  if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -ffree-line-length-none -fno-range-check -Wmissing-include-dirs")
  endif()

  if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable 5117 -fpp")
  message("sdfkjlslkfsdkfkljdsfjkds")
  message("0sdfds0f0f0000000000000000000000000000000")
  endif()
endif(UNIX)
#--------------------------------------
# Setup the Preprocessor to be applied to the Fortran source files
#       (Osiris uses the C-preprocessor on Fortran files rather then the standard Fortran preprocessor)
#
# For GCC (and intel's ICC):
#       FPP = gcc -C -E -x assembler-with-cpp
#               Where:
#                   -C: Do not discard comments. (needed because // is a valid operator in Fortran and we can't just nuke it)
#                   -E: Only so Preprocessing... then stop
#                   -x: treat input file as type "assembler-with-cpp" 
# For Windows Visual Studio Compiler:
#       FPP = cl -nologo /EP"
#               Where:
#                   /nologo:    stops annoying Microsoft compiler banner
#                   /EP:        Only do preprocessing... then stop
#                   /C:         Do not discard comments (needed because // is a valid operator in Fortran and we can't just nuke it)
if(MSVC)
  # Windows and  Visual Studio
  SET(FPP "cl")
  LIST(APPEND FPP_FLAGS "/E" "/C" "/nologo" "-I${OSIRIS_ROOT_DIR}/source" "-I${CMAKE_BINARY_DIR}")
else()
  # Everything else (unixy things)
  set(FPP "gcc")
	LIST(APPEND FPP_FLAGS -E -C -nostdinc -x assembler-with-cpp "-I${OSIRIS_ROOT_DIR}/source" "-I${CMAKE_BINARY_DIR}")
endif(MSVC)


# make directory for preprocessed files
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/preprocessed")


# Now that we have a list of source files, our preprocessor now join them together.
#   The following 'register' functions fix up files path (make them absolute)
#   and create 'custom rule' that assigns the C preprocessor to be applied to each of 
#   the Fortran source files. The compiler definitions/flags used by the preprocessor
#       to process the source files are in the variable 'OSIRIS_FPP_FLAGS' which
if(MSVC_IDE)
	register_fortran_source_for_IDE(OBJSALL OBJSALL)
  register_fortran_source_for_IDE(OBJSMAIN OBJSMAIN)
else()
    #create_compiler_defintions( OSIRIS_FPP_FLAGS OSIRIS_FPP_FLAGS)
    register_fortran_source_list(OBJSALL OBJSALL)
    register_fortran_source_list(OBJSMAIN OBJSMAIN)
    register_fortran_source_list(OBJS_CUDA_FORTRAN OBJS_CUDA_FORTRAN)
endif(MSVC_IDE)

if( OsirisOptions_ENABLE_Coroutines )
  foreach(f ${OBJSBOOST})
    list(APPEND OBJSALL "${OSIRIS_ROOT_DIR}/${f}")
  endforeach(f)
endif(OsirisOptions_ENABLE_Coroutines)

# Note: I have utility functions that might some in handy here:
#       remove_flag_from_list( CMAKE_C_FLAGS_DEBUG "D_DEBUG") <--- to specifically get rid of some annoying flag (like warnings)

# This is an annoyance in CMAKE. These variables are actually string when they should be a "list".. so convert to be safe.
TO_LIST_SPACES( "CMAKE_C_FLAGS"  "CMAKE_C_FLAGS")
TO_LIST_SPACES( "CMAKE_CXX_FLAGS"  "CMAKE_CXX_FLAGS")
#TO_LIST_SPACES( "CMAKE_Fortran_FLAGS"  "CMAKE_Fortran_FLAGS")
#TO_LIST_SPACES( "CMAKE_Fortran_FLAGS_DEBUG"  "CMAKE_Fortran_FLAGS_DEBUG")
#TO_LIST_SPACES( "CMAKE_Fortran_FLAGS_RELEASE"  "CMAKE_Fortran_FLAGS_RELEASE")

if( OsirisOptions_WRAPACIOUS )

    # ensure that we find a python that we can use
    if(UNIX)
      set(OsirisOptions_WRAPACIOUS_PYTHON      "python3"          CACHE String      "Python Command to use for Wrapacious")
    else()
      set(OsirisOptions_WRAPACIOUS_PYTHON      "python"           CACHE String      "Python Command to use for Wrapacious")
    endif()

    # TODO: error if we cannot find it.. I will also put in a python script to test...
    message("Using python3 found at: ${Python3_EXECUTABLE}")

    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/processed_wrapper")
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/processed_wrapper/cache")
    #file(TO_NATIVE_PATH  "${OSIRIS_ROOT_DIR}/tools/wrapacious/wrapicious_input_deck.py" TEMP_WRAPACIOUS_INPUT_DECK_PATH)
    file(TO_NATIVE_PATH  "${OSIRIS_ROOT_DIR}/tools/wrapacious/run.py" TEMP_WRAPACIOUS_INPUT_DECK_PATH)

    # assemble the names of the output wrapper files
    # TODO: make the outputlist for Wrapacious automatic.
    set(OBJS_WRAP "${OBJSALL}" ) 
    remove_matching_from_list(OBJS_WRAP "os-main")                  # This is an altered file. We add it explictly below.
    remove_matching_from_list(OBJS_WRAP "zdf-interface")            # This is broken. Come back to it just for curiosity..
                                                                    #  because it's strcitly not needed since iwould be just wrapping a wrapper.
                                                                    #  It has an interesting corner cas in wrapping ISO_C binding.
    remove_matching_from_list(OBJS_WRAP "memory")
    remove_matching_from_list(OBJS_WRAP "stringutil")
    remove_matching_from_list(OBJS_WRAP "fparser")
    register_wrappers(OBJS_WRAP OBJS_WRAP)

    set(EXTERNAL_DLL_DECLARES   "${CMAKE_BINARY_DIR}/processed_wrapper/__external_shared_lib_exports.c" )
    set(WRAPPER_UTILS           "${CMAKE_BINARY_DIR}/processed_wrapper/__wrapper__utils.f90" )
    set(ALTERED_FILES           "${CMAKE_BINARY_DIR}/processed_wrapper/__altered__os-main.f90" 
                                "${CMAKE_BINARY_DIR}/processed_wrapper/__wrapper____altered__os-main.f90" )

    set_source_files_properties(${EXTERNAL_DLL_DECLARES}  PROPERTIES GENERATED TRUE)
    set_source_files_properties(${WRAPPER_UTILS}          PROPERTIES GENERATED TRUE)
    set_source_files_properties(${ALTERED_FILES}          PROPERTIES GENERATED TRUE)

    set(OBJS_WRAP ${OBJS_WRAP} ${EXTERNAL_DLL_DECLARES} ${WRAPPER_UTILS} ${ALTERED_FILES})

    if(UNIX)
      # GNU and Intel (and all other, really - they just might rename the option) need to be compiled BUILD_WITH_INSTALL_RPATH
      #   PIC (Position Independent )
      #if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
      #endif()
    endif()
endif( OsirisOptions_WRAPACIOUS )


#-------------------------------------------------------------
#
#   Start building the files....
#   
#-------------------------------------------------------------

# setup the final include directories.
include_directories(    ${MPI_CXX_INCLUDE_PATH} 
                        ${MPI_Fortran_INCLUDE_PATH} 
                        ${HDF5_INCLUDE_DIR} 
                        ${CMAKE_BINARY_DIR} 
                       "${OSIRIS_ROOT_DIR}/source"
                       "${OSIRIS_ROOT_DIR}/tools/misc/windows"
					              ${FFTW_INCLUDE_DIRS}
                        ${CUDA_INCLUDE_DIRS}
                       ${CUDA_TOOLS_INCLUDE_DIRS}
                    ) 

set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)
set(output_name "osiris-${Osiris_P_X_DIM}D.e")

if(MSVC_IDE)
    #message("&&&&&&&&&&&&&& OBJSMAIN")
    #message("*** ${OBJSMAIN} ***")
	# final compile step.. assemble all source files
	add_library( osiris_fortran STATIC ${OBJSALL} ) #${OBJSLOG} ${OBJSBASE} )
	add_library( osiris_c STATIC ${OBJS_C} )
	add_executable(${output_name}  ${OBJSMAIN} )
	#
	# specify in the needed libraries for the linker.
	TARGET_LINK_LIBRARIES( ${output_name}  ${HDF5_LIBRARIES}   ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES} osiris_fortran osiris_c ${FFTW_LIBRARIES})
else()
    if( OsirisOptions_WRAPACIOUS )
        # Alertnate way to declare target: generate a 'fake' library and use it.
        #add_library(osiris_bulk OBJECT ${OBJSALL} ${OBJS_C}) #${OBJSLOG} ${OBJSBASE} ${OBJS_C})
        #add_executable(${output_name} $<TARGET_OBJECTS:osiris_bulk> ${OBJSMAIN} )
        #add_library(wrapacious SHARED $<TARGET_OBJECTS:osiris_bulk> ${OBJS_WRAP} ${EXTERNAL_DLL_DECLARES} ${WRAPPER_UTILS} ${ALTERED_FILES} )

        # make a real library
        add_library(osiris_bulk ${OBJSALL} ${OBJS_C}) #${OBJSLOG} ${OBJSBASE} ${OBJS_C})
        add_executable(${output_name} ${OBJSMAIN} )
        TARGET_LINK_LIBRARIES( ${output_name}  osiris_bulk ${HDF5_LIBRARIES}   ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES}  ${FFTW_LIBRARIES})

        
        add_library(pyosiris_lib SHARED ${OBJS_WRAP} ) #${EXTERNAL_DLL_DECLARES} ${WRAPPER_UTILS} ${ALTERED_FILES} )
        #add_library(pyosiris_lib SHARED ${OBJS_WRAP} ${EXTERNAL_DLL_DECLARES} ${WRAPPER_UTILS} ${ALTERED_FILES} )
        TARGET_LINK_LIBRARIES( pyosiris_lib osiris_bulk ${HDF5_LIBRARIES}   ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES}  ${FFTW_LIBRARIES})
        ADD_DEPENDENCIES ( pyosiris_lib ${output_name} )
        set_target_properties( pyosiris_lib PROPERTIES PREFIX "")
        set_target_properties( pyosiris_lib PROPERTIES OUTPUT_NAME  "pyosiris_lib")

        # One method is in 'register_wrappers' in OsirisBuildUtils'. It add a custom rule for every file in the wrapper.
        # It's the moethod I am currently using.. go see 'register_wrappers' in OsirisBuildUtils' for more details.

        # Note: Another way. to get things to build is to the following. It has the unwanted effect that 'custom targets' are alwayas run even if no changes were done to
        #       any source files... but wrapacious checks upfront for that condition and quickly exits. And running 'make' when there is nothing to do
        #       is not an important case to optimise/design for.. so this may be ok. What's happening here is that we are running the command every time
        #       the osiris-bulk target is run.. rather then hanging execution of file generation.
        #
        # add_custom_target( wrap_this_stuff
        #             python $ENV{WRAPACIOUS_HOME}/wrapacious.py "${TEMP_WRAPACIOUS_INPUT_DECK_PATH}"
        #             DEPENDS osiris_bulk
        #             COMMENT "Running Wrapacious to generate Python wrapper." 
        #             )
        # ADD_DEPENDENCIES( wrapacious wrap_this_stuff)

        # Yet Another way.. this one solves the issues of the previous way.. it is only run if '${output_name}' is updated and run!
        #      Note: we need to hang this rule off the osiris executable because building the executable triggers the preprocessing of 'os-main.f03'
        #       which wrapacious needs to properly make the wrapper (and we cant put os-main.f03 into the static library)
        #
        #      Downside is that you can't just build 'osiris' executable by itself without triggering the running of wrapacious. This is not a huge thing
        #       but is an inconvience.
        #
        # Two ways to run the command: 
        #   direct path:
        #     COMMAND python $ENV{WRAPACIOUS_HOME}/wrapacious.py "${TEMP_WRAPACIOUS_INPUT_DECK_PATH}"
        #   using Python module loading machinary (some extra string escaping of the path to the input config is needed)
        #     string(REPLACE "\\" "\\\\" TEMP_WRAPACIOUS_INPUT_DECK_PATH "${TEMP_WRAPACIOUS_INPUT_DECK_PATH}")
        #     COMMAND python -c\"import wrapacious.Wrapacious as w;w.startup("${TEMP_WRAPACIOUS_INPUT_DECK_PATH}"

        # COMMAND "${OsirisOptions_WRAPACIOUS_PYTHON}" ARGS -c"import wrapacious.Wrapacious as w;w.startup('${TEMP_WRAPACIOUS_INPUT_DECK_PATH}');"
        string(REPLACE "\\" "\\\\" TEMP_WRAPACIOUS_INPUT_DECK_PATH "${TEMP_WRAPACIOUS_INPUT_DECK_PATH}")
        add_custom_command( TARGET ${output_name} POST_BUILD
                          COMMAND "${OsirisOptions_WRAPACIOUS_PYTHON}" ARGS "${TEMP_WRAPACIOUS_INPUT_DECK_PATH}"
                          COMMENT "Generating python wrapper using Wrapacious."
                          VERBATIM )

    else()
       if(OsirisOptions_ENABLE_CUDA)
        cuda_add_executable(${output_name} ${OBJSALL} ${OBJSLOG} ${OBJSBASE} ${OBJS_C} ${OBJS_CUDA_FORTRAN} ${OBJS_CUDA_C} ${OBJSMAIN} )
       else()
    	  add_executable(${output_name} ${OBJSALL} ${OBJSMAIN} ${OBJS_C} ${OBJS_CUDA_FORTRAN} ${OBJS_CUDA_C})
       endif(OsirisOptions_ENABLE_CUDA)
      
    	# specify in the needed libraries for the linker.
    	TARGET_LINK_LIBRARIES( ${output_name}  ${HDF5_LIBRARIES}   ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES}  ${FFTW_LIBRARIES})
    endif( OsirisOptions_WRAPACIOUS )
endif(MSVC_IDE)

# Add a filthy windows library if needed.
if(MSVC)
        TARGET_LINK_LIBRARIES(${output_name} "Ws2_32.lib")
        if( OsirisOptions_WRAPACIOUS )
            TARGET_LINK_LIBRARIES(pyosiris_lib "Ws2_32.lib")
            SET_TARGET_PROPERTIES( pyosiris_lib PROPERTIES LINKER_LANGUAGE Fortran)
        endif( OsirisOptions_WRAPACIOUS )
endif(MSVC)

if(MSVC)
    SET_TARGET_PROPERTIES( ${output_name}  PROPERTIES LINKER_LANGUAGE Fortran)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /extfor:f90,f03  /free")
    #TARGET_LINK_LIBRARIES(${output_name} "Ws2_32.lib")
    if(MSVC_IDE)
        #target_compile_definitions(osiris_fortran PUBLIC "${OSIRIS_FPP_FLAGS}")
        #target_compile_definitions(${output_name} PUBLIC "${OSIRIS_FPP_FLAGS}" )

        target_compile_options(osiris_fortran PUBLIC $<$<CONFIG:Release>:/O3>)
        target_compile_options(osiris_fortran PUBLIC $<$<CONFIG:Debug>:/Od /debug:full >)

        target_compile_options( osiris_fortran PUBLIC -fpp)
        target_compile_options( osiris_fortran PUBLIC  "${OSIRIS_FPP_FLAGS}")
        target_compile_options( ${output_name} PUBLIC -fpp)
        target_compile_options( ${output_name} PUBLIC  "${OSIRIS_FPP_FLAGS}")
        SET_TARGET_PROPERTIES( osiris_fortran  PROPERTIES Fortran_FORMAT "Free")
    endif(MSVC_IDE)
endif(MSVC)
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /extfor:f90,f03")
#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /extfor:.f03")
#set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} /extfor:.f03")


#SET_TARGET_PROPERTIES( osiris_fortran  PROPERTIES Fortran_FORMAT "Free")


if(MSVC)
    #target_compile_options( osiris_fortran PUBLIC -fpp)
    #target_compile_options( ${output_name} PUBLIC fpp)
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /extfor:f90,f03")
endif(MSVC)

# Setup the Visual Studio IDE, if needed.
WINDOWS_VS_IDE_CONFIG("${output_name}")


# mark so that  we know if this file is run again in same directory so we won't force
#   defaults to write over user entered values
set( IS_INITAL_SETUP "0" CACHE INTERNAL "")
#message("${IS_INITAL_SETUP}")

#Uncomment the following line for detailed infomration about build flags
VERBOSE_INFO()

#-------------------------------------------------------------
#
#   END of main build script...
#   
#-------------------------------------------------------------

#SET_TARGET_PROPERTIES(${APP_NAME} PROPERTIES LINKER_LANGUAGE CXX)


# variables used in producing a configuariton file to tell what the build options were
if(Osiris_USE_DOUBLE)
  set(OSIRIS_PRECISION "DOUBLE")
else(Osiris_USE_DOUBLE)
  set(OSIRIS_PRECISION "SINGLE")
endif(Osiris_USE_DOUBLE)
configure_file( "build-info.template" "${CMAKE_BINARY_DIR}/__generated_config_info.f90" )
