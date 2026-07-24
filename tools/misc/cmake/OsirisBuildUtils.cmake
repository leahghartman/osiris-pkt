macro(VERBOSE_INFO)
	Message("----------------------------------------------------------")
	Message("      The build type is '${CMAKE_BUILD_TYPE}'")
	Message("----------------------------------------------------------")
	Message("")
	message("HDF5_DIR:                     ${HDF5_DIR}")
	message("HDF5 includes:                ${HDF5_INCLUDE_DIR}")
	message("HDF5 libraries:               ${HDF5_LIBRARIES}")
	message("-----------------------------------------------------")
	message("      Compiler Flags For this Osiris Build")
	message("-----------------------------------------------------")
	message(" ")

	message("-------------------")
	message("Fortran compilier:                ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER}")
	message("Fortran Flags:                   ${CMAKE_Fortran_FLAGS}")
	message("Fortran Flags Release:            ${CMAKE_Fortran_FLAGS_RELEASE}")
	message("Fortran Flags Debug:	          ${CMAKE_Fortran_FLAGS_DEBUG}")
	message("-------------------")
	message("C Preprocessor Flags For Fortran: ${OSIRIS_FPP_FLAGS}")
	message("-------------------")
	message("C FLAGS:                          ${CMAKE_C_FLAGS}")
	message("C FLAGS Release:                  ${CMAKE_C_FLAGS_RELEASE}")
	message("C FLAGS Debug  :                  ${CMAKE_C_FLAGS_DEBUG}")
	message("-------------------")
	message("C++ FLAGS:                        ${CMAKE_CXX_FLAGS}")
	message("C++ FLAGS Release:                ${CMAKE_CXX_FLAGS_RELEASE}")
	message("C++ FLAGS Debug  :                ${CMAKE_CXX_FLAGS_DEBUG}")
	message("-------------------")
	message("Linker Flags:                  ${LINKER_FLAGS}")
	message("Linker Flags (shared libs):      ${CMAKE_SHARED_LINKER_FLAGS}")
	message("Linker Flags Release:             ${LINKER_FLAGS_RELEASE}")
	message("Linker Flags Debug:               ${LINKER_FLAGS_DEBUG}")


	if(MSVC)
		message("")
		message("-----------------------------------------------------")
		message("     Basic Visual Studio Settings with IFORT")
		message("")
		message(" 0 underscores, _POSIX_TIMER_, and no SSE support")
		message("-----------------------------------------------------")
	endif(MSVC)
endmacro(VERBOSE_INFO)

# Setup specifically for the Visual Studio IDE, if needed.
function(WINDOWS_VS_IDE_CONFIG osiris_executable_name)
	# Some settings for The Visual Studio IDE (this needs work..)
	if(MSVC_IDE) 
		#add files into visual studio (files also have go into an exxecutable target)
		#SOURCE_GROUP(name [REGULAR_EXPRESSION regex] [FILES src1 src2 ...])
		message("Setting Microsoft Visual Studio Properties for target ${osiris_executable_name}.")
		set_target_properties(${osiris_executable_name} PROPERTIES LINK_FLAGS_DEBUG "/SUBSYSTEM:CONSOLE")
		set_target_properties(${osiris_executable_name} PROPERTIES COMPILE_DEFINITIONS_DEBUG "_CONSOLE")
		set_target_properties(${osiris_executable_name} PROPERTIES LINK_FLAGS_RELWITHDEBINFO "/SUBSYSTEM:CONSOLE")
		set_target_properties(${osiris_executable_name} PROPERTIES COMPILE_DEFINITIONS_RELWITHDEBINFO "_CONSOLE")
		set_target_properties(${osiris_executable_name} PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:WINDOWS")
		set_target_properties(${osiris_executable_name} PROPERTIES LINK_FLAGS_MINSIZEREL "/SUBSYSTEM:WINDOWS")
	endif(MSVC_IDE)
endfunction(WINDOWS_VS_IDE_CONFIG osiris_executable_name)

# C:/Program Files/HDF_Group/HDF5/1.8.15/cmake
macro(WINDOWS_VS_CONFIG)
if(MSVC)
	# On windows, use the more modern approach to find HDF5 (the one the HDF5.org suggests.. eventually, Unix will switch to this too)
	# 	This rapport relies on using an explicit " hdf5-config.cmake" found alongside the HDF5 files. The path to this the directory where this file resides should
	# be set in the HDF5_DIR. 
	#message("IS INTIAL ${IS_INITAL_SETUP}")
	# For ease, check if there is an HDF5_DIR environment variable and use this value to set the CMAKE 
	#    variable (so that you can just set this once in your environment var and never worry again)
	#if( IS_INITAL_SETUP)
	#message("CCCCCCCCCCCCCCCUUUUUUNNNT")
	#	set(HDF5_DIR "" CACHE STRING "Set this Variable to the directory where  the file 'hdf5-config.cmake' resides" FORCE)
	#	if(  NOT "$ENV{HDF5_DIR}" STREQUAL "" )
	#		set(HDF5_DIR "$ENV{HDF5_DIR}" CACHE STRING  "Set this Variable to the directory where  the file 'hdf5-config.cmake' resides" FORCE)
	#		message("--- Setting HDF_DIR from the 'HDF5_DIR' environment variable to ${HDF5_DIR}. If this not desired value, edit value in the CMAKE user interface.")
	#	endif(  NOT "$ENV{HDF5_DIR}" STREQUAL "" )
	#endif(IS_INITAL_SETUP)
	
	set(HDF5_DIR "C:/Program Files (x86)/HDF_Group/HDF5/1.10.8/cmake")
	# Under Windows, static link the HDF5 library (to avoid various DLL Hell issues with zlib.dll when using
	# python (expecially matploblic/hdf5 python libraries)
	#FIND_PACKAGE(HDF5 NAMES hdf5 COMPONENTS Fortran C static NO_MODULE  )
	find_package (HDF5)
	set(HDF5_DIR "C:/Program Files/HDF_Group/HDF5/1.10.8/cmake")

	if(NOT HDF5_FOUND)
		ERROR_MESSAGE("Fatal Error.. Could not find HDF5. You need to set the HDF5-DIR cmake variable to the directory that holds the file  'hdf5-config.cmake'. You can also set the HDF5_DIR environment variable (but this requires to you exit CMAKE, completely delete the contents of the build directory ( ${CMAKE_BINARY_DIR} ), then restart CMAKE.")
	endif(NOT HDF5_FOUND)
	
	set( HDF5_INCLUDE_DIR "${HDF5_INCLUDE_DIR}/static" )
	#set( HDF5_LIBRARIES "${HDF5_DIR}/../lib/libhdf5.lib" "${HDF5_DIR}/../lib/libhdf5_f90cstub.lib" "${HDF5_DIR}/../lib/libhdf5_fortran.lib" "${HDF5_DIR}/../lib/zlib.lib")
  set( HDF5_LIBRARIES "${HDF5_DIR}/../lib/hdf5.lib" "${HDF5_DIR}/../lib/hdf5_f90cstub.lib" "${HDF5_DIR}/../lib/libhdf5_fortran.lib")
	
endif(MSVC)

#
# Some windows specific cleans up... hopefully can be taken out one day when understanding of CMAKE get better
#	or windows gets better or if get better or the universergets better etc
#
# The HDF5 detection is a bit broken on windows or probably I am just misusing it.
if(MSVC)
	# Windows compiler uses "/" for compiler switches.
	set(COMP_SWITCH "/")
	
	# MPI detection on Windows has a weird behaviour where the Fortran libray gets set to the C one.
	#	look for that and fix it.
	if( (MPI_Fortran_LIBRARIES STREQUAL MPI_C_LIBRARIES) OR (NOT MPI_Fortran_LIBRARIES AND MPI_C_LIBRARIES) )
		unset(MPI_Fortran_LIBRARIES CACHE)
		# look in directory where we found the C library...
		get_filename_component( TRY_THIS_PATH "${MPI_C_LIBRARIES}" PATH)
		# look fir the fortran library we want.. we can be smater but this OK for now..
		find_library(MPI_Fortran_LIBRARIES NAMES fmpich2 HINTS ${TRY_THIS_PATH} )
	endif( (MPI_Fortran_LIBRARIES STREQUAL MPI_C_LIBRARIES) OR (NOT MPI_Fortran_LIBRARIES AND MPI_C_LIBRARIES) )
	
endif(MSVC)

# Windows compilers (Fortran in this case.) needs to be told which version of the C/C++ runtime to use. 
#	/MD means use DLL (shared library) thread-safe versions of the c runtime (this is de factor one in Windows)
#
if(MSVC)
	set(CMAKE_Fortran_FLAGS				"${CMAKE_Fortran_FLAGS} /MD")
endif(MSVC)

# Setup the build rules for Windows
#if(MSVC)
	#add_library(osiris_lib  STATIC  ${OBJSLOG_FINAL} ${OBJSBASE_FINAL} ${OBJSALL_FINAL}   )
	#TARGET_LINK_LIBRARIES(osiris_lib)
	
	#add_library(osiris_clib STATIC ${OBJS_C} )
	#TARGET_LINK_LIBRARIES(osiris_clib )
	
	#add_executable(osiris ${OBJSMAIN_FINAL})
	#TARGET_LINK_LIBRARIES(osiris osiris_clib osiris_lib  ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${HDF5_LIBRARIES} "Ws2_32.lib")
	
#	add_executable(osiris  ${OBJSMAIN_FINAL} ${OBJSLOG_FINAL} ${OBJSBASE_FINAL} ${OBJSALL_FINAL}  ${OBJS_C} )
#	TARGET_LINK_LIBRARIES(osiris ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${HDF5_LIBRARIES} "Ws2_32.lib")

#endif(MSVC)

endmacro(WINDOWS_VS_CONFIG)


## Notes
# Read properties from source files!
#list (GET OBJSMAIN_FINAL 0 PRO_FILE_NAME)
#message("PRO_FILE_NAME ${PRO_FILE_NAME}")
#get_source_file_property( PROP_VAL ${PRO_FILE_NAME} GENERATED )
#message("PRO_FILE_NAME GENERATED is ${PROP_VAL}")


# 
# This function is specifically for compiling code that has an initial
#	prepossessing step that is not built into the language. For example,
#	this routine was initially written to preprocess Fortran files using the
#	C preprocesser.
# 
# Cmake has a number of cross-cutting, confusing ways to accomplish this. This
#	implementation may not be the best way.. It's main virtue is that it works.
# 	And it can updated later when I grow some more patience with Cmake's documentation.
# 
# The prime issue is that if the original source file is edited, 
# we must make sure that Cmake automatically regenerates the intermediate file.
# (which it does.. but it checks every file every time.. I want to fix this to speed things up)
# TODO: implictly uses 'OSIRIS_FPP_FLAGS' and 'OSIRIS_ROOT_DIR'. Get rid of this output_required_files
#		or make it more explicit
function(register_fortran_source_for_IDE name_of_source_list name_of_output_list)
		string(REPLACE " " ";" var_as_list "${${name_of_source_list}}")
		set(temp_list "")
		GET_PROPERTY(OSIRIS_ROOT_DIR GLOBAL PROPERTY OSIRIS_ROOT_DIR)

		foreach(f ${var_as_list})  
			if(f MATCHES ".f")
				#file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/source/${f}" TEMP_SOURCE)	
				file(TO_NATIVE_PATH "${OSIRIS_ROOT_DIR}/source/${f}" TEMP_SOURCE)	
				#message( "Hi there: ${CMAKE_SOURCE_DIR}/source/${f}" )
				list(APPEND temp_list "${OSIRIS_ROOT_DIR}/source/${f}")
				
				set_source_files_properties("${OSIRIS_ROOT_DIR}/source/${f}" PROPERTIES LANGUAGE Fortran)
				set_source_files_properties("${OSIRIS_ROOT_DIR}/source/${f}" PROPERTIES Fortran_FORMAT "FREE")
				#get_source_file_property(VARME "${CMAKE_SOURCE_DIR}/source/${f}" LANGUAGE)
				#message("ffffff: ${f}    ${VARME}")
			endif(f MATCHES ".f")
		endforeach(f)
		
		set( ${name_of_output_list} "${temp_list}" PARENT_SCOPE)
		#message("ffffff:    ${VARME}")
		#set_source_files_properties(${name_of_output_list} PROPERTIES LANGUAGE Fortran)
endfunction(register_fortran_source_for_IDE name_of_source_list name_of_output_list)

# TODO: implictly uses 'OSIRIS_FPP_FLAGS' and 'OSIRIS_ROOT_DIR'. Get rid of this output_required_files
#		or make it more explicit
function(register_fortran_source_list name_of_source_list name_of_output_list)
		string(REPLACE " " ";" var_as_list "${${name_of_source_list}}")
		set(temp_list "")
		GET_PROPERTY(OSIRIS_ROOT_DIR GLOBAL PROPERTY OSIRIS_ROOT_DIR)

		foreach(f ${var_as_list})  
			if(f MATCHES ".f")
				
				get_filename_component(dest_filename "${f}" NAME_WE )
				#message("filename is: ${dest_filename}")
				#message( "${CMAKE_BINARY_DIR}/processed/${dest_filename}.f90")
				#file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/source/${f}" TEMP_SOURCE)
				file(TO_NATIVE_PATH "${OSIRIS_ROOT_DIR}/source/${f}" TEMP_SOURCE)
				file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/preprocessed/${dest_filename}.f90" TEMP_DEST)		
			add_custom_command(
			   OUTPUT "${CMAKE_BINARY_DIR}/preprocessed/${dest_filename}.f90" PRE_BUILD
			   COMMAND ${FPP}
			   ARGS ${FPP_FLAGS} ${OSIRIS_FPP_FLAGS} ${TEMP_SOURCE} > ${TEMP_DEST}
			   DEPENDS "${OSIRIS_ROOT_DIR}/source/${f}"
			   VERBATIM
			)
			
			list(APPEND temp_list "${CMAKE_BINARY_DIR}/preprocessed/${dest_filename}.f90")
			
			endif(f MATCHES ".f")
		endforeach(f)
		
		set( ${name_of_output_list} "${temp_list}" PARENT_SCOPE)
		set_source_files_properties(${name_of_output_list} PROPERTIES GENERATED TRUE)
endfunction(register_fortran_source_list name_of_source_list name_of_output_list)

#
# Testing for wrapacious
#

# TODO: implictly uses 'OSIRIS_FPP_FLAGS' and 'OSIRIS_ROOT_DIR'. Get rid of this output_required_files
#		or make it more explicit
function(register_wrappers name_of_source_list name_of_output_list)
		string(REPLACE " " ";" var_as_list "${${name_of_source_list}}")
		set(temp_list "")

		file(TO_NATIVE_PATH  "${OSIRIS_ROOT_DIR}/tools/wrapacious/wrapicious-input-deck.py" TEMP_WRAPACIOUS_INPUT_DECK_PATH)

		foreach(f ${var_as_list})  
			if(f MATCHES ".f")
				
				get_filename_component(dest_filename "${f}" NAME_WE )
				file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/preprocessed/${dest_filename}.f90" TEMP_SOURCE)
				file(TO_NATIVE_PATH "${CMAKE_BINARY_DIR}/processed_wrapper/__wrapper__${dest_filename}.f90" TEMP_DEST)

				# add_custom_command(
				# 	OUTPUT "${CMAKE_BINARY_DIR}/processed_wrapper/__wrapper__${dest_filename}.f90"
				# 	COMMAND python $ENV{WRAPACIOUS_HOME}/wrapacious.py "${TEMP_WRAPACIOUS_INPUT_DECK_PATH}"
				# 	DEPENDS "${CMAKE_BINARY_DIR}/processed/${dest_filename}.f90"
				# 	COMMENT "Running wrapacious to generate wrapper." )
			
				list(APPEND temp_list "${CMAKE_BINARY_DIR}/processed_wrapper/__wrapper__${dest_filename}.f90")
				set_source_files_properties(${TEMP_DEST} PROPERTIES GENERATED TRUE)
			endif(f MATCHES ".f")
		endforeach(f)
		
		set( ${name_of_output_list} "${temp_list}" PARENT_SCOPE)
		set_source_files_properties(${name_of_output_list} PROPERTIES GENERATED TRUE)
		set_source_files_properties(${TEMP_DEST} PROPERTIES GENERATED TRUE)

endfunction(register_wrappers name_of_source_list name_of_output_list)

