#	CMake Liquid Concealer
#	Adam Tableman 2012
# 
#	Cmake is good at lots of stuff but once you move beyond simple projects, 
#	strange issues pop up.. and some simple things are harder then they should be.
# 
#	This library tries to cover-up what it can and draw the attention away from 
# CMAKE's various blemishes and 'beauty marks' so you can just focus on 
# creating a decent build-system easily.
# 
# 




# Call in 2 ways:
# 	compile_flag(FLAG_NAME, VARIBALE_TO_WRITE_INTO)
#   compile_flag(FLAG_NAME, FLAG_VALUE, VARIBALE_TO_WRITE_INTO)
# 
# Example:
#	compile_flag("SSE_ENABLE", "DEST")		will set DEST to /DSSE_ENABLE or -DSSE_ENABLE
#	compile_flag("NUM_DIMS", 3, "DEST") 	will set DEST to /DNUM_DIMS=3 or -DNUM_DIMS=3
#	compile_flag("NUM_DIMS", "2", "DEST") 	will set DEST to /DNUM_DIMS=2 or -DNUM_DIMS=2
# 

function( compile_flag flag)
	if(MSVC)
		set(COMP_SWITCH "/")
	else()
		set(COMP_SWITCH "-")
	endif(MSVC)
	
	# one argument
	# output a flag like '-DFLAG' or '/DFLAG'
	if( ${ARGC} EQUAL 2)
		set( ${ARGV1} "${COMP_SWITCH}D${flag}" PARENT_SCOPE)
		#set( ${ARGV1} "${flag}" PARENT_SCOPE)
	endif ( ${ARGC} EQUAL 2)
	
	# two arguments
	# output a flag like '-DFLAG=value' or '/DFLAG=value'
	if( ${ARGC} EQUAL 3)
		set( ${ARGV2} "${COMP_SWITCH}D${flag}=${ARGV1}" PARENT_SCOPE)
		#set( ${ARGV2} "${flag}=${ARGV1}" PARENT_SCOPE)
	endif ( ${ARGC} EQUAL 3)
endfunction( compile_flag flag_name value)

# Call in 2 ways:
# 	append_compile_flag(FLAG_NAME, VARIBALE_TO_WRITE_INTO)
#   append_compile_flag(FLAG_NAME, FLAG_VALUE, VARIBALE_TO_WRITE_INTO)
#		Note: VARIBALE_TO_WRITE_INTO is a CMAKE LIST object.
# 
# Example:
#	append_compile_flag("SSE_ENABLE", "DEST")		will add /DSSE_ENABLE or -DSSE_ENABLE to the end of DEST
#	append_compile_flag("NUM_DIMS", 3, "DEST") 		will add /DNUM_DIMS=3 or -DNUM_DIMS=3 to the end of DEST
#	append_compile_flag("NUM_DIMS", "2", "DEST") 	will add /DNUM_DIMS=2 or -DNUM_DIMS=2 to the end of DEST
# 
# Also, if you wish to write the same thing in multiple lists (like CMAKE_C_FLAGS list and CMAKE_CXX_FLAGS) you
#		can just put multiple names in second argument, separated by semi-colons.
# Example:
#	append_compile_flag("SSE_ENABLE", "CMAKE_C_FLAGS;CMAKE_CXX_FLAGS;Cmake_FORTRAN_FLAGS")		will add /DSSE_ENABLE or -DSSE_ENABLE to the end CMAKE_C_FLAGS and CMAKE_CXX_FLAGS and Cmake_FORTRAN_FLAGS
# 
function( append_compile_flag flag)
	if(MSVC)
		set(COMP_SWITCH "/")
	else()
		set(COMP_SWITCH "-")
	endif(MSVC)
	
	foreach(target_var ${ARGV1})
		# one argument
		# output a flag like '-DFLAG' or '/DFLAG'
		if( ${ARGC} EQUAL 2)
			list( APPEND "${target_var}" "${COMP_SWITCH}D${flag}")
			#list( APPEND "${target_var}" "${flag}")
			set( ${target_var} "${${target_var}}" PARENT_SCOPE)
		endif ( ${ARGC} EQUAL 2)
		
		# one arguments
		# output a flag like '-DFLAG=value' or '/DFLAG=value'
		if( ${ARGC} EQUAL 3)
			list(APPEND "${ARGV2}" "${COMP_SWITCH}D${flag}=${target_var}")
			#list(APPEND "${ARGV2}" "${flag}=${target_var}")
			set( ${ARGV2} "${${ARGV2}}" PARENT_SCOPE)
		endif ( ${ARGC} EQUAL 3)
	endforeach()
endfunction( append_compile_flag flag_name value)


function( append_compile_flaga flag)
	if(MSVC)
		set(COMP_SWITCH "/")
	else()
		set(COMP_SWITCH "-")
	endif(MSVC)
	
	foreach(flag ${var_targets})
	# one argument
	# output a flag like '-DFLAG' or '/DFLAG'
	if( ${ARGC} EQUAL 2)
		list( APPEND "${ARGV1}" "${COMP_SWITCH}D${flag}")
		#list( APPEND "${ARGV1}" "${flag}")
		set( ${ARGV1} "${${ARGV1}}" PARENT_SCOPE)
	endif ( ${ARGC} EQUAL 2)
	
	# one arguments
	# output a flag like '-DFLAG=value' or '/DFLAG=value'
	if( ${ARGC} EQUAL 3)
		list(APPEND "${ARGV2}" "${COMP_SWITCH}D${flag}=${ARGV1}")
		#list(APPEND "${ARGV2}" "${flag}=${ARGV1}")
		set( ${ARGV2} "${${ARGV2}}" PARENT_SCOPE)
	endif ( ${ARGC} EQUAL 3)
endfunction( append_compile_flaga flag_name value)
# 
# Note for later life:
# 	dependencies for a custom_command can also be added using the syntax
#		add_custom_command( TARGET target forms which will execute the commands as soon as 'target' is built
# 
# the argument 'MAIN_DEPENDENCY' is for Visual Studio to help it make a custom build step (ingored on other platform)

# might not work.. might have to use function.. end function instead...
#macro(remove_flag_from_list list_to_change flag_to_remove)
#	foreach(f ${list_to_change})
#			message("remove: look at ${f}")
#			if (f MATCHES "${flag_to_remove}")
#				list(REMOVE_ITEM CMAKE_C_FLAGS_RELEASE "${f}")
#			endif (f MATCHES "${flag_to_remove}")
#	endforeach(f)
#endmacro(remove_flag_from_list)

function(create_compiler_defintions _LIST_NAME OUTPUT_VAR)
	if(MSVC)
		set(COMP_SWITCH "/")
	else()
		set(COMP_SWITCH "-")
	endif(MSVC)
	
  set(NEW_LIST)
  foreach(ITEM ${${_LIST_NAME}})
	list(APPEND NEW_LIST  "${COMP_SWITCH}D${ITEM}")
  endforeach()
  set( ${OUTPUT_VAR} "${NEW_LIST}" PARENT_SCOPE)
endfunction(create_compiler_defintions )

# 
# operates on Cmake LIST object 'list_to_change' and removes flags that MATCH 'flag_to_remove'
#	'list_to_change' can be a Cmake list OR a Cmake string with white space separated values.
#		Cmake uses both in list contexts, but they actually don't mix well. This code, however,
#		ensures that they do play well by scrubbing the input and converting into a proper Cmake 
#		list.. doing operations on this List.. and then returning a Cmake list object.
# 
function(remove_matching_from_list list_to_change flag_to_remove)
	# 'list_to_change' is the name of the list (in parent scope which we have read acces to)
	# "${${list_to_change}}" will give the values in the list
	
	# create a proper list object to edit converting the incomingl list name
	#	to a string then making it back into a list.
	string(REPLACE " " ";" var_as_list "${${list_to_change}}")								# this works
	#string(REGEX MATCHALL "[a-zA-Z]+ |[a-zA-Z]+$" var_as_list "${${list_to_change}}")		# this is a bit more robust
	
	foreach(f ${var_as_list})
			#message("look at ${f} and ${flag_to_remove}")
			if (f MATCHES "${flag_to_remove}")
				list(REMOVE_ITEM var_as_list "${f}")
				#message("     REMOVING ${f}")
			endif (f MATCHES "${flag_to_remove}")
	endforeach(f)
	set( ${list_to_change} "${var_as_list}" PARENT_SCOPE)
endfunction(remove_matching_from_list)


macro(LIST_REPLACE LIST INDEX NEWVALUE)
    list(INSERT ${LIST} ${INDEX} ${NEWVALUE})
    MATH(EXPR __INDEX "${INDEX} + 1")
    list (REMOVE_AT ${LIST} ${__INDEX})
endmacro(LIST_REPLACE)



# 
# Converts a CMake list to a string containing elements separated by spaces
# 
# Example:
#		TO_LIST_SPACES( "CMAKE_C_FLAGS", "NEW_VAR")
#			puts cleaned up "CMAKE_C_FLAGS" in "NEW_VAR"
#			Many important built-in variables, like CMAKE_C_FLAGS, are logically
#			lists but work better if you convert/keep them as a single string 
#			of white space delimited values.
#			I call these strings that contain white space delimited values "fake lists".
# 
function(TO_LIST_SPACES _LIST_NAME OUTPUT_VAR)
  set(NEW_LIST_SPACE)
  foreach(ITEM ${${_LIST_NAME}})
    set(NEW_LIST_SPACE "${NEW_LIST_SPACE} ${ITEM}")
  endforeach()
  string(STRIP ${NEW_LIST_SPACE} NEW_LIST_SPACE)
  #message(  "${NEW_LIST_SPACE}" )
  #message(  ${OUTPUT_VAR}  )
  set(${OUTPUT_VAR} "${NEW_LIST_SPACE}" PARENT_SCOPE)
endfunction()

# 
#	Appends a list of items to a "fake-list" (a "fake-list" is a
#		string with whitespace delimited values).
#	It will only append a value if that values does not exists with in "fake list"
# 
#	Example:
#		LIST_SPACES_APPEND_ONCE( "CMAKE_C_FLAGS" "-DO0" "-Dg" "-DG" )
# 
function(LIST_SPACES_APPEND_ONCE LIST_NAME)
	set(_LIST)
  string(REPLACE " " ";" _LIST ${${LIST_NAME}})
  list(APPEND _LIST ${ARGN})
  list(REMOVE_DUPLICATES _LIST)
  to_list_spaces(_LIST NEW_LIST_SPACE)
  set(${LIST_NAME} "${NEW_LIST_SPACE}" PARENT_SCOPE)
endfunction()

# 
# Convenience function that does the same as LIST(FIND ...) but with a TRUE/FALSE return value.
# Ex: IN_STR_LIST("CMAKE_C_FLAGS" "-DO3" "WAS_FOUND")  (variable WAS_FOUND will be True or False)
# 
function(IN_STR_LIST LIST_NAME ITEM_SEARCHED RETVAL)
  list(FIND ${LIST_NAME} ${ITEM_SEARCHED} FIND_POS)
  if(${FIND_POS} EQUAL -1)
    set(${RETVAL} FALSE PARENT_SCOPE)
  else()
    set(${RETVAL} TRUE PARENT_SCOPE)
  endif()
endfunction()


# 
#	Helper Functions for Configuring User Input
# 
function(bool_yes_no__flag  var_name var_targets  value_if_true  value_if_false)
	foreach(var_target ${var_targets})
		if(${var_name})
			# if the value is an empty string, don't make the compiler flag at all
			if (NOT "${value_if_true}" STREQUAL "")
				append_compile_flag( "${value_if_true}" "${var_target}")
			endif()
		else()
			# if the value is an empty string, don't make the compiler flag at all
			if (NOT "${value_if_false}" STREQUAL "")
				append_compile_flag( "${value_if_false}"  "${var_target}")
			endif()
		endif(${var_name})
		set( "${var_target}" "${${var_target}}" PARENT_SCOPE)
	endforeach()
endfunction(bool_yes_no__flag var_name var_target value_if_true value_if_false)

# 
# 
# Pass this:
#	input_var_name: a string literal that is the name of variable that holds the value to examine
#	output_var_name: a string literal telling what "input_var_name" should be called in the final output list
#	var_targets: A semi-colon string literal holding a semi-colon separated list of list variable names 
#					that will be appended with the results of this function.
#	possible_values: a semi-colon string literal holding the possible allowed values for the input value
#						held in variable name given by 'var_name'
#	err: A string error message to display if the value passed in is not found in the list of allowed values.
# 
#   If the output_var_name is "" then the compile options output will be /D{valid input value}
# 
function(add_flag_restrict_to_list_or_error input_var_name output_var_name var_targets  possible_values err )
	foreach(var_target ${var_targets})
		#message("DEBUG add_flag_restrict_to_list_or_error: var_target ${var_target}")
		#message("DEBUG add_flag_restrict_to_list_or_error: input_var_name ${input_var_name}")
		#message("DEBUG add_flag_restrict_to_list_or_error: possible values string ${possible_values} input_value ${${input_var_name}}")
		string(REGEX MATCH "${${input_var_name}}" IS_MATCH ${possible_values}) #"${${possible_values}}" )
		if( NOT "${IS_MATCH}" STREQUAL "" )
			if( NOT "${output_var_name}" STREQUAL "" )
				append_compile_flag(  "${output_var_name}" "${${input_var_name}}" "${var_target}")
			else()
				append_compile_flag(  "${${input_var_name}}" "${var_target}")
			endif()
			
		else()
			ERROR_MESSAGE(FATAL_ERROR "${err} ${${input_var_name}} is an invalid value! (valid values:  ${possible_values} )" )
		endif()
		set( "${var_target}" "${${var_target}}" PARENT_SCOPE)
	endforeach()
endfunction( add_flag_restrict_to_list_or_error input_var_name output_var_name var_target  possible_values err )

# 
# 
# Pass this:
#	input_var_name: a string literal that is the name of variable that holds the value to examine
#	output_var_name: a string literal telling what "input_var_name" should be called in the final output list
#	var_targets: A semi-colon string literal holding a semi-colon separated list of list variable names 
#					that will be appended with the results of this function.
#	possible_values: a semi-colon string literal holding the possible allowed values for the input value
#						held in variable name given by 'var_name'
#	final_mapped_values: a semi-colon string literal holding the actual string written out in the final compiler flag for entry in "possible_values" 
# 
#	If the output_var_name is "" then the compilw options will be /D{matching value from list final_mapped_values}
# 
#	err: A string error message to display if the value passed in is not found in the list of allowed values.
# 
function(add_flag_restrict_to_list__map_values__or_error input_var_name output_var_name var_targets  possible_values final_mapped_values err )
	foreach(var_target ${var_targets})
		#message("DEBUG add_flag_restrict_to_list__map_values__or_error: var_target ${var_target}")
		#message("DEBUG add_flag_restrict_to_list__map_values__or_error: input_var_name ${input_var_name}")
		#message("DEBUG add_flag_restrict_to_list__map_values__or_error: possible values string ${possible_values} (mapped to ${final_mapped_values} ) input_value ${${input_var_name}}")
		string(REGEX MATCH "${${input_var_name}}" IS_MATCH ${possible_values}) #"${${possible_values}}" )
		if( NOT "${IS_MATCH}" STREQUAL "" )
			# This is the formal way to convert a space delimited or ' demlinted string into a list... (TODO: be carefull about this conversion and put it in all these functions)
			set(temp_list ${possible_values})
			set(temp_final_values  ${final_mapped_values})
			list(FIND temp_list "${${input_var_name}}" matched_idx )
			#message("Found the index for values ${${input_var_name}}: ${matched_idx}")
			list(GET temp_final_values "${matched_idx}" final_value_string )
			#message("Found the mappinf got index ${matched_idx} is: ${final_value_string}")
			if( NOT "${output_var_name}" STREQUAL "" )
				# if output_var_name is specified then set the compile flag like "/D=output_name=final_value_string"
				append_compile_flag(  "${output_var_name}" "${final_value_string}" "${var_target}")
			else()
				# if output_var_name is NOT specified then set theh compile flag like "/Dfinal_value_string"
				append_compile_flag("${final_value_string}" "${var_target}")
			endif()
		else()
			ERROR_MESSAGE("${err} ${${input_var_name}} is an invalid value! (valid values:  ${possible_values} )" )
		endif()
		set( "${var_target}" "${${var_target}}" PARENT_SCOPE)
	endforeach()
endfunction( add_flag_restrict_to_list__map_values__or_error input_var_name output_var_name var_target  possible_values final_mapped_values err )


# 
# Macro to make a nicer Fatal Error for the CMAKE_USER_MAKE_RULES_OVERRIDE
# 
macro(ERROR_MESSAGE message_text)
	message("")
	message("******************************************")
	message("Fatal Error                            ")
	message("******************************************")
	message( FATAL_ERROR "${message_text}" )
endmacro(ERROR_MESSAGE)

# 
# Dump all variables in this CMAKE instance
#	(includes otherwise 'hidden' variables.. ie. variables defined with the various FIND packages
#	that the GUI won't show)
# 
macro(DUMP_VARIABLES)
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
endmacro(DUMP_VARIABLES)

# 
# Dump all variables with whose name starts with a certian prefix. 
#	(includes otherwise 'hidden' variables.. ie. variables defined with the various FIND packages
#	that the GUI won't show)
# 
macro(DUMP_VARIABLES_STARTING_WITH prefix)
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        string(FIND "${_variableName}" "${prefix}" output)
        if(output EQUAL 0)
            message(STATUS "${_variableName}=${${_variableName}}")
        endif(output EQUAL 0)
    endforeach()
endmacro(DUMP_VARIABLES_STARTING_WITH)

# 
# Dump all variables with whose name contians the given pattern
#	(includes otherwise 'hidden' variables.. ie. variables defined with the various FIND packages
#	that the GUI won't show)
# 
macro(DUMP_VARIABLES_CONTAINING pattern)
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        string(FIND "${_variableName}" "${pattern}" output)
        if(output GREATER -1)
            message(STATUS "${_variableName}=${${_variableName}}")
        endif(output GREATER -1)
    endforeach()
endmacro(DUMP_VARIABLES_CONTAINING)
