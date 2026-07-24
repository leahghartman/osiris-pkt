"""
  Too much of this is just direct copy paste from the Namelistulator.. I will make a helper package to collect all the Gnu Make-like
     and boilerplate stuff into it's own reusable package.

  
  'go':
    Function is where execution starts.

  'parse_and_get_polymorphic_types':
    Scans the source code... and returns a list of all the Polymorphic Derived Types in the code.
  
  'main_work':
    Takes the list of 'Polymorphic Derived Types' and checks all declaration in the code to make sure
       that the 'class' keyeord is used these these types (verusus the improper 'type' keyword)


"""


import os
import os.path
import sys
import traceback
import imp
import traceback
import fnmatch

import argparse
import pickle
import datetime
import re

import fortransfigure.Fortransfigure as Fortransfigure
import fortransfigure.HighLevelParseAbstractions as HLPA

import time

from namelistulator.Utils import WRITE_OUT_OUTPUT_TIMESTEP, GET_TEMPLATE_FILE_PATH, READ_CONFIGURATION_SETUP_FILES, PICKLE_LOAD, PICKLE_DUMP1, do_logging, log
from fortransfigure.Utils import  GET_PREVIOUS_OUTPUT_TIMESTAMP, MERGE_DICTIONARIES__UPDATE_BASE
from fortransfigure.SourceCodeEditUtils import FileBeingAltered, recursive_file_find, find_original_source_file


# TODO: take this out
FORTRANSFIGURE = None


"""
###########
    Primary bootstrap called by command line and programmatic interfaces to stat things going.
    If you reading this code, START HERE!
###########
"""
def go( config_path, preprocessed_source_path, source_path, 
        include_childless_classes):
  
  configurationOptions = setup( config_path )
  _polymorphic_types_ = parse_and_get_polymorphic_types( configurationOptions, include_childless_classes)

  print( " Here are all the polymorphic derived types seen in the code: ")
  _polymorphic_types_.sort(key=lambda x: x.getInheritancePathAsStr(), reverse=False)
  for p in _polymorphic_types_:
    print(p.getInheritancePathAsStr())

  print("")
  if( source_path[-1] == '\\'):
    source_path = source_path[0:-1]
  if( source_path[-1] == '/'):
    source_path = source_path[0:-1]
  source_path = os.path.normpath( source_path )

  main_work(_polymorphic_types_, source_path)

"""
###########
    The main code where we look trough the source files....
###########
"""
def main_work(  _polymorphic_types_ , source_path ):

  varManipulator = HLPA.GhettoVariableManiulator()
  # I have to finish up the global namespace scope object.. 
  #   that will fix this limitation.. but for now it doesn't matter since
  #   Osiris does not repeat type names.
  #
  #  So for now, we assume that there are no two
  #    derived types with the same name. And we check this.
  class_names = [x.name() for x in _polymorphic_types_ ]
  class_names_set = { x.name() for x in _polymorphic_types_ }
  if len(class_names_set) != len(class_names):
    raise ValueError("The name of a derived type is used multiple times (repeated). "
                     "This code cannot distinguish between these types with the same name.")

  for f in FORTRANSFIGURE.FileLevelParsedFileRancher.parsedFilelevelItems:
    out = []
    
    # The __altered main is leaking up from Wrapacious here.. not sure why.. weird..
    # TODO: fix this
    if  f._sourceFileName.find("__altered") > -1:
      continue
    pass

    #if  not f._sourceFileName.find("os-neutral.") > -1:
    #  continue


    alteredFile = FileBeingAltered( f._sourceFileName )

    ####################
    # Check all the variables within Subroutines/Function (i.e.Procedures) that ARE IN THE ROOT
    #   Namespace (that is, not within a Module another procedure.)
    ####################
    (funcs, subs, funcssubs) = FORTRANSFIGURE.FileLevelParsedFileRancher.getAllRootLevelSubsAndFuncs(fromFiles=[f])
    deriviedTypeDeclarations = FORTRANSFIGURE.FileLevelParsedFileRancher.getAllDerivedTypes(include_polymorphic=False,fromFiles=[f])
    polyDeriviedTypeDeclarations = FORTRANSFIGURE.FileLevelParsedFileRancher.getPolyDerivedTypes(fromFiles=[f])
    modulesDeclarations =  FORTRANSFIGURE.FileLevelParsedFileRancher.getAllModules(fromFiles=[f])
    for entity in funcssubs:
      alteredFile.startChangesProceedure(entity.name())
      for var in entity.getAllVars():
        if var._isDerivedType and ( var.getTypeName() in class_names):
          if not var.isDeclaredPolymorphic():
            if var.isDummy():
              out.append( "\t Dummy Variable '%s' "
                        "\n\t       in procedure '%s'" % (var.name(),entity.name()))
              alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(entity.name()), r'class\1')
            else:
              out.append( "\t Local Variable '%s' "
                        "\n\t       in procedure '%s'" % (var.name(), entity.name()))
              if(var.isPointer()):
                alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(entity.name()), r'class\1')
              else:
                # This is dangerous to replace... changing the local variable requires that the variable also be
                #   a pointer and thus 'allocated' explicitly.. this should be handled by a human so we 
                #   will simple place a comment
                as_pointer = varManipulator.typeDeclarationToClass(varManipulator.addPointerAttr( var.tofortran() ))
                #as_pointer = varManipulator.addPointerAttr( var.tofortran() )
                as_pointer = re.sub(r'type([^\n]*?%s)' %(var.name()), r'class\1', as_pointer)
                alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()),
                                             r'type\1\n'
                                             r'! CLASS UP THE CODE! THE ABOVE VARIABLE - WHEN CHANGED TO "class" '
                                               r'HAS TO BE MADE A POINTER AND ALOOCATED!:\n'
                                             r'!%s\n' %(as_pointer,))
                                           
          pass

      # Then any interfaces shoved within a procedure.
      for iface in entity._interfaceStatements:
        procInterface = iface._procedureInterface
        if procInterface is None:
          continue
        alteredFile.startChangesInterface( iface.name(), procedure_name = entity.name() )
        for var in procInterface.getAllVars():
          if var._isDerivedType and ( var.getTypeName() in class_names):
            if not var.isDeclaredPolymorphic():
              message = str("\t Dummy Variable '%s' "
                          "\n\t       in interface '%s' "
                          "\n\t       in procedure '%s'"% (var.name(), iface.name(), entity.name()))
              if var.isDummy():
                out.append( message)
                alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()), r'class\1')
              else:
                raise ValueError("A Error occurred.. the parser gave impossible results: "
                                 "an interface cannot have non-dummy variables.\n\n%s" % ( message,))
    pass
    
    ####################
    # Check all the components within all the various Derived Types
    ####################
    all_types = polyDeriviedTypeDeclarations + deriviedTypeDeclarations
    for entity in all_types:
      alteredFile.startChangesTypeDeclaration( entity.name() )
      for var in entity.type_bound_data:
        if var._isDerivedType and ( var.getTypeName() in class_names):
          if not var.isDeclaredPolymorphic():
            out.append( "\t Variable       '%s' "
                      "\n\t       part of type '%s') " % (var.name(), var.getTypeName()))
            alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()), r'class\1')
          pass
        pass
      pass
    pass

    ####################
    # Check everything within Modules
    ####################
    #print( alteredFile.full_region.text)

    for module in modulesDeclarations:

      #----------------------
      # Check the Subroutines/Functions (i.e. Procedures) within a Module
      #----------------------
      all_procs = module._subroutines + module._functions
      # First the variables
      for entity in all_procs:
        for var in entity.getAllVars():
          if var._isDerivedType and ( var.getTypeName() in class_names):
            if not var.isDeclaredPolymorphic():
              alteredFile.startChangesProceedure( entity.name() )
              if var.isDummy():
                out.append( "\t Dummy Variable '%s' "
                          "\n\t       in procedure '%s' "
                          "\n\t       in module    '%s'" % (var.name(),entity.name(),module.name()))
                alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()), r'class\1')
              else:
                out.append( "\t Local Variable '%s' "
                          "\n\t       in procedure '%s' "
                          "\n\t       in module    '%s'" % (var.name(), entity.name(),module.name()))
                alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()), 
                                             r'\n! CLASS UP THE CODE: THIS VARIABLE - WHEN CHANGED TO \"class\" '
                                                  'HAS TO BE MADE A POINTER AND ALOOCATED!:\n'
                                                  '! class\1')
            pass
        pass
        # Then any interfaces shoved within a procedure.
        for iface in entity._interfaceStatements:
          procInterface = iface._procedureInterface
          if procInterface is None:
            continue
          alteredFile.startChangesInterface( iface.name(), procedure_name = entity.name(), module_name = module.name() )
          for var in procInterface.getAllVars():
            if var._isDerivedType and ( var.getTypeName() in class_names):
              if not var.isDeclaredPolymorphic():
                message = str("\t Dummy Variable '%s'"
                            "\n\t       in interface '%s' "
                            "\n\t       in procedure '%s' "
                            "\n\t       in module    '%s'" % (var.name(), iface.name(), entity.name(). module.name()))
                if var.isDummy():
                  out.append( message)
                  alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()), r'class\1')
                else:
                  raise ValueError("A Error occurred.. the parser gave impossible results: "
                                   "an interface cannot have non-dummy variables.\n\n%s" % ( message,))

      #----------------------
      # Check the Specific Interface Blocks within a Module
      #----------------------
      for entity in module._interfaces:
        procInterface = entity._procedureInterface
        if procInterface is None:
          continue

        for var in procInterface.getAllVars():
          if var._isDerivedType and ( var.getTypeName() in class_names):
            alteredFile.startChangesInterface(  entity.name(), module_name = module.name() )
            if not var.isDeclaredPolymorphic():
              message = str("\t Dummy Variable '%s' is "
                          "\n\t       in interface '%s' "
                          "\n\t       in module    '%s'" % (var.name(), entity.name(), module.name()))
              if var.isDummy():
                out.append( message )
                alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()), r'class\1')
              else:
                raise ValueError("A Error occurred.. the parser gave impossible results: "
                                 "an interface cannot have non-dummy variables.\n\n%s" % ( message,))
            pass
        pass
      pass

      #----------------------
      # Check the Module level Variables
      #----------------------
      for var in module._declarations:
        alteredFile.startChangesModuleDeclarations( module.name() )
        if var._isDerivedType and ( var.getTypeName() in class_names):
          if not var.isDeclaredPolymorphic():
            out.append( "\t Variable       '%s' (module-level variable) "
                      "\n\t       in module    '%s'" % (var.name(), module.name()))
            alteredFile.searchAndReplace(r'type([^\n]*?%s)' %(var.name()), r'class\1')
          pass
        pass
      pass
      alteredFile.commitChanges()

      # I disabled output... uncomment the following line to write changes back to source code files.
      #alteredFile.saveToOriginalSourceFiles()

    # Print out our findings, if any, for this file.
    if len(out) > 0:
      print("----In file %s" % (f._sourceFileName,))
      for line in out:
        print(line)
    pass
  pass

def find_polymorphic_type_right_now( _polymorphic_types_, ignoreItems, include_childless_classes = False ):

  do_logging(True)
  print("")
  print("")
  print("--------------  Finding all polymorphic derived types -----------------------")
  print("")
  polys = FORTRANSFIGURE.FileLevelParsedFileRancher.getPolyDerivedTypes(include_private = True)
  polys_temp = []
  ignored = []
  for entity in polys:
    if(entity.name() in ignoreItems):
      ignored.append("      '%s'"  % (entity.name(),) )
      continue
    polys_temp.append( entity )
  do_logging(False)


  if len(ignored) > 0:
    print("   Note: ignoring the following types as requested by the 'ignore items' config option:")
    for ig in ignored:
      print(ig)
    pass


  # we also have the optional parameter to filter out any derived types that do no actually have any subclasses.
  #  We will just use path manipulations to calculate this for now..
  #  TODO: make a nicer API function for this.. is a common thing to do.
  excluded_poly_types_list = []
  if not include_childless_classes:
    full_class_paths = [x.getInheritancePathAsStr() 
                        for x in polys_temp]  
    all_base_classes = [x 
                        for x in polys_temp 
                        if x.getInheritancePathAsStr().count("/") == 1]
    valid_base_classes = []
    for base_class in all_base_classes:
      base_class_search_name = base_class.name() + "/"
      temp = list(map(lambda x: x.find(base_class_search_name), full_class_paths))
      number_subclasses = sum(i > -1 for i in temp)
      if number_subclasses > 0:
        valid_base_classes.append( base_class.name() )
    # now make a list that has all these bases classes and their children.
    final_poly_types_list =    [x for x in polys_temp if x.getInheritancePathAsStr().split('/')[1] in valid_base_classes]
    excluded_poly_types_list = [x for x in polys_temp if x.getInheritancePathAsStr().split('/')[1] not in valid_base_classes]
    polys_temp = final_poly_types_list

  if len(excluded_poly_types_list) > 0:
    print("")
    print("   Note: ignoring the following types because 'include-childless-classes' option is False:")
    print("         (to include these classes pass the --include-childless-classes option to run.py)")
    for ept in excluded_poly_types_list:
      print("      '%s'" % (ept.name(),))
  
  print("")
  print("-----------------------------------------------------------------------------")
  _polymorphic_types_.extend( polys_temp )



"""
###########
    Get down to work....
      Assemble list of input files, check if input is out of date etc..

    Construct list of input files to analyze. Then to 'make'-like chores to
      if any work needs to be done. We check
      - Did any input files change since list run? If not, no need to do anything
      - UNLESS: are any output file missing? If so, then we do need to run.
###########
"""
def parse_and_get_polymorphic_types (configurationOptions, exclude_childless_classes = False):

  # construct the list of input files we will process
  #   also returned is the list of files excluded from the list by
  #   the filepatterns in the 'SourceFilesExcludes' input setting.
  (sourcefiles, excludedfiles) = Fortransfigure.constructInputSourceFileList()
  
  # TODO: when I make output, put checks in here.
  any_missing = False
  #Fortransfigure.check_if_input_work_is_needed( any_output_files_missing = any_missing)


  # Do all the parsing of the input Fortran files.
  Fortransfigure.start()
  
  # similar timestamp, set at the end of the program, to stamp when the output was written.
  #   In production, this timestamp doesn't do much and we could just have the 'PARSE' timestamp.
  #   But this timestamp is usefull in debugging etc.. and doesn't hurt. If an exception occurs
  #   the timestamp isn't written.. so it keeps things consistent on reruns.
  # It can be disabled in the configuration.
  (lastDatetime_output, timestampFound_output) = GET_PREVIOUS_OUTPUT_TIMESTAMP(configurationOptions['regenerate_if_output_wrapper_changed'])
  

  _polymorphic_types_ = []

  find_polymorphic_type_right_now(_polymorphic_types_, configurationOptions['ignore items'], exclude_childless_classes )

  WRITE_OUT_OUTPUT_TIMESTEP()

  return _polymorphic_types_


"""
###########
    Reads configuration options and setup the final, top level
        exception handling code.
###########
"""
def setup(config_filename):
  from namelistulator import SetupConfiguration
  global FORTRANSFIGURE
  
  try:
    print ("")
    print ("Reading the options from %s " % config_filename)
    incomingConfigurationOptions = READ_CONFIGURATION_SETUP_FILES( config_filename )
  except:
    print ("Error reading configuration options. Aborting.. ")
    traceback.print_exc(file=sys.stdout)
    exit(-1)

  #
  # Create, init, and configure the embedded Fortransfigure object.
  #
  FORTRANSFIGURE = Fortransfigure.STATE
  fortansfigureOptions = Fortransfigure.init( incomingConfigurationOptions )

  #
  # Add to the input configuration
  #   by filling in default values for any unspecified parameters
  defaultConfigurationOptions = SetupConfiguration.NamelistilatorConfiguration()
  configurationOptions = MERGE_DICTIONARIES__UPDATE_BASE( incomingConfigurationOptions, 
                                                          defaultConfigurationOptions)

  #
  # Finish input configuration
  #   by merging Fortranfigure's input into that of the Namelistulator to get
  #   the final - with all defaults fill in - input.
  configurationOptions.update( FORTRANSFIGURE.configurationOptions )

  return configurationOptions




  

