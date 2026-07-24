"""
-----------------------------------------------------------------------------
                           class_up_that_code

          This is a small convience file that calls checks input config arguments 
          and then calls the main code in class_up_that_code.py
-----------------------------------------------------------------------------
"""


import argparse
import imp
import os
import os.path
import sys
import warnings


def check_input( config_path, ostools_path, preprocessed_source_path, source_path ):
  # the 'imp' module throws unimportant warnings.. suppress.
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # --------------------
    #  Make Sure we can read the configuration file
    # --------------------
    if config_path is None:
      config_path = os.path.dirname(__file__)
      config_path = os.path.join(config_path, "class_up_that_code_config.py")

    if not os.path.exists( config_path ):
      print ("")
      print ("Error: the config file location (%s) does not exist." % (config_path,))

    try:
      temp = imp.load_source('namelistulatorconfig',  config_path )
    except:
      print("")
      print("   Unable to load the configuration file at:")
      print("     %s" % (config_path,) )
      print("   please use the --config-path options to specify a path to a valid config file.")
      exit(-1)

    # --------------------
    #  Make Sure we can see and use the namelistulator package in 'os-tools'
    # --------------------
    if ostools_path is not None:
      sys.path.append(ostools_path)
    try:
      import namelistulator
    except:
      print("")
      print("   Error: Unable to import the Namelistualtor!")
      print("")
      print("   Please use the --ostools-path options to specify a path to a valid to the 'os-tools' repository")
      print("   Or add 'os-tools' to your PYTHONPATH.")
      print("")
      exit(-1)

    # --------------------
    #  Make Sure we can see the source directory.
    # --------------------
    if source_path is None:
      source_path = os.path.dirname(__file__)
      source_path = os.path.join(source_path, "../../../source")
      source_path = os.path.abspath( source_path )
      source_path = os.path.normpath( source_path )
    if not os.path.isdir(source_path):
      print("")
      print("   Error: Unable to find the Osiris source code directory")
      print("")
      print("   We looked for the Osiris source code in the directory: %s " %(source_path,))
      print("   Please use the --source-path option to specify a path to a valid path.")
      print("")
      exit(-1)

  # --------------------
  #  Make Sure we can see the preprocessed source directory.
  # --------------------
  if preprocessed_source_path is None:
    preprocessed_source_path = os.getcwd()
    preprocessed_source_path = os.path.join(preprocessed_source_path, "./processed")
    preprocessed_source_path = os.path.abspath( preprocessed_source_path )
    preprocessed_source_path = os.path.normpath( preprocessed_source_path )
  if not os.path.isdir(preprocessed_source_path):
    preprocessed_source_path = None
  if preprocessed_source_path is None:
    preprocessed_source_path = os.path.dirname(__file__)
    preprocessed_source_path = os.path.join(source_path, "../../../source/build/preprocessed")
    preprocessed_source_path = os.path.abspath( preprocessed_source_path )
  if not os.path.isdir(preprocessed_source_path):
    print("")
    print("   Error: Unable to find the Osiris source code directory")
    print("")
    print("   We looked for the Osiris source code in the directory: %s " %(preprocessed_source_path,))
    print("   Please use the --preprocessed-source-path option to specify a path to a valid path.")
    print("")
    exit(-1)

  return (preprocessed_source_path, source_path )

if __name__ == "__main__":
  parser1 = argparse.ArgumentParser(description='Scans Osiris looking for places where Polymorphic types are declared as "type" versus "class"')
  parser1.add_argument('--config-path', 
                       default=None, 
                       help='Filesystem path to the configuration file.')
  parser1.add_argument('--ostools-path', 
                       default=None, 
                       help="Filesystem path to the 'os-tools' distribution repository if it is not found in the Python library path. ")
  parser1.add_argument('--preprocessed-source-path', 
                       default=None, 
                       help="Filesystem path to the a version of the code run through the preprocesor (such code is generated as part of building Osiris)")
  parser1.add_argument('--source-path', 
                       default=None, 
                       help="Filesystem path to the /source directory of the Osiris repository.")
  parser1.add_argument('--include-childless-classes', 
                       action='store_true', 
                       help="Have the program make changes and save them to the output directory.")
  args = parser1.parse_args()

  # if the config_path is no specified, create a default one to check.
  config_path =  args.config_path
  if config_path is None:
    config_path = os.path.dirname(__file__)
    config_path = os.path.join(config_path, "class_up_that_code_config.py")

  # check the args
  (preprocessed_source_path, source_path) = check_input( config_path, 
                                               args.ostools_path,
                                               args.preprocessed_source_path, 
                                               args.source_path )
  
  # run the program!!!
  import class_up_that_code as CUTC

  CUTC.go( config_path, preprocessed_source_path, source_path, args.include_childless_classes )

  

















  """
  # # We use python's Module import machinery to import the configuration. As such, we:
  # #   1. get config file's directory
  # #   2. Add this directory to the PYTHONPATH
  # #   3. get the config file's filename.
  # #   4. strip out the .py extension and import the result.
  directory = os.path.dirname(  args.config_path )
  if directory == "":
    directory = '.'
  module_name = os.path.basename(  args.config_path )
  if( module_name.rfind(".") > -1):
    idx = module_name.rfind(".")
    module_name = module_name[:idx]

  print "ee%srr" % (directory,)
  # # add to the PYTHONPATH
  sys.path.append(directory)
  """
