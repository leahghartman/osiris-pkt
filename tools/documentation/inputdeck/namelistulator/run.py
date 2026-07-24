"""
-----------------------------------------------------------------------------
                           Namelistulator

          The Namelistulator is a tool that scans the Osiris source code
          to help generate documentation for Osiris's input parameters.
          
          This is a small convience file that calls the Namelistulator (the
          Namelistulator lives in the github repository 'os-tools')
-----------------------------------------------------------------------------
"""


import argparse
import imp
import os
import os.path
import sys
import warnings


def check_input( config_path, ostools_path, metadata_path, output_path ):
  # the 'imp' module throws unimportant warnings.. suppress.
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # --------------------
    #  Make Sure we can read the configuration file
    # --------------------
    if config_path is None:
      config_path = os.path.dirname(__file__)
      config_path = os.path.join(config_path, "namelistulator_config.py")

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
    #  Make Sure we can see and use the metadata we need.
    # --------------------
    if metadata_path is None:
      metadata_path = os.path.dirname(__file__)
      metadata_path = os.path.join(metadata_path, "..")
      metadata_path = os.path.abspath( metadata_path )
      metadata_path = os.path.join(metadata_path, "__namelist_ordering.py")
    try:
      temp = imp.load_source('namelistulator.input',  metadata_path )
    except:
      print("")
      print("   Error: Unable to import the need metadata!")
      print("")
      print("   We looked for the metadata in the directory: %s " %(metadata_path,))
      print("   Please use the --metadata-path option to specify a path to a valid metadata.")
      print("   Valid meta data consists: __namelist_ordering.py ")
      print("   This metadata should already be included in the Osiris repository at:")
      print("         tools/documentation/inputdeck")
      print("")
      exit(-1)

    # --------------------
    #  Make Sure we can see the output directory.
    # --------------------
    if output_path is None:
      output_path = os.path.dirname(__file__)
      output_path = os.path.join(output_path, "..")
      output_path = os.path.abspath( output_path )
      output_path = os.path.join(output_path, "generated")
    if not os.path.isdir(output_path):
      print("")
      print("   Error: Unable to find the output directory")
      print("")
      print("   We looked for the metadata in the directory: %s " %(output_path,))
      print("   Please use the --output-path option to specify a path to a valid output path.")
      print("")
      exit(-1)
  return (metadata_path, output_path )

if __name__ == "__main__":
  parser1 = argparse.ArgumentParser(description='Scans Osiris namelists to generate a list of input parameters')
  parser1.add_argument('--config-path', 
                       default=None, 
                       help='Filesystem path to the configuration file.')
  parser1.add_argument('--metadata-path', 
                       default=None, 
                       help="Filesystem path to the directory where the needed Osiris metadata python data files reside. (e.g. 'OsirisNamelistOrdering.py'")
  parser1.add_argument('--ostools-path', 
                       default=None, 
                       help="Filesystem path to the 'os-tools' distribution repository. ")
  parser1.add_argument('--output-path', 
                       default=None, 
                       help="Where the final generated html files whould be written.")
  parser1.add_argument('--use-cached', 
                       action='store_true', 
                       help="Used, if availiable, previously cached scanning results. Mainly used in development/debugging.")
  args = parser1.parse_args()

  # if the config_path is no specified, create a default one to check.
  config_path =  args.config_path
  if config_path is None:
    config_path = os.path.dirname(__file__)
    config_path = os.path.join(config_path, "namelistulator_config.py")

  # check the args
  (metadata_path, output_path ) = check_input( config_path, args.ostools_path, args.metadata_path, args.output_path )
  
  # run the program!!!
  import namelistulator.namelistulator as prog

  prog.go( config_path, metadata_path, output_path, args.use_cached )
  

















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
