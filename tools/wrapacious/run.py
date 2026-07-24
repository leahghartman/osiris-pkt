"""
-----------------------------------------------------------------------------
                           Wrapacious (python wrapper)
          
          This is a small convience file that calls the Wrapacaious
          (Wrapacious lives in the github repository 'os-tools')

          And checks that all the needed libraries are installed and helpful (hopefully)
          error messages are given.
-----------------------------------------------------------------------------
"""


import argparse
import imp
import os
import os.path
import sys
import warnings


def check_input( config_path ):
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
    except Exception as ex:
      print("")
      print("   The configuration file at:")
      print("     %s" % (config_path,) )
      print("   has found, but contains errors (is not a valid Python file). Please correct any errors.")
      print("")
      print("   ----------------------------------------------------------")
      print("   Python reported the following errors in that file:")
      print("   ----------------------------------------------------------")
      print("")
      print(ex)
      exit(-1)

    # --------------------
    #  Make Sure we can see and use the 'wrapacious' package in 'os-tools'
    # --------------------
    try:
      import wrapacious
    except:
      print("")
      print("   Error: Unable to importing 'wrapacious' from the 'osiris-tools' python package.")
      print("")
      print("   Please install the 'osiris-tools' python package and/or add 'os-tools' to your PYTHONPATH.")
      print("")
      print("   See https://github.com/GoLP-IST/osiris-tools for instructions on downloading 'osiris-tools'")
      print("")
      exit(-1)

    try:
      import mako
    except:
      print("")
      print("   Error: Unable to import the 'mako' REQUIRED python package.")
      print("   It's a common python package that can be installed in any of the usual ways (i.e using 'pip' or 'pip3'.")
      print("")

  return config_path

if __name__ == "__main__":
  parser1 = argparse.ArgumentParser(description='Simple launcher that checks that all libraries are installed when making python wrapper.')
  parser1.add_argument('--dry-run-for-gnu', 
                       default=None, 
                       help='Filesystem path to the configuration file.')
  parser1.add_argument('--config-path', 
                       default=None, 
                       help='Filesystem path to the configuration file.')
  args = parser1.parse_args()

  # if the config_path is no specified, create a default one to check.
  config_path =  args.config_path
  if config_path is None:
    config_path = os.path.dirname(__file__)
    config_path = os.path.join(config_path, "wrapicious_input_deck.py")

  # check the args
  print("")
  print("  ----------------------------Python Wrapper ---------------------------")
  print("  Using the python installation located at:")
  print("     '%s'" % (sys.executable,) )
  print("")
  print("  Checking that the Python installation has all the needed libraries:")

  config_path = check_input( config_path )
  print("                        Yes.. checks passed all ok.")
  print("  ----------------------------------------------------------------------")

  # run the program!!!
  import wrapacious.Wrapacious as w
  w.startup( config_path, args.dry_run_for_gnu)


















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
