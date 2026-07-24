"""
  Compile with syou app with the full boost availiable.. add in the -H flag (with GCC) and redirect the output to a file:
     make 2> boost_out.txt
  Then set INPUT_FILENAME to whatever filename you used.
  Set BOOST_BASE_DIR to whatever the base directory for Boost is.
  Set OUTPUT_DIR to wahtever dir you want to copy jsut the used boost header too.
  
  Note: In Windows with the Microsoft compiler, use /showIncludes rather then -H
        The -H or /showIncludes print out all the headers included.. this script takes that
        and jsut sifts though it, see which headers are below BOOST_BASE_DIR, and copies the relative path
        into the OUTPUT_DIR
"""
INPUT_FILENAME = "boost_out.txt"
BOOST_BASE_DIR = "/mnt/c/lab_code/slim_boost"
OUTPUT_DIR     = "/mnt/c/lab_code/slim_boost/minimal"

import os
from shutil import copyfile

def main():

  a = open(INPUT_FILENAME, 'r')
  lines = a.readlines()
  a.close()

  for line in lines:
    try:
      do_a_line(line)
    except:
      pass
  pass

def do_a_line(line):
  if line.find(  BOOST_BASE_DIR ) > -1:
    idx = line.find(  BOOST_BASE_DIR ) + len( BOOST_BASE_DIR )
    relative_dir = line[idx:]
    if not relative_dir.startswith("/"):
      relative_dir = "/" + relative_dir.strip()
    relative_dir = relative_dir.strip()

    source_path = BOOST_BASE_DIR + relative_dir
    output_path = OUTPUT_DIR + relative_dir #os.path.join(OUTPUT_DIR, relative_dir)
    print(source_path)
    print("    to %s" %(output_path,))

    try:
      os.makedirs(os.path.dirname(output_path))
    except:
      pass
    copyfile(source_path, output_path )
  pass

main()