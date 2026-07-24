"""
-----------------------------------------------------------------------------
              Input Deck (Configuration) for the Namelistulator

              We will base this off the input deck for wrapacious 
              as it uses the same underlying machinery. Then customize.
              (see ./toold/wrapacious/wrapicious_input_deck.py)
-----------------------------------------------------------------------------
"""

# -------------------------- Setup -------------------
import os.path
import sys

# find 'wrapacious' input configuration file. Then add the directory to the
#  python path
current_path = os.path.dirname(__file__)
wrapacious_dir_relative = os.path.join(current_path, "../../wrapacious")
wrapaciour_dir = os.path.abspath( wrapacious_dir_relative )

# add the path to the python path
sys.path.append(wrapaciour_dir)

# 'inherit' from the wrapacious input deck.
import wrapicious_input_deck as wid
deck = wid.WrapaciousConfiguration()


# This function is caLLed from namelistulator. Customize the inputdeck here.
def NamelistilatorConfiguration():
  """
  -----------------------------------------------------------------------------
                    LOOK HERE
                    Class Up The Code specific settings ARE HERE
                    LOOK HERE
  -----------------------------------------------------------------------------
  """
  
  deck['ignore items'].append('t_vdf')

  # get rid of the any altered source files that Wrapacious makes 
  #  (e.g. '__altered__os-main.f90' any other other files starting
  #   with '__altered__')
  deck['SourceFilesExcludes'].append("__altered*")
  deck["alterOriginalSourceCode"] =  False
  
  return deck

