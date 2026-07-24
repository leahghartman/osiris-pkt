"""
-----------------------------------------------------------------------------
              Input Deck (Configuration) for the Namelistulator

              We will base this off the input deck for wrapacious 
              as it uses the same underlying machinery. Then customize.
-----------------------------------------------------------------------------
"""

# -------------------------- Setup -------------------
import os.path
import sys

# find 'wrapacious' input configuration file. Then add the directory to the
#  python path
current_path = os.path.dirname(__file__)
wrapacious_dir_relative = os.path.join(current_path, "../../../wrapacious")
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
                    Namelistulator specific settings ARE HERE
                    LOOK HERE
  -----------------------------------------------------------------------------
  """
  #deck['SourceFilesGlobs'] = ["preprocessed/*.f90", "preprocessed/*.f03" ]
  deck['bin_dir'] = "./"
  deck['always_read_input'] = False
  
  return deck

