#!/bin/bash
#
# Calculates the dependencies of Osiris Modules
#
# Author: Axel Huebl, a.huebl@hzdr.de
# Date:   Wed, 24th of Jan 2013
#

# configure here
#
FEXT=".f90"
OEXT=".o"
MEXT=".mod"
#quiet="enabled"

MODNAME=$1
FPP=$2

if [ -z "$FPP" ]; then
  FPP="gcc -C -E -x assembler-with-cpp"
fi


# Module naming conventions
#
echo -n "Module file naming convention: "

if [ -z "$MODNAME" ]
then
  echo -n "(none) - "
  MODNAME="STD"
fi
if [ "$MODNAME" = "STD" ]
then
  echo $MODNAME
else
  echo "($MODNAME) is a unknown MOD_NAME."
  exit 1
fi


# Progress indicator
#
if [ -z "$quiet" ]
then
  echo -n "Calculating module dependencies "
fi


# Replace the extension .f90 with .o
#
# \param  filename.f90
# \return filename.o
#
function fsrc2obj
{
  echo $1 | sed -e 's/..\/source\///g' | sed -e "s/$FEXT$/$OEXT/g"
}

# Just playing around... Set $quiet to "enabled" to
# disable this progress bar.
#
# \param "n|d|.." no new line OR disable OR newline
#
function showprogress
{
  opt=$1$quiet
  if [ "$opt" == "" ];  then echo    "."; fi
  if [ "$opt" == "n" ]; then echo -n "."; fi
}

# remove old files
#
rm -rf f_deps.mk f_prov.txt

# progress bar
#
showprogress "n"

# Find all provided modules and list them in a file structured: ###############
#   source_file SPACE modules
# To force you to write clean code, just one module definition
# per file is supported :)
#
for f in `ls ../source/*$FEXT`
do
  # object file name
  f_obj=`fsrc2obj $f`
  # echo -n $f_obj":"

  # which modules provided?
  prov_mods=`grep -i "^[[:space:]]*module[[:space:]+]" $f | \
    grep -iv "procedure" | awk '{print $2}'`
  echo $f" "$prov_mods >> f_prov.txt
done

# progress bar
showprogress "n"

#
for f in `ls ../source/*$FEXT`
do
  # "Who am I?" (object and if provided module)
  f_obj=`fsrc2obj $f`
  f_mod=`grep "$f" f_prov.txt | awk '{print $2}'`
  if [ -n "$f_mod" ];
  then
    f_mod=`echo " "$f_mod"$MEXT"`
    echo "$f_mod: $f_obj" >> f_deps.mk
  fi

  echo -n "$f_obj$f_mod:" >> f_deps.mk

  # which mods used?
  dep_mods=`${FPP} $f | grep -i "^[[:space:]]*use[[:space:]+]" | \
    awk '{print $2}' | sed -e 's/,//g'`
  #echo $f ": "$dep_mods

  for mod in $dep_mods
  do
    # find corresponding files
    dep_files=`grep "$mod$" f_prov.txt | awk '{print $1}'`
    #echo $f" ("$mod") : " $dep_files

    # format output
    if [ -n "$dep_files" ]; then echo -n " $mod$MEXT" >> f_deps.mk; fi
  done

  # new line
  echo >> f_deps.mk
done

# progress bar
showprogress
