################################################################################
#!/bin/sh
# Filename:    installFunctions.sh
# Author:      Mark A. Hunt (CFDRC)
# Date:        2024-07-05
# Description: This file contains common bash functions that are needed by the
#              install scripts. This file is to be sourced, not executed.
################################################################################



################################################################################
# Function used to find executable programs.
################################################################################
find_exec() {
  RETURN_VALUE=0
  for i in ${PATH//:/ }; do
    if [ -e $i/$1 ]; then
      RETURN_VALUE=$i
      break ;
    fi
  done
}



################################################################################
# Function used to find library files.
################################################################################
find_lib() {
    RETURN_VALUE=0
    for i in ${LD_LIBRARY_PATH//:/ }; do
      if [ -e $i/lib$1.a ]; then
        RETURN_VALUE=$i
        break;
      fi
      if [ -e $i/lib$1.so ]; then
        RETURN_VALUE=$i
        break ;
      fi
      if [ -e $i/lib$1.dylib ]; then
        RETURN_VALUE=$i
        break ;
      fi
    done
}