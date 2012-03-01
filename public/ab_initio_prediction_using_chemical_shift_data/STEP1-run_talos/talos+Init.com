#!/bin/csh 

#
#
# Example usage:
#
#   if (-e /u/shenyang/com/talos+Init.com) then
#      source /u/shenyang/com/talos+Init.com
#   endif
#   

## PLEASE CHECK THE TALOS INSTALLATION DIRECTORY DEFINED NEXT LINE!!!
setenv talospDir  /work/shenyang/TALOS+

set path = (. $path ${talospDir})
