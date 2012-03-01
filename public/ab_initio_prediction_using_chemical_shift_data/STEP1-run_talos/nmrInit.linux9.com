#!/bin/csh -f

#
# NMRPipe System Environment script, Sun Aug 26 19:29:14 PDT 2007.

#
# Override software expiration.

setenv NMR_CONT CORRECT


#
# NMRPipe I/O Optimization.
# These settings may improve NMRPipe pipeline performance:

setenv NMR_IO_TIMEOUT 0
setenv NMR_IO_SELECT  0
setenv NMR_AUTOSWAP   1

alias u      'chmod a+rx *.com *.tcl'

#set   history = 256
#alias h      'history'
#alias rm     'rm -i'
#alias cd     'cd \!*; set prompt = "$cwd% "'
#alias dirs   'ls -l | egrep -e ^d'
#alias links  'ls -l | egrep -e ^l'

if ($?MANPATH) then
   setenv MANPATH \
    /work/shenyang/nmrPipe//man:/usr/share/man:/usr/local/man:$MANPATH
else
   setenv MANPATH \
    /work/shenyang/nmrPipe//man:/usr/share/man:/usr/local/man
endif

set path = \
 (. /work/shenyang/nmrPipe//nmrbin.linux9 /work/shenyang/nmrPipe//com $path)


if ($?LD_LIBRARY_PATH) then
   setenv LD_LIBRARY_PATH /work/shenyang/nmrPipe//xview/linux9/lib:/usr/openwin/lib:/work/shenyang/nmrPipe//nmrbin.linux9/lib:/usr/local/lib:$LD_LIBRARY_PATH
else
   setenv LD_LIBRARY_PATH /work/shenyang/nmrPipe//xview/linux9/lib:/usr/openwin/lib:/work/shenyang/nmrPipe//nmrbin.linux9/lib:/usr/local/lib
endif


if (!($?OPENWINHOME)) then
   if (-d /work/shenyang/nmrPipe//nmrbin.linux9/openwin) then
      setenv OPENWINHOME /work/shenyang/nmrPipe//nmrbin.linux9/openwin
   endif
endif

setenv NMRCHECK        ALL
setenv NMRBASE         /work/shenyang/nmrPipe/
setenv NMRTXT          /work/shenyang/nmrPipe//nmrtxt
setenv NMRBIN          /work/shenyang/nmrPipe//nmrbin.linux9
setenv TCLPATH         /work/shenyang/nmrPipe//com
setenv TALOS_DIR       /work/shenyang/nmrPipe//talos


setenv NMR_TCLTK8 TRUE

if ($?NMR_TCLTK8) then
   setenv TCL_LIBRARY     /work/shenyang/nmrPipe//nmrtcl/tcl8.4
   setenv TK_LIBRARY      /work/shenyang/nmrPipe//nmrtcl/tk8.4
   setenv BLT_LIBRARY     /work/shenyang/nmrPipe//nmrtcl/blt2.4z
   setenv NMRPIPE_TCL_LIB /work/shenyang/nmrPipe//nmrtcl/tcl8.4
   setenv NMRPIPE_TK_LIB  /work/shenyang/nmrPipe//nmrtcl/tk8.4
   setenv NMRPIPE_BLT_LIB /work/shenyang/nmrPipe//nmrtcl/blt2.4z
else
   setenv TCL_LIBRARY     /work/shenyang/nmrPipe//nmrtcl/tcl7.6
   setenv TK_LIBRARY      /work/shenyang/nmrPipe//nmrtcl/tk4.2
   setenv BLT_LIBRARY     /work/shenyang/nmrPipe//nmrtcl/blt2.4
   setenv NMRPIPE_TCL_LIB /work/shenyang/nmrPipe//nmrtcl/tcl7.6
   setenv NMRPIPE_TK_LIB  /work/shenyang/nmrPipe//nmrtcl/tk4.2
   setenv NMRPIPE_BLT_LIB /work/shenyang/nmrPipe//nmrtcl/blt2.4
endif

# setenv NMRPIPE_SMALLFONT  "-adobe-helvetica-medium-r-*-*-*-100-*-*-*-*-*-*"
# setenv NMRPIPE_BIGFONT    "-adobe-helvetica-medium-r-*-*-*-180-*-*-*-*-*-*"
# setenv NMRPIPE_STDFONT    "-adobe-helvetica-medium-r-*-*-*-120-*-*-*-*-*-*"
# setenv NMRPIPE_BOLDFONT   "-adobe-helvetica-bold-r-*-*-*-120-*-*-*-*-*-*"
# setenv NMRPIPE_FIXEDFONT  "-*-courier-medium-r-*-*-*-120-*-*-*-*-*-*"
# setenv NMRPIPE_TYPINGFONT "-*-courier-medium-r-*-*-*-100-*-*-*-*-*-*"
# setenv NMRPIPE_AXISFONT   "-*-lucida-bold-r-*-*-*-120-*-*-*-*-*-*"

