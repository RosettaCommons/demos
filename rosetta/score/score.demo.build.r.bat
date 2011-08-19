@echo off
rem Demo Release Build Script for MinGW GCC on 32-bit Windows

rem Adjust this source path for your system
rem set RosettaRoot=\Projects\Rosetta\mini\dev

rem Build the demo
g++ -pipe -std=c++98 -pedantic -Wall -Wextra -fmessage-length=0 -ffor-scope -fno-exceptions -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -O3 -DNDEBUG -s -malign-double -march=pentium4 score.demo.cc -o score.demo.exe -I. -I%RosettaRoot%\src -I%RosettaRoot%\src\platform\windows\xp\32\x86\gcc -L%RosettaRoot%\build\src\release\windows\xp\32\x86\gcc -L%RosettaRoot%\external\lib -lrosetta -lnumeric -lutility -lObjexxFCL -lz
