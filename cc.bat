echo on
prompt $g
rem
del %1.exe
gcc -march=x86-64 -O3 -fopenmp -o %1.exe %1.c -Wall -Wextra
REM

