ifort -c 1dsteady.f90
gcc -c realarray-c.cpp
ifort 1dsteady.o realarray-c.o
