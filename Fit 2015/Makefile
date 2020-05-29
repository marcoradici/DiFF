CORE = EvolveD1.cc
COREo = EvolveD1.o

FF = zd1u.f90 zd1d.f90 zd1s.f90 zd1c.f90
FFo = zd1u.o zd1d.o zd1s.o zd1c.o

MODULES = evo_grid.f90
MODULESo = evo_grid.o

#SOURCE = DiFFD1grid.f90
#SOURCEo = DiFFD1grid.o
SOURCE = test.f90
SOURCEo = test.o

#EXE = DiFFD1grid.exe
EXE = test.exe

F90 = gfortran
CXX = g++
FLAGS = -O3

compile:
	$(F90) $(FLAGS) -c $(FF) $(MODULES) $(SOURCE)
	$(CXX) $(FLAGS) -fPIC -std=c++11 `apfelxx-config --cppflags` -c $(CORE)

link:
	$(F90) -o $(EXE) $(COREo) $(FFo) $(MODULESo) $(SOURCEo) -lc++ -lstdc++ `apfelxx-config --ldflags`

cleanmod:
	rm *.mod
	
cleano:
	rm *.o
	
cleanall:
	rm *.o *.mod *.exe
