# define the C++ compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -Wall -g

# define any directories containing header files other than /usr/include
INCLUDES = -I../include -I/home/vultor/PROGRAMAS/dlib-18.16

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS =  -L/dlib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lm

# define the C++ source files
SRCS = CheckRessonanceType.cpp CheckTypeOfIntegralToBeCalculated.cpp CoreMatrixCalculations.cpp DiagonalizeDlib.cpp DiatomicOverlaps.cpp DiisProcedure.cpp DlibLinearSystem.cpp EularAngle.cpp FourCenterIntegrals.cpp FourCenterMatrix.cpp main.cpp MatrixDiagonalization.cpp Matrixop.cpp MatrixRotation.cpp Molecule.cpp OptimizeWithDlib.cpp Params.cpp PopleOverlaps.cpp PrintAll.cpp ReadMopFile.cpp ReadqInput.cpp ReadXyzFile.cpp RhfMethod.cpp ScfCycle.cpp ScfMethod.cpp ScfProcedure.cpp UhfMethod.cpp ZmatoCart.cpp FitnessOpt.cpp  


# define the C++ object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
MAIN = qRM1.exe

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  qrm1.exe foi compilado com sucesso:

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it



