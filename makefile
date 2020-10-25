
OBJECTS = constants.o \
		  io.o \
		  calculate.o \
		  main.o

MODULES = constants.mod \
		  io.mod \
		  calculate.mod

# FLAGS
# -c flag means compile to one file (very common if program is in many files)
# -o FILENAME means rename output from a.out to FILENAME.exe
# -g generates extra debugging information usable by GDB
# -03 level 3 optimisation for compling code

FFLAGS = -g  -Wall -fcheck=all -llapack -lblas

.PHONY: run clean

run: main
	./main

main: $(MODULES) $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o main

%.o : %.f95
	gfortran $(FFLAGS) -c $<

%.mod: %.f95
	gfortran $(FFLAGS) -c $<

clean:
	rm -f *.o main *.mod
