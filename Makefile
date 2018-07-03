COMPILER   = gfortran
COMPFLAGS  = -c 
OBJECTS    = bmo_model.o bmo_evolution.o parameters.o parameters_bmo.o splat.o drspln.o
EXECUTABLE = thermal_history_bmo


all : $(OBJECTS) $(EXECUTABLE).f90
	$(COMPILER) $(COMPFLAGS) $(EXECUTABLE).f90
	$(COMPILER) -o $(EXECUTABLE) $(EXECUTABLE).o $(OBJECTS)

bmo_evolution.o : bmo_evolution.f90 bmo_model.o parameters.o
	$(COMPILER) $(COMPFLAGS) bmo_evolution.f90

bmo_model.o : bmo_model.f90 parameters_bmo.o parameters.o 
	$(COMPILER) $(COMPFLAGS) bmo_model.f90

parameters_bmo.o : parameters_bmo.f90 
	$(COMPILER) $(COMPFLAGS) parameters_bmo.f90

parameters.o : parameters.f90
	$(COMPILER) $(COMPFLAGS) parameters.f90

splat.o : splat.f
	$(COMPILER) $(COMPFLAGS) splat.f
drspln.o : drspln.f
	$(COMPILER) $(COMPFLAGS) drspln.f

clean        :
	rm -rf *.o *.mod $(EXECUTABLE)

