## Generated by Migo's make.py
rm=/bin/rm -f
CC= g++
PROGNAME= ThermalLC
CFLAGS= -O3
SRCS = ThermalLC.cpp ../InertiaTimesCraterLagerros.cpp ../asteroid.cpp ../HybridModel.cpp ../ThermalInertiaOnlyConvex.cpp ../HemisphericCrater.cpp ../TriangulatedConvex.cpp ../ConvexFile.cpp ../vector3.cpp ParmfileLC.cpp ../ThermalInertia.cpp ../HemisphericNoInertia.cpp 
OBJS = ThermalLC.o ../InertiaTimesCraterLagerros.o ../asteroid.o ../HybridModel.o ../ThermalInertiaOnlyConvex.o ../HemisphericCrater.o ../TriangulatedConvex.o ../ConvexFile.o ../vector3.o ParmfileLC.o ../ThermalInertia.o ../HemisphericNoInertia.o 

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

all: $(SRCS) $(PROGNAME)

$(PROGNAME) : $(OBJS)
	$(CC) $(CFLAGS) -o $(PROGNAME) $(OBJS)

clean:
	$(rm) $(OBJS) $(PROGNAME) core *~
