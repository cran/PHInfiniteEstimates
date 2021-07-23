FFLAGS =  -g -ffpe-trap=zero,invalid,overflow,underflow  -pedantic -Wall -fbounds-check 
a.out: network.o drive.o expand.o dobnds.o cntdup.o copytox.o
	gfortran  $(FFLAGS) network.o drive.o dobnds.o expand.o cntdup.o copytox.o
dobnds.o: dobnds.f90
	gfortran $(FFLAGS) -c dobnds.f90
network.o: network.f90
	gfortran $(FFLAGS) -c network.f90
drive.o: drive.f90
	gfortran  $(FFLAGS) -c drive.f90
expand.o: expand.f90
	gfortran  $(FFLAGS) -c expand.f90
cntdup.o: cntdup.f90
	gfortran  $(FFLAGS) -c cntdup.f90
copytox.o: copytox.f90
	gfortran  $(FFLAGS) -c copytox.f90
clean:
	/bin/rm *.o
