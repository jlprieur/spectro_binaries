
##################################################################
# Makefile for c programs 
#
# JLP
# Version 05-01-93
##################################################################
include $(JLPSRC)/jlp_make.mk

.SUFFIXES:
.SUFFIXES: .o .c .exe $(SUFFIXES) 

.c.o:
	cc -c $(CFLAGS) $*.c

.o.exe:
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

.c.exe:
	cc -c $(CFLAGS) $*.c
	cc -o $(EXEC)/$*.exe $*.o \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

all: simu_mass.exe

absp1.exe:	absp1.c cel_meca1.c
	cc -c $(CFLAGS) absp1.c
	cc -o $(EXEC)/absp1.exe absp1.o \
	$(JLIB) $(MIDLIB) $(F77LIB) -lm

simu_mass.exe: simu_mass_set.c
#	cc -I$JLPSRC/jlp_incl simu_mass_set.c $JLPLIB/jlp/jlputil.a -lm
	cc -c $(CFLAGS) simu_mass_set.c 
	cc -g -o simu_mass.exe simu_mass_set.o $(JLIB) -lm

simu_mass_set.o : simu_mass_set.c
	cc -g -o $(EXEC)/$*.exe $*.o ../dcv/dcv_cgrad.o $(JLPLIB)/math/lbfgs.a \
	$(JLIB) $(MIDLIB) $(MATHLIB) $(F77LIB) $(XLIB) -lm

clear:
	rm -f *.o
