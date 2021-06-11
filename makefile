CMP     = gfortran
F132FORM = -ffixed-line-length-132
OPTFLAG = -O3
FREEFORM = -ffree-form
<<<<<<< HEAD
STND = #-std=gnu
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
LAPACK =   -framework accelerate
ARPACK =  -L/usr/local/lib/ -larpack
INCLUDE =  -I/usr/local/include
OBJS  = besselnew.o threejsixj.o Bsplines.o quadrature.o POTGENLI2.o LiScattering.o units.o matrix_stuff.o 

LiScattering.x: ${OBJS} 
	${CMP} ${STND} ${DEBUG} ${OBJS} ${INCLUDE} ${LAPACK}  ${ARPACK} ${OPTFLAG} ${FORCEDP} -o LiScattering.x

LiScattering.o: LiScattering.f90 units.o
	${CMP} ${STND} ${DEBUG} ${FORCEDP} ${FREEFORM} ${OPTFLAG} -c LiScattering.f90
=======
STND = #-fdec #-std=gnu
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-
LAPACK =   -framework accelerate
ARPACK =  -L/usr/local/lib/ -larpack
INCLUDE =  -I/usr/local/include
OBJS  = besselnew.o threejsixj.o Bsplines.o quadrature.o POTGENLI2.o AlkaliScattering.o units.o matrix_stuff.o 

AlkaliScattering.x: ${OBJS} 
	${CMP} ${STND} ${DEBUG} ${OBJS} ${INCLUDE} ${LAPACK}  ${ARPACK} ${OPTFLAG} ${FORCEDP} -o AlkaliScattering.x

AlkaliScattering.o: AlkaliScattering.f90 units.o
	${CMP} ${STND} ${DEBUG} ${FORCEDP} ${FREEFORM} ${OPTFLAG} -c AlkaliScattering.f90
>>>>>>> MQDT

#RMATPROP2016.o: RMATPROP2016.f90 quadrature.o
#	${CMP} ${STND} ${DEBUG} ${FORCEDP} ${FREEFORM} ${OPTFLAG} -c RMATPROP2016.f90

matrix_stuff.o: matrix_stuff.f 
	${CMP} ${STND} ${FORCEDP} ${F132FORM} ${INCLUDE} ${LAPACK} ${ARPACK} ${OPTFLAG} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${STND} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c Bsplines.f

nrtype.mod: modules_qd.o
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG}-c modules_qd.f90

besselnew.o:	besselnew.f
	${CMP} ${STND} ${DEBUG} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c besselnew.f

quadrature.mod: quadrature.o
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} quadrature.o

quadrature.o: quadrature.f90
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} -c quadrature.f90

POTGENLI2.o: POTGENLI2.f
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} -c POTGENLI2.f

threejsixj.o: threejsixj.f
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} -c threejsixj.f

units.o: units.f90
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} -c units.f90

clean:
	rm -f *.mod *.o *.x
