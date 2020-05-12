CMP     = gfortran
F132FORM = -ffixed-line-length-132
OPTFLAG = -O3
FREEFORM = -ffree-form
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
INCLUDE =  -I/opt/lapack/include
LAPACK =  -framework accelerate
ARPACK =  -I/usr/local/include -L /usr/local/lib/ -larpack
OBJS  = besselnew.o zgensub.o Bsplines.o matrix_stuff.o RMATPROP2016.o Quadrature.o POTGENLI2.o LiScattering.o

LiScattering.x: ${OBJS}
	${CMP} ${DEBUG} ${OBJS} ${INCLUDE} ${ARPACK} ${LAPACK} ${OPTFLAG} ${FORCEDP} -o LiScattering.x

LiScattering.o: LiScattering.f90
	${CMP} ${DEBUG} ${FORCEDP} ${FREEFORM} ${OPTFLAG} -c LiScattering.f90

RMATPROP2016.o: RMATPROP2016.f90
	${CMP} ${DEBUG} ${FORCEDP} ${FREEFORM} ${OPTFLAG} -c RMATPROP2016.f90

matrix_stuff.o: matrix_stuff.f
	${CMP} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c Bsplines.f

nrtype.mod: modules_qd.o
	${CMP} ${FORCEDP} ${OPTFLAG} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${FORCEDP} ${OPTFLAG}-c modules_qd.f90

besselnew.o:	besselnew.f
	${CMP} ${DEBUG} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c besselnew.f

Quadrature.mod: Quadrature.o
		${CMP} ${FORCEDP} ${OPTFLAG} Quadrature.o

Quadrature.o: Quadrature.f90
	${CMP} ${FORCEDP} ${OPTFLAG} -c Quadrature.f90

POTGENLI2.o: POTGENLI2.f
	${CMP} ${FORCEDP} ${OPTFLAG} -c POTGENLI2.f

zgensub.o: zgensub.f
		${CMP} ${FORCEDP} ${OPTFLAG} -c zgensub.f

clean:
	rm -f *.mod *.o *.x
