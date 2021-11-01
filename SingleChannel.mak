CMP     = gfortran
F132FORM = -ffixed-line-length-132
OPTFLAG = -O3
FREEFORM = -ffree-form
STND = #-std=gnu
LEGACY = -std=legacy
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
LAPACK =   -framework accelerate
ARPACK =  -L/usr/local/lib/ -larpack
INCLUDE =  -I/usr/local/include
OBJ = besselnew.o POTGENLI2.o units.o Bsplines.o matrix_stuff.o linfit.o EREdata.o minpack.o modules_qd.o bspline90_22.o SingleChannel.o

SingleChannel.x: ${OBJ}
	${CMP}  ${OBJ}  ${LAPACK} ${ARPACK} ${INCLUDE} -o SingleChannel.x

SingleChannel.o: bspline90_22.o EREdata.o units.o modules_qd.o SingleChannel.f90
	${CMP}  ${F132FORM}  -c SingleChannel.f90

matrix_stuff.o:	matrix_stuff.f
	${CMP} ${F132FORM} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${F132FORM} -c Bsplines.f	

linfit.o:	linfit.f
	${CMP}  ${LEGACY}  ${F132FORM}  -c linfit.f

EREdata.o: EREdata.f90
	${CMP} ${F132FORM}  -c EREdata.f90 

minpack.o: minpack.f90
	${CMP} -c minpack.f90

units.o: units.f90
	${CMP} -c units.f90

modules_qd.o:	modules_qd.f90
	${CMP}  ${F132FORM} -c modules_qd.f90

POTGENLI2.o: POTGENLI2.f
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} -c POTGENLI2.f

bspline90_22.o: bspline90_22.f90
	${CMP}    ${F132FORM} -c bspline90_22.f90

besselnew.o:	besselnew.f
	${CMP} ${STND} ${DEBUG} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c besselnew.f

