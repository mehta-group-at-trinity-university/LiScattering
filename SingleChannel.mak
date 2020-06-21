CMP     = gfortran
F132FORM = -ffixed-line-length-132
OPTFLAG = -O3
FREEFORM = -ffree-form
STND = #-std=gnu
DEBUG   = -fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-8
LAPACK =   -framework accelerate
ARPACK =  -L/usr/local/lib/ -larpack
INCLUDE =  -I/usr/local/include

SingleChannel.x: SingleChannel.f90 Bsplines.f matrix_stuff.f SingleChannel.o units.o Bsplines.o matrix_stuff.o linfit.o EREdata.o minpack.o
	${CMP}   SingleChannel.o matrix_stuff.o Bsplines.o linfit.o minpack.o units.o EREdata.o ${LAPACK} ${ARPACK} ${INCLUDE} -o SingleChannel.x

SingleChannel.o: bspline90_22.o EREdata.o units.o modules_qd.o SingleChannel.f90
	${CMP} bspline90_22.o modules_qd.o EREdata.o units.o  ${F132FORM}  -c SingleChannel.f90

matrix_stuff.o:	matrix_stuff.f
	${CMP}    ${F132FORM}  -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP}    ${F132FORM}  -c Bsplines.f	

linfit.o:	linfit.f
	${CMP}    ${F132FORM}  -c linfit.f

EREdata.o:   EREdata.f90
	${CMP}    ${F132FORM}  -c EREdata.f90 

minpack.o: minpack.f90
	${CMP} -c minpack.f90

units.o: units.f90
	${CMP} -c units.f90

modules_qd.o:	modules_qd.f90
	${CMP}    ${F132FORM} -c modules_qd.f90	

bspline90_22.o: bspline90_22.f90
	${CMP}    ${F132FORM} -c bspline90_22.f90

