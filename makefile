CMP     = gfortran
F132FORM = -ffixed-line-length-132
OPTFLAG = -O4 -floop-nest-optimize
FREEFORM = -ffree-form
STND = #-fdec #-std=gnu
LEGACY = -std=legacy
DEBUG   = -fcheck=bounds -fbounds-check # -Wall
FORCEDP = #-fdefault-real-8 -fdefault-double-
LAPACK =   -framework accelerate
PROMOTEQUAD = #-freal-8-real-16
ARPACK =  -L/usr/local/lib/ -larpack
INCLUDE =  -I/usr/local/include
OBJS  = bspline90_22.o linfit.o besselnew.o threejsixj.o Bsplines.o quadrature.o POTGENLI2.o AlkaliScattering.o units.o matrix_stuff.o 

AlkaliScattering.x: ${OBJS} 
	${CMP} ${PROMOTEQUAD} ${STND} ${DEBUG} ${OBJS} ${INCLUDE} ${LAPACK}  ${ARPACK} ${OPTFLAG} ${FORCEDP} -o AlkaliScattering.x

AlkaliScattering.o: AlkaliScattering.f90 units.o bspline90_22.o
	${CMP} ${PROMOTEQUAD}${STND} ${DEBUG} ${FORCEDP} ${FREEFORM} ${OPTFLAG} -c AlkaliScattering.f90

matrix_stuff.o: matrix_stuff.f 
	${CMP} ${LEGACY} ${PROMOTEQUAD}${STND} ${FORCEDP} ${F132FORM} ${INCLUDE} ${LAPACK} ${ARPACK} ${OPTFLAG} -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	${CMP} ${PROMOTEQUAD}${STND} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c Bsplines.f

nrtype.mod: modules_qd.o
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} -c modules_qd.f90

bspline90_22.o: bspline90_22.f90
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} -c bspline90_22.f90

besselnew.o:	besselnew.f
	${CMP} ${PROMOTEQUAD} ${STND} ${DEBUG} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c besselnew.f

quadrature.mod: quadrature.o
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} quadrature.o

quadrature.o: quadrature.f90
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} -c quadrature.f90

POTGENLI2.o: POTGENLI2.f
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} -c POTGENLI2.f

threejsixj.o: threejsixj.f
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} -c threejsixj.f

units.o: units.f90
	${CMP} ${PROMOTEQUAD} ${STND} ${FORCEDP} ${OPTFLAG} -c units.f90

linfit.o: linfit.f
	${CMP} ${PROMOTEQUAD} ${LEGACY} ${STND} ${FORCEDP} ${OPTFLAG} -c linfit.f

clean:
	rm -f *.mod *.o *.x
