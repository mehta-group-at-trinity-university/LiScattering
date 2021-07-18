CMP     = ifort
F132FORM = -extend-source #-ffixed-line-length-132
OPTFLAG = -O4  -i8  
FREEFORM = #-ffree-form
STND = #-fdec #-std=gnu
LEGACY = -std=legacy
DEBUG   = #-fcheck=all
FORCEDP = #-fdefault-real-8 -fdefault-double-
LAPACK =   ${MKLROOT}/lib/libmkl_blas95_ilp64.a ${MKLROOT}/lib/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl  -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
ARPACK =  -L/usr/local/lib/ -larpack
INCLUDE =   -I/usr/local/include  -I${MKLROOT}/include/intel64/ilp64 -i8  -I"${MKLROOT}/include"
OBJS  = bspline90_22.o besselnew.o threejsixj.o Bsplines.o quadrature.o POTGENLI2.o AlkaliScattering.o units.o matrix_stuff.o 

AlkaliScattering.x: ${OBJS} 
	${CMP} ${STND} ${DEBUG} ${OBJS} ${INCLUDE}  ${OPTFLAG} ${FORCEDP} -o AlkaliScattering.x ${LAPACK}  ${ARPACK}

AlkaliScattering.o: AlkaliScattering.f90 units.o bspline90_22.o
	${CMP} ${STND} ${DEBUG} ${FORCEDP} ${FREEFORM} ${OPTFLAG} -c AlkaliScattering.f90

matrix_stuff.o: matrix_stuff.f 
	${CMP} ${STND} ${FORCEDP} ${F132FORM} ${INCLUDE}  ${OPTFLAG} -c matrix_stuff.f ${LAPACK} ${ARPACK}

Bsplines.o:	Bsplines.f
	${CMP} ${STND} ${FORCEDP} ${F132FORM} ${OPTFLAG} -c Bsplines.f

nrtype.mod: modules_qd.o
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} modules_qd.o

modules_qd.o:	modules_qd.f90
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} -c modules_qd.f90

bspline90_22.o: bspline90_22.f90
	${CMP} ${STND} ${FORCEDP} ${OPTFLAG} -c bspline90_22.f90

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
