ARPACK =  -L/usr/local/lib/ -larpack

LoopOverCalcBound.x:	LoopOverCalcBound.f Bsplines.f matrix_stuff.f modules_qd.o besselnew.o LoopOverCalcBound.o Bsplines.o matrix_stuff.o bspline90_22.o
	gfortran LoopOverCalcBound.o matrix_stuff.o Bsplines.o modules_qd.o bspline90_22.o besselnew.o -framework Accelerate ${ARPACK} -C -o LoopOverCalcBound.x
#	gfortran LoopOverCalcBound.o matrix_stuff.o Bsplines.o 1Dpot.o modules_qd.o bspline90_22.o besselnew.o  -framework vecLib -larpack_OSX -lg2c -C -o LoopOverCalcBound.x

#1Dpot.o:	1Dpot.f
#	gfortran -ffixed-line-length-132 -c -C 1Dpot.f

LoopOverCalcBound.o:	LoopOverCalcBound.f
	gfortran  -ffixed-line-length-132 -c -C LoopOverCalcBound.f

matrix_stuff.o:	matrix_stuff.f
	gfortran -ffixed-line-length-132 -c -C matrix_stuff.f

Bsplines.o:	Bsplines.f
	gfortran -ffixed-line-length-132 -c -C Bsplines.f	

besselnew.o :	besselnew.f
	gfortran -ffixed-line-length-132 -c -C besselnew.f	

modules_qd.o :	modules_qd.f90
	gfortran -ffixed-line-length-132 -c -C modules_qd.f90

bspline90_22.o: bspline90_22.f90
	gfortran -c bspline90_22.f90
