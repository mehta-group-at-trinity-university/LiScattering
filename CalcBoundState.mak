LAPACK =   -framework accelerate
ARPACK =  -L/usr/local/lib/ -larpack


CalcBoundState.x:	Bsplines.f matrix_stuff.f CalcBoundState.o Bsplines.o matrix_stuff.o nrmat.o TwoBodydata.o CalcBoundState.f90 
	gfortran      CalcBoundState.o matrix_stuff.o Bsplines.o nrmat.o TwoBodydata.o ${LAPACK}  ${ARPACK} -o CalcBoundState.x

CalcBoundState.o:	CalcBoundState.f90
	gfortran    -ffixed-line-length-132  -c CalcBoundState.f90

matrix_stuff.o:	matrix_stuff.f
	gfortran    -ffixed-line-length-132  -c matrix_stuff.f

Bsplines.o:	Bsplines.f
	gfortran    -ffixed-line-length-132  -c Bsplines.f	

nrmat.o:	nrmat.f
	gfortran    -ffixed-line-length-132  -c nrmat.f

TwoBodydata.o:   TwoBodydata.f90
	gfortran    -ffixed-line-length-132  -c TwoBodydata.f90 



       
