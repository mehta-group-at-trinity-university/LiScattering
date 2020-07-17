! THIS PROGRAM CALCULATES THE BOUND STATES FOR A SINGLE-CHANNEL TWO-BODY POTENTIAL.  
! RIGHT NOW, IT'S SET UP FOR THE SECH**2 POTENTIAL.
program Calc2BodyBoundState
  use TwoBodydata
  implicit none
  double precision, external :: VSech, phirecon,Vsqrc6
  double precision, allocatable :: Energies(:), xgrid(:)
  integer xDimMin,iEnergy,NumEnergySteps,iE,ikeep,n,numgrid
  double precision mass,x,lwave,norm,dx,ediff,eexact
  double precision, allocatable :: evec(:,:),eval(:)
  character*64 LegendreFile
  Format1="(F18.15,F18.15,F18.15,F18.15)"

  ! read in legendre data
  read(5,*)
  read(5,1002) LegendreFile
  write(6,1002) LegendreFile
  read(5,*)
  read(5,*)
  read(5,*) LegPoints
  write(6,*) 'Using Legendre integration of order: ',LegPoints
  write(6,*) ' '
  !     read in boundary conditions
  read(5,*)
  read(5,*)
  read(5,*) Order,Left,Right
  write(6,*) 'Order, Left, Right'
  write(6,'(I5,I5,I5)') Order,Left,Right
  write(6,*) ' '
  !     read in potential parameters
  read(5,*)
  read(5,*)
  read(5,*) alpha, mass, DD, x0
  write(6,*) 'alpha,            mass,              DD,               x0'
!  DD = 10.d0/x0**6
  write(6,'(F15.10,F15.10,F15.10,F15.10)') alpha, mass, DD, x0
  write(6,*)

  mu=mass/2.0d0
  write(6,*) 'mu = ', mu
  write(6,*)
  read(5,*)
  read(5,*)
  read(5,*) xNumPoints,xMin,xMax
  write(6,*) 'xNumPoints,   xMin,      xMax'
  write(6,'(I8,F15.10,F15.10)') xNumPoints,xMin,xMax
  write(6,*)
  read(5,*)
  read(5,*)
  read(5,*) lwave
  write(6,*) 'lwave'
  write(6,*) lwave
  write(6,*)
  xDimMin=xNumPoints+Order-3
  xDim=xDimMin
  if (Left .eq. 2) xDim = xDim + 1
  if (Right .eq. 2) xDim = xDim + 1
  MatrixDim=xDim
  print*, 'xDim=',xDim
  allocate(xLeg(LegPoints),wLeg(LegPoints))
  write(6,*) 'Getting the Legendre points and weights'
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  allocate(xPoints(xNumPoints))
  allocate(xBounds(xNumPoints+2*Order))
  allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)) ! 

  allocate(S(MatrixDim,MatrixDim))
  allocate(H0(MatrixDim,MatrixDim))

  allocate(evec(MatrixDim,MatrixDim))
  allocate(eval(MatrixDim))

  write(6,*) 'Writing the potential plot to file fort.1000'
  numgrid=1000
  allocate(xgrid(numgrid))
  dx=(xMax-xMin)/numgrid
  do n=1,numgrid
     xgrid(n)=n*dx
     write(1000,*) xgrid(n), VSech(DD,x0,xgrid(n),lwave)
  end do

  write(6,*) 'Calling GridMaker.'
  call GridMaker(xNumPoints/2,xMin,x0,"linear",xPoints(1:xNumPoints/2))
  dx=xPoints(xNumPoints/2)-xPoints(xNumPoints/2-1)
  call GridMaker(xNumPoints/2,x0+dx,xMax,"quadratic",xPoints(xNumPoints/2+1:xNumPoints))
  call printmatrix(xPoints,xNumPoints,1,6)
!  call GridMaker(xNumPoints,xMin,xMax,"linear",xPoints)
  write(6,*) 'Calculating the basis functions.'
  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,0,u)

  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,2,uxx)
!  write (6,*) 'Calling CheckBasis...'
!  call CheckBasis(u,xDim,xNumPoints,LegPoints,xLeg,xPoints,700)

  H0=0.0d0
  S=0.0d0

  write(6,*) 'Calculating Hamiltonian'
  call CalcHamiltonian(lwave)
  call Mydggev(MatrixDim,H0,MatrixDim,S,MatrixDim,eval,evec)
  eval=eval*(-1.d0)
  call eigsrt(eval,evec,MatrixDim,MatrixDim)
  eval=eval*(-1.d0)
  do n=1,MatrixDim
     if(eval(n).lt.0d0) then
        eexact=-1d0/(8*mu*x0**2) * (-(1d0+2*(2*n-1)) + Sqrt(1d0 + 8*mu*DD*x0**2))**2
        write(6,*) 'Energy Eigenvalue (',n,') = ', eval(n), &
           eexact, (eval(n)-eexact)/abs(eexact)
     endif
  end do
  
  ikeep=28

  x=0.0d0
  xMax=5.d0
  write(6,*) 'Writing the wavefunction for the ground state to fort.777 ',ikeep
  do while(x.lt.xMax)
     write(777,10) x, -phirecon(x,ikeep,evec,Left,Right,xDim,MatrixDim,xNumPoints,xPoints,order)
     x=x+0.001d0
  enddo
  write(777,*)

  deallocate(H0,S)
  deallocate(xBounds,xPoints)
  deallocate(u,uxx,xLeg,wLeg)

10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)

end program Calc2BodyBoundState
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcHamiltonian(lwave)
  use TwoBodydata
  implicit none
  double precision, external :: VSech,Vsqrc6
  double precision ax, bx,x,xScaledZero,xIntScale,TempT,TempS,TempV,a,lwave,mcoeff
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  integer kx,ix,ixp,lx,n

  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))

  !c     seting up kxMin and kxMax
  do ix = 1,xDim
     do ixp = 1,xDim
        kxMin(ix,ixp) = max(xBounds(ix),xBounds(ixp))
        kxMax(ix,ixp) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1 ! 
     enddo
  enddo

  do ix = 1,xDim
     do ixp = max(1,ix-Order),min(xDim,ix+Order)
        H0(ix,ixp)=0.0d0
        S(ix,ixp)=0.0d0
!        L(ix,ixp)=0.0d0
!        V0(ix,ixp)=0.0d0
        do kx = kxMin(ix,ixp),kxMax(ix,ixp)
           ax = xPoints(kx)
           bx = xPoints(kx+1)
           xIntScale=0.5d0*(bx-ax)
           xScaledZero=0.5d0*(ax+bx)
           TempT = 0.0d0
           TempS = 0.0d0
           TempV = 0.0d0
           do lx = 1,LegPoints
              a = wLeg(lx)*xIntScale
              !x=xIntPoints(lx,kx)!
              x=xIntScale*xLeg(lx)+xScaledZero
              TempS = TempS + a*u(lx,kx,ix)*u(lx,kx,ixp)
              !TempV = TempV + a*u(lx,kx,ix)*(alpha*Vsqrc6(DD,x0,x,lwave))*u(lx,kx,ixp)
              TempV = TempV + a*u(lx,kx,ix)*(alpha*VSech(DD,x0,x,lwave))*u(lx,kx,ixp)
              !iRall = (kx-1)*LegPoints + lx
              !TempV = TempV + a*u(lx,kx,ix)*(alpha*VR(iRall) + lwave*(lwave+1d0)/(2d0*mu*xIntPoints(lx,kx)**2))*u(lx,kx,ixp)
              TempT = TempT + a*0.5d0/mu*(-u(lx,kx,ix)*uxx(lx,kx,ixp))
!              print*,a, TempS, TempT, TempV
           enddo
           S(ix,ixp) = S(ix,ixp) + TempS ! place values into overlap matrix
           H0(ix,ixp) = H0(ix,ixp) + TempV + TempT ! place values into Hamiltonian
!           write(6,*) ix, ixp, H0(ix,ixp), S(ix,ixp)
        enddo
     enddo
  enddo

  deallocate(kxMax,kxMin)
10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)
  
end subroutine CalcHamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GridMaker(numpts,E1,E2,scale,grid)
  implicit none
  DOUBLE PRECISION grid(numpts)
  DOUBLE PRECISION E1,E2,LE1,LE2,DE,LDE
  INTEGER numpts, iE
  CHARACTER(LEN=*), INTENT(IN) :: scale
  !--------------------------------------------
  ! Linear grid:
  !--------------------------------------------
  grid(1)=E1
  IF((scale.EQ."linear").and.(numpts.gt.1)) THEN
     DE=(E2-E1)/DBLE(numpts-1)
     DO iE=1,numpts
        grid(iE) = E1 + (iE-1)*DE
     ENDDO
  ENDIF
  !--------------------------------------------
  ! Log grid:
  !--------------------------------------------
  IF((scale.EQ."log").and.(numpts.gt.1)) THEN
     LE1=dlog(E1)
     LE2=dlog(E2)

     LDE=(LE2-LE1)/DBLE(numpts-1d0)
     DO iE=1,numpts
        grid(iE) = dexp(LE1 + (iE-1)*LDE)
!        write(6,*) LE1, LE2, LDE, grid(iE)
     ENDDO
  ENDIF
  !--------------------------------------------
  ! Quadratic grid:
  !--------------------------------------------
  IF((scale.EQ."quadratic").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**2*DE
     ENDDO
  ENDIF
END SUBROUTINE GridMaker
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerOld(xNumPoints,xMin,xMax,xPoints)
  implicit none
  integer xNumPoints
  double precision xMin,xMax,xPoints(xNumPoints)

  integer i,j,k
  double precision Pi,xmid
  double precision r0New
  double precision xRswitch
  double precision xDelt,x0,x1,x2

  Pi = 3.1415926535897932385d0

  x0 = xMin
  x1 = xMax
  k = 1
  xDelt = (x1-x0)/dfloat(xNumPoints-1)
  do i = 1,xNumPoints
     xPoints(k) = (i-1)*xDelt + x0
     k = k + 1
  enddo

15 format(6(1x,1pd12.5))


  return
end subroutine GridMakerOld

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine Mydggev(N,G,LDG,L,LDL,eval,evec)
!!$  implicit none
!!$  integer LDG,N,LDL,info
!!$  double precision G(LDG,N),L(LDL,N),eval(N),evec(N,N),b
!!$  double precision, allocatable :: alphar(:),alphai(:),beta(:),work(:),VL(:,:) ! 
!!$  integer lwork,i,im,in
!!$
!!$  allocate(alphar(N),alphai(N),beta(N))
!!$
!!$
!!$  info = 0
!!$  lwork = -1
!!$  allocate(work(1))
!!$  call dggev('N','V',N,G,LDG,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info) ! 
!!$  do im = 1,N
!!$     alphar(im)=0.0d0
!!$     alphai(im)=0.0d0
!!$     beta(im)=0.0d0
!!$     eval(im)=0.0d0
!!$     do in = 1,N
!!$        evec(im,in)=0.0d0
!!$     enddo
!!$  enddo
!!$
!!$  lwork=work(1)
!!$  deallocate(work)
!!$  allocate(work(lwork))
!!$  call dggev('N','V',N,G,LDG,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info)
!!$
!!$  do i = 1, N
!!$     if (abs(alphai(i)).ge.1d-14) then
!!$        print*, '#eigenvalue may be complex! alphai(',i,')=',alphai(i)
!!$     endif
!!$     if(abs(beta(i)).ge.1d-14) then
!!$        eval(i) = alphar(i)/beta(i)
!!$     endif
!!$     !c     if(abs(alphar(i)).ge.1d-14) then
!!$     !c     print*, alphar(i), alphai(i), beta(i)
!!$     !c     endif
!!$  enddo
!!$  deallocate(alphar,alphai,beta)
!!$  deallocate(work)
!!$
!!$end subroutine Mydggev
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function phirecon(R,beta,evec,left,right,RDim,MatrixDim,RNumPoints,RPoints,order)
  implicit none
  double precision, external :: BasisPhi
  integer MatrixDim,RDim,nch,beta,i,RNumPoints,left,right,order
  double precision R,evec(MatrixDim,MatrixDim),RPoints(RNumPoints)
  phirecon = 0.0d0
  do i = 1,RDim
     phirecon = phirecon + evec(i,beta)*BasisPhi(R,left,right,order,RDim,RPoints, &
          RNumPoints,0,i)
  enddo
  return
end function phirecon
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CheckBasis(u,RDim,RNumPoints,LegPoints,xLeg,RPoints,file)
  implicit none
  integer RNumPoints,NumChan,LegPoints,kx,lx,ix,RDim,file
  double precision u(LegPoints,RNumPoints,RDim),ax,bx,xLeg(LegPoints)
  double precision RPoints(RNumPoints),xScaledZero,xIntScale, x
  
  do ix=1,RDim
     do kx = 1,RNumPoints-1
        ax = RPoints(kx)
        bx = RPoints(kx+1)
        xIntScale = 0.5d0*(bx-ax)
        xScaledZero = 0.5d0*(bx+ax)
        do lx = 1,LegPoints
           x = xIntScale*xLeg(lx)+xScaledZero
           write(file,*) x, u(lx,kx,ix)
        enddo
     enddo
     write(file,*) ' '
  enddo
  
end subroutine CheckBasis
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printmatrix(M,nr,nc,file)
  implicit none
  integer nr,nc,file,j,k
  double precision M(nr,nc)

  do j = 1,nr
     write(file,20) (M(j,k), k = 1,nc)
  enddo
20 format(1P,100e16.8)
end subroutine printmatrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function VSech(DD,x0,x,lwave)
  implicit none
  double precision x,x0,DD,lwave
  VSech = -DD/dcosh(x/x0)**2.d0  + lwave*(lwave+1.d0)/(x*x)
  return
end function VSech
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function Vsqrc6(V0,x0,x,lwave)
  ! This returns the van der waals potential with a short-range square well
  ! It is in van der waals units where x is measured in units of (2*mu*C6/hbar**2)**0.25

  implicit none
  double precision x, x0, V0,lwave
  if (x.gt.x0) then
     VsqrC6 = -0.994596*x**(-6) - 0.00537041*x**(-8) - 0.000033133*x**(-10) + lwave*(lwave+1.d0)/(x*x)
  else
     vsqrC6 = -V0 + lwave*(lwave+1.d0)/(x*x)
  endif
end function Vsqrc6
