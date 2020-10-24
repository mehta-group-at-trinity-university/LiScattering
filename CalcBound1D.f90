! This program is a driver for Calc2BodyBoundState, which calculates the bound states for a single channel (1D) potential.
program BoundStateDriver
  implicit none
    character*64 LegendreFile
    integer LegPoints, Order, Left, Right, xNumPoints, xMin, xMax, NR, Ubcoef, kR, NumStates,iLam, PsiFlag
    double precision kleft, kright, alpha, mu, lwave
    double precision, allocatable :: Rknot(:), Energies(:)

    LegendreFile = 'Legendre.dat                          '
    LegPoints = 10
    
    Order = 5  ! order of the b-spline basis for the wave function
    Left = 0
    Right = 0
    alpha = 1.d0

    kR = 3  ! this is the order of the spline interpolation of the potential energy function
    NR = 100  ! Number of interpolated points for the potential energy function
    
    mu = 0.5d0
    NumStates = 20
    xNumPoints = 100
    xMin = 0.d0
    xMax = 20.d0
    
    allocate(Energies(NumStates),Rknot(kR+NR))

    PsiFlag = 0
    iLam = 0
    
    call Calc2BodyBoundState(LegendreFile, LegPoints, Order, Left, Right, kLeft, kRight, alpha, mu, xNumPoints, &
         xMin, xMax, lwave, NR, Ubcoef, kR, Rknot, NumStates, Energies, iLam, PsiFlag)

    
end program BoundStateDriver


subroutine Calc2BodyBoundState(LegendreFile,LegPoints,Order,Left,Right,kLeft,kRight,alpha,&
     mu,xNumPoints,xMin,xMax,lwave,NR,Ubcoef,kR,Rknot,NumStates,Energies,iLam,PsiFlag)
  use bspline
  implicit none
  integer MatrixDim,LegPoints,xDim,xNumPoints,Order,Left,Right,NumStates
  double precision, allocatable :: xLeg(:),wLeg(:)
  double precision mu,xMin,xMax,lam,r2b,DD,x0,alpha
  double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:),xPoints(:)
  double precision, allocatable :: S(:,:),L(:,:),H0(:,:),G0(:,:),V0(:,:),delta(:) 
  integer, allocatable :: xBounds(:)   
  character(LEN=30) :: Format1

  integer NR,kR,iLam,PsiFlag,i,NSteps
  double precision Ubcoef(NR),Rknot(kR+nR),Energies(NumStates)
  double precision, external :: VSech, phirecon
  double precision, allocatable :: xgrid(:)
  integer xDimMin,iEnergy,NumEnergySteps,iE,ikeep,n,numgrid
  double precision mass,x,lwave,norm,dx
  double precision, allocatable :: evec(:,:),eval(:)
  character*64 LegendreFile, iLamstr
  double precision kLeft,kRight
  double precision, external :: BasisPhi
  Format1="(F18.15,F18.15,F18.15,F18.15)"


  xDimMin=xNumPoints+Order-3
  xDim=xDimMin
  if (Left .eq. 2) xDim = xDim + 1
  if (Right .eq. 2) xDim = xDim + 1
  MatrixDim=xDim
!  print*, 'xDim=',xDim
  allocate(xLeg(LegPoints),wLeg(LegPoints))
!  write(6,*) 'Getting the Legendre points and weights'
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  allocate(xPoints(xNumPoints))
  allocate(xBounds(xNumPoints+2*Order))
  allocate(u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)) ! 

  allocate(S(MatrixDim,MatrixDim))
  allocate(H0(MatrixDim,MatrixDim))

  allocate(evec(MatrixDim,MatrixDim))
  allocate(eval(MatrixDim))

!  write(6,*) 'Writing the potential plot to file fort.1000'
!  numgrid=100
!  allocate(xgrid(numgrid))
!  dx=(xMax-xMin)/dble(numgrid)
!  do n=1,numgrid
!     xgrid(n)=n*dx
!!     write(1000,*) xgrid(n), VSech(DD,x0,xgrid(n)) + lwave*(lwave+1.d0)/(xgrid(n)**2.d0)
!     write(100000,*) xgrid(n), dbsval(xgrid(n),kR,Rknot,NR,Ubcoef)
!  end do

!  write(6,*) 'Calling GridMakerLinear.'
!  call GridMakerLinear(xNumPoints,xMin,xMax,xPoints)
  call GridMakerQuadratic(xNumPoints,xMin,xMax,xPoints)
!  write(6,*) 'Calculating the basis functions.'
!  write(6,*) left,right,aleft,aright
  call CalcBasisFuncs(Left,Right,kLeft,kRight,Order,xPoints,LegPoints,xLeg,&
       xDim,xBounds,xNumPoints,0,u)

  call CalcBasisFuncs(Left,Right,kLeft,kRight,Order,xPoints,LegPoints,xLeg,&
       xDim,xBounds,xNumPoints,1,ux)

  !print*,"uprime1/u=",ux(1,1,1)/u(1,1,1),'-1/a=',-1.d0/kLeft
  !print*,"uprime2/u=",BasisPhi(0.d0,left,right,kLeft,kRight,order,xDim,xPoints,xNumPoints,1,1)/&
  !     BasisPhi(0.d0,left,right,kLeft,kRight,order,xDim,xPoints,xNumPoints,1,1)

  call CalcBasisFuncs(Left,Right,kLeft,kRight,Order,xPoints,LegPoints,xLeg,&
       xDim,xBounds,xNumPoints,2,uxx)
!  write (6,*) 'Calling CheckBasis...'
  !call CheckBasis(u,xDim,xNumPoints,LegPoints,xLeg,xPoints,700)

  H0=0.0d0
  S=0.0d0

!  write(6,*) 'Calculating Hamiltonian'
  call CalcHamiltonian1D(mu,lwave,alpha,NR,Ubcoef,kR,Rknot,LegPoints,xDim,MatrixDim,H0,S,&
       xBounds,xNumPoints,Order,xPoints,wLeg,xLeg,u,ux,uxx,Left,Right,kLeft,kRight)
  call Mydggev(MatrixDim,H0,MatrixDim,S,MatrixDim,eval,evec)
  eval=eval*(-1.d0)
  call eigsrt(eval,evec,MatrixDim,MatrixDim)
  eval=eval*(-1.d0)
  do n=1,NumStates
     !write(6,*) 'Energy Eigenvalue (',n,') = ', eval(n)
     Energies(n) = eval(n)
  end do
  ikeep=11

  do ikeep = 1,1
     if (PsiFlag .ne. 0) then
        write(iLamstr,'(I1)') iLam
        open(unit=3,file='./psi'//trim(iLamstr)//'.dat')
        x=0.d0
        NSteps = 1000
        do i = 1,NSteps
           x= xMin + (dble(i) - 1.d0)/(dble(NSteps-1))*(xMax - xMin)
           write(3,20) x, phirecon(x,1,evec,Left,Right,kLeft,kRight,xDim,MatrixDim,xNumPoints,xPoints,Order,0)
        enddo
        close(unit=3)
     endif
  enddo

  deallocate(H0,S)
  deallocate(xBounds,xPoints)
  deallocate(u,ux,uxx,xLeg,wLeg)

10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)

end subroutine Calc2BodyBoundState
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcHamiltonian1D(mu,lwave,alpha,NR,Ubcoef,kR,Rknot,LegPoints,xDim,MatrixDim,&
     H0,S,xBounds,xNumPoints,Order,xPoints,wLeg,xLeg,u,ux,uxx,Left,Right,bcLeft,bcRight)
  !  use TwoBodydata
  use bspline
  implicit none
  integer NR,kR,LegPoints,xDim,MatrixDim,Order,xNumPoints,Left,Right
  integer xBounds(xNumPoints+2*Order)
  double precision Ubcoef(NR),Rknot(kR+nR), H0(MatrixDim,MatrixDim),S(MatrixDim,MatrixDim),xPoints(xNumPoints)
  double precision wLeg(LegPoints),xLeg(LegPoints)
  double precision u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim),alpha,mu
  double precision, external :: VSech
  double precision ax, bx,x,xScaledZero,xIntScale,TempT,TempS,TempV,a,lwave,mcoeff,vpot,DD,x0
  double precision bcLeft,bcRight
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

  DD = 2.d0
  x0 = 1.d0

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
              x=xIntScale*xLeg(lx)+xScaledZero
              TempS = TempS + a*u(lx,kx,ix)*u(lx,kx,ixp)
              vpot = dbsval(x,kR,Rknot,NR,Ubcoef)
              TempV = TempV + a*u(lx,kx,ix)*(alpha*VSech(DD,x0,x) + lwave*(lwave+1.d0)/(x*x))*u(lx,kx,ixp)
!              TempV = TempV + a*u(lx,kx,ix)*(alpha*vpot + lwave*(lwave+1.d0)/(x*x))*u(lx,kx,ixp)
              TempT = TempT + a*0.5d0/mu*(ux(lx,kx,ix)*ux(lx,kx,ixp))
!              print*,a, TempS, TempT, TempV
           enddo
           S(ix,ixp) = S(ix,ixp) + TempS ! place values into overlap matrix
           H0(ix,ixp) = H0(ix,ixp) + TempV + TempT ! place values into Hamiltonian
           !print*, TempV, TempT, lwave,vpot
           !write(6,*) ix, ixp, H0(ix,ixp), S(ix,ixp)
        enddo
        if ( (Left.eq.3).and.(ix.eq.1).and.(ixp.eq.1) ) then
           H0(ix,ixp) = H0(ix,ixp) + 0.5d0/mu*bcLeft
        endif
        if ( (Right.eq.3).and.(ix.eq.xDim).and.(ixp.eq.xDim) ) then
           H0(ix,ixp) = H0(ix,ixp) - 0.5d0/mu*bcRight
        endif
     enddo
  enddo

  do ix = 1,xDim
     do ixp = 1,xDim
        if(dabs(H0(ix,ixp) - H0(ixp,ix))/abs(H0(ix,ixp) + H0(ixp,ix)).gt.1d-8) then
             print*,"Hamiltonian not symmetric"
             print*,'H(ix,ixp) = ',H0(ix,ixp)
             print*,'H(ixp,ix) = ',H0(ixp,ix)
         endif
      enddo
   enddo

  deallocate(kxMax,kxMin)
10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)
  
end subroutine CalcHamiltonian1D
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerLinear(xNumPoints,xMin,xMax,xPoints)
  implicit none
  integer xNumPoints
  double precision xMin,xMax,xPoints(xNumPoints)

  integer i,j,k
  double precision Pi
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
end subroutine GridMakerLinear
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GridMakerQuadratic(xNumPoints,xMin,xMax,xPoints)
  implicit none
  integer xNumPoints
  double precision xMin,xMax,xPoints(xNumPoints)

  integer i,j,k
  double precision Pi
  double precision r0New
  double precision xRswitch
  double precision xDelt,x0,x1,x2

  Pi = 3.1415926535897932385d0

  x0 = xMin
  x1 = xMax
  k = 1
  xDelt = (x1-x0)/dfloat(xNumPoints-1)
  do i = 1,xNumPoints
     xPoints(k) = ((i-1)/dble(xNumPoints-1))**2*x1 + x0
     k = k + 1
  enddo

15 format(6(1x,1pd12.5))


  return
end subroutine GridMakerQuadratic

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Mydggev(N,G,LDG,L,LDL,eval,evec)
  implicit none
  integer LDG,N,LDL,info
  double precision G(LDG,N),L(LDL,N),eval(N),evec(N,N),b
  double precision, allocatable :: alphar(:),alphai(:),beta(:),work(:),VL(:,:) ! 
  integer lwork,i,im,in

  allocate(alphar(N),alphai(N),beta(N))


  info = 0
  lwork = -1
  allocate(work(1))
  call dggev('N','V',N,G,LDG,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info) ! 
  do im = 1,N
     alphar(im)=0.0d0
     alphai(im)=0.0d0
     beta(im)=0.0d0
     eval(im)=0.0d0
     do in = 1,N
        evec(im,in)=0.0d0
     enddo
  enddo

  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
  call dggev('N','V',N,G,LDG,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info)

  do i = 1, N
     if (abs(alphai(i)).ge.1d-14) then
        print*, '#eigenvalue may be complex! alphai(',i,')=',alphai(i)
     endif
     if(abs(beta(i)).ge.1d-14) then
        eval(i) = alphar(i)/beta(i)
     endif
     !c     if(abs(alphar(i)).ge.1d-14) then
     !c     print*, alphar(i), alphai(i), beta(i)
     !c     endif
  enddo
  deallocate(alphar,alphai,beta)
  deallocate(work)

end subroutine Mydggev
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
               write(file,*), x, u(lx,kx,ix)
            enddo
         enddo
         write(file,*), ' '
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
double precision function VSech(DD,x0,x)
  implicit none
  double precision x,x0,DD
  VSech = -DD/dcosh(x/x0)**2.d0
  return
end function VSech
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE eigsrt(d,v,n,np)
  INTEGER n,np
  DOUBLE PRECISION d(np),v(np,np)
  INTEGER i,j,k
  DOUBLE PRECISION p
  do i=1,n-1
     k=i
     p=d(i)
     do j=i+1,n
        if(d(j).ge.p)then
           k=j
           p=d(j)
        endif
     enddo
     if(k.ne.i)then
        d(k)=d(i)
        d(i)=p
        do j=1,n
           p=v(j,i)
           v(j,i)=v(j,k)
           v(j,k)=p
        enddo
     endif
  enddo
  return
END SUBROUTINE eigsrt

        


