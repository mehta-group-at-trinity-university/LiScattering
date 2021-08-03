program main
  use EREdata
  use units
  implicit none
  double precision, external :: VSech, Vsqrc6,phirecon
  double precision, allocatable :: Energies(:)
  integer xDimMin,iEnergy,NumEnergySteps,iE,iR,ISTATE,iatom,satom,ESTATE,ix
!  integer xNumPts1,xNumPts2,i
  double precision Ahfatom,giatom,Cvals(4),VLIM
  double precision mass,x,EnergyMin,EnergyMax,DeltaEnergy,deltatemp,lwave,norm,ascat
  character*64 LegendreFile
  CHARACTER(LEN=20)  scale
  Format1="(F18.15,F18.15,F18.15,F18.15)"

  ! read in legendre data
  read(5,*)
  read(5,1002) LegendreFile
  write(6,1002) LegendreFile
  read(5,*)
  read(5,*)
  read(5,*) LegPoints, ISTATE, ESTATE
  write(6,*) "LegPoints = ", LegPoints
  write(6,*) "ISTATE = ", ISTATE
  write(6,*) "ESTATE = ", ESTATE
  
  !     read in boundary conditions
  read(5,*)
  read(5,*)
  read(5,*) Order,Left,Right

  write(6,*) "order = ", order
  write(6,*) "Left = ", Left
  write(6,*) "Right = ", Right
  
  !     read in potential parameters
  read(5,*)
  read(5,*)
  read(5,*) alpha, PotInterpPoints, DD, x0

  write(6,*) "alpha = ", alpha
  write(6,*) "PotInterpPoints = ", PotInterpPoints
  write(6,*) "DD = ", DD
  write(6,*) "x0 = ", x0
!  mu=mass/2.0d0

  read(5,*)
  read(5,*)
  read(5,*) xNumPoints,xMin,xMax, scale
!  write(6,*) 'xNumPoints,xMin,xMax, scale'
 ! write(6,*) xNumPoints,xMin,xMax, scale
  write(6,*) "xNumPoints = ", xNumPoints
  write(6,*) "xMin = ", xMin
  write(6,*) "xMax = ", xMax
  write(6,*) "scale = ", scale

  read(5,*)
  read(5,*)
  read(5,*) lwave, EnergyMin, EnergyMax, NumEnergySteps
!  write(6,*) 'lwave, EnergyMin, EnergyMax, NumEnergySteps'
  !  write(6,*) lwave, EnergyMin, EnergyMax, NumEnergySteps
  write(6,*) "lwave = ", lwave
  write(6,*) "EnergyMin = ", EnergyMin
  write(6,*) "EnergyMax = ", EnergyMax
  write(6,*) "NumEnergySteps = ", NumEnergySteps
  EnergyMin = EnergyMin/HartreePermK
  EnergyMax = EnergyMax/HartreePermK
  VLIM=0d0
  call AtomData(ISTATE, AHfatom, iatom, satom, giatom, MU, MUREF, mass)
  !  mu = 0.5d0
  
  write(6,*) "mu = ", mu
  !DD = 10*4.899828818866035d0
  !x0=1d0
  
  xDimMin=xNumPoints+Order-3
  xDim=xDimMin
  if (Left .eq. 2) xDim = xDim + 1
  if (Right .eq. 2) xDim = xDim + 1
  MatrixDim=xDim
  write(6,*) 'xDim=',xDim, MatrixDim
  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  allocate(Energies(NumEnergySteps))
!  DeltaEnergy=(EnergyMax-EnergyMin)/dble(NumEnergySteps-1)
  call GridMaker(Energies,NumEnergySteps,EnergyMin,EnergyMax,'linear')
  
  allocate(delta(NumEnergySteps))

  allocate(xPoints(xNumPoints))
  allocate(xBounds(xNumPoints+2*Order))
  allocate(u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim)) ! 
  allocate(L(MatrixDim,MatrixDim))
  allocate(S(MatrixDim,MatrixDim))
  allocate(G0(MatrixDim,MatrixDim))
  allocate(V0(MatrixDim,MatrixDim))


  PotInterpOrder = 3
  allocate(R(PotInterpPoints),VS(PotInterpPoints),VT(PotInterpPoints))
  call GridMaker(R,PotInterpPoints,xMin,xMax,'linear')
!!$  do iR=1,PotInterpPoints
!!$     R(iR) = xmin + dble((iR-1)/(PotInterpPoints-1))**2d0 * xMax
!!$     write(6,*) iR, R(iR)
!!$  enddo
 allocate(Potknots(PotInterpPoints+PotInterpOrder),Vbcoef(PotInterpPoints))
 !  call GetSTPots(R,VS,VT,PotInterpPoints) ! use atomic units for R, VS, and VT

!  call SetupPotential(ISTATE, ESTATE, MU, MUREF, PotInterpPoints, VLIM, R, VS, Cvals)
!  call setup_Vinterp(R,PotInterpOrder,Potknots,PotInterpPoints,VS,Vbcoef)
!  do iR = 1, PotInterpPoints
!     write(200,*) R(iR), VS(iR)
!  enddo
!!$  do x=xMin,xMax,(xMax-xMin)/1000.0d0
!!$     !write(1000,*) x, VSech(DD,x0,x,lwave) 
!!$     write(1000,*) x, Vsqrc6(DD,x0,x,lwave) 
!!$  end do


!  call GridMaker(xPoints,xNumPoints,xMin,xMax,scale)
 call SplitGridMaker(xPoints,xNumPoints/2,xNumPoints/2,xMin,30d0,xMax)
!  do i=1,xNumPoints
!     write(6,*) i, xPoints(i)
!  enddo
  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,0,u)

  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,1,ux)

!  call CheckBasis(u,xDim,xNumPoints,LegPoints,xLeg,xPoints,700)

  G0=0.0d0
  write(6,*) "Calculating the Gamma matrix"
  call CalcGamma0(lwave,ISTATE,ESTATE)
  write(6,*) "Done"
  do iE=1,1!NumEnergySteps
     call CalcPhaseShift(Energies(iE),lwave,deltatemp,norm,ascat)
     delta(iE)=deltatemp
     write(6,*) Energies(iE), delta(iE), ascat
     call PlotCheckBessel(0d0, xMin, xMax, Energies(iE), mu, tan(delta(iE)), 10)
  enddo

  deallocate(G0,V0,S,L,delta)
  deallocate(xBounds,xPoints)
  deallocate(u,ux,xLeg,wLeg)

10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)

end program Main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$subroutine GetSTPots(R,VS,VT,NPP)
!!$  use units
!!$  implicit none
!!$  integer NPP,IMN1,IMN2,sym,lwave,ISTATE,iR
!!$  double precision R(NPP), VS(NPP), VT(NPP),XO(NPP)
!!$  double precision weights(NPP),RM2(NPP)
!!$  double precision VLIM
!!$  XO = R*BohrPerAngstrom !Convert to Angstrom to send to the Le Roy code
!!$  VLIM = 0d0
!!$
!!$  IMN1 = 7
!!$  IMN2 = 7
!!$  sym = 1 ! set to +1 for bosonic Li-7, -1 for fermionic Li-6, and 0 for mixture Li-6 - Li-7.
!!$  lwave = 0 ! s wave collisions
!!$
!!$  ISTATE = 1                ! Find the "singlet" X(^1\Sigma_g^+) curve
!!$  ! Call the Le Roy routine to generate the MLR potential
!!$  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VS)
!!$
!!$  ! print singlet potential to fort.10
!!$  do iR=1, NPP
!!$     VS(iR) = VS(iR)*InvcmPerHartree
!!$     write(10,*) R(iR), VS(iR)
!!$  enddo
!!$
!!$  ISTATE = 2                !Find the triplet potential
!!$  ! Call the Le Roy routine to generate the MLR potential
!!$  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VT)
!!$  do iR=1, NPP
!!$     VT(iR) = VT(iR)*InvcmPerHartree
!!$     write(30,*) R(iR), VT(iR)
!!$  enddo
!!$
!!$end subroutine GetSTPots

subroutine CalcPhaseShift(energy, lwave, deltatemp,norm,ascat)
  use EREdata
  implicit none
  double precision, external :: phirecon
  double precision, allocatable :: G(:,:), LL(:,:), ksqr(:), kcotandel(:),evec(:,:),eval(:),sig(:) ! 
  double precision ascat,energy, b, k, tandel, x, deltatemp,sj,sy,sjp,syp,lwave,norm!, ainv, halfr0, q, chi2,siga, sigb,a2inv ! 

  integer ix, ixp, n,ikeep

  allocate(G(MatrixDim,MatrixDim))
  allocate(LL(MatrixDim,MatrixDim))
  allocate(evec(MatrixDim,MatrixDim))
  allocate(eval(MatrixDim))

  b=0.0d0
  do ix=1,MatrixDim
     do ixp=1,MatrixDim
        G(ix,ixp)=0.0d0
        LL(ix,ixp)=0.0d0
        G(ix,ixp) = G0(ix,ixp) + 2.0d0*mu*energy*S(ix,ixp) - 2.0d0*mu*V0(ix,ixp) ! 
     enddo
  enddo
!!$  G(:,:) = 0d0
!!$  LL(:,:) = 0d0
!!$  G(:,:) = G0(:,:) + 2.0d0*mu*energy*S(:,:) - 2.0d0*mu*V0(:,:) ! 
  LL(MatrixDim,MatrixDim)=1.0d0
  write(6,*) "Calling Mydggev..."
  call Mydggev(MatrixDim,G,MatrixDim,LL,MatrixDim,eval,evec)
  write(6,*) "done."
  do n=1,MatrixDim
     if(dabs(eval(n)).gt.1d-14) then
        b=eval(n)
        ikeep=n
        write(6,*) 'eval(',n,') = ', b
     endif
  enddo

! b = negative log-derivative = -Y
  k=dsqrt(2.0d0*mu*energy)
  x=xMax
!  call sphbes(lwave,k*x,sj,sy,sjp,syp)
!  call sphbesjy(lwave,k*x,sj,sy,sjp,syp) ! change norm
!  tandel = ( (b*x + 1.d0) * sj + k*x*sjp)/( (b*x + 1.d0) * sy + k*x*syp)  
  tandel=(k*dcos(k*xMax) + b*dsin(k*xMax))/(-b*dcos(k*xMax) + k*dsin(k*xMax))
  deltatemp=datan(tandel)
!  norm = 2.d0*k*x*(sj-tandel*sy)/phirecon(xMax,ikeep,evec,Left,Right,xDim,MatrixDim,xNumPoints,xPoints,order)
  norm = (sin(k*x) + tandel*cos(k*x))/phirecon(xMax,ikeep,evec,Left,Right,xDim,MatrixDim,xNumPoints,xPoints,order)
  ascat = -tandel/k
!  write(6,*) 'b = ',b, 'freesol_b = ',-k*dcos(k*xMax)/dsin(k*xMax)
  !write(6,*) 'b = ',b, 'freesol_b = ',-k*sjp/sj
  x=0.0d0
  write(777,*) '# energy=',energy
  write(6,*) 'calling phirecon with ikeep = ',ikeep
  do while(x.lt.xMax)
     write(777,20) x, norm*phirecon(x,ikeep,evec,Left,Right,xDim,MatrixDim,xNumPoints,xPoints,order)
     x=x+0.01d0
  enddo
  write(777,*)
  
  deallocate(evec)
  deallocate(eval)
  deallocate(G, LL)

20 format(1P,100e16.8)      
end subroutine CalcPhaseShift
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine CalcGamma0(lwave,ISTATE,ESTATE)
  use EREdata
  use units
  use bspline
  implicit none
  double precision, external :: VSech,vsqrc6
  double precision ax, bx,x,xScaledZero,xIntScale,TempG,TempS,TempV,a,lwave,xdum(1),VV(1)
  double precision Cvals(4)
  integer, allocatable :: kxMin(:,:),kxMax(:,:)
  integer kx,ix,ixp,lx,n,ISTATE,ESTATE

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
        G0(ix,ixp)=0.0d0
        S(ix,ixp)=0.0d0
        L(ix,ixp)=0.0d0
        V0(ix,ixp)=0.0d0
        do kx = kxMin(ix,ixp),kxMax(ix,ixp)
           ax = xPoints(kx)
           bx = xPoints(kx+1)
           xIntScale=0.5d0*(bx-ax)
           xScaledZero=0.5d0*(ax+bx)
           TempG = 0.0d0
           TempS = 0.0d0
           TempV = 0.0d0
           do lx = 1,LegPoints
              a = wLeg(lx)*xIntScale
              x=xIntScale*xLeg(lx)+xScaledZero
              TempS = TempS + a*u(lx,kx,ix)*u(lx,kx,ixp)
              !VV(1) = alpha*VSech(DD,x0,x,lwave)
              !TempV = TempV + a*u(lx,kx,ix)*(alpha*VSech(DD,x0,x,lwave))*u(lx,kx,ixp)
              !TempV = TempV + a*u(lx,kx,ix)*(alpha*Vsqrc6(DD,x0,x,lwave))*u(lx,kx,ixp)

              xdum(1) = x*BohrPerAngstrom
              call SetupPotential(ISTATE, ESTATE, MU, MUREF, 1, 0d0, xdum, VV, Cvals)

              write(201,*) x, VV(1)
              !TempV = TempV + a*u(lx,kx,ix)* &
              !     (alpha*dbsval(x,PotInterpOrder,Potknots,PotInterpPoints,Vbcoef))* &
              !     u(lx,kx,ixp)
              TempV = TempV + a*u(lx,kx,ix)*(alpha*VV(1))*u(lx,kx,ixp)
              TempG = TempG + a*(-ux(lx,kx,ix)*ux(lx,kx,ixp))
              !c                  print*,a, TempS, TempG
           enddo
           S(ix,ixp) = S(ix,ixp) + TempS ! place values into overlap matrix
           G0(ix,ixp) = G0(ix,ixp) + TempG ! place values into Gamma
           V0(ix,ixp) = V0(ix,ixp) + TempV ! place values into Potential matrix
        enddo
!        write(6,*) ix,ixp,' G0 = ',G0(ix,ixp), ' S = ',S(ix,ixp), ' V0 = ',V0(ix,ixp)
     enddo
  enddo

  L(xDim,xDim)=1.0d0

  deallocate(kxMax,kxMin)
10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)
  
end subroutine CalcGamma0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$subroutine GridMaker(xNumPoints,xMin,xMax,xPoints)
!!$  implicit none
!!$  integer xNumPoints
!!$  double precision xMin,xMax,xPoints(xNumPoints)
!!$
!!$  integer i,j,k
!!$  double precision Pi
!!$  double precision r0New
!!$  double precision xRswitch
!!$  double precision xDelt,x0,x1,x2
!!$
!!$  Pi = 3.1415926535897932385d0
!!$
!!$  x0 = xMin
!!$  x1 = xMax
!!$  k = 1
!!$  xDelt = (x1-x0)/dfloat(xNumPoints-1)
!!$  do i = 1,xNumPoints
!!$     xPoints(k) = (i-1)*xDelt + x0
!!$     k = k + 1
!!$  enddo
!!$
!!$15 format(6(1x,1pd12.5))
!!$
!!$
!!$  return
!!$end subroutine GridMaker

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
     write(file,*) 
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setup_Vinterp(R,PotInterpOrder,Potknots,PotInterpPoints,VData,Vbcoef)
  use bspline

  implicit none
  integer PotInterpPoints,PotInterpOrder
  double precision VData(PotInterpPoints),R(PotInterpPoints) ! 
  double precision Potknots(PotInterpPoints+PotInterpOrder), Vbcoef(PotInterpPoints) ! 
  
  call dbsnak(PotInterpPoints,R,PotInterpOrder,Potknots)      
  call dbsint(PotInterpPoints,R,VData,PotInterpOrder,Potknots,Vbcoef) ! 
end subroutine setup_Vinterp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$double precision function VEffective(x)
!!$  use Vparams
!!$  double precision, external :: VSech
!!$  implicit none
!!$  VEffective = VSech(DD,x0,x) + (leff*(leff+1))/(2.0d0*mu*x*x)
!!$  return
!!$end function VEffective
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  The following subroutine will plot spherical bessel functions for the following inputs:
! lwave: the angular momentum partial wave (I allow this to be noninteger for certain hyperspherical calculations)
! energy: the energy
! (rMin, rMax) is the interval over which the functions are to be plotted
! file: an integer denoting that the x-y plot is to be written to a file named fort.file
! tandel: the tangent of the phase-shift delta_l.
! mu: the reduced mass
subroutine PlotCheckBessel(lwave,rMin,rMax,energy,mu,tandel,file)
  implicit none
  integer NumSteps,iStep,file
  double precision k,dr,sj,sy,sjp,syp,tandel,cosdel,sindel,delta,rMin,rMax,energy,mu,lwave,x
  double precision, allocatable :: r(:)
  delta=datan(tandel)
!  sindel=dsin(delta)
!  cosdel=dcos(delta)
  NumSteps=1000
  k=dsqrt(2.0d0*mu*energy)
  dr=(rMax-rMin)/(NumSteps-1)

  allocate(r(NumSteps))
  do iStep=1,NumSteps
     r(iStep)=rMin+(iStep-1)*dr
     x = r(iStep)*k     
!     call sphbesjy(lwave,x,sj,sy,sjp,syp) ! change norm
!     call sphbes(lwave,k*r(iStep),sj,sy,sjp,syp)
!     write(file,*) r(iStep), k*r(iStep)*sj, k*r(iStep)*sy, 2.d0*k*r(iStep)*(sj-tandel*sy)
!     write(file,*) r(iStep),  (sj-tandel*sy)!2.d0*k*r(iStep)*(sj-tandel*sy)
     write(file,*) r(iStep),  sin(x) + tandel*cos(x)
  enddo
  write(file, *) 

end subroutine PlotCheckBessel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      SUBROUTINE sphbes(xnu,x,sj,sy,sjp,syp)
!!$!      use nrtype, only : dbl
!!$!      INTEGER n
!!$!      double precision sj,sjp,sy,syp,x,xnu
!!$      double precision sj,sjp,sy,syp,x,xnu
!!$      !    USES bessjy
!!$      double precision factor,order,rj,rjp,ry,ryp,RTPIO2
!!$      PARAMETER (RTPIO2=1.25331413731550d0)
!!$      if(n.lt.0.d0.or.x.le.0.d0) then
!!$         write(6,*) 'Bad arguments in sphbes. Press ENTER to continue.'
!!$         read(5,*)
!!$      endif
!!$      order=xnu+0.5d0
!!$      call bessjy(x,order,rj,ry,rjp,ryp)
!!$      factor=RTPIO2/Dsqrt(x)
!!$      sj=factor*rj
!!$      sy=factor*ry
!!$      sjp=factor*rjp-sj/(2.d0*x)
!!$      syp=factor*ryp-sy/(2.d0*x)
!!$
!!$      return
!!$      END
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
!!$!      use nrtype, only : dbl
!!$      INTEGER MAXIT
!!$      double precision rj,rjp,ry,ryp,x,xnu,XMIN
!!$      DOUBLE PRECISION EPS,FPMIN,PI
!!$      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0, &
!!$           PI=3.141592653589793d0)
!!$!    USES beschb
!!$      INTEGER i,isign,l,nl
!!$      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e, &
!!$           f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q, &
!!$           r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1, &
!!$           temp,w,x2,xi,xi2,xmu,xmu2
!!$      !write(66,*)xnu  
!!$      if(x.le.0.d0.or.xnu.lt.0.d0) then
!!$         write(6,*) 'Bad arguments in bessjy. Press ENTER to continue.'
!!$         read(5,*)
!!$      endif
!!$
!!$      if(x.lt.XMIN)then
!!$        nl=int(xnu+.5d0)
!!$      else
!!$        nl=max(0,int(xnu-x+1.5d0))
!!$      endif
!!$      xmu=xnu-nl
!!$      xmu2=xmu*xmu
!!$      xi=1.d0/x
!!$      xi2=2.d0*xi
!!$      w=xi2/PI
!!$      isign=1
!!$      h=xnu*xi
!!$      if(h.lt.FPMIN)h=FPMIN
!!$      b=xi2*xnu
!!$      d=0.d0
!!$      c=h
!!$      do 11 i=1,MAXIT
!!$        b=b+xi2
!!$        d=b-d
!!$        if(abs(d).lt.FPMIN)d=FPMIN
!!$        c=b-1.d0/c
!!$        if(abs(c).lt.FPMIN)c=FPMIN
!!$        d=1.d0/d
!!$        del=c*d
!!$        h=del*h
!!$        if(d.lt.0.d0)isign=-isign
!!$        if(abs(del-1.d0).lt.EPS)goto 1
!!$11    continue
!!$!      pause 'x too large in bessjy; try asymptotic expansion'
!!$      goto 31
!!$1     continue
!!$      rjl=isign*FPMIN
!!$      rjpl=h*rjl
!!$      rjl1=rjl
!!$      rjp1=rjpl
!!$      fact=xnu*xi
!!$      do 12 l=nl,1,-1
!!$        rjtemp=fact*rjl+rjpl
!!$        fact=fact-xi
!!$        rjpl=fact*rjtemp-rjl
!!$        rjl=rjtemp
!!$12    continue
!!$      if(rjl.eq.0.d0)rjl=EPS
!!$      f=rjpl/rjl
!!$      if(x.lt.XMIN) then
!!$        x2=.5d0*x
!!$        pimu=PI*xmu
!!$        if(abs(pimu).lt.EPS)then
!!$          fact=1.d0
!!$        else
!!$          fact=pimu/sin(pimu)
!!$        endif
!!$        d=-log(x2)
!!$        e=xmu*d
!!$        if(abs(e).lt.EPS)then
!!$          fact2=1.d0
!!$        else
!!$          fact2=sinh(e)/e
!!$        endif
!!$        call beschb(xmu,gam1,gam2,gampl,gammi)
!!$        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
!!$        e=exp(e)
!!$        p=e/(gampl*PI)
!!$        q=1.d0/(e*PI*gammi)
!!$        pimu2=0.5d0*pimu
!!$        if(abs(pimu2).lt.EPS)then
!!$          fact3=1.d0
!!$        else
!!$          fact3=sin(pimu2)/pimu2
!!$        endif
!!$        r=PI*pimu2*fact3*fact3
!!$        c=1.d0
!!$        d=-x2*x2
!!$        sum=ff+r*q
!!$        sum1=p
!!$        do 13 i=1,MAXIT
!!$          ff=(i*ff+p+q)/(i*i-xmu2)
!!$          c=c*d/i
!!$          p=p/(i-xmu)
!!$          q=q/(i+xmu)
!!$          del=c*(ff+r*q)
!!$          sum=sum+del
!!$          del1=c*p-i*del
!!$          sum1=sum1+del1
!!$          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
!!$13      continue
!!$          write(6,*) 'bessy series failed to converge'
!!$          read(5,*)
!!$2       continue
!!$        rymu=-sum
!!$        ry1=-sum1*xi2
!!$        rymup=xmu*xi*rymu-ry1
!!$        rjmu=w/(rymup-f*rymu)
!!$      else
!!$        a=.25d0-xmu2
!!$        p=-.5d0*xi
!!$        q=1.d0
!!$        br=2.d0*x
!!$        bi=2.d0
!!$        fact=a*xi/(p*p+q*q)
!!$        cr=br+q*fact
!!$        ci=bi+p*fact
!!$        den=br*br+bi*bi
!!$        dr=br/den
!!$        di=-bi/den
!!$        dlr=cr*dr-ci*di
!!$        dli=cr*di+ci*dr
!!$        temp=p*dlr-q*dli
!!$        q=p*dli+q*dlr
!!$        p=temp
!!$        do 14 i=2,MAXIT
!!$          a=a+2*(i-1)
!!$          bi=bi+2.d0
!!$          dr=a*dr+br
!!$          di=a*di+bi
!!$          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
!!$          fact=a/(cr*cr+ci*ci)
!!$          cr=br+cr*fact
!!$          ci=bi-ci*fact
!!$          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
!!$          den=dr*dr+di*di
!!$          dr=dr/den
!!$          di=-di/den
!!$          dlr=cr*dr-ci*di
!!$          dli=cr*di+ci*dr
!!$          temp=p*dlr-q*dli
!!$          q=p*dli+q*dlr
!!$          p=temp
!!$          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
!!$14      continue
!!$          write(6,*) 'cf2 failed in bessjy. Press enter to continue.'
!!$          read(5,*)
!!$3       continue
!!$        gam=(p-f)/q
!!$        rjmu=sqrt(w/((p-f)*gam+q))
!!$        rjmu=sign(rjmu,rjl)
!!$        rymu=rjmu*gam
!!$        rymup=rymu*(p+q/gam)
!!$        ry1=xmu*xi*rymu-rymup
!!$      endif
!!$      fact=rjmu/rjl
!!$      rj=rjl1*fact
!!$      rjp=rjp1*fact
!!$      do 15 i=1,nl
!!$        rytemp=(xmu+i)*xi2*ry1-rymu
!!$        rymu=ry1
!!$        ry1=rytemp
!!$15    continue
!!$      ry=rymu
!!$      ryp=xnu*xi*rymu-ry1
!!$      return
!!$31    Continue
!!$!c  asymptotic expansion
!!$      rj = Dsqrt(2.d0/(pi*x))*Dcos(x-xnu*pi/2.d0-pi/4.d0)
!!$      ry = Dsqrt(2.d0/(pi*x))*Dsin(x-xnu*pi/2.d0-pi/4.d0)
!!$      rjp = -ry
!!$      ryp = rj
!!$      return
!!$      END
!!$
!!$
!!$      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
!!$!      use nrtype, only: dbl  
!!$      INTEGER NUSE1,NUSE2
!!$      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
!!$      PARAMETER (NUSE1=7,NUSE2=8)
!!$!CU    USES chebev
!!$      double precision xx,c1(7),c2(8),chebev
!!$      SAVE c1,c2
!!$      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4, &
!!$           -3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
!!$      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3, &
!!$           -4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
!!$      xx=8.d0*x*x-1.d0
!!$      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
!!$      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
!!$      gampl=gam2-x*gam1
!!$      gammi=gam2+x*gam1
!!$      return
!!$    END SUBROUTINE beschb
!!$
!!$    
!!$    FUNCTION chebev(a,b,c,m,x)
!!$!      use nrtype, only : dbl
!!$      INTEGER m
!!$      double precision chebev,a,b,x,c(m)
!!$      INTEGER j
!!$      double precision d,dd,sv,y,y2
!!$      if ((x-a)*(x-b).gt.0.d0) then
!!$         write(6,*) 'x not in range in chebev.  Press enter to continue'
!!$         read(5,*)
!!$      endif
!!$      d=0.d0
!!$      dd=0.d0
!!$      y=(2.d0*x-a-b)/(b-a)
!!$      y2=2.d0*y
!!$      do 11 j=m,2,-1
!!$        sv=d
!!$        d=y2*d-dd+c(j)
!!$        dd=sv
!!$11    continue
!!$      chebev=y*d-dd+0.5d0*c(1)
!!$      return
!!$      END         
!!$
!!$      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
!!$      INTEGER MAXIT
!!$      DOUBLE PRECISION ri,rip,rk,rkp,x,xnu,XMIN
!!$      DOUBLE PRECISION EPS,FPMIN,PI
!!$      PARAMETER (EPS=1.e-10,FPMIN=1.e-30,MAXIT=10000,XMIN=2., &
!!$           PI=3.141592653589793d0)
!!$!CU    USES beschb
!!$      INTEGER i,l,nl
!!$      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff, &
!!$           gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1, &
!!$           ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
!!$      if(x.le.0..or.xnu.lt.0.) then
!!$         write(6,*) 'bad arguments in bessik. Press enter to continue'
!!$         read(5,*)
!!$      endif
!!$      nl=int(xnu+.5d0)
!!$      xmu=xnu-nl
!!$      xmu2=xmu*xmu
!!$      xi=1.d0/x
!!$      xi2=2.d0*xi
!!$      h=xnu*xi
!!$      if(h.lt.FPMIN)h=FPMIN
!!$      b=xi2*xnu
!!$      d=0.d0
!!$      c=h
!!$      do 11 i=1,MAXIT
!!$        b=b+xi2
!!$        d=1.d0/(b+d)
!!$        c=b+1.d0/c
!!$        del=c*d
!!$        h=del*h
!!$        if(abs(del-1.d0).lt.EPS)goto 1
!!$11    continue
!!$        write(6,*)'x too large in bessik; try asymptotic expansion.  Press enter to continue'
!!$        read(5,*)
!!$        
!!$1     continue
!!$      ril=FPMIN
!!$      ripl=h*ril
!!$      ril1=ril
!!$      rip1=ripl
!!$      fact=xnu*xi
!!$      do 12 l=nl,1,-1
!!$        ritemp=fact*ril+ripl
!!$        fact=fact-xi
!!$        ripl=fact*ritemp+ril
!!$        ril=ritemp
!!$12    continue
!!$      f=ripl/ril
!!$      if(x.lt.XMIN) then
!!$        x2=.5d0*x
!!$        pimu=PI*xmu
!!$        if(abs(pimu).lt.EPS)then
!!$          fact=1.d0
!!$        else
!!$          fact=pimu/sin(pimu)
!!$        endif
!!$        d=-log(x2)
!!$        e=xmu*d
!!$        if(abs(e).lt.EPS)then
!!$          fact2=1.d0
!!$        else
!!$          fact2=sinh(e)/e
!!$        endif
!!$        call beschb(xmu,gam1,gam2,gampl,gammi)
!!$        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
!!$        sum=ff
!!$        e=exp(e)
!!$        p=0.5d0*e/gampl
!!$        q=0.5d0/(e*gammi)
!!$        c=1.d0
!!$        d=x2*x2
!!$        sum1=p
!!$        do 13 i=1,MAXIT
!!$          ff=(i*ff+p+q)/(i*i-xmu2)
!!$          c=c*d/i
!!$          p=p/(i-xmu)
!!$          q=q/(i+xmu)
!!$          del=c*ff
!!$          sum=sum+del
!!$          del1=c*(p-i*ff)
!!$          sum1=sum1+del1
!!$          if(abs(del).lt.abs(sum)*EPS)goto 2
!!$13      continue
!!$          write(6,*) 'bessk series failed to converge.  Press enter to continue.'
!!$          read(5,*)
!!$2       continue
!!$        rkmu=sum
!!$        rk1=sum1*xi2
!!$      else
!!$        b=2.d0*(1.d0+x)
!!$        d=1.d0/b
!!$        delh=d
!!$        h=delh
!!$        q1=0.d0
!!$        q2=1.d0
!!$        a1=.25d0-xmu2
!!$        c=a1
!!$        q=c
!!$        a=-a1
!!$        s=1.d0+q*delh
!!$        do 14 i=2,MAXIT
!!$          a=a-2*(i-1)
!!$          c=-a*c/i
!!$          qnew=(q1-b*q2)/a
!!$          q1=q2
!!$          q2=qnew
!!$          q=q+c*qnew
!!$          b=b+2.d0
!!$          d=1.d0/(b+a*d)
!!$          delh=(b*d-1.d0)*delh
!!$          h=h+delh
!!$          dels=q*delh
!!$          s=s+dels
!!$          if(abs(dels/s).lt.EPS)goto 3
!!$14      continue
!!$          write(6,*) 'bessik: failure to converge in cf2.  Press enter to continue.'
!!$          read(5,*)
!!$3       continue
!!$        h=a1*h
!!$        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
!!$        rk1=rkmu*(xmu+x+.5d0-h)*xi
!!$      endif
!!$      rkmup=xmu*xi*rkmu-rk1
!!$      rimu=xi/(f*rkmu-rkmup)
!!$      ri=(rimu*ril1)/ril
!!$      rip=(rimu*rip1)/ril
!!$      do 15 i=1,nl
!!$        rktemp=(xmu+i)*xi2*rk1+rkmu
!!$        rkmu=rk1
!!$        rk1=rktemp
!!$15    continue
!!$      rk=rkmu
!!$      rkp=xnu*xi*rkmu-rk1
!!$      return
!!$      END
!!$      !C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
      
Subroutine AtomData(ISTATE, AHf, i, s, gi, MU, MUREF, mass)
  implicit none
  integer ISTATE, i, s 
  double precision AHf, gi, MU, MUREF, mass
  double precision amuAU

  amuAU = 1.66053906660d-27/9.109383701528d-31  !mass conversion to atomic units
  s = 1                    !electronic spin

  ! Ahf is the hyperfine coupling in MHz.  gu
  select case (ISTATE)
  case (1) !Rb-87
     write(6,'(A)') "Setting up the potential for Rb-87"
     gi = -0.000995141410d0
     i = 3
     AHf = 3417.3413064215d0
     mass = 86.909180527*amuAU
     MU = mass/2
     MUREF = MU

  case (2) !Rb-85
     write(6,'(A)') "Setting up the potential for Rb-85"
     gi = -0.00029364006d0
     i = 5
     AHf = 1011.9108132d0
     mass = 84.911789738*amuAU
     MU = mass/2  
     MUREF = (86.909180527*amuAU)/2 !Reference reduced mass is 87-Rb2

  case (3) !K-39
     write(6,'(A)') "Setting up the potential for K-39"
     gi = -0.00014193489d0
     i = 3
     AHf = 230.85986013d0
     mass = 38.963708d0*amuAU
     MU = mass/2
     MUREF = MU

  case (4) !K-40
     write(6,'(A)') "Setting up the potential for K-40"
     gi = 0.00017649034d0
     i = 8
     AHf = -285.730824d0
     mass = 39.964008d0*amuAU
     MU = mass/2
     MUREF = (38.963708d0*amuAU)/2 !Reference reduced mass is 39-K2

  case (5) !Na-23
     write(6,'(A)') "Setting up the potential for Na-23"
     gi =  -0.00080461088d0
     i = 3
     AHf = 885.81306445d0
     mass = 22.98976928*amuAU
     MU = mass/2
     MUREF = MU

  case (6) !Li-6
     write(6,'(A)') "Setting up the potential for Li-6"
     gi = -0.00044765403d0
     i = 2
     AHf = 152.136840720d0
     mass = 6.0151223*amuAU
     MU = mass/2
     MUREF = (7.016004*amuAU)/2 !Reference reduced mass is 7-Li2

  case(7) !Li-7
     write(6,'(A)') "Setting up the potential for Li-7"
     gi =  -0.0011822136d0
     i = 3
     AHf = 401.75204335d0
     mass = (7.016004*amuAU)
     MU = mass/2
     MUREF = MU

  case(8)  !Cs-133
     write(6,'(A)') "Setting up the potential for Cs-133"
     gi = -0.0003988539552
     i = 7
     AHf = 2298.1579425
     mass = (132.905451933*amuAU)
     MU = mass/2
     MUREF = MU

  case default
     write(6,'(A)') "Invalid ISTATE value"
     stop

  end select



End Subroutine AtomData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine VdWLength(Cvals,beta,mu)
  use units
  implicit none

  double precision beta4,beta, C6g, C8g, C10g, Cvals(4), mu
  double precision OneThird
  C6g=Cvals(1)
  C8g=Cvals(2)
  C10g=Cvals(3)
  OneThird = 0.3333333333333333333333333d0
  if((C8g.gt.0d0).and.(C10g.gt.0d0)) then
     ! This nasty expression is from a Mathematica Solve[] command, exported to FortranForm
     beta4 = C6g*mu + Sqrt(6*C6g**2*mu**2 - 4*mu*(-C10g + C6g**2*mu) + &
          (4*mu**2*(4*C10g**2 + 4*C10g*C6g**2*mu + C6g*mu*(-3*C8g**2 + C6g**3*mu)))/ &
          (64*C10g**3*mu**3 + 96*C10g**2*C6g**2*mu**4 + 27*C8g**4*mu**4 - &
          36*C6g**3*C8g**2*mu**5 + 8*C6g**6*mu**6 + &
          24*C10g*C6g*mu**4*(-3*C8g**2 + 2*C6g**3*mu) + &
          3*sqrt(3d0)*Sqrt(C8g**4*mu**7* &
          (128*C10g**3 + 128*C10g**2*C6g**2*mu + &
          C8g**2*mu*(27*C8g**2 - 8*C6g**3*mu) + &
          16*C10g*C6g*mu*(-9*C8g**2 + 2*C6g**3*mu))))**OneThird + &
          (64*C10g**3*mu**3 + 96*C10g**2*C6g**2*mu**4 + 27*C8g**4*mu**4 - &
          36*C6g**3*C8g**2*mu**5 + 8*C6g**6*mu**6 + &
          24*C10g*C6g*mu**4*(-3*C8g**2 + 2*C6g**3*mu) + &
          3*sqrt(3d0)*Sqrt(C8g**4*mu**7* &
          (128*C10g**3 + 128*C10g**2*C6g**2*mu + C8g**2*mu*(27*C8g**2 - 8*C6g**3*mu) + &
          16*C10g*C6g*mu*(-9*C8g**2 + 2*C6g**3*mu))))**OneThird)/sqrt(6d0) + &
          Sqrt(8*C6g**2*mu**2 - (16*mu*(-C10g + C6g**2*mu))/3. - &
          (8*mu**2*(4*C10g**2 + 4*C10g*C6g**2*mu + C6g*mu*(-3*C8g**2 + C6g**3*mu)))/ &
          (3.*(64*C10g**3*mu**3 + 96*C10g**2*C6g**2*mu**4 + 27*C8g**4*mu**4 - &
          36*C6g**3*C8g**2*mu**5 + 8*C6g**6*mu**6 + &
          24*C10g*C6g*mu**4*(-3*C8g**2 + 2*C6g**3*mu) + &
          3*sqrt(3d0)*Sqrt(C8g**4*mu**7* &
          (128*C10g**3 + 128*C10g**2*C6g**2*mu + &
          C8g**2*mu*(27*C8g**2 - 8*C6g**3*mu) + &
          16*C10g*C6g*mu*(-9*C8g**2 + 2*C6g**3*mu))))**OneThird) - &
          (2*(64*C10g**3*mu**3 + 96*C10g**2*C6g**2*mu**4 + 27*C8g**4*mu**4 - &
          36*C6g**3*C8g**2*mu**5 + 8*C6g**6*mu**6 + &
          24*C10g*C6g*mu**4*(-3*C8g**2 + 2*C6g**3*mu) + &
          3*sqrt(3d0)*Sqrt(C8g**4*mu**7* &
          (128*C10g**3 + 128*C10g**2*C6g**2*mu + &
          C8g**2*mu*(27*C8g**2 - 8*C6g**3*mu) + &
          16*C10g*C6g*mu*(-9*C8g**2 + 2*C6g**3*mu))))**OneThird)/3. + &
          (4*sqrt(6d0)*C8g**2*mu**2)/ &
          Sqrt(6*C6g**2*mu**2 - 4*mu*(-C10g + C6g**2*mu) + &
          (4*mu**2*(4*C10g**2 + 4*C10g*C6g**2*mu + C6g*mu*(-3*C8g**2 + C6g**3*mu)))/ &
          (64*C10g**3*mu**3 + 96*C10g**2*C6g**2*mu**4 + 27*C8g**4*mu**4 - &
          36*C6g**3*C8g**2*mu**5 + 8*C6g**6*mu**6 + &
          24*C10g*C6g*mu**4*(-3*C8g**2 + 2*C6g**3*mu) + &
          3*sqrt(3d0)*Sqrt(C8g**4*mu**7* &
          (128*C10g**3 + 128*C10g**2*C6g**2*mu + &
          C8g**2*mu*(27*C8g**2 - 8*C6g**3*mu) + &
          16*C10g*C6g*mu*(-9*C8g**2 + 2*C6g**3*mu))))**OneThird + &
          (64*C10g**3*mu**3 + 96*C10g**2*C6g**2*mu**4 + 27*C8g**4*mu**4 - &
          36*C6g**3*C8g**2*mu**5 + 8*C6g**6*mu**6 + &
          24*C10g*C6g*mu**4*(-3*C8g**2 + 2*C6g**3*mu) + &
          3*sqrt(3d0)*Sqrt(C8g**4*mu**7* &
          (128*C10g**3 + 128*C10g**2*C6g**2*mu + &
          C8g**2*mu*(27*C8g**2 - 8*C6g**3*mu) + &
          16*C10g*C6g*mu*(-9*C8g**2 + 2*C6g**3*mu))))**OneThird))/2.
     

     beta = beta4**(0.25d0)
  else
     beta = (2d0*mu*C6g)**0.25d0
  endif
     
  
end subroutine VdWLength
SUBROUTINE SetupPotential(ISTATE, ESTATE, MU, MUREF, NPP, VLIM, XO, VV, Cvals)
  use units
  implicit none
  integer NPP, ISTATE, ESTATE, i, j, N, IMN1, IMN2, iphi
  double precision, intent(in) :: XO(NPP)
  double precision VV(NPP), RM2(NPP), Cvals(4)
  double precision MU, VLIM, RSR, RLR, NS, U1, U2, B, RM, xi, MUREF
  double precision C6, C8, C10, C26, Aex, gamma, beta, nu0, nu1
  double precision, allocatable :: A(:)
  double precision De, re, ref, alpha, p, q, s, bds, cds, rhoAB
  double precision yqref, ypref, ypeq, uLR, phiMLR, phiMLRtemp
  double precision phiinf, uLRre, DDS6, DDS8, DDS10, DDS6re, DDS8re, DDS10re
  double precision Scorr ! the strength of the short-range correction for Li
  
  if (ESTATE .EQ. 1) then  !Ground state singlet 
     select case (ISTATE)  

     case (1:2)  !Rubidium
        N = 26
        allocate(A(N))
        RSR = 3.126d0
        NS = 4.53389d0
        U1 = -0.638904880d4
        U2 = 0.112005361d7
        B = -0.13d0
        RLR = 11.00d0
        RM = 4.209912706d0
        C6 = 0.2270032d8
        C8 = 0.7782886d9
        C10 = 0.2868869d11
        C26 = 0.2819810d26
        Aex = 0.1317786d5
        gamma = 5.317689d0
        beta = 2.093816d0
        nu0 = 0d0
        nu1 = 0d0
        A = (/-3993.592873d0, 0.0d0, 0.282069372972346137d5, 0.560425000209256905d4, -0.423962138510562945d5, &
             -0.598558066508841584d5, &
             -0.162613532034769596d5,-0.405142102246254944d5, 0.195237415352729586d6, 0.413823663033582852d6, &
             -0.425543284828921501d7, 0.546674790157210198d6, 0.663194778861331940d8,-0.558341849704095051d8, &
             -0.573987344918535471d9, 0.102010964189156187d10, 0.300040150506311035d10,-0.893187252759830856d10, &
             -0.736002541483347511d10, 0.423130460980355225d11,-0.786351477693491840d10,-0.102470557344862152d12, &
             0.895155811349267578d11,0.830355322355692902d11,-0.150102297761234375d12,0.586778574293387070d11 /)
     case (3:4)  !Potassium
        N = 32
        allocate(A(N))
        RSR = 2.870d0
        NS = 12d0
        U1 = -0.263145571d4
        U2 = 0.813723194d9
        B = -0.400d0
        RLR = 12.00d0
        RM = 3.92436437d0
        C6 = 0.1892652670d8
        C8 =  0.5706799527d9
        C10 = 0.1853042722d11
        C26 = 0.0d0
        Aex = 0.90092159d4
        gamma = 5.19500d0
        beta = 2.13539d0
        A = (/-4450.899484d0, 0.30601009538111d-1, 0.13671217000518d5,0.10750910095361d5, &
             -0.20933401680991d4,-0.19385874804675d5,-0.49208915890513d5,0.11026639220148d6, &
             0.72867339500920d6,-0.29310679369135d7,-0.12407070106619d8,0.40333947198094d8, &
             0.13229848871390d9, -0.37617673798775d9,-0.95250413275787d9,0.24655585744641d10, &
             0.47848257695164d10,-0.11582132109947d11,-0.17022518297651d11,0.39469335034593d11, &
             0.43141949844339d11,-0.97616955325128d11,-0.77417530685917d11,0.17314133615879d12, &
             0.96118849114926d11,-0.21425463041449d12,-0.78513081754125d11, 0.17539493131251d12, &
             0.37939637008662d11,-0.85271868691526d11,-0.82123523240949d10,0.18626451751424d11/)
        nu0 = 0.13148609d0
        nu1 = 2.08523853d0
     case (5)   !Sodium
        N = 27
        allocate(A(N))
        RSR = 2.181d0
        U1 = -0.78531844d4
        U2 = 0.84258655d6
        RM = 3.0788576d0
        NS = 6
        RLR = 11.00d0
        B = -0.140d0
        nu0 = 0d0
        nu1 = 0d0
        C6 = 0.75186131d7
        C8 = 0.1686430d9
        C10 = 0.3081961d10
        C26 = 0d0
        Aex = 0.40485835d5
        gamma = 4.59105d0
        beta = 2.36594d0
        A = (/-6022.04193d0, -0.200727603516760356d1, 0.302710123527149044d5, 0.952678499004718833d4, &
             -0.263132712461278206d5,-0.414199125447689439d5, 0.100454724828577862d6, &
             0.950433282843468915d5, -0.502202855817934591d7, -0.112315449566019326d7, &
             0.105565865633448541d9,-0.626929930064849034d8, -0.134149332172454119d10,&
             0.182316049840707183d10, 0.101425117010240822d11, -0.220493424364290123d11, &
             -0.406817871927934494d11, 0.144231985783280396d12, 0.379491474653734665d11, &
             -0.514523137448139771d12, 0.342211848747264038d12, 0.839583017514805054d12, &
             -0.131052566070353687d13, -0.385189954553600769d11, 0.135760501276292969d13, &
             -0.108790546442390417d13, 0.282033835345282288d12/)

     case (6) !Lithium-6
        IMN1 = 6
        IMN2 = IMN1
        !C6 = 6.718500D+06
        !C8 = 1.126290D+08
        !C10 = 2.786830D+09

        C6 = 6.71527d6
        C8 = 1.125880D+08
        C10 = 2.786040D+09
        call POTGENLI2(1,IMN1,IMN2,NPP,VLIM,XO,RM2,VV)!,.FALSE.)
        
        re = 2.6729932d0
        Scorr = 1.337d-6 * HartreePerInvcm/(BohrPerAngstrom**2)
        do i=1,NPP
           if(XO(i).lt.re) then
              VV(i) = VV(i) + Scorr*(XO(i) - re)**2
           endif
        enddo
        
     case (7) !Lithium-7
        IMN1 = 7
        IMN2 = IMN1
        C6 = 6.71527d6
        C8 = 1.12588d8
        C10 = 2.78604d9
        call POTGENLI2(1,IMN1,IMN2,NPP,VLIM,XO,RM2,VV)!,.FALSE.)
        ! No short-range correction needed for Li-7 in the singlet state.  This already gives a scattering length of 34.331 bohr
        
     case (8) !Cesium
        N = 23
        allocate(A(N))
        p = 5d0
        q = p
        phiinf = 0.9469102524695
        re = 4.647967771d0
        De = 3650.041851d0
        C6 = 3.3164d7
        C8 = 1.38d9
        C10 = 6.01d10
        ref = 6.2d0
        A = (/0.0949905d0, -0.372698d0, -0.04090736d0, 0.129657d0, &
             0.1486696d0, 0.171395d0, 0.295883d0, 0.547375d0, &
             -1.14615d0, -2.7883d0, 9.98557d0, 1.69149d1, -4.17899d1, &
             -5.76544d1, 1.08881d2, 1.24037d2, -1.716d2, -1.6159d2, &
             1.5781d2, 1.175d2, -7.5d1, -3.67d1, 1.3d1/)

     case default
        write(6,*) "Invalid ISTATE value"
        stop
     end select

  else if (ESTATE .EQ. 3) then  !Ground state triplet
     select case (ISTATE)

     case (1:2) !Rubidium
        N = 13
        allocate(A(N))
        RSR = 5.07d0
        NS = 4.5338950d0
        U1 = -0.619088543d3
        U2 = 0.956231677d6
        B = -0.33d0
        RLR = 11.00d0
        RM = 6.093345d0
        C6 = 0.2270032d8
        C8 = 0.7782886d9
        C10 = 0.2868869d11
        C26 = 0.2819810d26
        Aex = 0.1317786d5
        gamma = 5.317689d0
        beta = 2.093816d0
        nu0 = 0d0
        nu1 = 0d0
        A = (/-241.503352d0, -0.672503402304666542d0, 0.195494577140503543d4, -0.141544168453406223d4,&
             -0.221166468149940465d4, 0.165443726445793004d4, -0.596412188910614259d4, &
             0.654481694231538040d4, 0.261413416681972012d5, -0.349701859112702878d5,&
             -0.328185277155018630d5,0.790208849885562522d5, -0.398783520249289213d5 /)

     case (3:4) !Potassium
        N = 21
        allocate(A(N))
        RSR = 4.750d0
        NS = 6d0
        U1 = -0.6948000684d3
        U2 = 0.7986755824d7
        B = -0.300d0
        RLR = 12.00d0
        RM = 5.73392370d0
        C6 = 0.1892652670d8
        C8 =  0.5706799527d9
        C10 = 0.1853042722d11
        C26 = 0.0d0
        Aex =0.90092159d4
        gamma = 5.19500d0
        beta = 2.13539d0
        A = (/-255.015289d0, -0.84057856111142d0, 0.20960112217307d4, -0.17090298954603d4, &
             -0.17873773359495d4, 0.29451253739583d4, -0.20200089247397d5, &
             -0.35699524005434d5, 0.59869055371895d6, -0.71054353363636d6,&
             -0.61711841390175d7,0.19365507566961d8, 0.67930587059121d7, &
             -0.12020061704172d9, 0.21603959986951d9,-0.63531969223760d8,-0.52391212820709d9, &
             0.15913304648629d10,-0.24792546567713d10,0.20326031881106d10,-0.68044508325774d9/)
        nu0 = 0.23803737d0
        nu1 = 0.0d0

     case (5) !Sodium
        N = 9
        allocate(A(N))
        RSR = 4.2780d0
        U1 = -0.2435819d3
        U2 = 0.1488425d7
        RM = 5.149085d0
        NS = 6
        RLR = 11.00d0
        B = -0.40d0
        nu0 = 0d0
        nu1 = 0d0
        C6 = 0.75186131d7
        C8 = 0.1686430d9
        C10 = 0.3081961d10
        C26 = 0d0
        Aex = 0.40485835d5
        gamma = 4.59105d0
        beta = 2.36594d0
        A = (/-172.90517d0, 0.355691862122135882d1, 0.910756126766199941d3, &
             -0.460619207631179620d3, 0.910227086296958532d3, -0.296064051187991117d4, &
             -0.496106499110302684d4, 0.147539144920038962d5, -0.819923776793683828d4/)

     case (6) !Lithium-6
        IMN1 = 6
        IMN2 = IMN1

        C6 = 6.7185d6
        C8 = 1.12629d8
        C10 = 2.78683d9
        call POTGENLI2(2,IMN1,IMN2,NPP,VLIM,XO,RM2,VV)

        ! The Le Roy potentials generated by the subroutine above do not produce the correct scattering lengths
        ! We therefore add a quadratic term a*(r-re)**2 for all r inside the minimum.
        re = 4.17005000d0
        Scorr =  1.5432d-6 * HartreePerInvcm/(BohrPerAngstrom**2)
        do i=1,NPP
           if(XO(i).lt.re) then
              VV(i) = VV(i) + Scorr*(XO(i) - re)**2
           endif
        enddo

     case (7) !Lithium-7
        IMN1 = 7
        IMN2 = IMN1

        call POTGENLI2(2,IMN1,IMN2,NPP,VLIM,XO,RM2,VV)
        
        C6 = 6.7185d6
        C8 = 1.12629d8
        C10 = 2.78683d9
        
        re = 4.17005000d0
        Scorr =  2.805d-7 * HartreePerInvcm/(BohrPerAngstrom**2)
        do i=1,NPP
           if(XO(i).lt.re) then
              VV(i) = VV(i) + Scorr*(XO(i) - re)**2
           endif
        enddo

     case (8) !Cesium 
        N = 5
        allocate(A(N))
        p = 5d0
        q = p
        phiinf = -3.971077022867d-1
        re = 6.226182057299d0
        De = 279.2222367314d0
        C6 = 3.3164d7
        C8 = 1.38d9
        C10 = 6.01d10
        ref = 8.7005d0
        A = (/-4.324429443667d-1, -9.206933982533d-2, -6.845846740405d-2, &
             -1.308218973148d-2, 3.457944786933d-1/)

     case default
        write(6,*) "Invalid ISTATE value"
        stop 
     end select

  else
     write(6,*) "Invalid ESTATE value"
  end if

  if((ISTATE.NE.6).AND.(ISTATE.NE.7))then  !Exclude Li-6 and Li-7 since those are generated by Le Roy routines
     if (ISTATE.EQ.8) then  ! Cesium only
        s = -1d0
        bds = 3.30d0
        cds = 0.423d0
        rhoAB = 0.434d0

        DDS6re = (1 - Exp(-(rhoAB*re)*((bds/6d0) + (cds*rhoAB*re)/Sqrt(6d0))))**(6d0+s)
        DDS8re = (1 - Exp(-(rhoAB*re)*((bds/8d0) + (cds*rhoAB*re)/Sqrt(8d0))))**(8d0+s)
        DDS10re = (1 - Exp(-(rhoAB*re)*((bds/10d0) + (cds*rhoAB*re)/Sqrt(10d0))))**(10d0+s)

        uLRre = DDS6re*(C6/(re**6.0d0)) + DDS8re*(C8/(re**8.0d0)) + DDS10re*(C10/(re**10.0d0))

        do i = 1, NPP
           DDS6 = (1 - Exp(-(rhoAB*XO(i))*(bds/6d0 + (cds*rhoAB*XO(i))/Sqrt(6d0))))**(6d0+s)
           DDS8 = (1 - Exp(-(rhoAB*XO(i))*(bds/8d0 + (cds*rhoAB*XO(i))/Sqrt(8d0))))**(8d0+s)
           DDS10 = (1 - Exp(-(rhoAB*XO(i))*(bds/10d0 + (cds*rhoAB*XO(i))/Sqrt(10d0))))**(10d0+s)
           uLR = DDS6*(C6/(XO(i)**6.0d0)) + DDS8*(C8/(XO(i)**8.0d0)) + DDS10*(C10/(XO(i)**10.0d0))
           ypref = (XO(i)**p - ref**p)/(XO(i)**p + ref**p)
           ypeq = (XO(i)**p - re**p)/(XO(i)**p + re**p)
           yqref = (XO(i)**q - ref**q)/(XO(i)**q + ref**q)
           phiMLR = 0d0
           phiMLRtemp = 0d0

           do iphi = 0, N-1
              phiMlRtemp = (1 - ypref)*(A(iphi+1))*yqref**iphi
              phiMLR = phiMLR + phiMLRtemp
           enddo

           phiMLR = phiMLR + ypref*phiinf

           VV(i) = De*(1-(uLR/uLRre)*Exp(-phiMLR*ypeq))**2 - De

        enddo

     else  ! All others
        do i = 1, NPP
           if(XO(i) .LE. RSR) then
              VV(i) = U1 + U2/(XO(i)**NS)
           else if ((XO(i) .GE. RSR) .AND. (XO(i) .LE. RLR)) then
              xi = (XO(i) - RM)/(XO(i) + B*RM)
              VV(i) = 0.0d0

              do j = 0, N-1
                 VV(i) = VV(i) + A(j+1)*xi**j
              enddo
           else
              VV(i) = -C6/(XO(i)**6.0d0) - C8/(XO(i)**8.0d0) - C10/(XO(i)**10.0d0)  - C26/(XO(i)**26.0d0)

              if(ESTATE.EQ.1) then
                 VV(i) = VV(i) - Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))
              elseif(ESTATE.EQ.3) then
                 VV(i) = VV(i) + Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))
              else
                 Write(6,*) "Invalid ESTATE value"
              endif

           endif

           VV(i) = VV(i) + (nu0 + nu1*((XO(i) - RM)/(XO(i) + B*RM)))*(1 - MU/MUREF)*((2*RM)/(XO(i)+RM))**6
        enddo

     endif

  endif
  Cvals(1) = C6 * AngstromPerBohr**6 * InvcmPerHartree
  Cvals(2) = C8 * AngstromPerBohr**8 * InvcmPerHartree
  Cvals(3) = C10 * AngstromPerBohr**10 * InvcmPerHartree
  Cvals(4) = C26 * AngstromPerBohr**26 * InvcmPerHartree
  VV = VV*InvcmPerHartree

END SUBROUTINE SetupPotential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SplitGridMaker(grid,numpts1,numpts2,E1,EM,E2)
  implicit none
  DOUBLE PRECISION grid(numpts1+numpts2)
  DOUBLE PRECISION E1,E2,EM,LE1,LE2,DE,LDE
  INTEGER numpts1,numpts2, iE
!  CHARACTER(LEN=*), INTENT(IN) :: scale

  if((numpts1+numpts2).ne.(size(grid))) then
     write(6,*) "SplitGridMaker must have numpts1+numpts2 = numpts. Stopping."
     stop
  endif
  !--------------------------------------------
  ! Linear grid:
  !--------------------------------------------
  grid(1)=E1
!  IF((scale.EQ."linear").and.(numpts.gt.1)) THEN
  DE=(EM-E1)/DBLE(numpts1-1)
  DO iE=1,numpts1
     grid(iE) = E1 + (iE-1)*DE
     write(6,*) iE, grid(iE)
  ENDDO
  DE=(E2-EM)/DBLE(numpts2)
  DO iE=1,numpts2
     grid(numpts1+iE) = EM + (iE)*DE
     write(6,*) numpts1+iE, grid(numpts1+iE)
  ENDDO
     
!  ENDIF

END SUBROUTINE SplitGridMaker

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GridMaker(grid,numpts,E1,E2,scale)
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
  ! NegLog grid:  (use this for a log spacing of negative numbers)
  !--------------------------------------------
  IF((scale.EQ."neglog").and.(numpts.gt.1)) THEN
     LE1=dlog(-E2)
     LE2=dlog(-E1)

     LDE=(LE2-LE1)/DBLE(numpts-1d0)
     DO iE=1,numpts
        grid(iE) = -dexp(LE2 - (iE-1)*LDE)
!        write(6,*) iE, grid(iE)
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
  !--------------------------------------------
  ! Cubic grid:
  !--------------------------------------------
  IF((scale.EQ."cubic").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**3*DE
     ENDDO
  ENDIF
  !--------------------------------------------
  ! quartic grid:
  !--------------------------------------------
  IF((scale.EQ."quartic").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**4*DE
     ENDDO
  ENDIF
  !--------------------------------------------
  ! sqrroot grid:
  !--------------------------------------------
  IF((scale.EQ."sqrroot").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**0.5d0*DE
     ENDDO
  ENDIF
  !--------------------------------------------
  ! cuberoot grid:
  !--------------------------------------------
  IF((scale.EQ."cuberoot").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**0.33333333333333333333333333333d0*DE
     ENDDO
  ENDIF

END SUBROUTINE GridMaker
