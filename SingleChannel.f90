program SingleChannel
  use EREdata
  implicit none
  double precision, external :: VSech, Vsqrc6,phirecon
  double precision, allocatable :: Energies(:)
  integer xDimMin,iEnergy,NumEnergySteps,iE,iR
  double precision mass,x,EnergyMin,EnergyMax,DeltaEnergy,deltatemp,lwave,norm
  character*64 LegendreFile
  Format1="(F18.15,F18.15,F18.15,F18.15)"

  ! read in legendre data
  read(5,*)
  read(5,1002) LegendreFile
  write(6,1002) LegendreFile
  read(5,*)
  read(5,*)
  read(5,*) LegPoints
  write(6,*) LegPoints,' LegPoints'

  !     read in boundary conditions
  read(5,*)
  read(5,*)
  read(5,*) Order,Left,Right
  print*,  'Order, Left, Right'
  print*, Order,Left,Right

  !     read in potential parameters
  read(5,*)
  read(5,*)
  read(5,*) alpha, mass, DD, x0
  write(6,*) 'alpha,            mass,              DD,               x0'
  write(6,Format1) alpha, mass, DD, x0

  mu=mass/2.0d0

  read(5,*)
  read(5,*)
  read(5,*) xNumPoints,xMin,xMax
  write(6,*) 'xNumPoints,xMin,xMax'
  write(6,*) xNumPoints,xMin,xMax

  read(5,*)
  read(5,*)
  read(5,*) lwave, EnergyMin, EnergyMax, NumEnergySteps
  write(6,*) 'lwave, EnergyMin, EnergyMax, NumEnergySteps'
  write(6,*) lwave, EnergyMin, EnergyMax, NumEnergySteps

  xDimMin=xNumPoints+Order-3
  xDim=xDimMin
  if (Left .eq. 2) xDim = xDim + 1
  if (Right .eq. 2) xDim = xDim + 1
  MatrixDim=xDim
  print*, 'xDim=',xDim
  allocate(xLeg(LegPoints),wLeg(LegPoints))
  call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

  allocate(Energies(NumEnergySteps))
  DeltaEnergy=(EnergyMax-EnergyMin)/dble(NumEnergySteps-1)
  do iE=1,NumEnergySteps
     Energies(iE)=EnergyMin+(iE-1)*DeltaEnergy
  end do
  allocate(delta(NumEnergySteps))

  allocate(xPoints(xNumPoints))
  allocate(xBounds(xNumPoints+2*Order))
  allocate(u(LegPoints,xNumPoints,xDim),ux(LegPoints,xNumPoints,xDim)) ! 
  allocate(L(MatrixDim,MatrixDim))
  allocate(S(MatrixDim,MatrixDim))
  allocate(G0(MatrixDim,MatrixDim))
  allocate(V0(MatrixDim,MatrixDim))

  PotInterpPoints = 100
  PotInterpOrder = 3
  allocate(R(PotInterpPoints),VS(PotInterpPoints),VT(PotInterpPoints))
  do iR=1,PotInterpPoints
     R(iR) = xmin + dble((iR-1)/(PotInterpPoints-1))**2d0 * xMax
     write(6,*) iR, R(iR)
  enddo
  allocate(Potknots(PotInterpPoints+PotInterpOrder),Vbcoef(PotInterpPoints))
  call GetSTPots(R,VS,VT,PotInterpPoints) ! use atomic units for R, VS, and VT
  call setup_Vinterp(R,PotInterpOrder,Potknots,PotInterpPoints,VData,Vbcoef)
!!$  do x=xMin,xMax,(xMax-xMin)/1000.0d0
!!$     !write(1000,*) x, VSech(DD,x0,x,lwave) 
!!$     write(1000,*) x, Vsqrc6(DD,x0,x,lwave) 
!!$  end do

  call GridMaker(xNumPoints,xMin,xMax,xPoints)

  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,0,u)

  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg, &
       xDim,xBounds,xNumPoints,1,ux)

  call CheckBasis(u,xDim,xNumPoints,LegPoints,xLeg,xPoints,700)

  G0=0.0d0
 
  call CalcGamma0(lwave)
  do iE=1,NumEnergySteps
     call CalcPhaseShift(Energies(iE),lwave,deltatemp,norm)
     delta(iE)=deltatemp
     write(6,*) Energies(iE), delta(iE)
     call PlotCheckBessel(0.0d0, 0.01d0, xMax, Energies(iE), mu, dtan(delta(iE)), 10)
  enddo

  deallocate(G0,V0,S,L,delta)
  deallocate(xBounds,xPoints)
  deallocate(u,ux,xLeg,wLeg)

10 format(1P,100e25.15)
20 format(1P,100e16.8)
1002 format(a64)

end program SingleChannel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetSTPots(R,VS,VT,NPP)
  use units
  implicit none
  integer NPP,IMN1,IMN2,sym,lwave,ISTATE
  double precision R(NPP), VS(NPP), VT(NPP),XO(NPP)
  double precision weights(NPP),RM2(NPP)
  double precision VLIM
  XO = R*BohrPerAngstrom !Convert to Angstrom to send to the Le Roy code
  VLIM = 0d0

  IMN1 = 7
  IMN2 = 7
  sym = 1 ! set to +1 for bosonic Li-7, -1 for fermionic Li-6, and 0 for mixture Li-6 - Li-7.
  lwave = 0 ! s wave collisions

  ISTATE = 1                ! Find the "singlet" X(^1\Sigma_g^+) curve
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VS)

  ! print singlet potential to fort.10
  do iR=1, NPP
     VS(iR) = VS(iR)*InvcmPerHartree
     write(10,*) R(iR), VS(iR)
  enddo

  ISTATE = 2                !Find the triplet potential
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VT)
  do iR=1, NPP
     VT(iR) = VT(iR)*InvcmPerHartree
     write(30,*) R(iR), VT(iR)
  enddo

end subroutine GetSTPots

subroutine CalcPhaseShift(energy, lwave, deltatemp,norm)
  use EREdata
  implicit none
  double precision, external :: phirecon
  double precision, allocatable :: G(:,:), LL(:,:), ksqr(:), kcotandel(:),evec(:,:),eval(:),sig(:) ! 
  double precision energy, b, k, tandel, x, deltatemp,sj,sy,sjp,syp,lwave,norm!, ainv, halfr0, q, chi2,siga, sigb,a2inv ! 

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
  LL(xDim,xDim)=1.0d0
  call Mydggev(MatrixDim,G,MatrixDim,LL,MatrixDim,eval,evec)
     
  do n=1,MatrixDim
     if(dabs(eval(n)).gt.1d-14) then
        b=eval(n)
        ikeep=n
        write(6,*) 'eval(',n,') = ', b
     endif
  enddo


  k=dsqrt(2.0d0*mu*energy)
  x=xMax
  call sphbes(lwave,k*x,sj,sy,sjp,syp)
  tandel = ( (b*x + 1.d0) * sj + k*x*sjp)/( (b*x + 1.d0) * sy + k*x*syp)  
!  tandel=(k*dcos(k*xMax) + b*dsin(k*xMax))/(-b*dcos(k*xMax) + k*dsin(k*xMax))
  deltatemp=datan(tandel)
  norm = 2.d0*k*x*(sj-tandel*sy)/phirecon(xMax,ikeep,evec,Left,Right,xDim,MatrixDim,xNumPoints,xPoints,order)

!  write(6,*) 'b = ',b, 'freesol_b = ',-k*dcos(k*xMax)/dsin(k*xMax)
  write(6,*) 'b = ',b, 'freesol_b = ',-k*sjp/sj
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
subroutine CalcGamma0(lwave)
  use EREdata
  implicit none
  double precision, external :: VSech,vsqrc6
  double precision ax, bx,x,xScaledZero,xIntScale,TempG,TempS,TempV,a,lwave
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
              !TempV = TempV + a*u(lx,kx,ix)*(alpha*VSech(DD,x0,x,lwave))*u(lx,kx,ixp)
              !TempV = TempV + a*u(lx,kx,ix)*(alpha*Vsqrc6(DD,x0,x,lwave))*u(lx,kx,ixp)
              TempV = TempV + a*u(lx,kx,ix)* &
                   (alpha*dbsval(x,PotInterpOrder,Potknots,PotInterpPoints,Vbcoef))* &
                   u(lx,kx,ixp)
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
subroutine GridMaker(xNumPoints,xMin,xMax,xPoints)
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
end subroutine GridMaker

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
  double precision k,dr,sj,sy,sjp,syp,tandel,cosdel,sindel,delta,rMin,rMax,energy,mu,lwave
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
     call sphbes(lwave,k*r(iStep),sj,sy,sjp,syp)
!     write(file,*) r(iStep), k*r(iStep)*sj, k*r(iStep)*sy, 2.d0*k*r(iStep)*(sj-tandel*sy)
     write(file,*) r(iStep),  2.d0*k*r(iStep)*(sj-tandel*sy)
  enddo
  write(file, *) ' ' 

end subroutine PlotCheckBessel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE sphbes(xnu,x,sj,sy,sjp,syp)
      use nrtype, only : dbl
!      INTEGER n
!      double precision sj,sjp,sy,syp,x,xnu
      double precision sj,sjp,sy,syp,x,xnu
      !    USES bessjy
      double precision factor,order,rj,rjp,ry,ryp,RTPIO2
      PARAMETER (RTPIO2=1.25331413731550d0)
      if(n.lt.0.d0.or.x.le.0.d0) then
         write(6,*) 'Bad arguments in sphbes. Press ENTER to continue.'
         read(5,*)
      endif
      order=xnu+0.5d0
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/Dsqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.d0*x)
      syp=factor*ryp-sy/(2.d0*x)

      return
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      use nrtype, only : dbl
      INTEGER MAXIT
      double precision rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0, &
           PI=3.141592653589793d0)
!    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e, &
           f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q, &
           r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1, &
           temp,w,x2,xi,xi2,xmu,xmu2
      !write(66,*)xnu  
      if(x.le.0.d0.or.xnu.lt.0.d0) then
         write(6,*) 'Bad arguments in bessjy. Press ENTER to continue.'
         read(5,*)
      endif

      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
!      pause 'x too large in bessjy; try asymptotic expansion'
      goto 31
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
          write(6,*) 'bessy series failed to converge'
          read(5,*)
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
          write(6,*) 'cf2 failed in bessjy. Press enter to continue.'
          read(5,*)
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
31    Continue
!c  asymptotic expansion
      rj = Dsqrt(2.d0/(pi*x))*Dcos(x-xnu*pi/2.d0-pi/4.d0)
      ry = Dsqrt(2.d0/(pi*x))*Dsin(x-xnu*pi/2.d0-pi/4.d0)
      rjp = -ry
      ryp = rj
      return
      END


      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      use nrtype, only: dbl  
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
!CU    USES chebev
      double precision xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4, &
           -3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3, &
           -4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
    END SUBROUTINE beschb

    
    FUNCTION chebev(a,b,c,m,x)
      use nrtype, only : dbl
      INTEGER m
      double precision chebev,a,b,x,c(m)
      INTEGER j
      double precision d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.d0) then
         write(6,*) 'x not in range in chebev.  Press enter to continue'
         read(5,*)
      endif
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END         

      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      INTEGER MAXIT
      DOUBLE PRECISION ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-10,FPMIN=1.e-30,MAXIT=10000,XMIN=2., &
           PI=3.141592653589793d0)
!CU    USES beschb
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff, &
           gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1, &
           ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) then
         write(6,*) 'bad arguments in bessik. Press enter to continue'
         read(5,*)
      endif
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=1.d0/(b+d)
        c=b+1.d0/c
        del=c*d
        h=del*h
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
        write(6,*)'x too large in bessik; try asymptotic expansion.  Press enter to continue'
        read(5,*)
        
1     continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
        ritemp=fact*ril+ripl
        fact=fact-xi
        ripl=fact*ritemp+ril
        ril=ritemp
12    continue
      f=ripl/ril
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS)goto 2
13      continue
          write(6,*) 'bessk series failed to converge.  Press enter to continue.'
          read(5,*)
2       continue
        rkmu=sum
        rk1=sum1*xi2
      else
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
        do 14 i=2,MAXIT
          a=a-2*(i-1)
          c=-a*c/i
          qnew=(q1-b*q2)/a
          q1=q2
          q2=qnew
          q=q+c*qnew
          b=b+2.d0
          d=1.d0/(b+a*d)
          delh=(b*d-1.d0)*delh
          h=h+delh
          dels=q*delh
          s=s+dels
          if(abs(dels/s).lt.EPS)goto 3
14      continue
          write(6,*) 'bessik: failure to converge in cf2.  Press enter to continue.'
          read(5,*)
3       continue
        h=a1*h
        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
      do 15 i=1,nl
        rktemp=(xmu+i)*xi2*rk1+rkmu
        rkmu=rk1
        rk1=rktemp
15    continue
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
