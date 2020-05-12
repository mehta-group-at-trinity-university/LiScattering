program main
  use units
  implicit none
  integer ISTATE, IMN1, IMN2, NPP, i
  double precision, allocatable :: XO(:), VSinglet(:), VTriplet(:), RM2(:)
  double precision VLIM,xmin,xmax,dx

  NPP = 100
  VLIM = 0d0


  allocate(XO(NPP),VSinglet(NPP),VTriplet(NPP),RM2(NPP))
  IMN1 = 6
  IMN2 = 6
  xmin = 1.0d0
  xmax = 15.0d0
  dx = (xmax-xmin)/(dble(NPP-1))
  do i=1, NPP
     XO(i) = xmin + (i-1)*dx
  enddo
  ISTATE = 1                ! Find the "singlet" X(^1\Sigma_g^+) curve
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VSinglet)

  ! print singlet potential to fort.10
  do i=1, NPP
     write(10,*) XO(i)*AngstrominBohr, VSinglet(i)*Hartreeininvcm
  enddo

  ISTATE = 2                !Find the triplet potential
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VTriplet)
  do i=1, NPP
     write(30,*) XO(i)*AngstrominBohr, VTriplet(i)*Hartreeininvcm
  enddo

end program 
