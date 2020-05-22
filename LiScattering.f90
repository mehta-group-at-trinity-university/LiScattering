module hyperfine
  use units
  implicit none
  double precision, parameter :: muB = 1.3996244936142    !Bohr magneton in units (MHz/Gauss)
  double precision, parameter :: muN = 7.62259328547d-4    !Nuclear Magneton in units (MHz/Gauss)
  integer, parameter :: sLi7 = 1 ! twice the electronic spin of Li-7
  integer, parameter :: iLi7 = 3 ! twice the nuclear spin of Li-7
  integer, parameter :: sLi6 = 1 ! twice the electronic spin of Li-6
  integer, parameter :: iLi6 = 2 ! twice the nuclear spin of Li-6
  double precision, parameter :: gs = 2.00231930436256d0 ! electron g factor
  double precision, parameter :: giLi7 = 2.170903 ! Nuclear g factor for LI-7
  double precision, parameter :: giLi6 = 0.822019 ! Nuclear g factor for LI-7
  
  !  Hyperfine from Arimondo et al Rev. Mod. Phys. 49, 31 (1977).  See table on pg. 67 
  double precision, parameter :: ALi7 = 401.7520435 ! MHz  
  double precision, parameter :: ALi6 = 152.136840720 ! MHz
  
  type hf1atom
     integer f, m
  end type hf1atom
  type hf2atom
     type(hf1atom) a1,a2
  end type hf2atom
  type symhf2atom
     type(hf2atom) state1, state2
     double precision norm, phase
  end type symhf2atom

  type(hf2atom), allocatable :: hf2TempGlobal(:)
  type(symhf2atom), allocatable :: hf2symTempGlobal(:)

  
contains
  !See Appendix B of Ian McAlexander's thesis (Rice University)
  ! All angular momentum quantum numbers are TWICE their actual
  ! values so that all arithemetic can be done with integers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MakeHF1Basis(i,s,size,hf1)
    implicit none
    ! Gives the size of the (unsymmetrized) hyperfine basis (dimension of the 1-atom hyperfine + zeeman Hamiltonian) for one atom
    integer i, s, size ! size is (in/out) if size = 0 on input -- determine the size.  Else, use it to calculate the basis stored in hf1
    integer f, m, count
    type(hf1atom) hf1(:)

    if(size.eq.0) then
       do f = iabs(i-s), i+s, 2
          do m = -f, f, 2
             size = size + 1
          enddo
       enddo
    else
       count = 0
       do f = iabs(i-s), i+s, 2
          do m = -f, f, 2
             count = count+1
             hf1(count)%f = f
             hf1(count)%m = m
          enddo
       enddo
    endif
    !    write(6,*) "size = ", size
  end subroutine MakeHF1Basis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function hf2eq(state1,state2)
    !determines if two-atom hyperfine states are equal
    implicit none
    logical hf2eq
    type(hf2atom) state1, state2

    hf2eq = .false.
    if(state1%a1%f.eq.state2%a1%f) then
       if(state1%a1%m.eq.state2%a1%m) then
          if(state1%a2%f.eq.state2%a2%f) then
             if(state1%a2%m.eq.state2%a2%m) then
                hf2eq = .true.
             endif
          endif
       endif
    endif

    return
  end function hf2eq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setzero2atom(state)
    implicit none
    type(hf2atom) state
    state%a1%f=0
    state%a1%m=0
    state%a2%f=0
    state%a2%m=0
  end subroutine setzero2atom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MakeHF2Basis(i1, s1, i2, s2, sym, lwave, mtot, size)
    implicit none
    ! Gives the size of the (unsymmetrized) hyperfine basis for 2 atoms with total M_F = mtot
    ! lwave denotes the partial wave and sym < 0 for fermions, sym>0 for bosons and sym=0 for distinguishable.
    integer i1, s1, i2, s2, size, count, keepflag, lwave, sym
    integer f1, m1, f2, m2, mtot, symsize, n1, n2
    type(hf2atom) zerostate
    type(hf2atom), allocatable :: tempstates(:)
    type(symhf2atom), allocatable :: symstates(:)
    double precision kd1, kd2
    double precision, external :: dkdelta
    !    zerostate=0
    call setzero2atom(zerostate)
    write(6,*) "zerostate = ",zerostate
    
    size = 0
    do f1 = iabs(i1-s1), i1+s1, 2
       do m1 = -f1, f1, 2
          do f2 = iabs(i2-s2), i2+s2, 2
             do m2 = -f2, f2, 2
                if((m1+m2).eq.mtot) then
                   size = size + 1
                endif
             enddo
          enddo
       enddo
    enddo
    
    allocate(hf2TempGlobal(size))
    do n1=1,size
       call setzero2atom(hf2TempGlobal(n1))
    enddo
    
    write(6,'(A,I3)') "Generating the unsymmetrized basis with total MF = ",mtot
    write(6,'(5A5)') "#","f1","m1","f2","m2"
    write(6,'(5A5)') "---","---","---","---","---"
    
    count=0
    do f1 = iabs(i1-s1), i1+s1, 2
       do m1 = -f1, f1, 2
          do f2 = iabs(i2-s2), i2+s2, 2
             do m2 = -f2, f2, 2
                if((m1+m2).eq.mtot) then
                   count = count + 1
                   hf2TempGlobal(count)%a1%f=f1
                   hf2TempGlobal(count)%a1%m=m1
                   hf2TempGlobal(count)%a2%f=f2
                   hf2TempGlobal(count)%a2%m=m2
                   write(6,'(5I5)') count, hf2TempGlobal(count)
                endif
             enddo
          enddo
       enddo
    enddo

    allocate(tempstates(size),symstates(size))
    count = 0
    do n1 = 1, size
       !construct a temporary state with the q-#s of atom 1 and 2 interchanged
       tempstates(n1)%a1 = hf2TempGlobal(n1)%a2
       tempstates(n1)%a2 = hf2TempGlobal(n1)%a1
       !write(6,'(4I3,A,4I3)') hf2TempGlobal(n1),"  ||", tempstates(n1)
       ! keep a tally of unique states
       !consider below the case of bosonic atoms with l=0 only.  Check if sum of two states is unique
       keepflag = 1
       if((sym.lt.0).and.hf2eq(tempstates(n1),hf2TempGlobal(n1))) keepflag = 0
       do n2 = 1, n1-1
          if(hf2eq(tempstates(n1),hf2TempGlobal(n2))) keepflag = 0
       enddo
       if(keepflag.eq.1) then
          count=count+1
          symstates(count)%state1 = hf2TempGlobal(n1)
          symstates(count)%state2 = tempstates(n1)
       endif
    enddo
    symsize = count
!    count = 0
    allocate(hf2symTempGlobal(symsize))
    write(6,'(A,I2)') "(+1/0/-1) = (boson/dist/fermion) Symmetry case: ", sym
    write(6,'(A,I2)') "partial wave = ", lwave
    do n1 = 1, symsize
       hf2symTempGlobal(n1)=symstates(n1)
       ! record the normalization and relative phase into the state
       kd1 = dkdelta(hf2symTempGlobal(n1)%state1%a1%f,hf2symTempGlobal(n1)%state1%a2%f)
       kd2 = dkdelta(hf2symTempGlobal(n1)%state1%a1%m,hf2symTempGlobal(n1)%state1%a2%m)
       
       hf2symTempGlobal(n1)%norm = 1d0/dsqrt(2d0*(1d0+kd1*kd2))
       hf2symTempGlobal(n1)%phase = 1d0
       if(sym.lt.0) hf2symTempGlobal(n1)%phase = -1d0
       hf2symTempGlobal(n1)%phase = hf2symTempGlobal(n1)%phase*(-1)**lwave
       
       write(6,'(I4,A,f5.2,A,4I2,A,f5.2,A,4I2,A)') n1,"  symmetrized basis:", hf2symTempGlobal(n1)%norm, &
            " ( |", hf2symTempGlobal(n1)%state1,"> +", hf2symTempGlobal(n1)%phase,"|", hf2symTempGlobal(n1)%state2,"> )"
    enddo
    size=symsize
    deallocate(symstates,tempstates)
    
  end subroutine MakeHF2Basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MakeHHZ1(HHZ,B,size,hf1,gs,gi,Ahf,s,i)
    
     implicit none
     integer f, m, fp, mp, i, s, bra, ket, ms, mi, size
     double precision, external :: THRJ
     double precision, intent(in) :: B
     double precision HHZ(size,size)
     double precision Zsum1, Zsum2, c1, c2, tj1, tj2, gs, gi, Ahf
     type(hf1atom) hf1(size)


     ! Calcualte the Hyperfine/Zeeman matrix in units of MHz with B measured in Gauss
     do bra = 1, size
        do ket = 1, size
           f = hf1(ket)%f
           m = hf1(ket)%m
           fp = hf1(bra)%f
           mp = hf1(bra)%m
           c1 = gs*muB*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))
           c2 = gi*muN*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))
           
           Zsum1 = 0d0
           Zsum2 = 0d0
           
           if((f.eq.fp).and.(m.eq.mp)) then
              HHZ(bra, ket) = HHZ(bra,ket) + &
                   0.5d0*Ahf*( 0.5d0*f*(0.5d0*f + 1d0) - 0.5d0*i*(0.5d0*i + 1d0) - 0.5d0*s*(0.5d0*s + 1d0) )
           endif
           if(m.eq.mp) then
              do ms = -s, s, 2
                 do mi = -i, i, 2
                    tj1 = THRJ(s,i,fp,ms,mi,-mp)
                    tj2 = THRJ(s,i,f,ms,mi,-m)
                    Zsum1 = Zsum1 + 0.5d0*ms*tj1*tj2
                    Zsum2 = Zsum2 + 0.5d0*mi*tj1*tj2
                 enddo
              enddo
              HHZ(bra,ket) = HHZ(bra, ket) + c1*Zsum1 + c2*Zsum2
           endif
        enddo
     enddo
     
   end subroutine MakeHHZ1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine OverlapSpinHF(i1,i2,s1,s2,S,mS,I,mI,f1,m1,f2,m2,me1)
     implicit none
     ! This subroutine calcualtes the overlap between the hyperfine state |f1 m1 f2 m2> and the
     ! spin state |(s1 s2) S mS (i1 i2) I mI >
     ! it returns in me1 = <(s1 s2) S mS (i1 i2) I mI | (i1 s1) f1 m1 (i2 s2) f2 m2 >
     ! I haven't checked this, but it should be equal to the 9-J coefficient corresponding to the above overal.
     integer S,mS,I,mI,f1,m1,f2,m2,i1,i2,s1,s2
     integer mi1,mi2,ms1,ms2
     double precision phase,prefact,tj1,tj2,tj3,tj4,me1
     double precision, external :: THRJ
     
     phase = (-1)**(2d0*i2-2d0*s1-m1-m2-mS-mI)/2d0
     prefact = dble((f1+1)*(f2+1)*(S+1)*(I+1))
     prefact = dsqrt(prefact)
     
     me1 = 0d0
     do ms1 = -s1, s1, 2
        do ms2 = -s2, s2, 2
           do mi1 = -i1, i1, 2
              do mi2 = -i2, i2, 2
                 tj1=THRJ(s1,i1,f1,ms1,mi1,-m1)
                 tj2=THRJ(s2,i2,f2,ms2,mi2,-m2)
                 tj3=THRJ(s1,s2,S,ms1,ms2,-mS)
                 tj4=THRJ(i1,i2,I,mi1,mi2,-mI)
                 me1 = me1 + tj1*tj2*tj3*tj4
              enddo
           enddo
        enddo
     enddo
     
     me1 = me1*phase*prefact
  end subroutine OverlapSpinHF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine STproj(S,i1,i2,s1,s2,hf2bra,hf2ket,res)
    implicit none
    ! This subroutine calculates the unsymmetrized matrix element for the singlet/triplet projection operator
    ! in the hyperfine basis |f1 m1 f2 m2>
    ! i1 = nuclear spin of atom 1, mi1 = z-projection of i1
    ! i2 = nuclear spin of atom 2, mi2 = z-projection of i2
    ! I = total nuclear spin
    ! s1 = electronic spin of atom 1
    ! s2 = electronic spin of atom 2
    ! S = total electronic spin
    ! mS, mI = z-projections of total nuclear and electronic spins
    ! hf2bra and hf2ket contain f1p m1p f2p m2p and f1 m1 f2 m2m respectively. 
    ! f = S+I = total spin of atom
    ! S = 0 for singlet, 2 for triplet (since this is really 2S)
    integer S,I,mS,mI,mi1,ms1,mi2,ms2
    integer i1,i2,s1,s2!,f1p,m1p,f2p,m2p,f1,m1,f2,m2
    double precision res,me1,me2
    type(hf2atom) hf2bra, hf2ket

!    write(6,*) "For i1,i2,s1,s2 = ",i1,i2,s1,s2
    res=0d0
    if((S.eq.0).or.(S.eq.2)) then
       if ((S.le.(s1+s2)) .and. (S.ge.iabs(s1-s2))) then
          do I = iabs(i1-i2), i1+i2, 2
             do mS = -S, S, 2
                do mI = -I, I, 2
                   call OverlapSpinHF(i1,i2,s1,s2,S,mS,I,mI,hf2ket%a1%f,hf2ket%a1%m,hf2ket%a2%f,hf2ket%a2%m,me1)
                   call OverlapSpinHF(i1,i2,s1,s2,S,mS,I,mI,hf2bra%a1%f,hf2bra%a1%m,hf2bra%a2%f,hf2bra%a2%m,me2)
                   res = res + me1*me2
                   !write(6,*) I,mI,S,mS, res
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine STproj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MakeSTHFProj(S, i1, i2, s1, s2, hfsym, size, PST)
    implicit none
    integer S, i1, i2, s1, s2
    integer bra, ket, size
    type(symhf2atom), intent(in) :: hfsym(size)
    double precision, intent(out) :: PST(size,size)
    double precision me1, me2, me3, me4

    PST = 0d0
    do bra = 1, size
       do ket = 1, size
          call STproj(S,i1,i2,s1,s2,hfsym(bra)%state1,hfsym(ket)%state1,me1)
          me1 = me1*hfsym(bra)%norm*hfsym(ket)%norm
          call STproj(S,i1,i2,s1,s2,hfsym(bra)%state1,hfsym(ket)%state2,me2)
          me2 = me2*hfsym(bra)%norm*hfsym(ket)%norm
          me2 = me2*hfsym(ket)%phase
          call STproj(S,i1,i2,s1,s2,hfsym(bra)%state2,hfsym(ket)%state1,me3)
          me3 = me3*hfsym(bra)%norm*hfsym(ket)%norm
          me3 = me3*hfsym(bra)%phase
          call STproj(S,i1,i2,s1,s2,hfsym(bra)%state2,hfsym(ket)%state2,me4)
          me4 = me4*hfsym(bra)%norm*hfsym(ket)%norm
          PST(bra, ket) = me1 + me2 + me3 + me4
       enddo
    enddo
    
  end subroutine MakeSTHFProj
end module hyperfine

!****************************************************************************************************
program main
  use units
  use hyperfine
  implicit none
  integer ISTATE, IMN1, IMN2, NPP, iR, iB, ihf
  integer i1,i2,s1,s2,S,sym,size1,size2,NBgrid,mtot, lwave
  double precision, allocatable :: XO(:), VSinglet(:), VTriplet(:), RM2(:)
  double precision, allocatable :: HHZ(:,:), Bgrid(:), EVAL(:), EVEC(:,:)
  double precision, allocatable :: SPmat(:,:), TPmat(:,:), VHZ(:,:,:)
  double precision VLIM,xmin,xmax,dx
  double precision B
  type(hf1atom) a1, a2
  type(hf2atom) mstate1, mstate2
  type(hf1atom), allocatable :: Li6hf(:), Li7hf(:)  
  ! initialize the angular momentum routines
  call setupam
  ! determine the size of the one-atom hyperfine/zeeman hamiltonian

  NBgrid = 40
  !make the magnetic field grid
  allocate(Bgrid(NBgrid))
  call makeEgrid(Bgrid,NBgrid,0.0d0,1000d0)
  !------------------------------------------------------------------
  ! solve the Li-7 1-atom Zeeman problem
  !call once with size1 = 0 to determine the size of the basis.
  size1 = 0
  call MakeHF1Basis(iLi7,sLi7,size1,Li7hf)
  write(6,*) "size of the Li-7 1-atom hyperfine basis = ",size1
  !allocate memory for the basis
  allocate(Li7hf(size1),EVEC(size1,size1),EVAL(size1),HHZ(size1,size1))
  call MakeHF1Basis(iLi7,sLi7,size1,Li7hf)


  ! construct the HZ Hamiltonian
  do iB = 1, NBgrid
     EVEC=0d0
     EVAL=0d0
     HHZ=0d0
     call MakeHHZ1(HHZ,Bgrid(iB),size1,Li7hf,gs,giLi7,ALi7,sLi7,iLi7)
!     call printmatrix(HHZ,size1,size1,6)
     call MyDSYEV(HHZ,size1,EVAL,EVEC)
     write(90,*) Bgrid(iB), EVAL
  enddo
  deallocate(EVEC,EVAL,HHZ)
  !____________________________________________________________________
  !------------------------------------------------------------------
  ! solve the Li-6 1-atom Zeeman problem
  size1=0
  call MakeHF1Basis(iLi6,sLi6,size1,Li6hf)
  write(6,*) "size of the Li-6 1-atom hyperfine basis = ",size1
  allocate(Li6hf(size1))
  call MakeHF1Basis(iLi6,sLi6,size1,Li6hf)
  allocate(EVEC(size1,size1),EVAL(size1),HHZ(size1,size1))
  ! construct the HZ Hamiltonian
  do iB = 1, NBgrid
     EVEC=0d0
     EVAL=0d0
     HHZ=0d0
     call MakeHHZ1(HHZ,Bgrid(iB),size1,Li6hf,gs,giLi6,ALi6,sLi6,iLi6)
!     call printmatrix(HHZ,size1,size1,6)
     call MyDSYEV(HHZ,size1,EVAL,EVEC)
     write(91,*) Bgrid(iB), EVAL
  enddo
  deallocate(EVAL,EVEC,HHZ)
  !____________________________________________________________________
  
  NPP = 1000
  VLIM = 0d0

  allocate(XO(NPP),VSinglet(NPP),VTriplet(NPP),RM2(NPP))

  IMN1 = 7
  IMN2 = 7
  sym = 1 ! set to +1 for bosonic Li-7, -1 for fermionic Li-6, and 0 for mixture Li-6 - Li-7.
  lwave = 0 ! s wave collisions
  xmin = 1.0d0
  xmax = 30.0d0
  dx = (xmax-xmin)/(dble(NPP-1))
  do iR=1, NPP
     XO(iR) = xmin + (iR-1)*dx
  enddo
  ISTATE = 1                ! Find the "singlet" X(^1\Sigma_g^+) curve
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VSinglet)

  ! print singlet potential to fort.10
  do iR=1, NPP
     write(10,*) XO(iR)*AngstromPerBohr, VSinglet(iR)*InvcmPerHartree*HartreePerTHz
  enddo

  ISTATE = 2                !Find the triplet potential
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VTriplet)
  do iR=1, NPP
     write(30,*) XO(iR)*AngstromPerBohr, VTriplet(iR)*InvcmPerHartree*HartreePerTHz
  enddo
  mtot = 4
  call MakeHF2Basis(iLi7, sLi7, iLi7, sLi7, sym, lwave, mtot, size2)
  write(6,*) "size of the symmetrized 2-atom hyperfine basis = ", size2
  allocate(SPmat(size2,size2), TPmat(size2,size2))

  S = 0
  call MakeSTHFProj(S, iLi7, iLi7, sLi7, sLi7, hf2symTempGlobal, size2, SPmat)
  write(6,*) "The singlet projection matrix:"
  write(6,*) "-------------------------------"
  call printmatrix(SPmat,size2,size2,6)

  S = 2
  call MakeSTHFProj(S, iLi7, iLi7, sLi7, sLi7, hf2symTempGlobal, size2, TPmat)
  write(6,*) "The triplet projection matrix:"
  write(6,*) "-------------------------------"
  call printmatrix(TPmat,size2,size2,6)

  !allocate(VHZ(NPP,size2,size2))

  
!  deallocate(SPmat,TPmat,VHZ)
  deallocate(Bgrid,XO,VSinglet,VTriplet,RM2)

10 format(100F10.4)
end program 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function dkdelta(a,b)
  integer a, b
  double precision dkdelta
  if(a.eq.b) then
     dkdelta = 1d0
  else
     dkdelta = 0d0
  endif
  return
end function dkdelta
