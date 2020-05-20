module hyperfine
  use units
  implicit none
  double precision, parameter :: muB = 1.3996244936142    !Bohr magneton in units (MHz/Gauss)
  double precision, parameter :: muN = 7.62259328547d-4    !Nuclear Magneton in units (MHz/Gauss)
  integer, parameter :: sLi7 = 1 ! twice the electronic spin of Li-7
  integer, parameter :: iLi7 = 3 ! twice the nuclear spin of Li-7
  double precision, parameter :: gsLi7 = 2.00231930436256d0 ! electron g factor
  double precision, parameter :: giLi7 = 2.170903 ! Nuclear g factor for LI-7
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
  end type symhf2atom
  
  type(hf2atom), allocatable :: hfstates(:)
  type(symhf2atom), allocatable :: symhfbasis(:)

  
contains
  !See Appendix B of Ian McAlexander's thesis (Rice University)
  ! All angular momentum quantum numbers are TWICE their actual
  ! values so that all arithemetic can be done with integers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine HHZSize1(i,s,size)
    ! Gives the size of the (unsymmetrized) hyperfine basis (dimension of the 1-atom hyperfine + zeeman Hamiltonian) for one atom
    integer i, s, size
    integer f, m
    size = 0
    do f = iabs(i-s), i+s, 2
       do m = -f, f, 2
          size = size + 1
       enddo
    enddo
!    write(6,*) "size = ", size
  end subroutine HHZSize1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function hf2eq(state1,state2)
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
  subroutine HHZSize2(i1, s1, i2, s2, mtot, size)
    ! Gives the size of the (unsymmetrized) hyperfine basis for 2 atoms with total M_F = mtot
    integer i1, s1, i2, s2, size, count, keepflag
    integer f1, m1, f2, m2, mtot, symsize, n1, n2
    type(hf2atom) zerostate
    type(hf2atom), allocatable :: tempstates(:)
    type(symhf2atom), allocatable :: symstates(:)
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
    
    allocate(hfstates(size))
    do n1=1,size
       call setzero2atom(hfstates(n1))
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
                   hfstates(count)%a1%f=f1
                   hfstates(count)%a1%m=m1
                   hfstates(count)%a2%f=f2
                   hfstates(count)%a2%m=m2
                   write(6,'(5I5)') count, hfstates(count)
                endif
             enddo
          enddo
       enddo
    enddo

    allocate(tempstates(size),symstates(size))
!    tempstates = 0
!    symstates = 0
    count = 0
    do n1 = 1, size
       !construct a temporary state with the q-#s of atom 1 and 2 interchanged
       tempstates(n1)%a1 = hfstates(n1)%a2
       tempstates(n1)%a2 = hfstates(n1)%a1
       !write(6,'(4I3,A,4I3)') hfstates(n1),"  ||", tempstates(n1)
       ! keep a tally of unique states
       !consider below the case of bosonic atoms with l=0 only.  Check if sum of two states is unique
       keepflag = 1
       do n2 = 1, n1-1
          if(hf2eq(tempstates(n1),hfstates(n2))) keepflag = 0
       enddo
       if(keepflag.eq.1) then
          count=count+1
          symstates(n1)%state1 = hfstates(n1)
          symstates(n1)%state2 = tempstates(n1)
       endif

    enddo
    symsize = count
    count = 0
    allocate(symhfbasis(symsize))
    do n1 = 1, size
       if(hf2eq(symstates(n1)%state1,zerostate).eqv..false.) then
          count = count+1
          symhfbasis(count)=symstates(n1)
          write(6,'(I4,A,A,4I1,A,4I1,A)') count,"  symmetrized basis", &
               "  |", symhfbasis(count)%state1,"> +/- (-1)^l |", symhfbasis(count)%state2,">"
       endif

    enddo

       

    deallocate(symstates,tempstates)
  end subroutine HHZSize2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine makeHFBasisLi7(mstate,fullsize)
    implicit none
    integer i1,i2,s1,s2,f1,f2,m1,m2,count,fullsize
    type(hf2atom) mstate(fullsize)
    
    s1=sLi7
    s2=sLi7
    i1=iLi7
    i2=iLi7
    
    count = 0
    write(6,*) iabs(s1-i1)
    write(6,*) s1+i1

    do f1 = iabs(s1-i1),s1+i1,2
       do m1 = -f1, f1, 2
           do f2 = iabs(s2-i2),s2+i2,2
              do m2 = -f2, f2, 2
                 
                 if((m1+m2).eq.8) then
                    count = count+1
                    write(6,10) count,"f1=",f1,"m1=",m1,"f2=",f2,"m2=",m2
                 endif
              enddo
           enddo
        enddo
     enddo
 10 format(I10,A5,I4,A5,I4,A5,I4,A5,I4)    
   end subroutine makeHFBasisLi7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine makeLi7HHZ(HHZ,B,size)
     implicit none
     integer f, m, fp, mp, i, s, bra, ket, ms, mi, size, fmin, fmax
     double precision, external :: THRJ
     double precision, intent(in) :: B
     double precision HHZ(size,size)
     double precision Zsum1, Zsum2, c1, c2, tj1, tj2

     i = iLi7
     s = sLi7
     bra = 0
     ket = 0
     fmin = iabs(i-s)
     fmax = i+s
     write(6,*) "fmin, fmax = ", fmin, fmax
     ! Calcualte the Hyperfine/Zeeman matrix in units of MHz with B measured in Gauss

     do f = fmin, fmax, 2
        do m = -f, f, 2
           bra = bra+1
           ket = 0
           do fp = fmin, fmax, 2
              c1 = gsLi7*muB*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))
              c2 = giLi7*muN*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))
              do mp = -fp, fp, 2
                 Zsum1 = 0d0
                 Zsum2 = 0d0
                 ket = ket+1
                 if((f.eq.fp).and.(m.eq.mp)) then
                    HHZ(bra, ket) = HHZ(bra,ket) + &
                         0.5d0*ALi7*( 0.5d0*f*(0.5d0*f + 1d0) - 0.5d0*i*(0.5d0*i + 1d0) - 0.5d0*s*(0.5d0*s + 1d0) )
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
        enddo
     enddo
     
   end subroutine makeLi7HHZ
   
   subroutine OverlapSpinHF(i1,i2,s1,s2,S,mS,I,mI,f1,m1,f2,m2,me1)
     implicit none
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
  
  subroutine STproj(S,i1,i2,s1,s2,f1p,m1p,f2p,m2p,f1,m1,f2,m2,res)
    implicit none
    ! i1 = nuclear spin of atom 1, mi1 = z-projection of i1
    ! i2 = nuclear spin of atom 2, mi2 = z-projection of i2
    ! I = total nuclear spin
    ! s1 = electronic spin of atom 1
    ! s2 = electronic spin of atom 2
    ! S = total electronic spin
    ! mS, mI = z-projections of total nuclear and electronic spins
    ! f = S+I = total spin of atom
    ! S = 0 for singlet, 2 for triplet (since this is really 2S)
    integer S,I,mS,mI,mi1,ms1,mi2,ms2
    integer i1,i2,s1,s2,f1p,m1p,f2p,m2p,f1,m1,f2,m2
    double precision res,me1,me2

!    write(6,*) "For i1,i2,s1,s2 = ",i1,i2,s1,s2
    res=0d0
    if((S.eq.0).or.(S.eq.2)) then
       if ((S.le.(s1+s2)) .and. (S.ge.iabs(s1-s2))) then
          do I = iabs(i1-i2), i1+i2, 2
             do mS = -S, S, 2
                do mI = -I, I, 2
                   call OverlapSpinHF(i1,i2,s1,s2,S,mS,I,mI,f1,m1,f2,m2,me1)
                   call OverlapSpinHF(i1,i2,s1,s2,S,mS,I,mI,f1p,m1p,f2p,m2p,me2)
                   res = res + me1*me2
                   !write(6,*) I,mI,S,mS, res
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine STproj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine makeProjMatricesLi7(SPmat,TPmat,mtot,size1)
    implicit none
    double precision SPmat(:,:), TPmat(:,:)
    integer size1, s1, i1, s2, i2, Stot, bra, ket,mtot
    integer f1, f2, m1, m2, f1p, f2p, m1p, m2p
    integer fmin, fmax
    double precision SME1, SME2, TME1, TME2
    
    s1 = sLi7
    i1 = iLi7
    s2 = sLi7
    i2 = iLi7
    fmin = iabs(s1-i1)
    fmax = s1+i1
   
    do bra = 1, 5
       do ket = 1, 5
          Stot = 0
          call STproj(Stot,i1,i2,s1,s2,f1p,m1p,f2p,m2p,f1,m1,f2,m2,SME1)
          Stot = 1
          call STproj(Stot,i1,i2,s1,s2,f1p,m1p,f2p,m2p,f1,m1,f2,m2,TME1)
          SPmat(bra, ket) = SME1
          TPmat(bra, ket) = TME1
       enddo
    enddo
    
  end subroutine makeProjMatricesLi7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SymSTproj(sym,l,S,i1,i2,s1,s2,f1p,m1p,f2p,m2p,f1,m1,f2,m2,res)
    integer i1,i2,s1,s2,f1p,m1p,f2p,m2p,f1,m1,f2,m2,S
    integer l,sym
    double precision res,res1,res2,denom

    ! sym > 0 for bosons
    ! sym < 0 for fermions
    ! sym = 0 for distinguishible atoms

    call STproj(S,i1,i2,s1,s2,f1p,m1p,f2p,m2p,f1,m1,f2,m2,res1)
    if (sym.eq.0) then
       res = res1
    else if (sym.lt.0) then
       denom = 2d0
       call STproj(S,i1,i2,s1,s2,f2p,m2p,f1p,m1p,f2,m2,f1,m1,res2)
       res = (res1 - (-1)**l * res2)/denom    !Note 'l' is NOT really 2*l here -- it's the real 'l'
    else 
       denom = 2d0
       call STproj(S,i1,i2,s1,s2,f2p,m2p,f1p,m1p,f2,m2,f1,m1,res2)
       res = (res1 + (-1)**l * res2)/denom    !Note 'l' is NOT really 2*l here -- it's the real 'l'
    endif
  end subroutine SymSTproj

end module hyperfine

!****************************************************************************************************
program main
  use units
  use hyperfine
  implicit none
  integer ISTATE, IMN1, IMN2, NPP, iR, iB, ihf
  integer i1,i2,s1,s2,S,sym,size1,size2,NBgrid,mtot
  double precision, allocatable :: XO(:), VSinglet(:), VTriplet(:), RM2(:)
  double precision, allocatable :: HHZ(:,:), Bgrid(:), EVAL(:), EVEC(:,:)
  double precision, allocatable :: SPmat(:,:), TPmat(:,:), VHZ(:,:,:)
  double precision VLIM,xmin,xmax,dx
  double precision B
  type(hf1atom) a1, a2
  type(hf2atom) mstate1, mstate2
  
  ! initialize the angular momentum routines
  call setupam
  ! determine the size of the one-atom hyperfine/zeeman hamiltonian
  call HHZSize1(iLi7,sLi7,size1)
  write(6,*) "size1 = ",size1
  allocate(HHZ(size1,size1),EVAL(size1),EVEC(size1,size1))

  ! construct the HZ Hamiltonian
  NBgrid = 40
  allocate(Bgrid(NBgrid))

  call makeEgrid(Bgrid,NBgrid,0.0d0,1000d0)
  do iB = 1, NBgrid
     EVEC=0d0
     EVAL=0d0
     HHZ=0d0
     call makeLi7HHZ(HHZ,Bgrid(iB),size1)
     call printmatrix(HHZ,size1,size1,6)
     call MyDSYEV(HHZ,size1,EVAL,EVEC)
     write(20,*) Bgrid(iB), EVAL
  enddo
  !call makeHFBasisLi7
10 format(100F10.4)

  NPP = 100
  VLIM = 0d0
  
  allocate(XO(NPP),VSinglet(NPP),VTriplet(NPP),RM2(NPP))

  IMN1 = 7
  IMN2 = 7
  sym = 1 ! set to +1 for bosonic Li-7, -1 for fermionic Li-6, and 0 for mixture Li-6 - Li-7.
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
     write(10,*) XO(iR)*AngstromPerBohr, VSinglet(iR)*HartreePerinvcm
  enddo

  ISTATE = 2                !Find the triplet potential
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VTriplet)
  do iR=1, NPP
     write(30,*) XO(iR)*AngstromPerBohr, VTriplet(iR)*HartreePerinvcm
  enddo
  mtot = 4
  call HHZSize2(iLi7, sLi7, iLi7, sLi7, mtot, size2)
!  write(6,*) "size2 = ", size2
!  size2=5
!  allocate(SPmat(size2,size2), TPmat(size2,size2))

  !  call makeProjMatricesLi7(SPmat,TPmat,mtot,size2)
  
!  allocate(VHZ(NPP,size2,size2))

  
!  deallocate(SPmat,TPmat,VHZ)
  deallocate(Bgrid,XO,VSinglet,VTriplet,RM2)
  deallocate(HHZ,EVAL,EVEC)

end program 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
