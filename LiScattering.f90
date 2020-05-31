!****************************************************************************************************
MODULE DataStructures
  IMPLICIT NONE

  !****************************************************************************************************
  TYPE ScatData

    DOUBLE PRECISION, ALLOCATABLE :: K(:,:), R(:,:), sigma(:,:), sigmatot(:)
    complex*16, ALLOCATABLE :: f(:,:), S(:,:), T(:,:)
    double precision delta, tandel, sindel,sin2del
  END TYPE ScatData

CONTAINS
  !****************************************************************************************************
  SUBROUTINE AllocateScat(SD,N)
    IMPLICIT NONE
    TYPE(ScatData) SD
    INTEGER N
    ALLOCATE(SD%K(N,N),SD%R(N,N),SD%sigma(N,N))
    ALLOCATE(SD%f(N,N),SD%S(N,N),SD%T(N,N))
    ALLOCATE(SD%sigmatot(0:2*N))
    SD%K=0d0
    SD%R=0d0
    SD%f=0d0
    SD%sigma=0d0
    SD%S=(0d0,0d0)
    SD%T=(0d0,0d0)
  END SUBROUTINE AllocateScat
  !****************************************************************************************************
  SUBROUTINE DeAllocateScat(SD)
    IMPLICIT NONE
    TYPE(ScatData) SD

    DEALLOCATE(SD%K,SD%R,SD%f,SD%sigma,SD%S,SD%T,SD%sigmatot)

  END SUBROUTINE DeAllocateScat

END MODULE DataStructures
module scattering
  use datastructures
  !****************************************************************************************************

CONTAINS
 SUBROUTINE CalcK(Y,rm,SD,mu,d,alpha,EE,Eth,NumChannels,NumOpen)
!   use DipoleDipole
   IMPLICIT NONE
   TYPE(ScatData) :: SD

   DOUBLE PRECISION mu, EE, rm, d, alpha,Y(NumChannels,NumChannels)
   DOUBLE PRECISION, ALLOCATABLE :: JJ(:),NN(:),JJp(:),NNp(:)
   double precision, allocatable :: Ktemp1(:,:),Ktemp2(:,:)
   DOUBLE PRECISION rhypj,rhypy,rhypjp,rhypyp,Pi,rhypi,rhypk,rhypip,rhypkp,ldrhk,ldrhi
   DOUBLE PRECISION k(NumChannels),Eth(NumChannels)
   complex*16, allocatable :: tmp(:,:),Identity(:,:)
   complex*16  II
   INTEGER i,j,no,nw,nc,beta,NumChannels,NumOpen

   
   II=(0d0,1d0)
   Pi=dacos(-1d0)

   no=0
   nw=0
   nc=0

   DO i = 1,NumChannels
      IF (EE.GE.Eth(i)) THEN
         k(i) = dsqrt(2d0*mu*(EE-Eth(i))) ! k is real
         no=no+1
      ELSE
         k(i) = dsqrt(2d0*mu*(Eth(i)-EE)) ! k->kappa, kappa is real
         IF( (k(i)*rm).LT.10d0) nw = nw+1
         IF( (k(i)*rm).GE.10d0) nc = nc+1
      ENDIF
   ENDDO
   !      write(6,*) "no = ", no
   IF((no+nw+nc).NE.NumChannels) THEN
      WRITE(6,*) "Channel miscount in calcK"
      STOP
   ENDIF
!   write(6,*) "no = ", no
   deallocate(SD%S,SD%T,SD%sigma)
   
   allocate(SD%S(no,no),SD%T(no,no),SD%sigma(no,no))
   ALLOCATE(JJ(NumChannels),NN(NumChannels),tmp(no,no))
   ALLOCATE(JJp(NumChannels),NNp(NumChannels))
   allocate(Ktemp1(NumChannels,NumChannels))
   allocate(Ktemp2(NumChannels,NumChannels))
   allocate(Identity(NumChannels,NumChannels))
   Identity = 0d0;

   DO i = 1,no
      Identity(i,i) = 1d0
      !write(6,*) k(i), rm
      CALL hyperrjry(INT(d),alpha,0d0,k(i)*rm,rhypj,rhypy,rhypjp,rhypyp)
      JJ(i) = rhypj/dsqrt(Pi*k(i))
      NN(i) = -rhypy/dsqrt(Pi*k(i))
      JJp(i) = dsqrt(k(i)/Pi)*rhypjp
      NNp(i) = -dsqrt(k(i)/Pi)*rhypyp
   ENDDO
   do i=no+1,NumChannels
      Identity(i,i) = 1d0
      CALL hyperrirk(INT(d),alpha,0d0,k(i)*rm,rhypi,rhypk,rhypip,rhypkp,ldrhi,ldrhk)
      JJ(i) = 1d0
      NN(i) = -1d0
      JJp(i) = ldrhi
      NNp(i) = ldrhk

   ENDDO

   Ktemp1=0d0
   Ktemp2=0d0
   SD%K=0d0
   do i=1,NumChannels
      Ktemp1(i,i) = NNp(i)
      Ktemp2(i,i) = JJp(i)
      do j=1,NumChannels
         Ktemp1(i,j) = Ktemp1(i,j) - Y(i,j)*NN(j)
         Ktemp2(i,j) = Ktemp2(i,j) - Y(i,j)*JJ(j)
      enddo
   enddo
   call sqrmatinv(Ktemp1,NumChannels)
   SD%K = -MATMUL(Ktemp1,Ktemp2)
   

   tmp = Identity(1:no,1:no) - II*SD%K(1:no,1:no)
 
   SD%S = Identity(1:no,1:no) + II*SD%K(1:no,1:no)
   call CompSqrMatInv(tmp,no)
   SD%S = MATMUL(SD%S,tmp)
   SD%T = -II*0.5d0*(SD%S-Identity(1:no,1:no))
   SD%sigma = conjg(SD%T)*SD%T*Pi/(2d0*mu*EE)

   SD%tandel = SD%K(1,1)
   SD%delta = atan(SD%tandel)
   SD%sindel = sin(SD%delta)
   SD%sin2del = SD%sindel**2

 END SUBROUTINE CalcK
end module scattering
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

  ! masses
  double precision, parameter :: mLi7 = 7.01600455000*amuAU
  double precision, parameter :: mLi6 = 6.0151228874*amuAU

  
  type hf1atom
     integer f, m
  end type hf1atom
  type hf2atom
     type(hf1atom) a1,a2
  end type hf2atom
  type symhf2atom
     type(hf2atom) state1, state2
     double precision norm
     integer phase
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
       hf2symTempGlobal(n1)%phase = 1
       if(sym.lt.0) hf2symTempGlobal(n1)%phase = -1
       hf2symTempGlobal(n1)%phase = hf2symTempGlobal(n1)%phase*(-1)**lwave
       
       write(6,'(I4,A,f5.2,A,4I2,A,I3,A,4I2,A)') n1,"  symmetrized basis:", hf2symTempGlobal(n1)%norm, &
            " ( |", hf2symTempGlobal(n1)%state1,"> +", hf2symTempGlobal(n1)%phase,"|", hf2symTempGlobal(n1)%state2,"> )"
    enddo
    size=symsize
    deallocate(symstates,tempstates)
    
  end subroutine MakeHF2Basis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine HHZ1ME(B,hf1bra,hf1ket,gs,gi,Ahf,s,i,res)
    implicit none
    double precision, intent(in) :: B
    integer f, m, fp, mp, i, s, bra, ket, ms, mi, size
    double precision Zsum1, Zsum2, c1, c2, tj1, tj2, gs, gi, Ahf, res
    double precision, external :: THRJ
    type(hf1atom) hf1bra, hf1ket

    f = hf1ket%f
    m = hf1ket%m
    fp = hf1bra%f
    mp = hf1bra%m
    c1 = gs*muB*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))
    c2 = gi*muN*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))
    
    Zsum1 = 0d0
    Zsum2 = 0d0
    res = 0d0
    if((f.eq.fp).and.(m.eq.mp)) then
       res = 0.5d0*Ahf*( 0.5d0*f*(0.5d0*f + 1d0) - 0.5d0*i*(0.5d0*i + 1d0) - 0.5d0*s*(0.5d0*s + 1d0) )
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
       res = res + c1*Zsum1 + c2*Zsum2
    endif
  end subroutine HHZ1ME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MakeHHZ1(HHZ,B,size,hf1,gs,gi,Ahf,s,i)
    implicit none
    integer  i, s, bra, ket, size
    double precision gs, gi, Ahf
    double precision, intent(in) :: B
    double precision HHZ(size,size)
    double precision res
    type(hf1atom) hf1(size)
    
    ! Calcualte the Hyperfine/Zeeman matrix in units of MHz with B measured in Gauss
    do bra = 1, size
       do ket = 1, size
          call HHZ1ME(B,hf1(bra),hf1(ket),gs,gi,Ahf,s,i,res)
          HHZ(bra,ket) = res
       enddo
    enddo
    
  end subroutine MakeHHZ1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MakeHHZ2(B,Ahf1,Ahf2,gs,gi1,gi2,i1,s1,i2,s2,hf2sym,size,HHZ2)
     implicit none
     !ah1b = <alphaprime | h1 | beta>
     !ah2b = <alphaprime | h2 | beta>, etc...
     !primed states are always "bra" states: <primed | hi | unprimed>
     double precision B, Ahf1, Ahf2, gs, gi1, gi2, ah1a, ah1b, bh1a, bh1b, bh2b, bh2a, ah2b, ah2a
     double precision kdaa, kdab, kdba, kdbb, prefact
     integer i1, s1, i2, s2, size, bra, ket, p, pprime
     type(symhf2atom) hf2sym(size)
     type(hf1atom) alpha, alphaprime, beta, betaprime
     double precision HHZ2(size,size)
     double precision, external :: dkdelta
     
     do bra = 1, size
        do ket = 1, size
           ! total of 8 terms when we don't know the symmetry of the atoms.
           ! the 1 particle states can be gleaned purlely from state1 of
           ! the symmetrized state, namely |alpha beta>, since state2 is just |beta alpha>
           prefact = hf2sym(bra)%norm*hf2sym(ket)%norm
           alpha = hf2sym(ket)%state1%a1
           beta = hf2sym(ket)%state1%a2
           alphaprime = hf2sym(bra)%state1%a1
           betaprime = hf2sym(bra)%state1%a2
           pprime = hf2sym(bra)%phase
           p = hf2sym(ket)%phase
           kdaa = dkdelta(alphaprime%f,alpha%f)*dkdelta(alphaprime%m,alpha%m)
           kdab = dkdelta(alphaprime%f,beta%f)*dkdelta(alphaprime%m,beta%m)
           kdbb = dkdelta(betaprime%f,beta%f)*dkdelta(betaprime%m,beta%m)
           kdba = dkdelta(betaprime%f,alpha%f)*dkdelta(betaprime%m,alpha%m)
           ! use a different nuclear g-factor, gi, and hyperfine coupling Ahf, for different 1-atom Zeeman Hamiltonians
           call HHZ1ME(B,alphaprime,alpha,gs,gi1,Ahf1,s1,i1,ah1a)
           call HHZ1ME(B,alphaprime,beta,gs,gi1,Ahf1,s1,i1,ah1b)
           call HHZ1ME(B,betaprime,alpha,gs,gi1,Ahf1,s1,i1,bh1a)
           call HHZ1ME(B,betaprime,beta,gs,gi1,Ahf1,s1,i1,bh1b)
           call HHZ1ME(B,alphaprime,alpha,gs,gi2,Ahf2,s2,i2,ah2a)
           call HHZ1ME(B,alphaprime,beta,gs,gi2,Ahf2,s2,i2,ah2b)
           call HHZ1ME(B,betaprime,alpha,gs,gi2,Ahf2,s2,i2,bh2a)
           call HHZ1ME(B,betaprime,beta,gs,gi2,Ahf2,s2,i2,bh2b)
           HHZ2(bra, ket) = prefact*(ah1a*kdbb + p*ah1b*kdba + pprime*bh1a*kdab + p*pprime*bh1b*kdaa + &
                bh2b*kdaa + p*bh2a*kdab + pprime*ah2b*kdba + p*pprime*ah2a*kdbb)
           
        enddo
     enddo

   end subroutine MakeHHZ2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine OverlapSpinHF(i1,i2,s1,s2,S,mS,I,mI,f1,m1,f2,m2,me1)
     implicit none
     ! This subroutine calcualtes the overlap between the hyperfine state |f1 m1 f2 m2> and the
     ! spin state |(s1 s2) S mS (i1 i2) I mI >
     ! it returns in me1 = <(s1 s2) S mS (i1 i2) I mI | (i1 s1) f1 m1 (i2 s2) f2 m2 >
     ! I haven't checked this, but it should be equal to the 9-J coefficient corresponding to the above overal.
     integer S,mS,I,mI,f1,m1,f2,m2,i1,i2,s1,s2,phaseexp
     integer mi1,mi2,ms1,ms2,phase
     double precision prefact,tj1,tj2,tj3,tj4,me1
     double precision, external :: THRJ

     phaseexp = 2*i2-2*s1-m1-m2-mS-mI
     phase = 1
     if(mod(phaseexp,2).eq.0) then
        phaseexp = phaseexp/2
        phase = (-1)**phaseexp
     else
        write(6,*) "Warning: phase factor will be complex, but you're using only real variables.  Phase is not set."
     endif
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
  use datastructures
  use scattering
  implicit none
  integer ISTATE, IMN1, IMN2, NPP, iR, iB, ihf, n, i, j
  integer i1,i2,s1,s2,S,sym,size1,size2,NBgrid,NEgrid,mtot, lwave,NumOpen
  double precision, allocatable :: XO(:), VSinglet(:), VTriplet(:), RM2(:)
  double precision, allocatable :: weights(:),ystart(:,:),yi(:,:),yf(:,:),Kmat(:,:),identity(:,:)
  double precision, allocatable :: HHZ(:,:), Bgrid(:),Egrid(:), EVAL(:), EVEC(:,:),RotatedVHZ(:,:,:),TEMP(:,:)
  double precision, allocatable :: SPmat(:,:), TPmat(:,:), VHZ(:,:,:), HHZ2(:,:),AsymChannels(:,:),Eth(:)
  double precision VLIM,xmin,xmax,dx,mured
  double precision Bmin, Bmax, CGhf,Energy
  type(hf1atom) a1, a2
  type(hf2atom) mstate1, mstate2
  type(hf1atom), allocatable :: Li6hf(:), Li7hf(:)  
  TYPE(ScatData) :: SD
  ! set the reduced mass (for Li-7 collisions)
  mured = 0.5d0*mLi7

  ! initialize the angular momentum routines
  call setupam
  
  ! determine the size of the one-atom hyperfine/zeeman hamiltonian
  NBgrid = 500
  NEgrid = 500
  Bmin = 0d0
  Bmax = 1000d0
  !make the magnetic field grid and energy grid
  allocate(Bgrid(NBgrid),Egrid(NEgrid))
  call GridMaker(Bgrid,NBgrid,Bmin,Bmax,'linear')
  call GridMaker(Egrid,NEgrid,0.0d0,0.001d0,'linear')
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
  EVEC=0d0
  EVAL=0d0
  HHZ=0d0
  call MakeHHZ1(HHZ,0d0,size1,Li7hf,gs,giLi7,ALi7,sLi7,iLi7)
  call MyDSYEV(HHZ,size1,EVAL,EVEC)
  CGhf = SUM(EVAL)!/dble(size1)
  write(6,*) "CGhf = ", CGhf
  do iB = 1, NBgrid
     EVEC=0d0
     EVAL=0d0
     HHZ=0d0
     call MakeHHZ1(HHZ,Bgrid(iB),size1,Li7hf,gs,giLi7,ALi7,sLi7,iLi7)
     !     call printmatrix(HHZ,size1,size1,6)
     call MyDSYEV(HHZ,size1,EVAL,EVEC)
     write(90,*) Bgrid(iB), EVAL - CGhf
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
  
  NPP = 100000
  VLIM = 0d0

  allocate(XO(NPP),weights(NPP),VSinglet(NPP),VTriplet(NPP),RM2(NPP))

  IMN1 = 7
  IMN2 = 7
  sym = 1 ! set to +1 for bosonic Li-7, -1 for fermionic Li-6, and 0 for mixture Li-6 - Li-7.
  lwave = 0 ! s wave collisions
  xmin = 0.05d0
  xmax = 400.0d0
  dx = (xmax-xmin)/(dble(NPP-1))
  do iR=1, NPP
     XO(iR) = xmin + (iR-1)*dx
  enddo
  ISTATE = 1                ! Find the "singlet" X(^1\Sigma_g^+) curve
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VSinglet)

  ! print singlet potential to fort.10
  do iR=1, NPP
     VSinglet(iR) = VSinglet(iR)*InvcmPerHartree
     write(10,*) XO(iR)*AngstromPerBohr, VSinglet(iR)
  enddo

  ISTATE = 2                !Find the triplet potential
  ! Call the Le Roy routine to generate the MLR potential
  call POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VTriplet)
  do iR=1, NPP
     VTriplet(iR) = VTriplet(iR)*InvcmPerHartree
     write(30,*) XO(iR)*AngstromPerBohr, VTriplet(iR)
  enddo
  XO(:) = XO(:)*AngstromPerBohr
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

  allocate(VHZ(size2,size2,NPP),RotatedVHZ(size2,size2,NPP),TEMP(size2,size2),AsymChannels(size2,size2))
  allocate(HHZ2(size2,size2))
  allocate(EVEC(size2,size2),EVAL(size2))

  ! prepare some arrays for the log-derivative propagator
  do iR = 1, NPP-1
     if (mod(iR,2).eq.1)  weights(iR)=4d0
     if (mod(iR,2).eq.0)  weights(iR)=2d0
  enddo
  weights(NPP)=1d0
  call AllocateScat(SD,size2)
  allocate(identity(size2,size2),ystart(size2,size2),yi(size2,size2),yf(size2,size2),Eth(size2))
  identity(:,:)=0d0
  ystart(:,:)=0d0
  yf(:,:)=0d0
  do i=1,size2
     identity(i,i)=1d0
     ystart(i,i)=1d20
  enddo


  do iB = 1, NBgrid ! consider only B=0 for now
     yi(:,:)=ystart(:,:)
     call MakeHHZ2(Bgrid(iB),ALi7,ALi7,gs,giLi7,giLi7,iLi7,sLi7,iLi7,sLi7,hf2symTempGlobal,size2,HHZ2)
     HHZ2(:,:) = HHZ2(:,:)*MHzPerHartree
     !Find the asymptotic channel states
     VHZ(:,:,NPP) = VSinglet(NPP)*SPmat(:,:) + VTriplet(NPP)*TPmat(:,:) + HHZ2(:,:)
     call MyDSYEV(VHZ(:,:,NPP),size2,EVAL,AsymChannels)
     Eth(:)=EVAL(:)
     
     NumOpen=0
     Energy = Eth(1)+1d-15
     do j = 1, size2
        if(energy.gt.Eth(j)) NumOpen = NumOpen+1
     enddo

     do iR = 1, NPP
        VHZ(:,:,iR) = VSinglet(iR)*SPmat(:,:) + VTriplet(iR)*TPmat(:,:) + HHZ2(:,:)
        call MyDSYEV(VHZ(:,:,iR),size2,EVAL,EVEC)
        call dgemm('T','N',size2,size2,size2,1d0,AsymChannels,size2,VHZ(:,:,iR),size2,0d0,TEMP,size2)
        call dgemm('N','N',size2,size2,size2,1d0,TEMP,size2,AsymChannels,size2,0d0,RotatedVHZ(:,:,iR),size2)
        !write(1000+iB,*) XO(iR), (RotatedVHZ(n,n,iR), n=1,NumOpen+1)!EVAL(:)!VHZ(:,:,iR)
     enddo
!     write(6,*) "yi = "
!     call printmatrix(yi,size2,size2,6)
     call logderprop(mured,Energy,identity,weights,NPP,yi,yf,XO,RotatedVHZ,size2)
!     write(6,*) "yf = "
!     call printmatrix(yf,size2,size2,6)
     call CalcK(yf,XO(NPP),SD,mured,3d0,1d0,Energy,Eth,size2,NumOpen)
!     write(2000+iB,*) Bgrid(iB), (SD%K(j,j), j=1,size2)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
     write(2000,'(100f20.6)') Bgrid(iB), (-SD%K(j,j)/dsqrt(2d0*mured*(Energy-Eth(j))), j=1,NumOpen)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
     write(6,'(100d16.6)') Bgrid(iB), (-SD%K(j,j)/dsqrt(2d0*mured*(Energy-Eth(j))), j=1,NumOpen)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
     
  enddo
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine logderprop(mu,Energy,identity,weights,NPP,yi,yf,XO,Pot,size)
  implicit none
  integer i,j,step,NPP,size
  double precision h,Energy,mu
  double precision xx(NPP),weights(NPP),XO(NPP)
  double precision yi(size,size),yf(size,size)
  double precision Pot(size,size,NPP) !make sure Pot includes the threshold offsets
  double precision tempy(size,size),ycurrent(size,size),yprevious(size,size),identity(size,size)
  double precision vtemp1(size,size), vtemp2(size,size), un(size,size)
  double precision, parameter :: onesixth = 0.166666666666666666666d0
  double precision, parameter :: onethird = 0.333333333333333333333d0
  
  h=XO(2)-XO(1)
  
  yprevious = yi
  do step = 1, NPP
     !write(6,*) weights(step)
     vtemp1 = 2d0*mu*(identity*Energy-Pot(:,:,step))
!     write(6,*) "vtemp1 = "
!     call printmatrix(vtemp1,size,size,6)
     if (mod(step,2).eq.0) then
        un = vtemp1
!        write(6,*) step, "un = "
!        call printmatrix(un,size,size,6)
     else
        vtemp2 = identity + h*h*onesixth*vtemp1
        call sqrmatinv(vtemp2,size)
        un = matmul(vtemp2,vtemp1)
!        write(6,*) step, "un = "
!        call printmatrix(un,size,size,6)
     endif     
     tempy = identity + h*yprevious
     call sqrmatinv(tempy,size)
     ycurrent = MATMUL(tempy,yprevious) - onethird*h*weights(step)*un
!     write(6,*) "ycurrent = "
!     call printmatrix(ycurrent,size,size,6)
     yprevious = ycurrent
!     write(6,*) "ycurrent = ", ycurrent
  enddo
  !stop
  yf(:,:) = ycurrent(:,:)
  
end subroutine logderprop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)
  
  DO j = 1,nr
     WRITE(file,30) (M(j,k), k = 1,nc)
  ENDDO
  
20 FORMAT(1P,100D16.8)
30 FORMAT(100D14.4)
END SUBROUTINE printmatrix
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
  ! Quadratic grid:
  !--------------------------------------------
  IF((scale.EQ."quadratic").and.(numpts.gt.1)) THEN
     DE=(E2-E1)
     DO iE=1,numpts
        grid(iE) = E1 + ((iE-1)/DBLE(numpts-1))**2*DE
     ENDDO
  ENDIF
END SUBROUTINE GridMaker
