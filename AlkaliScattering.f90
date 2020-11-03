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
  integer, parameter :: sRb87 = 1 ! twice the electronic spin of Rb-87
  integer, parameter :: iRb87 = 3 ! twice the nuclear spin of Rb-87
  integer, parameter :: sRb85 = 1 ! twice the electronic spin of Rb-85
  integer, parameter :: iRb85 = 5 ! twice the nuclear spin of Rb-85
  integer, parameter :: sK39 = 1 ! twice the electronic spin of K-39
  integer, parameter :: iK39 = 3 ! twice the nuclear spin of K-39
  integer, parameter :: sK40 = 1 ! twice the electronic spin of K-40
  integer, parameter :: iK40 = 8 ! twice the nuclear spin of K-40
  double precision, parameter :: gs = 2.00231930436256d0 ! electron g factor
  double precision, parameter :: giRb87 = -0.000995141410 ! Nuclear g factor for Rb-87
  double precision, parameter :: giRb85 = -0.00029364006; ! Nuclear g factor for Rb-85
  double precision, parameter :: giK39 = -0.00014193489 ! Nuclear g factor for K-39
  double precision, parameter :: giK40 = 0.00017649034; ! Nuclear g factor for K-40
  ! multiply nuclear g-factors by Bohr magneton 
  !  Hyperfine from Arimondo et al Rev. Mod. Phys. 49, 31 (1977).  See table on pg. 67 
  double precision, parameter :: ARb87 = 3417.3413064215 ! MHz  
  double precision, parameter :: ARb85 = 1011.9108132 ! MHz
  double precision, parameter :: AK39 = 230.85986013 ! MHz  
  double precision, parameter :: AK40 = -285.730824 ! MHz

  ! masses
  double precision, parameter :: mRb87 = 86.909180527*amuAU
  double precision, parameter :: mRb85 = 84.911789738*amuAU
  double precision, parameter :: mK39 = 38.963708d0*amuAU
  double precision, parameter :: mK40 = 39.964008d0*amuAU

  
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
    c2 = (-1)*gi*muB*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))
    
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
  integer NPP, iR, iB, ihf, n, i, j, A
  integer i1,i2,s1,s2,S,sym,size1,size2,NBgrid,NEgrid,mtot, lwave,NumOpen
  double precision, allocatable :: XO(:), VSinglet(:), VTriplet(:), RM2(:), R(:)
  double precision, allocatable :: weights(:),ystart(:,:),yi(:,:),yf(:,:),Kmat(:,:),identity(:,:)
  double precision, allocatable :: HHZ(:,:), Bgrid(:),Egrid(:), EVAL(:), EVEC(:,:),RotatedVHZ(:,:,:),TEMP(:,:)
  double precision, allocatable :: SPmat(:,:), TPmat(:,:), VHZ(:,:,:), HHZ2(:,:),AsymChannels(:,:),Eth(:)
  double precision VLIM,xmin,xmax,dx,mured
  double precision Bmin, Bmax, CGhf,Energy,h
  integer NA,iE
  double precision RX, RF
  double precision, allocatable :: alpha(:,:)
  double precision, external :: VLRLi, rint
  type(hf1atom) a1, a2
  type(hf2atom) mstate1, mstate2
  type(hf1atom), allocatable :: K39hf(:), K40hf(:)  
  TYPE(ScatData) :: SD
  ! set the reduced mass 
  mured = 0.5d0*mK39

  ! initialize the angular momentum routines
  call setupam

  ! determine the size of the one-atom hyperfine/zeeman hamiltonian
  NBgrid = 500
  NEgrid = 20
  Bmin = 100d0
  Bmax = 300d0
  !make the magnetic field grid and energy grid
  allocate(Bgrid(NBgrid),Egrid(NEgrid))
  call GridMaker(Bgrid,NBgrid,Bmin,Bmax,'linear')
  call GridMaker(Egrid,NEgrid,1d-15,1d-10,'linear') ! measure the collision energy in Hartree
  !------------------------------------------------------------------
  ! solve the K-39 1-atom Zeeman problem
  !call once with size1 = 0 to determine the size of the basis.
  size1 = 0
  call MakeHF1Basis(iK39,sK39,size1,K39hf)
  write(6,*) "size of the Li-7 1-atom hyperfine basis = ",size1
  !allocate memory for the basis
  allocate(K39hf(size1),EVEC(size1,size1),EVAL(size1),HHZ(size1,size1))
  call MakeHF1Basis(iK39,sK39,size1,K39hf)
  
  write(6,*) "1-atom K-39 states:"
  do i = 1, size1
     write(6,'(I4,I4)') K39hf(i)
  enddo


  ! construct the HZ Hamiltonian
  EVEC=0d0
  EVAL=0d0
  HHZ=0d0
  call MakeHHZ1(HHZ,0d0,size1,K39hf,gs,giK39,AK39,sK39,iK39)
  call MyDSYEV(HHZ,size1,EVAL,EVEC)
  CGhf = SUM(EVAL)!/dble(size1)

  write(6,*) "CGhf = ", CGhf
  do iB = 1, NBgrid
     EVEC=0d0
     EVAL=0d0
     HHZ=0d0
     call MakeHHZ1(HHZ,Bgrid(iB),size1,K39hf,gs,giK39,AK39,sK39,iK39)
     !     call printmatrix(HHZ,size1,size1,6)
     call MyDSYEV(HHZ,size1,EVAL,EVEC)
     write(90,*) Bgrid(iB), EVAL - CGhf
  enddo
  
  deallocate(EVEC,EVAL,HHZ)
  !____________________________________________________________________
  !------------------------------------------------------------------
  ! solve the K-40 1-atom Zeeman problem
  size1=0
  call MakeHF1Basis(iK40,sK40,size1,K40hf)
  write(6,*) "size of the K-40 1-atom hyperfine basis = ",size1
  allocate(K40hf(size1))
  call MakeHF1Basis(iK40,sK40,size1,K40hf)
  write(6,*) "1-atom K-40 states:"
  do i = 1, size1
     write(6,'(I4,I4)') K40hf(i)
  enddo

  allocate(EVEC(size1,size1),EVAL(size1),HHZ(size1,size1))

  ! construct the HZ Hamiltonian
  
  do iB = 1, NBgrid
     EVEC=0d0
     EVAL=0d0
     HHZ=0d0
     call MakeHHZ1(HHZ,Bgrid(iB),size1,K40hf,gs,giK40,AK40,sK40,iK40)
!     call printmatrix(HHZ,size1,size1,6)
     call MyDSYEV(HHZ,size1,EVAL,EVEC)
     write(91,*) Bgrid(iB), EVAL
  enddo
  deallocate(EVAL,EVEC,HHZ)
 ! stop
 !____________________________________________________________________
  
  NPP = 100000
  VLIM = 0d0
  A = 40  !A = 39 for K-39, 40 for K-40
  allocate(XO(NPP),R(NPP),weights(NPP),VSinglet(NPP),VTriplet(NPP),RM2(NPP))

  sym = -1 ! set to +1 for bosonic K-39, Rb-85, and Rb-87; -1 for fermionic K-40; 0 for any mixed collision
  lwave = 0 ! s wave collisions
  xmin = 0.03d0 ! use atomic units here (bohr)
  xmax = 300.0d0 ! use atomic units here (bohr)
  dx = (xmax-xmin)/(dble(NPP-1))
  do iR=1, NPP
     R(iR) = xmin + (iR-1)*dx
  enddo

  write(6,*) 'setting up potentials' 
  XO(:) = R(:)*BohrPerAngstrom
  ! Find the "singlet" X(^1\Sigma_g^+) curve
   call PotassiumSinglet(A,XO,NPP,VLIM,VSinglet)
!  call Rubidium87Singlet(XO,NPP,VLIM,VSinglet)

  ! print singlet potential to fort.10
  do iR=1, NPP
     VSinglet(iR) = VSinglet(iR)*InvcmPerHartree
     write(10,*) R(iR), VSinglet(iR)     
  enddo
  
 !Find the triplet potential
   call PotassiumTriplet(A,XO,NPP,VLIM,VTriplet)
 ! call Rubidium87Triplet(XO,NPP,VLIM,VTriplet)
 
  do iR=1, NPP
     VTriplet(iR) = VTriplet(iR)*InvcmPerHartree
     write(30,*) R(iR), VTriplet(iR)
  enddo

  
  mtot = -16
  
  call MakeHF2Basis(iK40, sK40, iK40, sK40, sym, lwave, mtot, size2)
  write(6,*) "size of the symmetrized 2-atom hyperfine basis = ", size2
  allocate(SPmat(size2,size2), TPmat(size2,size2))

  S = 0
  call MakeSTHFProj(S, iK40, iK40, sK40, sK40, hf2symTempGlobal, size2, SPmat)
  write(6,*) "The singlet projection matrix:"
  write(6,*) "-------------------------------"
  call printmatrix(SPmat,size2,size2,6)

  S = 2
  call MakeSTHFProj(S, iK40, iK40, sK40, sK40, hf2symTempGlobal, size2, TPmat)
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
  open(unit = 50, file = "K40FR.dat")

  do iB = 1, NBgrid 
     yi(:,:)=ystart(:,:)
     call MakeHHZ2(Bgrid(iB),AK40,AK40,gs,giK40,giK40,iK40,sK40,iK40,sK40,hf2symTempGlobal,size2,HHZ2)
     HHZ2(:,:) = HHZ2(:,:)*MHzPerHartree
     !Find the asymptotic channel states
     VHZ(:,:,NPP) = VSinglet(NPP)*SPmat(:,:) + VTriplet(NPP)*TPmat(:,:) + HHZ2(:,:)
     call MyDSYEV(VHZ(:,:,NPP),size2,EVAL,AsymChannels)
     Eth(:)=EVAL(:)
     

     do iE = 1, 1
        NumOpen=0        
        Energy = Eth(1) + Egrid(iE)!*nKPerHartree
        !write(6,*) "energy = ", energy
        do j = 1, size2
           if(energy.gt.Eth(j)) NumOpen = NumOpen+1
        enddo

        do iR = 1, NPP
           VHZ(:,:,iR) = VSinglet(iR)*SPmat(:,:) + VTriplet(iR)*TPmat(:,:) + HHZ2(:,:)
           !call MyDSYEV(VHZ(:,:,iR),size2,EVAL,EVEC)  ! Calculate the adiabatic curves

           !----- Rotate into the asymptotic eigenstates -------!
           call dgemm('T','N',size2,size2,size2,1d0,AsymChannels,size2,VHZ(:,:,iR),size2,0d0,TEMP,size2)
           call dgemm('N','N',size2,size2,size2,1d0,TEMP,size2,AsymChannels,size2,0d0,RotatedVHZ(:,:,iR),size2)
           !-----------------------------------------------------!
           !write(1000+iB,*) R(iR), (RotatedVHZ(n,n,iR), n=1,NumOpen+1)!EVAL(:)!VHZ(:,:,iR)
        enddo
        !     write(6,*) "yi = "
        !     call printmatrix(yi,size2,size2,6) 
        call logderprop(mured,Energy,identity,weights,NPP,yi,yf,R,RotatedVHZ,size2)
        !     write(6,*) "yf = "
        !     call printmatrix(yf,size2,size2,6)
        call CalcK(yf,R(NPP),SD,mured,3d0,1d0,Energy,Eth,size2,NumOpen)
        !     write(2000+iB,*) Bgrid(iB), (SD%K(j,j), j=1,size2)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
        write(50,'(100f20.10)') Bgrid(iB), (-SD%K(j,j)/dsqrt(2d0*mured*(Energy-Eth(j))), j=1,NumOpen)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
        !        write(50,'(100d20.10)') Egrid(iE), (-SD%K(j,j)/dsqrt(2d0*mured*(Energy-Eth(j))), j=1,NumOpen)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
        !        write(6,'(100d20.10)') Egrid(iE), (-SD%K(j,j)/dsqrt(2d0*mured*(Energy-Eth(j))), j=1,NumOpen)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
        write(6,'(100d16.6)') Bgrid(iB), (-SD%K(j,j)/dsqrt(2d0*mured*(Energy-Eth(j))), j=1,NumOpen)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
     enddo
  enddo
  close(50)
!  deallocate(SPmat,TPmat,VHZ)
  deallocate(Bgrid,XO,R,VSinglet,VTriplet,RM2)

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RK4StepMilne(y,mu,lwave,energy,h,R)
  implicit none
  double precision h,y(2),mu,energy,f(2),k1(2),k2(2),k3(2),k4(2),R
  integer lwave

  call dydR_Milne(R,y,mu,lwave,energy,f)
  k1 = h * f
  call dydR_Milne(R + 0.5d0*h,y,mu,lwave,energy,f)
  k2 = h * f
  call dydR_Milne(R + 0.5d0*h,y,mu,lwave,energy,f)
  k3 = h * f
  call dydR_Milne(R + h,y,mu,lwave,energy,f)
  k4 = h * f

  y = y + k1/6.0d0 + k2/3.0d0 + k3/3.0d0 + k4/6.0d0

end subroutine RK4StepMilne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function VLRLi(mu,lwave,R)
  implicit none
  integer lwave
  double precision R, mu
  ! in atomic units
!  double precision, parameter :: C6=1393.39D0, C8=83425.5D0, C10=7.37211D6
  ! in van der Waals units
  !double precision, parameter :: C6=0.985829, C8=0.0138825, C10=0.000288538516
  double precision, parameter :: C6=1d0, C8=0d0, C10=0d0
  
  VLRLi = -C6/R**6 - C8/R**8 - C10/R**10 + lwave*(lwave+1)/(2*mu*R*R)
  ! van der Waals units (mu->1/2)
  !VLRLi = -C6/R**6 - C8/R**8 - C10/R**10 + lwave*(lwave+1)/(R*R)
  
end function VLRLi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function VLRLiPrime(mu,lwave,R)
  implicit none
  integer lwave
  double precision R, mu
  ! in atomic units
!  double precision, parameter :: C6=1393.39D0, C8=83425.5D0, C10=7.37211D6
  ! in van der Waals units
  !double precision, parameter :: C6=0.985829, C8=0.0138825, C10=0.000288538516
  double precision, parameter :: C6=1d0, C8=0d0, C10=0d0
  
  VLRLiPrime = +6*C6/R**7 + 8*C8/R**9 + 10*C10/R**11 - 2*lwave*(lwave+1)/(2*mu*R*R*R)
  ! van der Waals units (mu->1/2)
  !VLRLi = -C6/R**6 - C8/R**8 - C10/R**10 + lwave*(lwave+1)/(R*R)
  
end function VLRLiPrime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function ksqrLRLi(energy,mu,lwave,R)
  implicit none
  double precision mu, R, energy
  double precision, external :: VLRLi
  integer lwave
  ksqrLRLi = 2d0*mu*(energy - VLRLi(mu,lwave,R))
end function ksqrLRLi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dydR_Milne(R,y,mu,lwave,energy,f)
  implicit none
  integer lwave
  double precision R,y(2),mu,energy,f(2)
  double precision, external :: ksqrLRLi
  
  f(1) = y(2)
  f(2) = y(1)**(-3) - ksqrLRLi(energy,mu,lwave,R)*y(1)
  
end subroutine dydR_Milne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcMilne(R,alpha,NA,energy,lwave,mu)
  implicit none
  double precision RX, RF, h, energy, y(2), mu
  integer NA, iR, lwave
  double precision alpha(NA,2), R(NA)
  double precision, external :: ksqrLRLi, VLRLi, VLRLiPrime
  ! make a radial grid

  h = R(2)-R(1)
  RX=R(1)
  RF=R(NA)
  ! set the initial conditions (WKB-like boundary conditions at RX)
  alpha(1,1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLRLi(mu,lwave,RX))/RX**2))**(-0.25d0)
  
  !ksqrLRLi(energy,mu,lwave,RX)**(-0.25d0)
  alpha(1,2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRLiPrime(mu, lwave, RX))) &
       /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLRLi(mu,lwave,RX))**1.25d0)
     
     !-0.5d0*ksqrLRLi(energy,mu,lwave,RX)**(-0.75d0)
  y = alpha(1,:)
  do iR = 1, NA
     alpha(iR,:) = y
     call RK4StepMilne(y,mu,lwave,energy,h,R(iR))
  enddo
  alpha(NA,:)=y
  
end subroutine CalcMilne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Rubidium87Singlet(XO,NPP,VLIM,VSinglet)
  implicit none
  integer NPP, i, j
  double precision VLIM, VSinglet(NPP), XO(NPP)
  double precision Rsr, Rlr, Ns, u1, u2, b, Rm, xi
  double precision C6, C8, C10, C26, Aex, gamma, beta
  double precision, DIMENSION(26) :: a

  Rsr = 3.126d0
  Ns = 4.53389d0
  u1 = -0.638904880d4
  u2 = 0.112005361d7
  b = -0.13d0
  Rlr = 11.00d0
  Rm = 4.209912706d0
  C6 = 0.2270032d8
  C8 = 0.7782886d9
  C10 = 0.2868869d11
  C26 = 0.2819810d26
  Aex = 0.1317786d5
  gamma = 5.317689d0
  beta = 2.093816d0
  a = (/-3993.592873d0, 0.0d0, 0.282069372972346137d5, 0.560425000209256905d4, -0.423962138510562945d5, &
       -0.598558066508841584d5, &
       -0.162613532034769596d5,-0.405142102246254944d5, 0.195237415352729586d6, 0.413823663033582852d6, &
       -0.425543284828921501d7, 0.546674790157210198d6, 0.663194778861331940d8,-0.558341849704095051d8, &
       -0.573987344918535471d9, 0.102010964189156187d10, 0.300040150506311035d10,-0.893187252759830856d10, &
       -0.736002541483347511d10, 0.423130460980355225d11,-0.786351477693491840d10,-0.102470557344862152d12, &
       0.895155811349267578d11,0.830355322355692902d11,-0.150102297761234375d12,0.586778574293387070d11 /)

  do i = 1, NPP
     if (XO(i) .LE. Rsr) then
        VSinglet(i) = u1 + u2/(XO(i)**Ns)
     else if ((XO(i) .GE. Rsr) .AND. (XO(i) .LE. Rlr)) then
        xi = (XO(i) - Rm)/(XO(i) + b*Rm)
        VSinglet(i) = 0d0

        do j = 0, 25
           VSinglet(i) = VSinglet(i) + a(j+1)*xi**j
        enddo

     else 
        VSinglet(i) = -C6/(XO(i)**6.0d0) - C8/(XO(i)**8.0d0) - C10/(XO(i)**10.0d0) - C26/(XO(i)**26.0d0) &
             - Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))

     endif

  enddo
   
end subroutine Rubidium87Singlet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Rubidium87Triplet(XO,NPP,VLIM,VTriplet)
  implicit none
  integer NPP, i, j
  double precision VLIM, VTriplet(NPP), XO(NPP)
  double precision Rsr, Rlr, Ns, u1, u2, b, Rm, xi
  double precision C6, C8, C10, C26, Aex, gamma, beta
  double precision, DIMENSION(13) :: a

  Rsr = 5.07d0
  Ns = 4.5338950d0
  u1 = -0.619088543d3
  u2 = 0.956231677d6
  b = -0.33d0
  Rlr = 11.00d0
  Rm = 6.093345d0
  C6 = 0.2270032d8
  C8 = 0.7782886d9
  C10 = 0.2868869d11
  C26 = 0.2819810d26
  Aex = 0.1317786d5
  gamma = 5.317689d0
  beta = 2.093816d0
  a = (/-241.503352d0, -0.672503402304666542d0, 0.195494577140503543d4, -0.141544168453406223d4,&
       -0.221166468149940465d4, 0.165443726445793004d4, -0.596412188910614259d4, &
       0.654481694231538040d4, 0.261413416681972012d5, -0.349701859112702878d5,&
       -0.328185277155018630d5,0.790208849885562522d5, -0.398783520249289213d5 /)

  do i = 1, NPP
     if (XO(i) .LE. Rsr) then
        VTriplet(i) = u1 + u2/(XO(i)**Ns)
     else if ((XO(i) .GE. Rsr) .AND. (XO(i) .LE. Rlr)) then
        xi = (XO(i) - Rm)/(XO(i) + b*Rm)
        VTriplet(i) = 0d0

        do j = 0, 12
           VTriplet(i) = VTriplet(i) + a(j+1)*xi**j
        enddo
        
     else 
        VTriplet(i) = -C6/(XO(i)**6.0d0) - C8/(XO(i)**8.0d0) - C10/(XO(i)**10.0d0) - C26/(XO(i)**26.0d0) &
             + Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))
     endif

  enddo
endsubroutine Rubidium87Triplet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PotassiumTriplet(A,XO,NPP,VLIM,VTriplet) !A = 39 for K-39, A = 40 for K-40
  implicit none
  integer NPP, i, j, A
  double precision VLIM, VTriplet(NPP), XO(NPP)
  double precision Rsr, Rlr, Ns, u1, u2, b, Rm, xi, nu
  double precision C6, C8, C10, Aex, gamma, beta
  double precision, DIMENSION(21) :: acf
  double precision mu39, mu40

  Rsr = 4.750d0
  Ns = 6d0
  u1 = -0.6948000684d3
  u2 = 0.7986755824d7
  b = -0.300d0
  Rlr = 12.00d0
  Rm = 5.73392370d0
  C6 = 0.1892652670d8
  C8 =  0.5706799527d9
  C10 = 0.1853042722d11
  Aex =0.90092159d4
  gamma = 5.19500d0
  beta = 2.13539d0
  acf = (/-255.015289d0, -0.84057856111142d0, 0.20960112217307d4, -0.17090298954603d4, &
       -0.17873773359495d4, 0.29451253739583d4, -0.20200089247397d5, &
       -0.35699524005434d5, 0.59869055371895d6, -0.71054353363636d6,&
       -0.61711841390175d7,0.19365507566961d8, 0.67930587059121d7, &
       -0.12020061704172d9, 0.21603959986951d9,-0.63531969223760d8,-0.52391212820709d9, &
       0.15913304648629d10,-0.24792546567713d10,0.20326031881106d10,-0.68044508325774d9/)
  nu = 0.23803737d0
  mu39 = 35513.3d0
  mu40 = 36425d0
  
if (A.EQ. 39) then    !K-39 Triplet
  do i = 1, NPP
     if (XO(i) .LE. Rsr) then
        VTriplet(i) = u1 + u2/(XO(i)**Ns)
     else if ((XO(i) .GE. Rsr) .AND. (XO(i) .LE. Rlr)) then
        xi = (XO(i) - Rm)/(XO(i) + b*Rm)
        VTriplet(i) = 0d0

        do j = 0, 20
           VTriplet(i) = VTriplet(i) + acf(j+1)*xi**j
        enddo
        
     else 
        VTriplet(i) = -C6/(XO(i)**6.0d0) - C8/(XO(i)**8.0d0) - C10/(XO(i)**10.0d0) &
             + Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))
     endif

  enddo
else if (A.EQ.40) then  !K-40 Triplet
     do i = 1, NPP
     if (XO(i) .LE. Rsr) then
        VTriplet(i) = u1 + u2/(XO(i)**Ns)
     else if ((XO(i) .GE. Rsr) .AND. (XO(i) .LE. Rlr)) then
        xi = (XO(i) - Rm)/(XO(i) + b*Rm)
        VTriplet(i) = 0d0

        do j = 0, 20
           VTriplet(i) = VTriplet(i) + acf(j+1)*xi**j
        enddo
        
     else 
        VTriplet(i) = -C6/(XO(i)**6.0d0) - C8/(XO(i)**8.0d0) - C10/(XO(i)**10.0d0) &
             + Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))
     endif

     VTriplet(i) = VTriplet(i) + nu*(1 - mu39/mu40)*((2*Rm)/(XO(i)+Rm))**6
  enddo
else
   WRITE(6,*) "Invalid A value"
   stop
end if

endsubroutine PotassiumTriplet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PotassiumSinglet(A,XO,NPP,VLIM,VSinglet) !A = 39 for K-39, A = 40 for K-40
  implicit none
  integer NPP, i, j, A
  double precision VLIM, VSinglet(NPP), XO(NPP)
  double precision Rsr, Rlr, Ns, u1, u2, b, Rm, xi, yi
  double precision C6, C8, C10, Aex, gamma, beta
  double precision, DIMENSION(32) :: acf
  double precision mu39, mu40
  double precision nu0, nu1

  Rsr = 2.870d0
  Ns = 12d0
  u1 = -0.263145571d4
  u2 = 0.813723194d9
  b = -0.400d0
  Rlr = 12.00d0
  Rm = 3.92436437d0
  C6 = 0.1892652670d8
  C8 =  0.5706799527d9
  C10 = 0.1853042722d11
  Aex = 0.90092159d4
  gamma = 5.19500d0
  beta = 2.13539d0
  acf = (/-4450.899484d0, 0.30601009538111d-1, 0.13671217000518d5,0.10750910095361d5, &
       -0.20933401680991d4,-0.19385874804675d5,-0.49208915890513d5,0.11026639220148d6, &
       0.72867339500920d6,-0.29310679369135d7,-0.12407070106619d8,0.40333947198094d8, &
       0.13229848871390d9, -0.37617673798775d9,-0.95250413275787d9,0.24655585744641d10, &
       0.47848257695164d10,-0.11582132109947d11,-0.17022518297651d11,0.39469335034593d11, &
       0.43141949844339d11,-0.97616955325128d11,-0.77417530685917d11,0.17314133615879d12, &
       0.96118849114926d11,-0.21425463041449d12,-0.78513081754125d11, 0.17539493131251d12, &
       0.37939637008662d11,-0.85271868691526d11,-0.82123523240949d10,0.18626451751424d11/)
  nu0 = 0.13148609d0
  nu1 = 2.08523853d0
  mu39 = 35513.3d0
  mu40 = 36425d0
  
if (A.EQ.39) then
  do i = 1, NPP
     if (XO(i) .LE. Rsr) then
        VSinglet(i) = u1 + u2/(XO(i)**Ns)
     else if ((XO(i) .GE. Rsr) .AND. (XO(i) .LE. Rlr)) then
        xi = (XO(i) - Rm)/(XO(i) + b*Rm)
        VSinglet(i) = 0d0

        do j = 0, 31
           VSinglet(i) = VSinglet(i) + acf(j+1)*xi**j
        enddo
        
     else 
        VSinglet(i) = -C6/(XO(i)**6.0d0) - C8/(XO(i)**8.0d0) - C10/(XO(i)**10.0d0) &
             + Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))
     endif

  enddo
else if (A.EQ.40) then
     do i = 1, NPP
     if (XO(i) .LE. Rsr) then
        VSinglet(i) = u1 + u2/(XO(i)**Ns)
     else if ((XO(i) .GE. Rsr) .AND. (XO(i) .LE. Rlr)) then
        xi = (XO(i) - Rm)/(XO(i) + b*Rm)
        VSinglet(i) = 0d0

        do j = 0, 31
           VSinglet(i) = VSinglet(i) + acf(j+1)*xi**j
        enddo
        
     else 
        VSinglet(i) = -C6/(XO(i)**6.0d0) - C8/(XO(i)**8.0d0) - C10/(XO(i)**10.0d0) &
             + Aex*(XO(i)**gamma)*EXP((-1)*beta*XO(i))
     endif
     
     yi = (XO(i) - Rm)/(XO(i) + b*Rm)
     VSinglet(i) = VSinglet(i) + (nu0 + nu1*yi)*(1 - mu39/mu40)*((2*Rm)/(XO(i)+Rm))**6
  enddo
  
else
   WRITE(6,*) "Invalid A value"
   stop
end if

endsubroutine PotassiumSinglet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
