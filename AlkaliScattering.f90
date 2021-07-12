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
    !SD%sigma = conjg(SD%T)*SD%T*Pi/(2d0*mu*EE)

    SD%tandel = SD%K(1,1)
    SD%delta = atan(SD%tandel)
    SD%sindel = sin(SD%delta)
    SD%sin2del = SD%sindel**2
    SD%sigma = (4*pi/(2d0*mu*(EE-Eth(1))))*SD%sin2del

  END SUBROUTINE CalcK
end module scattering
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module hyperfine
  use units
  implicit none
  double precision, parameter :: muB = 1.3996244936142    !Bohr magneton in units (MHz/Gauss)
  double precision, parameter :: gs = 2.00231930436256d0 ! electron g factor
  double precision, parameter :: muN = 0.000762267d0 !Nuclear magneton in units (MHz/Gauss)


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MakeHF1Basis(i,s,size,hf1)
    implicit none
    ! Gives the size of the (unsymmetrized) hyperfine basis (dimension of the 1-atom hyperfine + zeeman Hamiltonian) for one atom
    integer i, s, size ! size is (in/out) if size = 0 on input -- determine the size.  Else, use it to calculate the basis stored in hf1
    integer f, m, count, n
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
       write(6,'(A)') "   #   |  f   mf  >"
       write(6,'(A)') "-------------------" 
       do n = 1, size
          write(6,'(I4,A4,I3,A1,I3,A4)') n, '   |', hf1(n)%f,' ', hf1(n)%m, '  >'
       enddo
    endif



    
    !    stop
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

    call setzero2atom(zerostate)

    size = 0
    do f1 = iabs(i1-s1), i1+s1, 2
       do m1 = -f1, f1, 2
          do f2 = iabs(i2-s2), i2+s2, 2
             do m2 = -f2, f2, 2
                if((m1+m2).eq.mtot) then   
                   size = size + 1      !Determines the size of the unsym basis
                endif
             enddo
          enddo
       enddo
    enddo

    allocate(hf2TempGlobal(size))
    do n1=1,size
       call setzero2atom(hf2TempGlobal(n1))
    enddo

    write(6,'(A,I3,A1)') "Generating the unsymmetrized basis with total MF = ", mtot/2,":"
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
    write(6,*)
    
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
    write(6,'(A,I3)') "(+1/0/-1) = (boson/dist/fermion) Symmetry case: ", sym
    write(6,'(A,I2)') "partial wave = ", lwave
    write(6,'(A)') "The symmetrized basis in term of unsymmetrized states | f1 m1 f2 m2 >:"
    do n1 = 1, symsize
       hf2symTempGlobal(n1)=symstates(n1)
       ! record the normalization and relative phase into the state
       kd1 = dkdelta(hf2symTempGlobal(n1)%state1%a1%f,hf2symTempGlobal(n1)%state1%a2%f)
       kd2 = dkdelta(hf2symTempGlobal(n1)%state1%a1%m,hf2symTempGlobal(n1)%state1%a2%m)

       hf2symTempGlobal(n1)%norm = 1d0/dsqrt(2d0*(1d0+kd1*kd2))
       hf2symTempGlobal(n1)%phase = 1
       if(sym.lt.0) hf2symTempGlobal(n1)%phase = -1
       hf2symTempGlobal(n1)%phase = hf2symTempGlobal(n1)%phase*(-1)**lwave

       write(6,'(I4,A,f5.2,A,4I3,A,I3,A,4I3,A)') n1," : ", hf2symTempGlobal(n1)%norm, &
            " ( |", hf2symTempGlobal(n1)%state1," > +", hf2symTempGlobal(n1)%phase,"|", hf2symTempGlobal(n1)%state2," > )"
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
    c2 = gi*muB*B*dsqrt((dble(fp)+1d0)*(dble(f)+1d0))

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
  use bspline
  use InterpType
  use units
  use hyperfine
  use datastructures
  use scattering
  use quadrature

  implicit none
  integer NPP, iR, iB, ihf, n, i, j, A, nc, nr, ESTATE, ISTATE, IMN1
  integer nspin1, nspin2, espin1, espin2 !nspin is the nuclear spin and espin is the electronic spin
  integer i1,i2,s1,s2,S,sym,size1,size2,NBgrid,NEgrid,mtot, lwave,NumOpen
  integer Nsr, Nlr, CALCTYPE
  logical writepot
  double precision, allocatable :: Rsr(:), Rlr(:),VsrSinglet(:),VsrTriplet(:),VlrSinglet(:),VlrTriplet(:)
!  double precision, allocatable :: RotatedVsrHZ(:,:,:),RotatedVlrHZ(:,:,:)
  double precision, allocatable :: wSR(:),wLR(:),Vdum(:),Rdum(:)
!  double precision, allocatable :: VSinglet(:), VTriplet(:), R(:)
  double precision, allocatable :: weights(:),ystart(:,:),yi(:,:),ym(:,:),yf(:,:),Kmat(:,:),identity(:,:)
  double precision, allocatable :: HHZ(:,:), Bgrid(:),Egrid(:), EVAL(:), EVEC(:,:)
  double precision, allocatable :: EThreshMat(:,:),TEMP(:,:)
  double precision, allocatable :: cotgamma(:,:,:),Ktemp1(:,:),Ktemp2(:,:),KPQ(:,:),KQP(:,:)
  double precision, allocatable :: SPmat(:,:), Sdressed(:,:), Tdressed(:,:), TPmat(:,:)
  double precision, allocatable :: VHZ(:,:), HHZ2(:,:),AsymChannels(:,:),Eth(:),Ksr(:,:),RmidArray(:)
  double precision VLIM,Rmin,dR,phiL,stepsize,Vmin,ascat
  double precision gi1,gi2,Ahf1,Ahf2,MU,MUREF,mass1,mass2
  double precision Bmin, Bmax, Emin, Emax, CGhf,Energy,h, betavdw,wavek
  integer iE, iRmid,NRmid,Ndum,NXM,NXF,multnsr, multnlr,ith
  double precision RX, Rmid, rmidmax, RF, Cvals(3), Ktilde
  double precision SingletQD, TripletQD, x
  double precision, external :: VLR, rint, abar
  character(len=20), external :: str
  CHARACTER(LEN=20) scale
  type(hf1atom) a1, a2
  type(hf2atom) mstate1, mstate2
  type(hf1atom), allocatable :: hf1(:) 
  TYPE(ScatData) :: SD
  type(InterpolatingFunction) :: InterpCotGamma

!  call testinterpolation
 
  
  read(5,*)
  read(5,*)
  read(5,*) ISTATE, sym, CALCTYPE
  read(5,*)
  read(5,*) lwave, mtot, writepot
  read(5,*)
  read(5,*) Rmin, Rmidmax, RX, RF
  read(5,*)
  read(5,*) NXM, NXF, Nsr, Nlr
  read(5,*)
  read(5,*) NEgrid, Emin, Emax  ! enter in mK
  read(5,*)
  read(5,*) NBgrid, Bmin, Bmax  ! enter in Gauss
  read(5,*)
  read(5,*) scale

  Emin = Emin/HartreePermK
  Emax = Emax/HartreePermK

  !make the magnetic field grid and energy grid
  allocate(Egrid(NEgrid),Bgrid(NBgrid))
  call GridMaker(Bgrid,NBgrid,Bmin,Bmax,'linear')
  call GridMaker(Egrid,NEgrid,Emin,Emax,'linear') ! measure the collision energy in Hartree
  !--------------------------------------------------!
  call AtomData(ISTATE, AHf1, nspin1, espin1, gi1, MU, MUREF, mass1)  !atom 1 data (and atom 2 for identical particles)
  !write(6,*) AHf1, nspin1, espin1, gi1, MU, MUREF, mass1
  ! initialize the angular momentum routines
  call setupam
  ! initialize the Cash-Karp RK parameters (this may not even be used)
  call initCashKarp

  !------------------------------------------------------------------
  ! determine the size of the one-atom hyperfine/zeeman hamiltonian
  !call once with size1 = 0 to determine the size of the basis.
  size1 = 0
  call MakeHF1Basis(nspin1,espin1,size1,hf1)
  write(6,'(A,I3)') "size of the 1-atom hyperfine basis = ",size1
  !allocate memory for the basis
  allocate(hf1(size1),EVEC(size1,size1),EVAL(size1),HHZ(size1,size1))
  call MakeHF1Basis(nspin1,espin1,size1,hf1)

  ! construct the HZ Hamiltonian
  EVEC=0d0
  EVAL=0d0
  HHZ=0d0
  call MakeHHZ1(HHZ,0d0,size1,hf1,gs,gi1,AHf1,espin1,nspin1)
  call MyDSYEV(HHZ,size1,EVAL,EVEC)
  CGhf = SUM(EVAL)!/dble(size1) ! hyperfine "center of gravity"

  ! Calculate and print the hyperfine/Zeeman one-atom levels
  write(6,*) "See fort.90 for the one-atom hyperfine/zeeman levels in MHz vs. Gauss..."
  write(6,*)
  write(90,*) "# One-atom hyperfine/zeeman levels in MHz vs. Gauss..."
  do iB = 1, NBgrid
     EVEC=0d0
     EVAL=0d0
     HHZ=0d0
     call MakeHHZ1(HHZ,Bgrid(iB),size1,hf1,gs,gi1,AHf1,espin1,nspin1)
     call MyDSYEV(HHZ,size1,EVAL,EVEC)
     write(90,*) Bgrid(iB), EVAL - CGhf
  enddo
  write(6,'(A)') '-----------------------------------'
  write(6,'(A,f6.1,A)') "Bmax = ", Bgrid(NBgrid), " gauss"
  write(6,'(A,f10.3,A)') "The one atom ground state has E = ", EVAL(1), " MHz"
  write(6,'(A)') "In the | f mf > basis, the ground state Eigenvector is: "
  do i = 1, size1
     write(6,'(f12.6, A4,I3,A1,I3,A4)') EVEC(i,1),'   |', hf1(i)%f,' ', hf1(i)%m, '  >'
  enddo
  write(6,'(A)') '-----------------------------------'

  deallocate(EVEC,EVAL,HHZ)
  !____________________________________________________________________
  call MakeHF2Basis(nspin1, espin1, nspin1, espin1, sym, lwave, mtot, size2)
  write(6,'(A,I3)') "size of the symmetrized 2-atom hyperfine basis = ", size2
  write(6,*)
!  stop
  allocate(SPmat(size2,size2), TPmat(size2,size2))
  allocate(Sdressed(size2,size2), Tdressed(size2,size2))
  allocate(EThreshMat(size2,size2),Ksr(size2,size2))
  
  S = 0
  call MakeSTHFProj(S, nspin1, nspin1, espin1, espin1, hf2symTempGlobal, size2, SPmat)
  write(6,*) "The singlet projection matrix:"
  write(6,*) "-------------------------------"
  call printmatrix(SPmat,size2,size2,6)
  write(6,*)
  S = 2
  call MakeSTHFProj(S, nspin1, nspin1, espin1, espin1, hf2symTempGlobal, size2, TPmat)
  write(6,*) "The triplet projection matrix:"
  write(6,*) "-------------------------------"
  call printmatrix(TPmat,size2,size2,6)
  write(6,*)
  allocate(TEMP(size2,size2),AsymChannels(size2,size2))
  allocate(HHZ2(size2,size2))
  allocate(EVEC(size2,size2),EVAL(size2))
  call AllocateScat(SD,size2)
  allocate(identity(size2,size2),ystart(size2,size2),yi(size2,size2),ym(size2,size2),yf(size2,size2),Eth(size2))
  allocate(cotgamma(size2-1,size2-1,NEgrid))
  allocate(Ktemp1(size2-1,size2-1),Ktemp2(1,1))
  allocate(KPQ(1,size2-1),KQP(size2-1,1))
  identity(:,:)=0d0
  ystart(:,:)=0d0
  yf(:,:)=0d0
  do i=1,size2
     identity(i,i)=1d0
     ystart(i,i)=1d20
  enddo

  open(unit = 20, file = "CollisionThresholds-"//trim(str(ISTATE))//"-"//trim(str(CALCTYPE))//".dat")
  open(unit = 50, file = "SigmaEnergyDependence-"//trim(str(ISTATE))//"-"//trim(str(CALCTYPE))//".dat")
  open(unit = 51, file = "ScatLenFieldDependence-"//trim(str(ISTATE))//"-"//trim(str(CALCTYPE))//".dat")
  open(unit = 52, file = "QDFieldDependence-"//trim(str(ISTATE))//"-"//trim(str(CALCTYPE))//".dat")
  write(6,'(A)') "See file CollisionThresholds.dat for the field dependence of the scattering thresholds."
  write(6,'(A)') "See file SigmaEnergyDependence.dat for the energy-dependent cross section"
  write(6,'(A)') "See file ScatLenFieldDependence.dat for the field-dependent scattering length at threshold"
  write(6,'(A)') "See file QDFieldDependence.dat for the field-dependent scattering length at threshold"

  ! None of the code above this line depends on the radial grid.
  !----------------------------------------------------------------------------------------------------
  NRmid = 10
  allocate(RmidArray(NRmid))
  call GridMaker(RmidArray,NRmid,15d0,Rmidmax,'linear')

  ! This is just a "dummy" call to SetupPotential to obtain the discpersion C6/C8/C10 coefficients
  ! and to calculate the van der waals length.
  VLIM = 0d0
  Ndum = 100
  allocate(Vdum(Ndum),Rdum(Ndum))
  call GridMaker(Rdum,Ndum,3d0,20d0,'linear')
  call SetupPotential(ISTATE,1,muref,muref,Ndum,VLIM,Rdum*BohrPerAngstrom,Vdum,Cvals)
  Vmin = MINVAL(Vdum,1)
  write(6,'(A,f12.5)') "Vmin = ", Vmin
  !stepsize = 2d0*pi*0.1d0*(2d0*mu*(Emax - Vmin))**(-0.5d0)
  !write(6,'(A,f12.5)') "The suggested maximum step size is ", stepsize
!  write(6,'(A,I10)') "The minimum value of Nsr = ", int((Rmidmax - Rmin)/stepsize)
  call VdWLength(Cvals,betavdw,mu)
  write(6,'(A,f12.4)') "The van der Waals length is rvdw = ", betavdw

  RX = RX*betavdw
  RF = RF*betavdw
  !Nsr = multnsr!*int(Rmidmax/stepsize)
  !Nlr = multnlr!Nsr*(int(RF/Rmidmax))
  write(6,*) "Nsr and Nrl = ", Nsr, Nlr
  
  ! Diagonalize the HF-Zeeman Hamiltonian for the largest field values so we know the range
  ! of energies needed for the closed channel MQDT parameter cot(gamma)
  call MakeHHZ2(Bgrid(NBgrid),AHf1,AHf1,gs,gi1,gi1,nspin1,espin1,nspin1,espin1,hf2symTempGlobal,size2,HHZ2)
  HHZ2(:,:) = HHZ2(:,:)*MHzPerHartree
  call MyDSYEV(HHZ2,size2,Eth,AsymChannels)
  !Now that we know Eth at the largest field value, we can compute an interpolated function for cot(gamma)

  write(6,*) "Rmin = ", Rmin
  write(6,*) "RX = ", RX
  write(6,*) "RMid = ", Rmidmax
  write(6,*) "RF = ", RF
    
  allocate(VHZ(size2,size2))
!  allocate(RotatedVsrHZ(size2,size2,Nsr),RotatedVlrHZ(size2,size2,Nlr))
  allocate(Rsr(Nsr),Rlr(Nlr),wSR(Nsr),wLR(Nlr),VsrSinglet(Nsr),VlrSinglet(Nlr),VsrTriplet(Nsr),VlrTriplet(Nlr))

  
  if(CALCTYPE.eq.1) goto 11
  do iRmid = NRmid, NRmid

     Rmid = RmidArray(iRmid)

     ! prepare some arrays for the log-derivative propagator
     call GridMaker(Rsr,Nsr,Rmin,Rmid,'linear')
     call GridMaker(Rlr,Nlr,Rmid,RF,'linear')

     call SetLogderWeights(wSR,Nsr)
     call SetLogderWeights(wLR,Nlr)
     
     ESTATE = 1
     call SetupPotential(ISTATE,ESTATE,muref,muref,Nsr,VLIM,Rsr*BohrPerAngstrom,VsrSinglet,Cvals)
     call SetupPotential(ISTATE,ESTATE,muref,muref,Nlr,VLIM,Rlr*BohrPerAngstrom,VlrSinglet,Cvals)
     ESTATE = 3
     call SetupPotential(ISTATE,ESTATE,muref,muref,Nsr,VLIM,Rsr*BohrPerAngstrom,VsrTriplet,Cvals)
     call SetupPotential(ISTATE,ESTATE,muref,muref,Nlr,VLIM,Rlr*BohrPerAngstrom,VlrTriplet,Cvals)
     
     if((iRmid.eq.NRmid).and.writepot) then
        write(6,'(A)') "See fort.10 and fort.30 for the singlet/triplet potentials"
        do iR=1, Nsr
           write(10,*) Rsr(iR), abs((VsrSinglet(iR) - VLR(mu,lwave,Rsr(iR),Cvals))/VsrSinglet(iR))
           write(30,*) Rsr(iR), abs((VsrTriplet(iR) - VLR(mu,lwave,Rsr(iR),Cvals))/VsrTriplet(iR))
           !           write(100,*) Rsr(iR), (VLR(mu,lwave,Rsr(iR),Cvals)+Eth(ith), ith=1,size2)
           write(100,*) Rsr(iR), (VsrSinglet(iR)+Eth(ith), ith=1,size2)
           write(101,*) Rsr(iR), (VsrTriplet(iR)+Eth(ith), ith=1,size2)
        enddo
     endif

     call CalcPhaseStandard(RX,RF,NXF,lwave,mu,betavdw,Cvals,phiL,scale) ! calculate the phase standardization for lwave = 0
     !call CalcCotGammaFunction(RX,RF,NXF,size2,lwave,mu,betavdw,Cvals,phiL,Eth,Emin, Emax,InterpCotGamma)
     !call CalcNewCotGammaFunction(RX,RF,NXF,size2,lwave,mu,betavdw,Cvals,phiL,Eth,Emin,Emax,InterpCotGamma,scale)
!  x= - (Eth(size2) - Eth(1)) + 1d0*nkPerHartree
!  do while (x.lt.-1d0*nkPerHartree)
!     write(12,*) x, Interpolated(x,InterpCotGamma)
!     x = x + (Eth(size2) - Eth(1))/300d0
!  enddo

     EThreshMat(:,:) = 0d0
     energy = 0d0 ! Calculate the zero-energy quantum defects
     call logderQD(lwave,mu,energy,NXM,Nsr,wsr,VsrSinglet,VsrTriplet,Rsr,TripletQD,SingletQD,phiL,betavdw,RX,Cvals,scale)
!     SingletQD = -0.12d0
!     TripletQD = -0.08d0
     do iB = 1, NBgrid
!        write(6,*) "phiL = ", phiL
        yi = ystart
        call MakeHHZ2(Bgrid(iB),AHf1,AHf1,gs,gi1,gi1,nspin1,espin1,nspin1,espin1,hf2symTempGlobal,size2,HHZ2)
        HHZ2(:,:) = HHZ2(:,:)*MHzPerHartree
        !Find the asymptotic channel states
        VHZ(:,:) = VlrSinglet(Nlr)*SPmat(:,:) + VlrTriplet(Nlr)*TPmat(:,:) + HHZ2(:,:) 
        !     call MyDSYEV(VHZ(:,:),size2,EVAL,AsymChannels)
        call MyDSYEV(VHZ,size2,Eth,AsymChannels)
        !     Eth(:)=EVAL(:)
        do i=1,size2
           EThreshMat(i,i) = Eth(i)
        enddo
        
        Write(20,*) Bgrid(iB), Eth
        !Calculate the tan(gamma) matrix

        !----- Rotate into the asymptotic eigenstates (Use the B-field dressed states) -------!
        !Rotate the singlet projection operator
        call dgemm('T','N',size2,size2,size2,1d0,AsymChannels,size2,SPmat,size2,0d0,TEMP,size2)
        call dgemm('N','N',size2,size2,size2,1d0,TEMP,size2,AsymChannels,size2,0d0,Sdressed,size2)
        
        !Rotate the triplet projection operator
        call dgemm('T','N',size2,size2,size2,1d0,AsymChannels,size2,TPmat,size2,0d0,TEMP,size2)
        call dgemm('N','N',size2,size2,size2,1d0,TEMP,size2,AsymChannels,size2,0d0,Tdressed,size2)
        
!        do iR = 1, Nsr
!           RotatedVsrHZ(:,:,iR) = VsrSinglet(iR)*Sdressed(:,:) + VsrTriplet(iR)*Tdressed(:,:) + EThreshMat(:,:)
!        enddo
!        do iR = 1, Nlr
!           RotatedVlrHZ(:,:,iR) = VlrSinglet(iR)*Sdressed(:,:) + VlrTriplet(iR)*Tdressed(:,:) + EThreshMat(:,:)
!        enddo
        
        ! Start the energy loop
        do iE = 1, 1!NEgrid
           NumOpen=0        
           Energy = Eth(1) + Egrid(iE)!*nKPerHartree
           wavek = sqrt(2*mu*Egrid(iE))
           !write(6,*) "energy = ", energy
           do j = 1, size2
              if(energy.gt.Eth(j)) NumOpen = NumOpen+1
           enddo
           if(NumOpen.gt.1) then
              write(6,*) "More than one open channel -- must modify the algorithm.  Stopping"
              stop
           endif

           call CalcTanGamma(RX,RF,NXF,size2,lwave,mu,betavdw,Cvals,phiL,Egrid,NEgrid,Eth,iE,cotgamma,InterpCotGamma,scale)
!           write(6,*) "done."
           if((CALCTYPE.eq.2).and.(iE.eq.1)) then
              !call logderprop(mu,Energy,identity,wSR,Nsr,yi,ym,Rsr,RotatedVsrHZ,size2)
              call logderpropB(mu,Energy,identity,wSR,Nsr,yi,ym,Rsr,Sdressed,Tdressed,EthreshMat,VsrSinglet,VsrTriplet,size2)
              call CalcKsr(ym, Ksr, size2, NXM,RX, Rsr(Nsr), energy, Eth, lwave, mu, Cvals, betavdw, phiL,scale)
           else if (CALCTYPE.eq.3) then
              Ksr(:,:) = tan(pi*TripletQD) * Tdressed(:,:) + tan(pi*SingletQD) * Sdressed(:,:)
           endif
           
           call MyDSYEV(Ksr,size2,EVAL,EVEC)


!           write(6,*) Bgrid(iB), Ktilde!(1d0-Ktilde)*abar(lwave) !Rmid, (atan(EVAL(i))/Pi, i=1,size2)
           write(52,*) Bgrid(iB), (atan(EVAL(i))/Pi, i=1,size2)!EVAL!,Ksr

           Ktemp1 = Ksr(2:size2,2:size2) + cotgamma(:,:,iE)
           KQP = Ksr(2:size2,1:1)
           KPQ = Ksr(1:1,2:size2)
           
           call sqrmatinv(Ktemp1,size2-1)
           Ktemp2 = MATMUL(KPQ,MATMUL(Ktemp1,KQP))
!           write(6,*) "Ktemp2:"
           !           call printmatrix(Ktemp2,1,1,6)

           Ktilde = Ksr(1,1) - Ktemp2(1,1)
           ascat = (1d0-Ktilde)*abar(lwave)*betavdw
           write(6,*) Bgrid(iB), ascat , 8d0*pi*ascat**2
           write(51,*) Bgrid(iB), ascat , 8d0*pi*ascat**2
!           write(6,*) Bgrid(iB),  (1d0-Ktilde)*abar(lwave)*betavdw/(1d0 + (abar(lwave)*betavdw*wavek)**2 * Ktilde)
!           write(51,*) Bgrid(iB), (1d0-Ktilde)*abar(lwave)*betavdw/(1d0 + (abar(lwave)*betavdw*wavek)**2 * Ktilde)
        enddo
     enddo
  enddo
  stop
11 continue  
  !----------------------------------------------------------------------------------------------------
  !This next loop is does the full log-derivative calculation

  ! prepare some arrays for the log-derivative propagator
  Rmid = RmidArray(NRmid)     
  call GridMaker(Rsr,Nsr,Rmin,Rmid,'linear')  ! Make the radial grids
  call GridMaker(Rlr,Nlr,Rmid,RF,'linear')
  call SetLogderWeights(wSR,Nsr)
  call SetLogderWeights(wLR,Nlr)

  
  ESTATE = 1
  call SetupPotential(ISTATE,ESTATE,muref,muref,Nsr,VLIM,Rsr*BohrPerAngstrom,VsrSinglet,Cvals)
  call SetupPotential(ISTATE,ESTATE,muref,muref,Nlr,VLIM,Rlr*BohrPerAngstrom,VlrSinglet,Cvals)
  ESTATE = 3
  call SetupPotential(ISTATE,ESTATE,muref,muref,Nsr,VLIM,Rsr*BohrPerAngstrom,VsrTriplet,Cvals)
  call SetupPotential(ISTATE,ESTATE,muref,muref,Nlr,VLIM,Rlr*BohrPerAngstrom,VlrTriplet,Cvals)
  
  write(51,*)

  do iB = 1, NBgrid
     yi = ystart

     call MakeHHZ2(Bgrid(iB),AHf1,AHf1,gs,gi1,gi1,nspin1,espin1,nspin1,espin1,hf2symTempGlobal,size2,HHZ2)
     HHZ2 = HHZ2*MHzPerHartree

     !Find the asymptotic channel states
     VHZ(:,:) = VlrSinglet(Nlr)*SPmat(:,:) + VlrTriplet(Nlr)*TPmat(:,:) + HHZ2(:,:) 
     call MyDSYEV(VHZ(:,:),size2,Eth,AsymChannels)

     do i=1,size2
        EThreshMat(i,i) = Eth(i)
     enddo

     Write(20,*) Bgrid(iB), Eth
     !----- Rotate into the asymptotic eigenstates (Use the B-field dressed states) -------!
     !Rotate the singlet projection operator
     call dgemm('T','N',size2,size2,size2,1d0,AsymChannels,size2,SPmat,size2,0d0,TEMP,size2)
     call dgemm('N','N',size2,size2,size2,1d0,TEMP,size2,AsymChannels,size2,0d0,Sdressed,size2)

     !Rotate the triplet projection operator
     call dgemm('T','N',size2,size2,size2,1d0,AsymChannels,size2,TPmat,size2,0d0,TEMP,size2)
     call dgemm('N','N',size2,size2,size2,1d0,TEMP,size2,AsymChannels,size2,0d0,Tdressed,size2)

!     do iR = 1, Nsr
!        RotatedVsrHZ(:,:,iR) = VsrSinglet(iR)*Sdressed(:,:) + VsrTriplet(iR)*Tdressed(:,:) + EThreshMat(:,:)
!     enddo
!     do iR = 1, Nlr
!        RotatedVlrHZ(:,:,iR) = VlrSinglet(iR)*Sdressed(:,:) + VlrTriplet(iR)*Tdressed(:,:) + EThreshMat(:,:)
!     enddo

     ! Start the energy loop
     do iE = 1, 1!NEgrid
        NumOpen=0        
        Energy = Eth(1) + Egrid(iE)!*nKPerHartree
!        write(6,*) "energy = ", energy
        do j = 1, size2
           if(energy.gt.Eth(j)) NumOpen = NumOpen+1
        enddo

        !call logderpropA(mu,Energy,identity,wSR,Nsr,yi,ym,Rsr,RotatedVsrHZ,size2)
        call logderpropB(mu,Energy,identity,wSR,Nsr,yi,ym,Rsr,Sdressed,Tdressed,EthreshMat,VsrSinglet,VsrTriplet, size2) 
        !call logderpropA(mu,Energy,identity,wLR,Nlr,ym,yf,Rlr,RotatedVlrHZ,size2)
        call logderpropB(mu,Energy,identity,wLR,Nlr,ym,yf,Rlr,Sdressed,Tdressed,EthreshMat,VlrSinglet, VlrTriplet, size2)
        call CalcK(yf,RF,SD,mu,3d0,1d0,Energy,Eth,size2,NumOpen)


        ! Various output statements:
        if(iE.eq.1) then ! Write the field dependence at the lowest energy (close to threshold)
           !write(51,'(100f20.10)') Bgrid(iB), (-SD%K(j,j)/dsqrt(2d0*muref*(Energy-Eth(j))), j=1,NumOpen)!, SD%K(1,2), SD%K(2,1), SD%K(2,2)
           ascat = -SD%K(1,1)/dsqrt(2d0*mu*(Energy-Eth(1)))
           write(51,*) Bgrid(iB), ascat, 8d0*pi*ascat**2 
           write(6,*) Bgrid(iB),  ascat, 8d0*pi*ascat**2,  2d0*SD%sigma
        endif

        write(50,'(100d20.10)') Egrid(iE)*HartreePermK, 2d0*SD%sigma*(Bohrpercm**2)




     enddo
  enddo
  
  close(50)
  close(51)

10 format(100F10.4)
end program main
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
subroutine SetLogderWeights(weights,Npts)
  implicit none
  integer Npts,iR
  double precision weights(Npts)


  do iR = 1, Npts-1
     if (mod(iR,2).eq.1)  weights(iR)=4d0
     if (mod(iR,2).eq.0)  weights(iR)=2d0
  enddo
  weights(Npts)=1d0

end subroutine SetLogderWeights
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! This integrates F''(r) = U(r) F(r) where F and U are matrices
!!$! U((r) = 2*mu (V(r) - E)
!!$! Let Z(r) = ( 1 - h^2 U(r)/12 ) F(r)
!!$! Let W(r) = h^2 U(r) + h^4 U(r)U(r)/ 12 = h^2 U(r) (1 + h^2 U(r)/12)
!!$! It also follows that F(r) = (1 + W(r)/12) Z(r)
!!$! The recursion relation is Z(r+h) = (W(r)+2)Z(r) - Z(r-h)
!!$! So we will propogate Z(r) and then find F(r) from that.  F' and the log-derivative
!!$subroutine NumerovBoundaryConditions(F0,F1,Z0,Z1,h,V0,V1,size)
!!$  implicit none
!!$  integer size,i,j
!!$  double precision F0(size,size),F1(size,size),Z0(size,size),Z1(size,size)
!!$
!!$  F0(:,:) = 0d0
!!$  F1(:,:) = 0d0
!!$  Z0(:,:) = 0d0
!!$  Z1(:,:) = 0d0
!!$  do i = 1,size
!!$   end subroutine NumerovBoundaryConditions
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$subroutine NumerovMC(mu,Energy,identity,weights,NPP,yi,yf,XO,Pot,size)  
!!$  implicit none
!!$  integer i,j,step,NPP,size
!!$  double precision h,Energy,mu
!!$  double precision xx(NPP),weights(NPP),XO(NPP)
!!$  double precision yi(size,size),yf(size,size)
!!$  double precision Pot(size,size,NPP) !make sure Pot includes the threshold offsets
!!$  double precision tempy(size,size),ycurrent(size,size),yprevious(size,size),identity(size,size)
!!$  double precision vtemp1(size,size), vtemp2(size,size), un(size,size)
!!$  double precision, parameter :: onesixth = 0.166666666666666666666d0
!!$  double precision, parameter :: onethird = 0.333333333333333333333d0
!!$  double precision, parameter :: OneTwelfth = 0.083333333333333333333d0
!!$
!!$end subroutine NumerovMC
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine logderprop(mu,Energy,identity,weights,NPP,yi,yf,XO,Pot,size)
subroutine logderpropB(mu,Energy,identity,weights,NPP,yi,yf,XO,Sdressed,Tdressed,EthreshMat,Vsinglet,Vtriplet,size)
  implicit none
  integer i,j,step,NPP,size
  double precision h,Energy,mu
  double precision xx(NPP),weights(NPP),XO(NPP)
  double precision yi(size,size),yf(size,size)
  double precision Tdressed(size,size), Sdressed(size,size), Ethreshmat(size,size)
  double precision Vtriplet(NPP), Vsinglet(NPP)
  double precision pottemp(size,size)
!  double precision Pot(size,size,NPP) !make sure Pot includes the threshold offsets
  double precision tempy(size,size),ycurrent(size,size),yprevious(size,size),identity(size,size)
  double precision vtemp1(size,size), vtemp2(size,size), un(size,size)
  double precision, parameter :: onesixth = 0.166666666666666666666d0
  double precision, parameter :: onethird = 0.333333333333333333333d0

  h=XO(2)-XO(1)
  yprevious = yi
  do step = 1, NPP
     pottemp = Vsinglet(step)*Sdressed + Vtriplet(step)*Tdressed + EthreshMat
     !vtemp1 = 2d0*mu*(identity*Energy-Pot(:,:,step))
     vtemp1 = 2d0*mu*(identity*Energy-pottemp)
     if (mod(step,2).eq.0) then
        un = vtemp1
     else
        vtemp2 = identity + h*h*onesixth*vtemp1
        call sqrmatinv(vtemp2,size)
        un = matmul(vtemp2,vtemp1)
     endif
     tempy = identity + h*yprevious
     call sqrmatinv(tempy,size)
     ycurrent = MATMUL(tempy,yprevious) - onethird*h*weights(step)*un
     yprevious = ycurrent
  enddo

  yf(:,:) = ycurrent(:,:)

end subroutine logderpropB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine logderpropA(mu,Energy,identity,weights,NPP,yi,yf,XO,Pot,size)
!subroutine logderprop(mu,Energy,identity,weights,NPP,yi,yf,XO,Sdressed,Tdressed,EthreshMat,size)
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
     vtemp1 = 2d0*mu*(identity*Energy-Pot(:,:,step))
     if (mod(step,2).eq.0) then
        un = vtemp1
     else
        vtemp2 = identity + h*h*onesixth*vtemp1
        call sqrmatinv(vtemp2,size)
        un = matmul(vtemp2,vtemp1)
     endif
     tempy = identity + h*yprevious
     call sqrmatinv(tempy,size)
     ycurrent = MATMUL(tempy,yprevious) - onethird*h*weights(step)*un
     yprevious = ycurrent
  enddo

  yf(:,:) = ycurrent(:,:)

end subroutine logderpropA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)

  DO j = 1,nr
     WRITE(file,40) (M(j,k), k = 1,nc)
  ENDDO

20 FORMAT(1P,100D16.8)
30 FORMAT(100D14.4)
40 FORMAT(100F20.10)

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

END SUBROUTINE GridMaker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RK4StepMilne(y,mu,lwave,energy,h,R,Cvals)
  implicit none
  double precision h,y(2),mu,energy,f(2),k1(2),k2(2),k3(2),k4(2),R,Cvals(3)
  integer lwave

  call dydR_Milne(R,y,mu,lwave,energy,f,Cvals)
  k1 = h * f
  call dydR_Milne(R + 0.5d0*h,y,mu,lwave,energy,f,Cvals)
  k2 = h * f
  call dydR_Milne(R + 0.5d0*h,y,mu,lwave,energy,f,Cvals)
  k3 = h * f
  call dydR_Milne(R + h,y,mu,lwave,energy,f,Cvals)
  k4 = h * f

  y = y + k1/6.0d0 + k2/3.0d0 + k3/3.0d0 + k4/6.0d0
end subroutine RK4StepMilne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RK4StepNewMilne(y,mu,lwave,energy,h,R,Cvals)
  implicit none
  double precision h,y(2),mu,energy,f(2),k1(2),k2(2),k3(2),k4(2),R,Cvals(3)
  integer lwave

  call dydR_NewMilne(R,y,mu,lwave,energy,f,Cvals)
  k1 = h * f
  call dydR_NewMilne(R + 0.5d0*h,y,mu,lwave,energy,f,Cvals)
  k2 = h * f
  call dydR_NewMilne(R + 0.5d0*h,y,mu,lwave,energy,f,Cvals)
  k3 = h * f
  call dydR_NewMilne(R + h,y,mu,lwave,energy,f,Cvals)
  k4 = h * f

  y = y + k1/6.0d0 + k2/3.0d0 + k3/3.0d0 + k4/6.0d0
end subroutine RK4StepNewMilne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RK5CashKarpStepMilne(y,mu,lwave,energy,h,R,Cvals,Delta)
  use quadrature
  implicit none
  double precision h,y(2),mu,energy,f(2),k1(2),k2(2),k3(2),k4(2)
  double precision k5(2),k6(2),R,Cvals(3),Delta(2)
  integer lwave

  call dydR_Milne(R,y,mu,lwave,energy,f,Cvals)
  k1 = h * f
  call dydR_Milne(R + cka2*h, y + b21*k1, mu,lwave,energy,f,Cvals)
  k2 = h * f
  call dydR_Milne(R + cka3*h, y + b31*k1 + b32*k2, mu, lwave,energy,f,Cvals)
  k3 = h * f
  call dydR_Milne(R + cka4*h, y + b41*k1 + b42*k2 + b43*k3, mu,lwave,energy,f,Cvals)
  k4 = h * f
  call dydR_Milne(R + cka5*h, y + b51*k1 + b52*k2 + b53*k3 + b54*k4, mu,lwave,energy,f,Cvals)
  k5 = h * f
  call dydR_Milne(R + cka6*h, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5, mu,lwave,energy,f,Cvals)
  k6 = h * f

  y = y + c1*k1 + c3*k3 + c4*k4 + c6*k6

  Delta = (c1-c1s)*k1 + (c3-c3s)*k3 + (c4-c4s)*k4 + (-c5s)*k5 + (c6-c6s)*k6  

end subroutine RK5CashKarpStepMilne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function VLR(mu,lwave,R,Cvals)
  implicit none
  integer lwave
  double precision R, mu
  double precision Cvals(3)
  ! in atomic units
  !  double precision, parameter :: C6=1393.39D0, C8=83425.5D0, C10=7.37211D6
  ! in van der Waals units
  !double precision, parameter :: C6=0.985829, C8=0.0138825, C10=0.000288538516
  !  double precision, parameter :: C6=1d0, C8=0d0, C10=0d0

  VLR = -Cvals(1)/R**6 - Cvals(2)/R**8 - Cvals(3)/R**10 + lwave*(lwave+1)/(2*mu*R*R)
  ! van der Waals units (mu->1/2)
  !VLRLi = -C6/R**6 - C8/R**8 - C10/R**10 + lwave*(lwave+1)/(R*R)

end function VLR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function VLRPrime(mu,lwave,R,Cvals)
  implicit none
  integer lwave
  double precision R, mu, Cvals(3)
  ! in atomic units
  !  double precision, parameter :: C6=1393.39D0, C8=83425.5D0, C10=7.37211D6
  ! in van der Waals units
  !double precision, parameter :: C6=0.985829, C8=0.0138825, C10=0.000288538516
  !  double precision, parameter :: C6=1d0, C8=0d0, C10=0d0

  VLRPrime = 6*Cvals(1)/R**7 + 8*Cvals(2)/R**9 + 10*Cvals(3)/R**11 - 2*lwave*(lwave+1)/(2*mu*R*R*R)
  ! van der Waals units (mu->1/2)
  !VLR = -C6/R**6 - C8/R**8 - C10/R**10 + lwave*(lwave+1)/(R*R)

end function VLRPrime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function ksqrLR(energy,mu,lwave,R,Cvals)
  implicit none
  double precision mu, R, energy, Cvals(3)
  double precision, external :: VLR
  integer lwave
  ksqrLR = 2d0*mu*(energy - VLR(mu,lwave,R,Cvals))
end function ksqrLR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dydR_Milne(R,y,mu,lwave,energy,f,Cvals)
  implicit none
  integer lwave
  double precision R,y(2),mu,energy,f(2),Cvals(3)
  double precision, external :: ksqrLR

  f(1) = y(2)
  f(2) = y(1)**(-3) - ksqrLR(energy,mu,lwave,R,Cvals)*y(1)

end subroutine dydR_Milne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dydR_NewMilne(R,y,mu,lwave,energy,f,Cvals)
  ! For this routine the Milne phase alpha = exp(x) and x is propogated by the RK-4 routine
  implicit none
  integer lwave
  double precision R,y(2),mu,energy,f(2),Cvals(3)
  double precision, external :: ksqrLR

  f(1) = y(2)
  f(2) = exp(-4*y(1)) - ksqrLR(energy,mu,lwave,R,Cvals) - y(2)**2

end subroutine dydR_NewMilne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcMilne(R,alpha,NA,energy,lwave,mu,Cvals)
  implicit none
  double precision h, energy, y(2), mu,Cvals(3)
  integer NA, iR, lwave
  double precision alpha(NA,2), R(NA)

  h = R(2)-R(1)
  y = alpha(1,:)
  do iR = 1, NA-1
     call RK4StepMilne(y,mu,lwave,energy,h,R(iR),Cvals)
     alpha(iR+1,:) = y
  enddo


end subroutine CalcMilne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcMilne2step(R,alpha,NA,energy,lwave,mu,Cvals,phaseint)
  implicit none
  double precision h, energy, y(2), mu,Cvals(3),y1(2),y2(2),y3(2),hhalf
  integer NA, iR, lwave
  double precision alpha(NA,2), R(NA),Delta(2),Delta0(2),eps, phaseint(NA),Rmid,OneThird
  OneThird = 0.333333333333333333333333d0
  eps = 1d-10
  ! make a radial grid

  y1 = alpha(1,:)
  Delta0= eps*y1
  phaseint(1) = 0d0

  do iR = 1, NA-1
     h = R(iR+1)-R(iR)
     hhalf = 0.5d0*h
     ! do two half-steps and use the 3-point Simpson integration rule to calculate the phase integral
     !step from R(i) to R(i)+h/2
     y2=y1
     !call RK4StepMilne(y2,mu,lwave,energy,hhalf,R(iR),Cvals)
     call RK5CashKarpStepMilne(y2,mu,lwave,energy,hhalf,R(iR),Cvals,Delta)
     y3=y2
     !step from R(i)+h/2 to R(i+1)
     Rmid = R(iR) + hhalf
     !call RK4StepMilne(y3,mu,lwave,energy,hhalf, Rmid, Cvals)
     call RK5CashKarpStepMilne(y3,mu,lwave,energy,hhalf,Rmid,Cvals,Delta)
     alpha(iR+1,:) = y3
     phaseint(iR+1) = phaseint(iR) + OneThird*hhalf*(y1(1)**(-2) + 4d0*y2(1)**(-2) + y3(1)**(-2))
     y1 = y3
  enddo


end subroutine CalcMilne2step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcNewMilne2step(R,alpha,y,NA,energy,lwave,mu,Cvals,phaseint)
  implicit none
  double precision h, energy, y(2), mu,Cvals(3),y1(2),y2(2),y3(2),hhalf
  integer NA, iR, lwave
  double precision alpha(NA,2), R(NA),Delta(2),Delta0(2),eps, phaseint(NA),Rmid,OneThird
  OneThird = 0.333333333333333333333333d0
  eps = 1d-10
  ! make a radial grid

  y1(1) = log(alpha(1,1))
  y1(2) = alpha(1,2)/alpha(1,1)
  
!  write(6,*) "y1 = ", y1
  Delta0= eps*y1
  phaseint(1) = 0d0
  !write(6,*) R(1), phaseint(iR+1),y1(1), y1(2)
  do iR = 1, NA-1
     h = R(iR+1)-R(iR)
     hhalf = 0.5d0*h
     ! do two half-steps and use the 3-point Simpson integration rule to calculate the phase integral
     !step from R(i) to R(i)+h/2
     y2=y1
     call RK4StepNewMilne(y2,mu,lwave,energy,hhalf,R(iR),Cvals)
     !call RK5CashKarpStepMilne(y2,mu,lwave,energy,hhalf,R(iR),Cvals,Delta)
     y3=y2
     !step from R(i)+h/2 to R(i+1)
     Rmid = R(iR) + hhalf
     call RK4StepNewMilne(y3,mu,lwave,energy,hhalf, Rmid, Cvals)
     !call RK5CashKarpStepMilne(y3,mu,lwave,energy,hhalf,Rmid,Cvals,Delta)

     phaseint(iR+1) = phaseint(iR) + OneThird*hhalf*(exp(-2*y1(1)) + 4d0*exp(-2*y2(1)) + exp(-2*y3(1)))
     y1 = y3
     y=y3
     alpha(iR+1,1) = exp(y(1))
     alpha(iR+1,2) = y(2)*exp(y(1))
     !write(6,*) R(iR+1), phaseint(iR+1),y(1), y(2)
     !write(102,*) R(iR+1), phaseint(iR+1)
  enddo
!stop

end subroutine CalcNewMilne2step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcNewMilne4step(R,alpha,y,NA,energy,lwave,mu,Cvals,phaseint)
  implicit none
  double precision h, energy, y(2), mu,Cvals(3),y1(2),y2(2),y3(2),y4(2),y5(2),h4,R1,R2,R3,R4
  integer NA, iR, lwave
  double precision alpha(NA,2), R(NA),Delta(2),Delta0(2),eps, phaseint(NA)
  double precision TwoOnFortyFive

  TwoOnFortyFive = 0.04444444444444444444444444d0
  
  y1(1) = log(alpha(1,1))
  y1(2) = alpha(1,2)/alpha(1,1)

  phaseint = 0d0
  phaseint(1) = 0d0
  !write(6,*) R(1), phaseint(1), y1(1), y1(2)
  do iR = 1, NA-1
     h = R(iR+1)-R(iR)
     h4 = 0.25d0*h
     R1=R(iR)
     R2=R1+h4
     R3=R2+h4
     R4=R3+h4

     ! do 4 quarter-steps and use Bode's 5-point integration rule to calculate the phase integral
     y2=y1
     call RK4StepNewMilne(y2,mu,lwave,energy,h4,R1,Cvals)
     y3=y2
     call RK4StepNewMilne(y3,mu,lwave,energy,h4,R2,Cvals)
     y4=y3
     call RK4StepNewMilne(y4,mu,lwave,energy,h4,R3,Cvals)
     y5=y4
     call RK4StepNewMilne(y5,mu,lwave,energy,h4,R4,Cvals)

     phaseint(iR+1) = phaseint(iR) + TwoOnFortyFive*h4* &
          (7d0*exp(-2*y1(1)) + 32d0*exp(-2*y2(1)) + 12d0*exp(-2*y3(1)) + 32d0*exp(-2*y4(1)) + 7d0*exp(-2*y5(1)))
     y1=y5
     y=y5
     alpha(iR+1,1) = exp(y(1))
     alpha(iR+1,2) = y(2)*exp(y(1))
     !write(102,*) R(iR+1), phaseint(iR+1)
     !write(6,*) R(iR+1), phaseint(iR+1), y(1), y(2)
  enddo
!stop

end subroutine CalcNewMilne4step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MQDTfunctions(R1, R2, NA, scale, Cvals, mu, lwave, energy, betavdw, phiL, alpha0, alphaf,f0, g0, f0p, g0p)
  use units
  implicit none
  integer NA, i, j, lwave
  double precision Cvals(3), mu, betavdw, energy, f0, g0, f0p, g0p, phiL
  double precision alpha0(2), alphaf(2)
  double precision, allocatable :: phaseint(:)
  double precision, allocatable :: alpha(:,:),R(:)
  double precision R1, R2
  CHARACTER(LEN=*), INTENT(IN) :: scale

  allocate(R(NA),alpha(NA,2),phaseint(NA))
  call GridMaker(R,NA,R1,R2,scale)
  alpha(1,:) = alpha0
  
  call CalcMilne2step(R,alpha,NA,energy,lwave,mu,Cvals,phaseint)
  alphaf = alpha(NA,:)
  do i = NA, NA
     f0 = Pi**(-0.5d0)*alpha(i,1)*sin(phaseint(i)+phiL)
     g0 = -Pi**(-0.5d0)*alpha(i,1)*cos(phaseint(i)+phiL)
     f0p = alpha(i,2)/alpha(i,1) * f0 - g0/alpha(i,1)**2
     g0p = alpha(i,2)/alpha(i,1) * g0 + f0/alpha(i,1)**2

  enddo
!  stop

end subroutine MQDTfunctions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MQDTNewfunctions(R1, R2, NA, scale, Cvals, mu, lwave, energy, &
     betavdw, phiL, phir, alpha0, alphaf,f0, g0, f0p, g0p, ldf0,ldg0)
  use units
  implicit none
  integer NA, i, j, lwave
  double precision Cvals(3), y(2), mu, betavdw, energy, f0, g0, f0p, g0p, phiL,phir
  double precision alpha0(2), alphaf(2), ldf0, ldg0
  double precision, allocatable :: phaseint(:)
  double precision, allocatable :: alpha(:,:),R(:)
  double precision R1, R2
  CHARACTER(LEN=*), INTENT(IN) :: scale

  allocate(R(NA),alpha(NA,2),phaseint(NA))
  call GridMaker(R,NA,R1,R2,scale)
  alpha(1,:) = alpha0
  
  !call CalcNewMilne2step(R,alpha,y,NA,energy,lwave,mu,Cvals,phaseint)
  call CalcNewMilne4step(R,alpha,y,NA,energy,lwave,mu,Cvals,phaseint)
!  write(6,*) "y = ", y
  alphaf = alpha(NA,:)
  phir = phaseint(NA)
  
  f0 = Pi**(-0.5d0)*alpha(NA,1)*sin(phir+phiL)
  g0 = -Pi**(-0.5d0)*alpha(NA,1)*cos(phir+phiL)
  f0p = alpha(NA,2)/alpha(NA,1) * f0 - g0/alpha(NA,1)**2
  g0p = alpha(NA,2)/alpha(NA,1) * g0 + f0/alpha(NA,1)**2

  ldf0 = y(2) + exp(-2*y(1))/tan(phir+phiL)
  ldg0 = y(2) - exp(-2*y(1))*tan(phir+phiL)
  
end subroutine MQDTNewfunctions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the phase standardization according to the Ruzic scheme.
subroutine CalcPhaseStandard(RX,RF,NXF,lwave,mu,betavdw,Cvals,phiL,scale)
  use units
  implicit none
!  double precision rj, ry, rjp, ryp
  double precision, intent(in) :: RX, RF, betavdw, mu, Cvals(3)
  double precision, intent(out) :: phiL
  double precision chim, chimp, energy, phir, ldg0, ldf0, ldchim
  double precision   f0, g0, f0p, g0p, tanphi, alpha0(2), alphaf(2)
  integer, intent(in) :: lwave
  integer i, NXF
  CHARACTER(LEN=*), INTENT(IN) :: scale
  double precision, external :: ksqrLR, VLR, VLRPrime

  !-------------------------------------------------------
  !uncomment this block to test the phase standardization
  !agreement with the analytical values for the C6 potential
!!$  Cvals(1) = 1d0
!!$  Cvals(2) = 0d0
!!$  Cvals(3) = 0d0
!!$  betavdw = 1d0
!!$  mu = 0.5d0
!!$  RX=0.1d0*betavdw
!!$  RF=4.0d0*betavdw
  !--------------------------------------------------------
  energy = 0d0

  alpha0(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
  alpha0(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
       /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)

  call MQDTNewfunctions(RX, RF, NXF, scale, Cvals, mu, lwave, energy, &
     betavdw, 0d0, phir, alpha0, alphaf,f0, g0, f0p, g0p, ldf0,ldg0)

  call chiminus(lwave,RF/betavdw,chim,chimp) !this requires van der Waals units, so divide by betavdw since x(i) is in bohr
  chimp = chimp/betavdw ! d(chi)/dR comes out in van der Waals units so this converts length to bohr
  ldchim = chimp/chim
  tanphi = (ldg0-ldchim)/(ldf0-ldchim)/tan(phir)
  !     write(200,*) X(i), tanphi- tan(-1d0/(2d0*RX**2) + dble(lwave)*Pi/4d0 -5d0*Pi/8d0)!, atan(tanphi)  

  write(6,*) "lwave = ", lwave, "tanphi = ", tanphi, "analytical result = ", &
       tan(-1d0/(2d0*RX**2) + dble(lwave)*Pi/4d0 -5d0*Pi/8d0), &
       "difference = ", tanphi- tan(-1d0/(2d0*RX**2) + dble(lwave)*Pi/4d0 -5d0*Pi/8d0)
  phiL = atan(tanphi)

end subroutine CalcPhaseStandard
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcKsr(ym, Ksr, size2,NXM, RX, RM, energy, Eth, lwave, mu, Cvals, betavdw, phiL,scale)
  implicit none
  integer size2, lwave, i, NXM
  double precision mu, betavdw, phiL, RX, RM, energy,alpha0(2),alpham(2),alphaf(2),f0,g0,f0p,g0p
  double precision ym(size2,size2), Eth(size2), Cvals(3),phir,ldg0,ldf0
  double precision f0mat(size2,size2), g0mat(size2,size2), f0pmat(size2,size2),g0pmat(size2,size2), Ksr(size2,size2)
  double precision temp1(size2,size2), temp2(size2,size2)
  double precision, external :: ksqrLR, VLR, VLRPrime
  CHARACTER(LEN=*), INTENT(IN) :: scale
  
  ! Ksr = (Y g - g')^(-1) (Y f - f')
  alpha0(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
  alpha0(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
       /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)
  f0mat = 0d0
  g0mat = 0d0
  f0pmat = 0d0
  g0pmat = 0d0
  !NXM = 1000000
  ! Construct the diagonal reference function matrices
  do i = 1, size2
!     call MQDTfunctions(RX, RM, NXM,'quadratic', Cvals, mu, lwave, energy-Eth(i), betavdw, phiL, alpha0, alphaf, f0, g0, f0p, g0p)
     call MQDTNewfunctions(RX, RM, NXM, scale, Cvals, mu, lwave, energy-Eth(i), &
          betavdw,phiL,phir, alpha0, alphaf,f0, g0, f0p, g0p,ldf0,ldg0)
     f0mat(i,i) = f0
     f0pmat(i,i) = f0p
     g0mat(i,i) = g0
     g0pmat(i,i) = g0p
  enddo
  
  temp1 = MATMUL(ym,g0mat) - g0pmat
  temp2 = MATMUL(ym,f0mat) - f0pmat
  call sqrmatinv(temp1,size2)
  Ksr = MATMUL(temp1,temp2)
  
end subroutine CalcKsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcTanGamma(RX,RF,NXF,size,lwave,mu,betavdw,Cvals,phiL,Egrid,NEgrid,Eth,iE,cotgamma,InterpCotGamma,scale)
  !This calculates the tan(gamma) MQDT parameter.  See Ruzic's PRA for details
  ! While in principle this can be calculated once as a function of energy, interpolated and
  ! evaluated at the relative energy for each channel, here I'm doing it the "brute force" way
  ! I simply re-calculate tan(gamma) and setup the cot(gamma) matrix for each channel at each needed energy
  use InterpType
  implicit none
  integer lwave, NEgrid, iE, ix, Nx,size,ith,NXF
  double precision RX, RF, mu, betavdw, Cvals(3), phiL, Egrid(NEgrid), cotgamma(size-1,size-1,NEgrid),EvdW
  double precision x, xscale, si, sk, sip, skp, ldk, ldi, kappa, f0, f0p, g0, g0p,Eth(size)
  double precision alpha0(2),alphaf(2),energy, tangamma, gam, phir, ldg0, ldf0 ,tp
  double precision, allocatable :: xgrid(:), Rgrid(:)
  double precision, allocatable :: phaseint(:)
  double precision, external :: ksqrLR, VLR, VLRPrime
  type(InterpolatingFunction) :: InterpCotGamma
  CHARACTER(LEN=*), INTENT(IN) :: scale
  
  EvdW = 1d0/(2d0*mu*betavdw**2)
  alpha0(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
  alpha0(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
       /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)
  
  xscale=0d0
  Nx = 200000
  cotgamma = 0d0

  do ith = 2, size
     ! For gamma, the energy represents a binding, so     
     energy = Egrid(iE) - (Eth(ith) - Eth(1))
     !------------------------------------------
     ! This part does the calculation of cot(gamma) from scratch
     kappa = sqrt(2*mu*abs(energy))
     alpha0(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
     alpha0(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
          /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)
     call MQDTNewfunctions(RX, RF, NXF, scale, Cvals, mu, lwave, energy, &
          betavdw,phiL,phir, alpha0, alphaf,f0, g0, f0p, g0p,ldf0,ldg0)
     x = kappa*RF
     
     tp = tan(phir+phiL)
     call Mysphbesik(lwave,x,xscale,si,sk,sip,skp,ldk,ldi) ! change norm

     cotgamma(ith-1,ith-1,iE) =tan(phir+phiL) * (ldf0 - kappa*ldk)/(ldg0 - kappa*ldk)! 1d0/tangamma

     !-----------------------------------------
     ! This part uses the interpolated function
!     cotgamma(ith-1,ith-1,iE) = Interpolated(energy, InterpCotGamma)
!     write(6,*) "Compare full: ", energy, cotgamma(ith-1,ith-1,iE)!,
!     write(6,*) "Compare interp: ", ith, Interpolated(energy, InterpCotGamma)
!     write(6,*) "Compare: ", ith, energy, cotgamma(ith-1, ith-1, iE), Interpolated(energy, InterpCotGamma)
     
     !-----------------------------------------
     !     write(400,*) xgrid(ix), ldk
     !write(400,*) energy/EvdW, gam(iE)
  enddo
 

  
end subroutine CalcTanGamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcCotGammaFunction(RX,RF,NXF,size,lwave,mu,betavdw,Cvals,phiL,Eth,Emin,Emax,InterpCotGamma,scale)
  !This calculates the tan(gamma) MQDT parameter.  See Ruzic's PRA for details
  ! While in principle this can be calculated once as a function of energy, interpolated and
  ! evaluated at the relative energy for each channel, here I'm doing it the "brute force" way
  ! I simply re-calculate tan(gamma) and setup the cot(gamma) matrix for each channel at each needed energy
  use InterpType
  use bspline
  use units
  implicit none
  integer lwave, iE, NE,ix, NXF,size, kx
  double precision RX, RF, mu, betavdw, Cvals(3), phiL,Emin,Emax,EvdW,E1,E2
  double precision x, xscale, si, sk, sip, skp, ldk, ldi, kappa, f0, f0p, g0, g0p,Eth(size)
  double precision alpha0(2),alphaf(2),energy, tangamma, gam,phir,ldg0,ldf0
  double precision, allocatable :: xgrid(:), Rgrid(:)
  double precision, allocatable :: phaseint(:)
  double precision, external :: ksqrLR, VLR, VLRPrime
  CHARACTER(LEN=*), INTENT(IN) :: scale
  type(InterpolatingFunction) :: InterpCotGamma

  EvdW = 1d0/(2d0*mu*betavdw**2)
  alpha0(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
  alpha0(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
       /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)

!  write(6,*) "alpha0 = ", alpha0
  xscale=0d0
  NE = 100
  kx=2

  E1 =  - (Eth(size)-Eth(1))
  E2 = -1d-12

  write(6,*) "E1, E2 = ",E1, E2
  call AllocateInterpolatingFunction(NE,kx,InterpCotGamma)
  call GridMaker(InterpCotGamma%x, nE, E1, E2, "linear")
  write(6,'(A)') "Calculating the MQDT parameter cot(gamma) as a function of energy"

  do iE = 1, NE

     energy = InterpCotGamma%x(iE)

     kappa = sqrt(2*mu*abs(energy))
     alpha0(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
     alpha0(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
          /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)
     call MQDTfunctions(RX, RF, NXF,scale, Cvals, mu, lwave, energy, betavdw, phiL, alpha0, alphaf, f0, g0, f0p, g0p)
     
     x = kappa*RF
!     write(6,*) "alphaf = ", alphaf
     call Mysphbesik(lwave,x,xscale,si,sk,sip,skp,ldk,ldi) ! change norm
     tangamma = -(g0p - kappa*ldk*g0)/(f0p - kappa*ldk*f0)
     InterpCotGamma%y(iE) = 1d0/tangamma
     write(6,'(A)',ADVANCE='NO') "."
 !    write(6,'(10f12.5)') energy, ldk, g0, g0p, f0, f0p, kappa, atan(tangamma)
  !   write(13,*) energy, atan(tangamma)
  enddo
  write(6,*) "done."
  call SetupInterpolatingFunction(InterpCotGamma)  
!  stop  
end subroutine CalcCotGammaFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcNewCotGammaFunction(RX,RF,NXF,size,lwave,mu,betavdw,Cvals,phiL,Eth,Emin,Emax,InterpCotGamma,scale)
  !This calculates the tan(gamma) MQDT parameter.  See Ruzic's PRA for details
  ! While in principle this can be calculated once as a function of energy, interpolated and
  ! evaluated at the relative energy for each channel, here I'm doing it the "brute force" way
  ! I simply re-calculate tan(gamma) and setup the cot(gamma) matrix for each channel at each needed energy
  use InterpType
  use bspline
  use units
  implicit none
  integer lwave, iE, NE,ix, NXF,size, kx
  double precision RX, RF, mu, betavdw, Cvals(3), phiL,Emin,Emax,EvdW,E1,E2,ldg0,ldf0,phir,tp
  double precision x, xscale, si, sk, sip, skp, ldk, ldi, kappa, f0, f0p, g0, g0p,Eth(size)
  double precision alpha0(2),alphaf(2),energy, tangamma, gam
  double precision, allocatable :: xgrid(:), Rgrid(:)
  CHARACTER(LEN=*), INTENT(IN) :: scale
  double precision, external :: ksqrLR, VLR, VLRPrime
  type(InterpolatingFunction) :: InterpCotGamma

  EvdW = 1d0/(2d0*mu*betavdw**2)

  xscale=0d0
  NE = 500
  kx=4

  E1 =  - (Eth(size)-Eth(1))
  E2 = -1d-12

  write(6,*) "E1, E2 = ",E1, E2
  call AllocateInterpolatingFunction(NE,kx,InterpCotGamma)
  call GridMaker(InterpCotGamma%x, nE, E1, E2, "linear")
  write(6,'(A)') "Calculating the MQDT parameter cot(gamma) as a function of energy"

  do iE = 1, NE
     energy = InterpCotGamma%x(iE)
     kappa = sqrt(2*mu*abs(energy))
     alpha0(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
     alpha0(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
          /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)
     call MQDTNewfunctions(RX, RF, NXF, scale, Cvals, mu, lwave, energy, &
          betavdw,phiL,phir, alpha0, alphaf,f0, g0, f0p, g0p,ldf0,ldg0)
!     call MQDTfunctions(RX, RF, NXF,'quadratic', Cvals, mu, lwave, energy, betavdw, phiL, alpha0, alphaf, f0, g0, f0p, g0p)     
     x = kappa*RF
!     write(6,*) "alphaf = ", alphaf
     call Mysphbesik(lwave,x,xscale,si,sk,sip,skp,ldk,ldi) ! change norm
!     tangamma = -(g0p - kappa*ldk*g0)/(f0p - kappa*ldk*f0)
!     InterpCotGamma%y(iE) = 1d0/tangamma
     tp = tan(phir+phiL)
     InterpCotGamma%y(iE) = tan(phir+phiL) * (ldf0 - kappa*ldk)/(ldg0 - kappa*ldk)
     write(6,'(A)',ADVANCE='NO') "."
     !write(6,'(10d15.5)') energy, ldk, ldf0, ldg0, phir
  !   write(13,*) energy, atan(tangamma)
  enddo
  write(6,*) "done."
  call SetupInterpolatingFunction(InterpCotGamma)  
  !stop  
end subroutine CalcNewCotGammaFunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine to setup the effect adiabatic potential energy function for any
!one of the electronic states for the alkali homonuclear molecules:
! for ground state singlet:
!   87-Rb2 [ISTATE = 1; ESTATE = 1]   from [Strauss et al.,] 
!   85-Rb2 [ISTATE = 2; ESTATE = 1]   from [Strauss et al.,]
!   39-K2 [ISTATE = 3; ESTATE = 1]    from [Falke et al., Phys. Rev. A 78, 012503 (2008)]
!   40-K2 [ISTATE = 4; ESTATE = 1]    from [Falke et al., Phys. Rev. A 78, 012503 (2008)]
!   23-Na2 [ISTATE = 5; ESTATE = 1]   from [Knoop et al., Phys. Rev. A 83, 042704 (2011)]
!   6-Li2 [ISTATE = 6; ESTATE = 1]    from LeRoy POTGENLI2.f
!   7-Li2 [ISTATE = 7; ESTATE = 1]    from LeRoy POTGENLI2.f
! 133-Cs2 [ISTATE = 8; ESTATE = 1]    from [J Baldwin, MSc. Thesis from University of Waterloo (2012)]
!
! for ground state triplet:
!   87-Rb2 [ISTATE = 1; ESTATE = 3]    from [Strauss et al.,]
!   85-Rb2 [ISTATE = 2; ESTATE = 3]    from [Strauss et al., ]
!   39-K2 [ISTATE = 3; ESTATE = 3]     from [Falke et al., Phys. Rev. a 78, 012503 (2008)]
!   40-K2 [ISTATE = 4; ESTATE = 3]     from [Falke et al., Phys. Rev. a 78, 012503 (2008)]
!   23-Na2 [ISTATE = 5; ESTATE = 3]    from [Knoop et al., Phys. Rev. A 83, 042704 (2011)]
!   6-Li2 [ISTATE = 6; ESTATE = 3]    from LeRoy POTGENLI2.f
!   7-Li2 [ISTATE = 7; ESTATE = 3]    from LeRoy POTGENLI2.f
! 133-Cs2 [ISTATE = 8; ESTATE = 3]    from [J Baldwin, MSc. Thesis from University of Waterloo (2012)]
!
! ON INPUT:
!  *integer ISTATE specifies the choice of the atomic species
!  *integer ESTATE specifies the choice of the electronic state (1 = Singlet, 3 = Triplet)
!  *MU is the reduced mass
!  *MUREF is the reduced mass of the reference isotopologue (used for BOB correction terms)
!     for Rb2, 87 is the reference
!     for K2, 39 is the reference
!  *integer NPP is the number of radial points
!  *X0(i) is an array of dimension NPP with the radial distances (in Angstrom) at which
!   the potential energy function is to be calculated
!  *VLIM is the (externally specified) asymptotic value of the potential energy (in cm-1)
!
! ON OUTPUT:
!  *VV(i) is an array of the potential function values (in cm-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetupPotential(ISTATE, ESTATE, MU, MUREF, NPP, VLIM, XO, VV, Cvals)
  use units
  implicit none
  integer NPP, ISTATE, ESTATE, i, j, N, IMN1, IMN2, iphi
  double precision, intent(in) :: XO(NPP)
  double precision VV(NPP), RM2(NPP), Cvals(3)
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

  if((ISTATE.NE.6).AND.(ISTATE.NE.7))then
     if (ISTATE.EQ.8) then
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

     else  
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

  VV = VV*InvcmPerHartree

END SUBROUTINE SetupPotential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for retrieving the nuclear spin, electronic spin, hyperfine constant (in MHz),
! nuclear gyromagnetic ratio, electronic gyromagnetic ratio, mass (in atomic units),
! for the selected atomic species:
!   87-Rb [ISTATE = 1]
!   85-Rb [ISTATE = 2]
!   39-K [ISTATE = 3]   
!   40-K [ISTATE = 4]    
!   23-Na [ISTATE = 5]
!    6-Li [ISTATE = 6]
!    7-Li [ISTATE = 7]
!  133-Cs [ISTATE = 8]
!
! The reduced mass (MU) and reduced mass for the reference isotopologue (MUREF) are for the
! homonuclear molecules 
!
! Electronic (s) and nuclear (i) spins are multiplied by two
! Nuclear g-factors need to be multiplied by the Bohr magneton
! Nuclear g-factors and hyperfine constants are from [Arimondo et al. Rev. Mod. Phys. 49, 31 (1977)]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine AtomData (ISTATE, AHf, i, s, gi, MU, MUREF, mass)
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

  double precision beta4,beta, C6g, C8g, C10g, Cvals(3), mu
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE polint(xa,ya,n,x,y,dy)
  INTEGER n,NMAX
  double precision dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER i,m,ns
  double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=abs(x-xa(1))
  do  i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do  i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) then
           write(6,*) "failure in polint.  Stopping."
           stop
        endif
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  return
END SUBROUTINE polint
SUBROUTINE ratint(xa,ya,n,x,y,dy)
  INTEGER n,NMAX
  double precision dy,x,y,xa(n),ya(n),TINY
  PARAMETER (NMAX=10,TINY=1.d-25)
  INTEGER i,m,ns
  double precision dd,h,hh,t,w,c(NMAX),d(NMAX)
  ns=1
  hh=abs(x-xa(1))
  do i=1,n
     h=abs(x-xa(i))
     if (h.eq.0.)then
        y=ya(i)
        dy=0.0
        return
     else if (h.lt.hh) then
        ns=i
        hh=h
     endif
     c(i)=ya(i)
     d(i)=ya(i)+TINY
  enddo
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        w=c(i+1)-d(i)
        h=xa(i+m)-x
        t=(xa(i)-x)*d(i)/h
        dd=t-c(i+1)
        if(dd.eq.0.) then
           write(6,*) "failure in ratint.  Stopping."
           stop
        endif
        dd=w/dd
        d(i)=c(i+1)*dd
        c(i)=t*dd
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  return
END SUBROUTINE ratint

double precision function abar(L)
  integer L
  double precision, parameter :: Pi = 3.14159265358979323846264338328

  abar=Pi**(2d0/(1d0 + 2d0*L))*(2**(-1.5d0 - 2d0*L)/ &
       (Gamma(1.25d0 + L/2.d0)*Gamma(0.5d0 + L)))**(2d0/(1d0 + 2d0*L))

end function abar

subroutine logderQD(lwave,mu,energy,NXM,Nsr,wsr,VSINGLET,VTRIPLET,R,TripletQD,SingletQD,phiL,betavdw,RX,Cvals,scale)
    use units
    implicit none
    integer n, n1, n2, Nsr, NXM,lwave
    double precision norm, k, energy,td1,td2,TripletQD,SingletQD, Ustart,psistart,s,mu,RX,ldf1,ldg1
    double precision R(Nsr),wsr(Nsr),VSINGLET(Nsr), VTRIPLET(Nsr),ystart(2,2),VMAT(2,2,Nsr)
    double precision alphax(2), alpham1(2), alpham2(2), phiL,betavdw,Cvals(3),identity(2,2),yf(2,2)!
    double precision, allocatable :: U(:)
    double precision f1,g1,f1p,g1p,f2,g2,f2p,g2p,RM,psi1,psi2,phir
    double precision, external :: VLR, VLRPrime
    CHARACTER(LEN=*), INTENT(IN) :: scale
    
  alphax(1) = (-((lwave + lwave**2 - 2*energy*mu*RX**2 + 2*mu*RX**2*VLR(mu,lwave,RX,Cvals))/RX**2))**(-0.25d0)
  alphax(2) = -(mu*((lwave*(1 + lwave))/(mu*RX**3) - VLRPrime(mu, lwave, RX,Cvals))) &
       /(4.d0*2**0.25d0*(energy*mu - (lwave*(1 + lwave))/(2.d0*RX**2) - mu*VLR(mu,lwave,RX,Cvals))**1.25d0)

  RM = R(Nsr)

  !call MQDTfunctions(RX, RM, NXM,'quadratic', Cvals, mu, lwave, energy, betavdw, phiL, alphax, alpham1, f1, g1, f1p, g1p)
  call MQDTNewfunctions(RX, RM, NXM, scale, Cvals, mu, lwave, energy, betavdw, phiL, phir, alphax, alpham1, &
       f1, g1, f1p, g1p, ldf1, ldg1)
  write(6,*) f1, g1, f1p, g1p

  identity = 0d0
  identity(1,1) = 1d0
  identity(2,2) = 1d0
  ystart = 0d0
  ystart(1,1) = 1d20
  ystart(2,2) = 1d20
  VMAT(:,:,:) = 0d0
  do n = 1, Nsr
     VMAT(1,1,n) = VSINGLET(n)
     VMAT(2,2,n) = VTRIPLET(n)
  enddo

  call logderpropA(mu,energy,identity,wSR,Nsr,ystart,yf,R,VMAT,2)

  td1 = (f1p - yf(1,1)*f1) / (g1p - yf(1,1)*g1)
  td2 = (f1p - yf(2,2)*f1) / (g1p - yf(2,2)*g1)

  SingletQD = atan(td1)/pi
  TripletQD = atan(td2)/pi
  write(6,*) "Quantum Defects:"
  write(6,*) "----------------"
  write(6,*) "Singlet:", SingletQD
  write(6,*) "Triplet:", TripletQD
  write(6,*)
!  stop
end subroutine LogderQD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testinterpolation
  use InterpType
  use bspline
  implicit none
  integer nx, kx
  double precision x
  type(InterpolatingFunction) :: finterpolant

  nx = 100
  kx = 3
  call AllocateInterpolatingFunction(nx,kx,finterpolant)
  call GridMaker(finterpolant%x, nx, 0.0d0, 100d0, "quadratic")
  finterpolant%y = finterpolant%x**2
  call SetupInterpolatingFunction(finterpolant)

  x=0d0
  do while (x.lt.100d0)
     write(12,*) x, Interpolated(x,finterpolant)
     x = x + 0.1d0
  enddo
  
end subroutine testinterpolation

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
