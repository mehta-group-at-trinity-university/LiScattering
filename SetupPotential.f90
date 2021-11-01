      
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
