MODULE Quadrature
  IMPLICIT NONE
  INTEGER LegPoints
  DOUBLE PRECISION, ALLOCATABLE :: xLeg(:),wLeg(:)
  double precision wbode10(10)
  CHARACTER*64 LegendreFile
  double precision cka2, cka3, cka4, cka5, cka6
  double precision b21, b31, b32, b41, b42, b43, b51, b52, b53, b54, b61, b62, b63, b64, b65
  double precision c1, c3, c4, c6, c1s, c3s, c4s, c5s, c6s

CONTAINS
  SUBROUTINE GetGaussFactors(File,Points,x,w)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c     This subroutine retrieves the points and weights for
    !c      Gaussian Quadrature from a given file
    !c
    !c     Variables:
    !c      File		name of file containing points and
    !c			 weights
    !c      Points		number of points to be used for
    !c			quadrature
    !c      x		array of nodes
    !c      w		array of weights
    !c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    IMPLICIT NONE
    INTEGER Points
    DOUBLE PRECISION x(Points),w(Points)
    CHARACTER*64 File
    INTEGER i,tempPoints

    OPEN(unit=7,file=File(1:INDEX(File,' ')-1))
    DO i = 1,18
       READ(7,*)
    ENDDO
    READ(7,*) tempPoints
    DO WHILE (tempPoints .NE. Points)
       DO i = 1,tempPoints
          READ(7,*)
       ENDDO
       READ(7,*) tempPoints
    ENDDO
    DO i = 1,Points
       READ(7,*) x(i),w(i)
    ENDDO
    CLOSE(unit=7)
    RETURN
  END SUBROUTINE GetGaussFactors

  subroutine initCashKarp
  ! These are the coefficients needed for the imbedded RK method of Cash and Karp. See Numerical Recipes RKCK.f
    cka2 = 0.2d0
    cka3 = 0.3d0
    cka4 = 0.6d0
    cka5 = 1d0
    cka6 = 0.875d0

    b21 = 0.2d0
    b31 = 0.075d0
    b32 = 0.225d0
    b41 = 0.3d0
    b42 = -0.9d0
    b43 = 1.2d0
    b51 = -11d0/54d0 !-0.203703703703703703703703703704d0
    b52 = 2.5d0
    b53 = -70d0/27d0 !-2.59259259259259259259259259259d0
    b54 = 35d0/27d0 !1.29629629629629629629629629630d0
    b61 = 1631d0/55296d0 !0.0294958043981481481481481481481d0
    b62 = 175d0/512d0 !0.341796875d0
    b63 = 575d0/13828d0 !0.0415943287037037037037037037037d0
    b64 = 44275d0/110529d0 !0.400345413773148148148148148148
    b65 = 253d0/4096d0 !0.061767578125d0

    c1 = 37d0/378d0 !0.0978835978835978835978835978836d0
    c3 = 250d0/621d0 !0.402576489533011272141706924316d0
    c4 = 125d0/594d0 !0.210437710437710437710437710438d0
    c6 = 512d0/1771d0 !0.289102202145680406549971767363

    c1s = 2825d0/27648d0 !0.102177372685185185185185185185d0
    c3s = 18575d0/48384d0 !0.383907903439153439153439153439d0
    c4s = 13525d0/55296d0!0.244592737268518518518518518519d0
    c5s = 277d0/14366d0!0.0192816371989419462620075177502d0
    c6s = 0.25d0
    
  end subroutine initCashKarp
  subroutine initbode10

    wbode10(1) = 2857d0
    wbode10(10) = 2857d0
    wbode10(2) = 15741d0
    wbode10(9) = 15741d0
    wbode10(3) = 1080d0
    wbode10(8) = 1080d0
    wbode10(4) = 19344d0
    wbode10(7) = 19344d0
    wbode10(5) = 5778d0
    wbode10(6) = 5778d0

  end subroutine initbode10

  !***********************************************************************
END MODULE Quadrature
