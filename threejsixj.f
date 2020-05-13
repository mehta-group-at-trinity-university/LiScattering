********************************************************************************
C      gensub.f  -- standard subroutines, to calculate wigner coeffs,
c			interpolation, root-finding, numerical integration
c			Coulomb function programs also.
C**********************************************************C
C
C WIGNER 3J, 6J, AND 9J FUNCTION SUBPROGRAMS.
C NOTE:  ALL ANGULAR MOMENTUM QUANTUM NUMBERS SUPPLIED TO THESE
C        FUNCTIONS ARE INTEGERS WHICH ARE TWICE THE VALUE OF THE
C        ACTUAL ANGULAR MOMENTUM DESIRED.  (THIS ALLOWS FOR HALF-
C        INTEGER VALUES CONVENIENTLY.)  ALSO YOU MUST CALL SETUPAM
C        ONE TIME BEFORE CALLING ANY OF THESE SO THAT THE RELE-
C        VANT FACTORIALS CAN BE CALCULATED ONCE AND FOR ALL AND
C        STORED, AS IN THE ABOVE EXAMPLE.
C
C**********************************************************C
C

      function coeff(l1,l2,l1p,l2p,l,k)
	    implicit real*8(a-h,o-z)
			front=(2*l1+1)*(2*l2+1)*(2*l1p+1)*(2*l2p+1)
			front=dsqrt(front)*(-1)**(l1+l1p+l)
			l1d=2*l1
			l2d=2*l2
			l1pd=2*l1p
			l2pd=2*l2p
			ld=2*l
			kd=2*k
			iz=0
			t1=thrj(l1d,kd,l1pd,iz,iz,iz)
			t2=thrj(l2d,kd,l2pd,iz,iz,iz)
			s1=sixj(l1d,l2d,ld,l2pd,l1pd,kd)
			coeff=front*t1*t2*s1
			return
   		end
c
      FUNCTION XNINEJ(J11,J12,J13,J21,J22,J23,J31,J32,J33)
      IMPLICIT REAL*8(A-H,O-Z)
      KMIN1 = IABS(J11-J33)
      KMIN2 = IABS(J32-J21)
      KMIN3 = IABS(J23-J12)
      KMAX1 = J11+J33
      KMAX2 = J32+J21
      KMAX3 = J23+J12
C
      IF(KMIN2.GT.KMIN1) KMIN1=KMIN2
       IF(KMIN3.GT.KMIN1) KMIN1=KMIN3
      IF(KMAX2.LT.KMAX1) KMAX1=KMAX2
      IF(KMAX3.LT.KMAX1) KMAX1=KMAX3
C
      KMIN1 = KMIN1 + 1
      KMAX1 = KMAX1 + 1
       XNINEJ = 0.D0
      IF(KMIN1.GT.KMAX1) GO TO 1000
      DO 100 K1 = KMIN1,KMAX1,2
      K = K1 - 1
      S1 = SIXJ(J11,J21,J31,J32,J33,K)
      S2 = SIXJ(J12,J22,J32,J21,K,J23)
      S3 = SIXJ(J13,J23,J33,K,J11,J12)
      P = (K+1)*((-1)**K)
      XNINEJ = XNINEJ + P*S1*S2*S3
  100 CONTINUE
 1000 CONTINUE
      RETURN
      END
C
      FUNCTION THRJ(J1D,J2D,J3D,M1D,M2D,M3D)
C     Gives the Wigner 3-j symbol:
C     (j1 j2 j3)
C     (m1 m2 m3)
      IMPLICIT REAL*8(A-H,O-Z)
      X1 = J1D/2.D0
      X2 = J2D/2.D0
      X3 = J3D/2.D0
      Y1 = M1D/2.D0
      Y2 = M2D/2.D0
      Y3 = M3D/2.D0
C
C -- NEXT COME THE TRIANGULARITY CHECKS:
C
      IF(J1D+J2D-J3D.LT.0) GO TO 9998
      IF(J2D+J3D-J1D.LT.0) GO TO 9998
      IF(J3D+J1D-J2D.LT.0) GO TO 9998
      IF(J3D.LT.IABS(J1D-J2D)) GO TO 9998
      IF(M1D+M2D+M3D.NE.0) GO TO 9998
      LLL = J1D+J2D+J3D
      IF(2*(LLL/2) .NE. LLL) GO TO 9998
C
      KMIN = (J3D-J1D-M2D)/2
      KMIN1 = KMIN
      KMIN2 = (J3D-J2D+M1D)/2
      IF(KMIN2.LT.KMIN) KMIN=KMIN2
      KMIN = (-1)*KMIN
      KMAX = X1+X2-X3 +0.1D0
      KMAX1 = KMAX
      KMAX2 = X1-Y1
      KMAX3 = X2+Y2
      IF(KMAX2.LT.KMAX) KMAX=KMAX2
      IF(KMAX3.LT.KMAX) KMAX=KMAX3
      IF(KMIN.LT.0) KMIN = 0
      IF(KMIN.GT.KMAX) GO TO 9998
C
      JMIN = KMIN+1
      JMAX = KMAX+1
      TERM1 = FRONTL(X1,X2,X3,Y1,Y2,Y3)
	iphase=iabs((j1d-j2d-m3d)/2)
	msign=(-1)**iphase
cg     MSIGN = (-1)**((J1D-J2D-M3D)/2)
      SUM = 0.D0
      DO 10 I1 = JMIN,JMAX
      I = I1 - 1
      TERM2 = FL(I1)+FL(KMIN1+I1)+FL(KMIN2+I1)
      TERM2 = TERM2+FL(KMAX1-I+1)+FL(KMAX2-I+1)+FL(KMAX3-I+1)
      TERM= DEXP(TERM1-TERM2)
      TERM = TERM*MSIGN*((-1)**I)
      SUM = SUM + TERM
  10  CONTINUE
      THRJ = SUM
      GO TO 9999
 9998 THRJ = 0.D0
 9999 CONTINUE
      RETURN
      END
C
      FUNCTION FL(I)
       IMPLICIT REAL*8(A-H,O-Z)
CCC   DIMENSION FACL(60)
      COMMON/FACTOR/FACL(200)
      FL = FACL(I)
      RETURN
      END
C
C** ** **
C-- THIS SUBROUTINE INITIALIZES BY FINDING THE LOGARITHM
C---OF THE FIRST 199 FACTORIALS AND STORING THEM.
C
      SUBROUTINE SETUPAM
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/FACTOR/FACL(200)
      N = 199
      FACL(1) = 0.D0
      DO 100 I = 2,N
      I1 = I - 1
      FACL(I) = FACL(I1) + DLOG(I1*1.D0)
  100 CONTINUE
      RETURN
      END
C** ** **
C
      FUNCTION FRONTL(X1,X2,X3,Y1,Y2,Y3)
      IMPLICIT REAL*8(A-H,O-Z)
      L1 = X1+X2-X3 +1.1D0
      L2 = X2+X3-X1 +1.1D0
      L3 = X3+X1-X2 +1.1D0
      L4 = X1+X2+X3+1+1.1D0
      L5 = X1+Y1+1.1D0
      L6 = X1-Y1+1.1D0
      L7 = X2+Y2+1.1D0
      L8 = X2-Y2+1.1D0
      L9 = X3+Y3+1.1D0
      L10 = X3-Y3+1.1D0
      FRONTL = FL(L1)+FL(L2)+FL(L3)-FL(L4)+FL(L5)+FL(L6)
      FRONTL = FRONTL +FL(L7)+FL(L8)+FL(L9)+FL(L10)
      FRONTL = FRONTL/2.D0
      RETURN
      END
C
      FUNCTION SIXJ(J1D,J2D,J3D,J4D,J5D,J6D)
      IMPLICIT REAL*8(A-H,O-Z)
C
C -- CHECK THAT TRIANGULARITY CONDITIONS ARE SATISFIED.
C
      IF(J3D.LT.IABS(J1D-J2D) .OR. J3D.GT.J1D+J2D) GO TO 9998
      IF(J6D.LT.IABS(J4D-J2D) .OR. J6D.GT.J4D+J2D) GO TO 9998
      IF(J6D.LT.IABS(J1D-J5D) .OR. J6D.GT.J1D+J5D) GO TO 9998
      IF(J3D.LT.IABS(J4D-J5D) .OR. J3D.GT.J4D+J5D) GO TO 9998
      K1=J1D+J2D+J3D
      K2=J4D+J2D+J6D
      K3=J6D+J1D+J5D
      K4=J3D+J4D+J5D
      IF(2*(K1/2).NE.K1 .OR. 2*(K2/2).NE.K2) GO TO 9998
      IF(2*(K3/2).NE.K3 .OR. 2*(K4/2).NE.K4) GO TO 9998
C
C -- NOW GO AHEAD AND CALCULATE THE SIXJ.
C
      JM1 = (J1D+J2D+J3D)/2
      JM2 = (J1D+J5D+J6D)/2
      JM3 = (J4D+J2D+J6D)/2
      JM4 = (J4D+J5D+J3D)/2
      JX1 = (J1D+J2D+J4D+J5D)/2
      JX2 = (J2D+J3D+J5D+J6D)/2
      JX3 = (J3D+J1D+J6D+J4D)/2
C
      JM = JM1
      IF(JM2.GT.JM) JM = JM2
      IF(JM3.GT.JM) JM = JM3
      IF(JM4.GT.JM) JM = JM4
      JX = JX1
      IF(JX2.LT.JX) JX = JX2
      IF(JX3.LT.JX) JX = JX3
      KM = JM+1
      KX = JX+1
      IF(KM.GT.KX) GO TO 9998
      TERM1 = FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)
      SIXJ = 0.D0
      DO 10 I1 = KM,KX
      I = I1 - 1
      TERM2 = FL(I+2)-FL(I+1-JM1)-FL(I+1-JM2)-FL(I+1-JM3)
      TERM2 = TERM2-FL(I+1-JM4)-FL(JX1-I+1)-FL(JX2-I+1)
      TERM2 = TERM2 -FL(JX3-I+1)
      TERM = DEXP(TERM1+TERM2) * ((-1)**I)
      SIXJ = SIXJ + TERM
   10 CONTINUE
      GO TO 9999
 9998 CONTINUE
      SIXJ = 0.D0
 9999 CONTINUE
      RETURN
      END
C
      FUNCTION FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)
      IMPLICIT REAL*8(A-H,O-Z)
      FRTSJL = DL(J1D,J2D,J3D) + DL(J1D,J5D,J6D)
      FRTSJL = FRTSJL + DL(J4D,J2D,J6D) + DL(J4D,J5D,J3D)
      RETURN
      END
C
      FUNCTION DL(J1D,J2D,J3D)
      IMPLICIT REAL*8(A-H,O-Z)
      L1 = (J1D+J2D-J3D)/2
      L2 = (J2D+J3D-J1D)/2
      L3 = (J3D+J1D-J2D)/2
      L4 = (J1D+J2D+J3D)/2 + 1
      DL = FL(L1+1)+FL(L2+1)+FL(L3+1)-FL(L4+1)
      DL = DL/2.D0
      RETURN
      END
C
      FUNCTION CLEBSH(J1D,J2D,JD,M1D,M2D,MD)
      IMPLICIT REAL*8(A-H,O-Z)
      CLEBSH = 0.D0
      IF(M1D+M2D.NE.MD) GO TO 100
      MMD = -MD
      Q = THRJ(J1D,J2D,JD,M1D,M2D,MMD)
      PHASE = ((-1)**(20 + (J2D-J1D-MD)/2))
      CLEBSH = Q*PHASE*DSQRT(JD+1.D0)
  100 CONTINUE
      RETURN
      END
************************************************************************
        function yfix(x,vx,v,ndata,igo)
        implicit real*8(a-h,o-z)
c
c --- interpolation routine yfix.
c --- vx and v are input data arrays for the x and y data points,
c -----  ndata is the number of points in these arrays, and x is the
c -----  final x-value at which y=yfix is desired.  Note that the
c -----  vx-values must be strictly increasing.  The parameter igo
c -----  should be set to 0 on the first call at a given x, and
c -----  different from 0 on subsequent calls at the same x.
c
        dimension vx(ndata),v(ndata)
        nx = ndata - 1
        if(x.ge.vx(nx)) go to 75
        if(x.le.vx(2)) go to 65
        if(igo.ne.0) go to 20
        n1 = 2
        n2 = ndata - 1
        i = 0
   10   continue
        i =   i + 1
        j = (n1+n2)/2
        if(vx(j).le.x) go to 15
        if(n2-n1 .le. 1) go to 20
        n2 = j
        go  to 10
   15   continue
        n1 = j
        if(n2-n1 .le. 1) go to 20
        go to 10
   20   continue
        n = n2-1
        np = n2
        np1 = np + 1
        n1 = n - 1
        x0 = vx(n1)
        y0 = v(n1)
        x1 = vx(n)
        y1 = v(n)
        x2 = vx(np)
        y2 = v(np)
        x3 = vx(np1)
        y3 = v(np1)
        q0 = (x-x1)*(x-x2)*(x-x3)/( (x0-x1)*(x0-x2)*(x0-x3) )
        q1 = (x-x0)*(x-x2)*(x-x3)/( (x1-x0)*(x1-x2)*(x1-x3) )
        q2 = (x-x0)*(x-x1)*(x-x3)/( (x2-x0)*(x2-x1)*(x2-x3) )
        q3 = (x-x0)*(x-x1)*(x-x2)/( (x3-x0)*(x3-x1)*(x3-x2) )
        y = q0*y0 + q1*y1 + q2*y2 + q3*y3
        go to 100
   65   n1 = 2
        n2 = 3
        go to 20
   75   nx = ndata - 1
        slope = (v(ndata)-v(nx))/(vx(ndata)-vx(nx))
        y = v(ndata) - slope*(vx(ndata)-x)
  100   yfix = y
        return
        end
***********************************************************************
**********************************************************************
	function rint(f,na,nb,nq,h)
c
c -- function rint performs a numerical integration over the function f,
c -----  assuming an equally spaced x-axis grid with step size h and a
c -----  total of nb mesh points.  Use na=1, nq=10, and nb>20 always.
c
	implicit real*8(a-h,o-z)
	dimension c(55),d(10),f(nb)
	data c/1.d0,
     1  2.d0,1.d0,
     2  23.d0, 28.d0, 9.d0,
     3  25.d0, 20.d0, 31.d0, 8.d0,
     4  1413.d0, 1586.d0, 1104.d0, 1902.d0, 475.d0,
     5  1456.d0, 1333.d0, 1746.d0, 944.d0,1982.d0,459.d0,
     6  119585.d0, 130936.d0, 89437.d0, 177984.d0, 54851.d0,
     7  176648.d0, 36799.d0,
     8  122175.d0, 111080.d0, 156451.d0,46912.d0,220509.d0,
     9  29336.d0, 185153.d0, 35584.d0,
     a  7200319.d0, 7783754.d0, 5095890.d0,12489922.d0,-1020160.d0,
     b  16263486.d0, 261166.d0, 11532470.d0, 2082753.d0,
     c  7305728.d0, 6767167.d0, 9516362.d0,1053138.d0, 18554050.d0,
     d   -7084288.d0, 20306238.d0, -1471442.d0, 11965622.d0,
     e  2034625.d0/
	data d/  2.d0,    2.d0, 24.d0, 24.d0,1440.d0,1440.d0,120960.d0,
     &     120960.d0, 7257600.d0, 7257600.d0/
	a=0.d0
	l=na
	m=nb
	i=nq*(nq+1)/2
	do j=1,nq
           a=a+c(i)*(f(l)+f(m))
           l=l+1
           m=m-1
           i=i-1
        enddo
	a=a/d(nq)
	do n=l,m
           a=a+f(n)
        enddo
	rint=a*h
	return
	end
c***********************************************************************
c***********************************************************************
	function zrint(zf,na,nb,nq,h)
c
c -- function rint performs a numerical integration over the function f,
c -----  assuming an equally spaced x-axis grid with step size h and a
c -----  total of nb mesh points.  Use na=1, nq=10, and nb>20 always.
c
	implicit real*8(a-h,o-y),complex*16(z)
	dimension c(55),d(10),zf(nb)
	data c/1.d0,
     1  2.d0,1.d0,
     2  23.d0, 28.d0, 9.d0,
     3  25.d0, 20.d0, 31.d0, 8.d0,
     4  1413.d0, 1586.d0, 1104.d0, 1902.d0, 475.d0,
     5  1456.d0, 1333.d0, 1746.d0, 944.d0,1982.d0,459.d0,
     6  119585.d0, 130936.d0, 89437.d0, 177984.d0, 54851.d0,
     7  176648.d0, 36799.d0,
     8  122175.d0, 111080.d0, 156451.d0,46912.d0,220509.d0,
     9  29336.d0, 185153.d0, 35584.d0,
     a  7200319.d0, 7783754.d0, 5095890.d0,12489922.d0,-1020160.d0,
     b  16263486.d0, 261166.d0, 11532470.d0, 2082753.d0,
     c  7305728.d0, 6767167.d0, 9516362.d0,1053138.d0, 18554050.d0,
     d   -7084288.d0, 20306238.d0, -1471442.d0, 11965622.d0,
     e  2034625.d0/
	data d/  2.d0,    2.d0, 24.d0, 24.d0,1440.d0,1440.d0,120960.d0,
     &     120960.d0, 7257600.d0, 7257600.d0/
	za=dcmplx(0.d0,0.d0)
	l=na
	m=nb
	i=nq*(nq+1)/2
	do j=1,nq
           za=za+c(i)*(zf(l)+zf(m))
           l=l+1
           m=m-1
           i=i-1
        enddo
	za=za/d(nq)
	do n=l,m
           za=za+zf(n)
        enddo
	zrint=za*h
	return
	end
c***********************************************************************
c*****************************************************************************
        double precision function delt(i,j)
        implicit double precision (a-h,o-z)
        delt=0.d0
        if(i.eq.j) delt=1.d0
        return
        end
!NPM REMOVED FROM THIS FILE ROUTINES FOR CALCULATING QDT FUNCTIONS.  SEE ORIGINAL ZGENSUB FOR THOSE ROUTINES.
