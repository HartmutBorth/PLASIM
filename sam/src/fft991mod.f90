!     =============
!     MODULE FFTMOD
!     =============

!     alternate module with factors 2, 3, 4, 5, 6, 8
!     this is a FORTRAN-90 version of the old FFT991 package
!     use this for resolutions not supported by the <fftmod>,
!     e.g. T63 (NLAT=96), T106 (NLAT=160)

      module fftmod
      parameter(NRES = 9)
      integer :: nallowed(NRES) = (/ 64, 96, 128, 192, 256, 320, 384, 512, 1024 /)
!     T21  - N64  : 8-4-2
!     T31  - N96  : 8-4-3
!     T42  - N128 : 8-4-4
!     T63  - N192 : 8-6-4
!     T85  - N256 : 8-4-4-2
!     T106 - N320 : 8-5-4-2
!     T127 - N384 : 8-4-4-3
!     T170 - N512 : 8-4-4-4
!     T341 - N1024 : 8-4-4-4-2  ! Added by T.Frisius
      integer :: lastn = 0
      integer :: ifax(10)
      real,allocatable :: trigs(:)
      end module fftmod


!     ================
!     SUBROUTINE GP2FC
!     ================

      subroutine gp2fc(a,n,lot)
      use fftmod
      real :: a(n,lot)
      real :: worka(n+2,lot)
      real :: workb(n+2,lot)

      if (n /= lastn) then
         if (allocated(trigs)) deallocate(trigs)
         allocate(trigs(n))
         lastn = n
         call fftini(n)
      endif

      worka(1:n,:) = a(:,:)
      worka(n+1:n+2,:) = 0.0

      call fft991(worka,workb,trigs,ifax,1,n+2,n,lot,-1)

      a(:,:) = worka(1:n,:)

      return
      end subroutine gp2fc


!     ================
!     SUBROUTINE FC2GP
!     ================

      subroutine fc2gp(a,n,lot)
      use fftmod
      real :: a(n,lot)
      real :: worka(n+2,lot)
      real :: workb(n+2,lot)

      if (n /= lastn) then
         if (allocated(trigs)) deallocate(trigs)
         allocate(trigs(n))
         lastn = n
         call fftini(n)
      endif

      worka(1:n,:) = a(:,:)
      worka(n+1:n+2,:) = 0.0

      call fft991(worka,workb,trigs,ifax,1,n+2,n,lot,1)

      a(:,:) = worka(1:n,:)

      end subroutine fc2gp


!     =================
!     SUBROUTINE FFTINI
!     =================

      subroutine fftini(n)
      use fftmod
      logical labort

!     check for allowed values of n

      labort = .true.
      do j = 1 , NRES
         if (n == nallowed(j)) labort = .false.
      enddo

      if (labort) then
         write (*,*) '*** FFT does not support n = ',n,' ***'
         write (*,*) 'Following resolutions may be used:'
         write (*,*) '----------------------------------'
         do j = 1 , NRES
            write (*,1000) nallowed(j), nallowed(j)/2, nallowed(j)/3
         enddo
         stop
      endif
 1000 format(' NLON=',I5,'  NLAT=',I5,'  NTRU=',I5)

      call set99(trigs,ifax,n) ! Factorization of n and sine/cosine values

      return
      end subroutine fftini




      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!     SUBROUTINE FFT991 - MULTIPLE FAST REAL PERIODIC TRANSFORM
!     SUPERSEDES PREVIOUS ROUTINE 'FFT991'
!
!     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
!     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
!
!     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*MIN(LOT,64)
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS
!     IN PARALLEL
!
!     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
      DIMENSION A(*),WORK(*),TRIGS(N),IFAX(10)
!
      NFAX=IFAX(1)
      NX=N+1
      IF (MOD(N,2).EQ.1) NX=N
      NBLOX=1+(LOT-1)/64
      NVEX=LOT-(NBLOX-1)*64
      IF (ISIGN.EQ.-1) GO TO 300
!
!     ISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
!     -----------------------------------------
100   CONTINUE
      ISTART=1
      DO 220 NB=1,NBLOX
        IA=ISTART
        I=ISTART
        DO 110 J=1,NVEX
          A(I+INC)=0.5*A(I)
          I=I+JUMP
110       CONTINUE
        IF (MOD(N,2).EQ.1) GO TO 130
        I=ISTART+N*INC
        DO 120 J=1,NVEX
          A(I)=0.5*A(I)
          I=I+JUMP
120       CONTINUE
130     CONTINUE
        IA=ISTART+INC
        LA=1
        IGO=+1
!
        DO 160 K=1,NFAX
          IFAC=IFAX(K+1)
          IERR=-1
          IF (IGO.EQ.-1) GO TO 140
          CALL RPASSM(A(IA),A(IA+LA*INC),WORK(1),WORK(IFAC*LA+1),TRIGS, &
                   INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
          GO TO 150
140       CONTINUE
          CALL RPASSM(WORK(1),WORK(LA+1),A(IA),A(IA+IFAC*LA*INC),TRIGS, &
                   1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
150       CONTINUE
          IF (IERR.NE.0) GO TO 500
          LA=IFAC*LA
          IGO=-IGO
          IA=ISTART
160       CONTINUE
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (MOD(NFAX,2).EQ.0) GO TO 190
        IBASE=1
        JBASE=IA
        DO 180 JJ=1,NVEX
          I=IBASE
          J=JBASE
          DO 170 II=1,N
            A(J)=WORK(I)
            I=I+1
            J=J+INC
170         CONTINUE
          IBASE=IBASE+NX
          JBASE=JBASE+JUMP
180       CONTINUE
190     CONTINUE
!
!     FILL IN ZEROS AT END
!     --------------------
        IX=ISTART+N*INC
        DO 210 J=1,NVEX
          A(IX)=0.0
          A(IX+INC)=0.0
          IX=IX+JUMP
210       CONTINUE
!
        ISTART=ISTART+NVEX*JUMP
        NVEX=64
220     CONTINUE
      RETURN
!
!     ISIGN=-1, GRIDPOINT TO SPECTRAL TRANSFORM
!     -----------------------------------------
300   CONTINUE
      ISTART=1
      DO 410 NB=1,NBLOX
        IA=ISTART
        LA=N
        IGO=+1
!
        DO 340 K=1,NFAX
          IFAC=IFAX(NFAX+2-K)
          LA=LA/IFAC
          IERR=-1
          IF (IGO.EQ.-1) GO TO 320
          CALL QPASSM(A(IA),A(IA+IFAC*LA*INC),WORK(1),WORK(LA+1),TRIGS, &
                   INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
          GO TO 330
320       CONTINUE
          CALL QPASSM(WORK(1),WORK(IFAC*LA+1),A(IA),A(IA+LA*INC),TRIGS, &
                   1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
330       CONTINUE
          IF (IERR.NE.0) GO TO 500
          IGO=-IGO
          IA=ISTART+INC
340       CONTINUE
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (MOD(NFAX,2).EQ.0) GO TO 370
        IBASE=1
        JBASE=IA
        DO 360 JJ=1,NVEX
          I=IBASE
          J=JBASE
          DO 350 II=1,N
            A(J)=WORK(I)
            I=I+1
            J=J+INC
350         CONTINUE
          IBASE=IBASE+NX
          JBASE=JBASE+JUMP
360       CONTINUE
370     CONTINUE
!
!     SHIFT A(0) & FILL IN ZERO IMAG PARTS
!     ------------------------------------
        IX=ISTART
        DO 380 J=1,NVEX
          A(IX)=A(IX+INC)
          A(IX+INC)=0.0
          IX=IX+JUMP
380       CONTINUE
        IF (MOD(N,2).EQ.1) GO TO 400
        IZ=ISTART+(N+1)*INC
        DO 390 J=1,NVEX
          A(IZ)=0.0
          IZ=IZ+JUMP
390       CONTINUE
400     CONTINUE
!
        ISTART=ISTART+NVEX*JUMP
        NVEX=64
410     CONTINUE
      RETURN
!
!     ERROR MESSAGES
!     --------------
500   CONTINUE
      GO TO (510,530,550) IERR
510   CONTINUE
      WRITE(6,520) NVEX
520   FORMAT(16H1VECTOR LENGTH =,I4,17H, GREATER THAN 64)
      GO TO 570
530   CONTINUE
      WRITE(6,540) IFAC
540   FORMAT( 9H1FACTOR =,I3,17H, NOT CATERED FOR)
      GO TO 570
550   CONTINUE
      WRITE(6,560) IFAC
560   FORMAT(9H1FACTOR =,I3,31H, ONLY CATERED FOR IF LA*IFAC=N)
570   CONTINUE
      RETURN
      END



      SUBROUTINE SET99(TRIGS,IFAX,N)
!     SUBROUTINE SET99 - COMPUTES FACTORS OF N & TRIGONOMETRIC
!     FUNCTINS REQUIRED BY FFT99 & FFT991
!
      DIMENSION TRIGS(N),IFAX(10),JFAX(10),LFAX(7)
      DATA LFAX/6,8,5,4,3,2,1/
      IXXX=1
!
      DEL=4.0*ASIN(1.0)/FLOAT(N)
      NIL=0
      NHL=(N/2)-1
      DO 10 K=NIL,NHL
        ANGLE=FLOAT(K)*DEL
        TRIGS(2*K+1)=COS(ANGLE)
        TRIGS(2*K+2)=SIN(ANGLE)
10      CONTINUE
!
!     FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
!     LOOK FOR SIXES FIRST, STORE FACTORS IN DESCENDING ORDER
      NU=N
      IFAC=6
      K=0
      L=1
20    CONTINUE
      IF (MOD(NU,IFAC).NE.0) GO TO 30
      K=K+1
      JFAX(K)=IFAC
      IF (IFAC.NE.8) GO TO 25
      IF (K.EQ.1) GO TO 25
      JFAX(1)=8
      JFAX(K)=6
25    CONTINUE
      NU=NU/IFAC
      IF (NU.EQ.1) GO TO 50
      IF (IFAC.NE.8) GO TO 20
30    CONTINUE
      L=L+1
      IFAC=LFAX(L)
      IF (IFAC.GT.1) GO TO 20
!
      WRITE(6,40) N
40    FORMAT(4H1N =,I4,27H - CONTAINS ILLEGAL FACTORS)
      RETURN
!
!     NOW REVERSE ORDER OF FACTORS
50    CONTINUE
      NFAX=K
      IFAX(1)=NFAX
      DO 60 I=1,NFAX
        IFAX(NFAX+2-I)=JFAX(I)
60      CONTINUE
      IFAX(10)=N
      RETURN
      END




      SUBROUTINE QPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA,IERR)
!     SUBROUTINE QPASSM - PERFORMS ONE PASS THROUGH DATA AS PART
!     OF MULTIPLE REAL FFT (FOURIER ANALYSIS) ROUTINE
!
!     A IS FIRST REAL INPUT VECTOR
!         EQUIVALENCE B(1) WITH A(IFAC*LA*INC1+1)
!     C IS FIRST REAL OUTPUT VECTOR
!         EQUIVALENCE D(1) WITH C(LA*INC2+1)
!     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
!     INC1 IS THE ADDRESSING INCREMENT FOR A
!     INC2 IS THE ADDRESSING INCREMENT FOR C
!     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
!     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
!     LOT IS THE NUMBER OF VECTORS
!     N IS THE LENGTH OF THE VECTORS
!     IFAC IS THE CURRENT FACTOR OF N
!     LA = N/(PRODUCT OF FACTORS USED SO FAR)
!     IERR IS AN ERROR INDICATOR:
!              0 - PASS COMPLETED WITHOUT ERROR
!              1 - LOT GREATER THAN 64
!              2 - IFAC NOT CATERED FOR
!              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
!
!-----------------------------------------------------------------------
!
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(N)
!
      parameter(SIN36 = 0.587785252292473)
      parameter(SIN72 = 0.951056516295154)
      parameter(QRT5  = 0.559016994374947)
      parameter(SIN60 = 0.866025403784437)
!
      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      IJUMP=(IFAC-1)*IINK
      KSTOP=(N-IFAC)/(2*IFAC)
!
      IBAD=1
      IF (LOT.GT.64) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO
!
!     CODING FOR FACTOR 2
!     -------------------
200   CONTINUE
      IA=1
      IB=IA+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
!
      IF (LA.EQ.M) GO TO 290
!
      DO 220 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 210 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=A(IA+I)-A(IB+I)
          I=I+INC3
          J=J+INC4
210       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
220     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JB) GO TO 260
      DO 250 K=LA,KSTOP,LA
        KB=K+K
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        JBASE=0
        DO 240 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 230 IJK=1,LOT
            C(JA+J)=A(IA+I)+(C1*A(IB+I)+S1*B(IB+I))
            C(JB+J)=A(IA+I)-(C1*A(IB+I)+S1*B(IB+I))
            D(JA+J)=(C1*B(IB+I)-S1*A(IB+I))+B(IA+I)
            D(JB+J)=(C1*B(IB+I)-S1*A(IB+I))-B(IA+I)
            I=I+INC3
            J=J+INC4
230         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
240       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB-JINK
250     CONTINUE
      IF (JA.GT.JB) GO TO 900
260   CONTINUE
      JBASE=0
      DO 280 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 270 IJK=1,LOT
          C(JA+J)=A(IA+I)
          D(JA+J)=-A(IB+I)
          I=I+INC3
          J=J+INC4
270       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
280     CONTINUE
      GO TO 900
!
290   CONTINUE
      Z=1.0/FLOAT(N)
      DO 294 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 292 IJK=1,LOT
          C(JA+J)=Z*(A(IA+I)+A(IB+I))
          C(JB+J)=Z*(A(IA+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
292       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
294     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 3
!     -------------------
300   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB
!
      IF (LA.EQ.M) GO TO 390
!
      DO 320 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 310 IJK=1,LOT
          C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
          C(JB+J)=A(IA+I)-0.5*(A(IB+I)+A(IC+I))
          D(JB+J)=SIN60*(A(IC+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
310       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
320     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JC) GO TO 360
      DO 350 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        JBASE=0
        DO 340 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 330 IJK=1,LOT
            A1=(C1*A(IB+I)+S1*B(IB+I))+(C2*A(IC+I)+S2*B(IC+I))
            B1=(C1*B(IB+I)-S1*A(IB+I))+(C2*B(IC+I)-S2*A(IC+I))
            A2=A(IA+I)-0.5*A1
            B2=B(IA+I)-0.5*B1
            A3=SIN60*((C1*A(IB+I)+S1*B(IB+I))-(C2*A(IC+I)+S2*B(IC+I)))
            B3=SIN60*((C1*B(IB+I)-S1*A(IB+I))-(C2*B(IC+I)-S2*A(IC+I)))
            C(JA+J)=A(IA+I)+A1
            D(JA+J)=B(IA+I)+B1
            C(JB+J)=A2+B3
            D(JB+J)=B2-A3
            C(JC+J)=A2-B3
            D(JC+J)=-(B2+A3)
            I=I+INC3
            J=J+INC4
330         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
340       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC-JINK
350     CONTINUE
      IF (JA.GT.JC) GO TO 900
360   CONTINUE
      JBASE=0
      DO 380 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 370 IJK=1,LOT
          C(JA+J)=A(IA+I)+0.5*(A(IB+I)-A(IC+I))
          D(JA+J)=-SIN60*(A(IB+I)+A(IC+I))
          C(JB+J)=A(IA+I)-(A(IB+I)-A(IC+I))
          I=I+INC3
          J=J+INC4
370       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
380     CONTINUE
      GO TO 900
!
390   CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 394 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 392 IJK=1,LOT
          C(JA+J)=Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
          C(JB+J)=Z*(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))
          D(JB+J)=ZSIN60*(A(IC+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
392       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
394     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 4
!     -------------------
400   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JB
!
      IF (LA.EQ.M) GO TO 490
!
      DO 420 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 410 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
          C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
          C(JB+J)=A(IA+I)-A(IC+I)
          D(JB+J)=A(ID+I)-A(IB+I)
          I=I+INC3
          J=J+INC4
410       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
420     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JC) GO TO 460
      DO 450 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        JBASE=0
        DO 440 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 430 IJK=1,LOT
            A0=A(IA+I)+(C2*A(IC+I)+S2*B(IC+I))
            A2=A(IA+I)-(C2*A(IC+I)+S2*B(IC+I))
            A1=(C1*A(IB+I)+S1*B(IB+I))+(C3*A(ID+I)+S3*B(ID+I))
            A3=(C1*A(IB+I)+S1*B(IB+I))-(C3*A(ID+I)+S3*B(ID+I))
            B0=B(IA+I)+(C2*B(IC+I)-S2*A(IC+I))
            B2=B(IA+I)-(C2*B(IC+I)-S2*A(IC+I))
            B1=(C1*B(IB+I)-S1*A(IB+I))+(C3*B(ID+I)-S3*A(ID+I))
            B3=(C1*B(IB+I)-S1*A(IB+I))-(C3*B(ID+I)-S3*A(ID+I))
            C(JA+J)=A0+A1
            C(JC+J)=A0-A1
            D(JA+J)=B0+B1
            D(JC+J)=B1-B0
            C(JB+J)=A2+B3
            C(JD+J)=A2-B3
            D(JB+J)=B2-A3
            D(JD+J)=-(B2+A3)
            I=I+INC3
            J=J+INC4
430         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
440       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC-JINK
        JD=JD-JINK
450     CONTINUE
      IF (JB.GT.JC) GO TO 900
460   CONTINUE
      SIN45=SQRT(0.5)
      JBASE=0
      DO 480 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 470 IJK=1,LOT
          C(JA+J)=A(IA+I)+SIN45*(A(IB+I)-A(ID+I))
          C(JB+J)=A(IA+I)-SIN45*(A(IB+I)-A(ID+I))
          D(JA+J)=-A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
          D(JB+J)=A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
          I=I+INC3
          J=J+INC4
470       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
480     CONTINUE
      GO TO 900
!
490   CONTINUE
      Z=1.0/FLOAT(N)
      DO 494 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 492 IJK=1,LOT
          C(JA+J)=Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
          C(JC+J)=Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
          C(JB+J)=Z*(A(IA+I)-A(IC+I))
          D(JB+J)=Z*(A(ID+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
492       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
494     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 5
!     -------------------
500   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC
      JE=JB
!
      IF (LA.EQ.M) GO TO 590
!
      DO 520 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 510 IJK=1,LOT
          A1=A(IB+I)+A(IE+I)
          A3=A(IB+I)-A(IE+I)
          A2=A(IC+I)+A(ID+I)
          A4=A(IC+I)-A(ID+I)
          A5=A(IA+I)-0.25*(A1+A2)
          A6=QRT5*(A1-A2)
          C(JA+J)=A(IA+I)+(A1+A2)
          C(JB+J)=A5+A6
          C(JC+J)=A5-A6
          D(JB+J)=-SIN72*A3-SIN36*A4
          D(JC+J)=-SIN36*A3+SIN72*A4
          I=I+INC3
          J=J+INC4
510       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
520     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JD) GO TO 560
      DO 550 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        JBASE=0
        DO 540 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 530 IJK=1,LOT
            A1=(C1*A(IB+I)+S1*B(IB+I))+(C4*A(IE+I)+S4*B(IE+I))
            A3=(C1*A(IB+I)+S1*B(IB+I))-(C4*A(IE+I)+S4*B(IE+I))
            A2=(C2*A(IC+I)+S2*B(IC+I))+(C3*A(ID+I)+S3*B(ID+I))
            A4=(C2*A(IC+I)+S2*B(IC+I))-(C3*A(ID+I)+S3*B(ID+I))
            B1=(C1*B(IB+I)-S1*A(IB+I))+(C4*B(IE+I)-S4*A(IE+I))
            B3=(C1*B(IB+I)-S1*A(IB+I))-(C4*B(IE+I)-S4*A(IE+I))
            B2=(C2*B(IC+I)-S2*A(IC+I))+(C3*B(ID+I)-S3*A(ID+I))
            B4=(C2*B(IC+I)-S2*A(IC+I))-(C3*B(ID+I)-S3*A(ID+I))
            A5=A(IA+I)-0.25*(A1+A2)
            A6=QRT5*(A1-A2)
            B5=B(IA+I)-0.25*(B1+B2)
            B6=QRT5*(B1-B2)
            A10=A5+A6
            A20=A5-A6
            B10=B5+B6
            B20=B5-B6
            A11=SIN72*B3+SIN36*B4
            A21=SIN36*B3-SIN72*B4
            B11=SIN72*A3+SIN36*A4
            B21=SIN36*A3-SIN72*A4
            C(JA+J)=A(IA+I)+(A1+A2)
            C(JB+J)=A10+A11
            C(JE+J)=A10-A11
            C(JC+J)=A20+A21
            C(JD+J)=A20-A21
            D(JA+J)=B(IA+I)+(B1+B2)
            D(JB+J)=B10-B11
            D(JE+J)=-(B10+B11)
            D(JC+J)=B20-B21
            D(JD+J)=-(B20+B21)
            I=I+INC3
            J=J+INC4
530         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
540       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC+JINK
        JD=JD-JINK
        JE=JE-JINK
550     CONTINUE
      IF (JB.GT.JD) GO TO 900
560   CONTINUE
      JBASE=0
      DO 580 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 570 IJK=1,LOT
          A1=A(IB+I)+A(IE+I)
          A3=A(IB+I)-A(IE+I)
          A2=A(IC+I)+A(ID+I)
          A4=A(IC+I)-A(ID+I)
          A5=A(IA+I)+0.25*(A3-A4)
          A6=QRT5*(A3+A4)
          C(JA+J)=A5+A6
          C(JB+J)=A5-A6
          C(JC+J)=A(IA+I)-(A3-A4)
          D(JA+J)=-SIN36*A1-SIN72*A2
          D(JB+J)=-SIN72*A1+SIN36*A2
          I=I+INC3
          J=J+INC4
570       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
580     CONTINUE
      GO TO 900
!
590   CONTINUE
      Z=1.0/FLOAT(N)
      ZQRT5=Z*QRT5
      ZSIN36=Z*SIN36
      ZSIN72=Z*SIN72
      DO 594 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 592 IJK=1,LOT
          A1=A(IB+I)+A(IE+I)
          A3=A(IB+I)-A(IE+I)
          A2=A(IC+I)+A(ID+I)
          A4=A(IC+I)-A(ID+I)
          A5=Z*(A(IA+I)-0.25*(A1+A2))
          A6=ZQRT5*(A1-A2)
          C(JA+J)=Z*(A(IA+I)+(A1+A2))
          C(JB+J)=A5+A6
          C(JC+J)=A5-A6
          D(JB+J)=-ZSIN72*A3-ZSIN36*A4
          D(JC+J)=-ZSIN36*A3+ZSIN72*A4
          I=I+INC3
          J=J+INC4
592       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
594     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 6
!     -------------------
600   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JC
      JF=JB
!
      IF (LA.EQ.M) GO TO 690
!
      DO 620 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 610 IJK=1,LOT
          A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
          C(JA+J)=(A(IA+I)+A(ID+I))+A11
          C(JC+J)=(A(IA+I)+A(ID+I)-0.5*A11)
          D(JC+J)=SIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
          A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
          C(JB+J)=(A(IA+I)-A(ID+I))-0.5*A11
          D(JB+J)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
          C(JD+J)=(A(IA+I)-A(ID+I))+A11
          I=I+INC3
          J=J+INC4
610       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
620     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JC.EQ.JD) GO TO 660
      DO 650 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        KF=KE+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        C5=TRIGS(KF+1)
        S5=TRIGS(KF+2)
        JBASE=0
        DO 640 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 630 IJK=1,LOT
            A1=C1*A(IB+I)+S1*B(IB+I)
            B1=C1*B(IB+I)-S1*A(IB+I)
            A2=C2*A(IC+I)+S2*B(IC+I)
            B2=C2*B(IC+I)-S2*A(IC+I)
            A3=C3*A(ID+I)+S3*B(ID+I)
            B3=C3*B(ID+I)-S3*A(ID+I)
            A4=C4*A(IE+I)+S4*B(IE+I)
            B4=C4*B(IE+I)-S4*A(IE+I)
            A5=C5*A(IF+I)+S5*B(IF+I)
            B5=C5*B(IF+I)-S5*A(IF+I)
            A11=(A2+A5)+(A1+A4)
            A20=(A(IA+I)+A3)-0.5*A11
            A21=SIN60*((A2+A5)-(A1+A4))
            B11=(B2+B5)+(B1+B4)
            B20=(B(IA+I)+B3)-0.5*B11
            B21=SIN60*((B2+B5)-(B1+B4))
            C(JA+J)=(A(IA+I)+A3)+A11
            D(JA+J)=(B(IA+I)+B3)+B11
            C(JC+J)=A20-B21
            D(JC+J)=A21+B20
            C(JE+J)=A20+B21
            D(JE+J)=A21-B20
            A11=(A2-A5)+(A4-A1)
            A20=(A(IA+I)-A3)-0.5*A11
            A21=SIN60*((A4-A1)-(A2-A5))
            B11=(B5-B2)-(B4-B1)
            B20=(B3-B(IA+I))-0.5*B11
            B21=SIN60*((B5-B2)+(B4-B1))
            C(JB+J)=A20-B21
            D(JB+J)=A21-B20
            C(JD+J)=A11+(A(IA+I)-A3)
            D(JD+J)=B11+(B3-B(IA+I))
            C(JF+J)=A20+B21
            D(JF+J)=A21+B20
            I=I+INC3
            J=J+INC4
630         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
640       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC+JINK
        JD=JD-JINK
        JE=JE-JINK
        JF=JF-JINK
650     CONTINUE
      IF (JC.GT.JD) GO TO 900
660   CONTINUE
      JBASE=0
      DO 680 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 670 IJK=1,LOT
          C(JA+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I)))+ SIN60*(A(IB+I)-A(IF+I))
          D(JA+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))-SIN60*(A(IC+I)+A(IE+I))
          C(JB+J)=A(IA+I)-(A(IC+I)-A(IE+I))
          D(JB+J)=A(ID+I)-(A(IB+I)+A(IF+I))
          C(JC+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I)))-SIN60*(A(IB+I)-A(IF+I))
          D(JC+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))+SIN60*(A(IC+I)+A(IE+I))
          I=I+INC3
          J=J+INC4
670       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
680     CONTINUE
      GO TO 900
!
690   CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 694 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 692 IJK=1,LOT
          A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
          C(JA+J)=Z*((A(IA+I)+A(ID+I))+A11)
          C(JC+J)=Z*((A(IA+I)+A(ID+I))-0.5*A11)
          D(JC+J)=ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
          A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
          C(JB+J)=Z*((A(IA+I)-A(ID+I))-0.5*A11)
          D(JB+J)=ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
          C(JD+J)=Z*((A(IA+I)-A(ID+I))+A11)
          I=I+INC3
          J=J+INC4
692       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
694     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 8
!     -------------------
800   CONTINUE
      IBAD=3
      IF (LA.NE.M) GO TO 910
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      IG=IF+IINK
      IH=IG+IINK
      JA=1
      JB=JA+LA*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JD+2*M*INC2
      Z=1.0/FLOAT(N)
      ZSIN45=Z*SQRT(0.5)
!
      DO 820 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 810 IJK=1,LOT
          C(JA+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+ &
                   ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
          C(JE+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))- &
                   ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
          C(JC+J)=Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
          D(JC+J)=Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
          C(JB+J)=Z*(A(IA+I)-A(IE+I)) &
                   +ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
          C(JD+J)=Z*(A(IA+I)-A(IE+I)) &
                   -ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
          D(JB+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I))) &
                   +Z*(A(IG+I)-A(IC+I))
          D(JD+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I))) &
                   -Z*(A(IG+I)-A(IC+I))
          I=I+INC3
          J=J+INC4
810       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
820     CONTINUE
!
!     RETURN
!     ------
900   CONTINUE
      IBAD=0
910   CONTINUE
      IERR=IBAD
      RETURN
      END



      SUBROUTINE RPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA,IERR)
!     SUBROUTINE 'RPASSM' - PERFORMS ONE PASS THROUGH DATA AS PART
!     OF MULTIPLE REAL FFT (FOURIER SYNTHESIS) ROUTINE
!
!     A IS FIRST REAL INPUT VECTOR
!         EQUIVALENCE B(1) WITH A (LA*INC1+1)
!     C IS FIRST REAL OUTPUT VECTOR
!         EQUIVALENCE D(1) WITH C(IFAC*LA*INC2+1)
!     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
!     INC1 IS THE ADDRESSING INCREMENT FOR A
!     INC2 IS THE ADDRESSING INCREMENT FOR C
!     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
!     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
!     LOT IS THE NUMBER OF VECTORS
!     N IS THE LENGTH OF THE VECTORS
!     IFAC IS THE CURRENT FACTOR OF N
!     LA IS THE PRODUCT OF PREVIOUS FACTORS
!     IERR IS AN ERROR INDICATOR:
!              0 - PASS COMPLETED WITHOUT ERROR
!              1 - LOT GREATER THAN 64
!              2 - IFAC NOT CATERED FOR
!              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
!
!-----------------------------------------------------------------------
!
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(N)
!
      DIMENSION A10(64),A11(64),A20(64),A21(64),B10(64),B11(64),B20(64),B21(64)
      parameter(SIN36 = 0.587785252292473)
      parameter(SIN72 = 0.951056516295154)
      parameter(QRT5  = 0.559016994374947)
      parameter(SIN60 = 0.866025403784437)
!
      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      KSTOP=(N-IFAC)/(2*IFAC)
!
      IBAD=1
      IF (LOT.GT.64) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO
!
!     CODING FOR FACTOR 2
!     -------------------
200   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      JA=1
      JB=JA+JINK
!
      IF (LA.EQ.M) GO TO 290
!
      DO 220 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 210 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=A(IA+I)-A(IB+I)
          I=I+INC3
          J=J+INC4
210       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
220     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB-IINK
      IBASE=0
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IB) GO TO 260
      DO 250 K=LA,KSTOP,LA
        KB=K+K
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        IBASE=0
        DO 240 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 230 IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            D(JA+J)=B(IA+I)-B(IB+I)
            C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))
            D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))
            I=I+INC3
            J=J+INC4
230         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
240       CONTINUE
        IA=IA+IINK
        IB=IB-IINK
        JBASE=JBASE+JUMP
250     CONTINUE
      IF (IA.GT.IB) GO TO 900
260   CONTINUE
      IBASE=0
      DO 280 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 270 IJK=1,LOT
          C(JA+J)=A(IA+I)
          C(JB+J)=-B(IA+I)
          I=I+INC3
          J=J+INC4
270       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
280     CONTINUE
      GO TO 900
!
290   CONTINUE
      DO 294 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 292 IJK=1,LOT
          C(JA+J)=2.0*(A(IA+I)+A(IB+I))
          C(JB+J)=2.0*(A(IA+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
292       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
294     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 3
!     -------------------
300   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
!
      IF (LA.EQ.M) GO TO 390
!
      DO 320 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 310 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=(A(IA+I)-0.5*A(IB+I))-(SIN60*(B(IB+I)))
          C(JC+J)=(A(IA+I)-0.5*A(IB+I))+(SIN60*(B(IB+I)))
          I=I+INC3
          J=J+INC4
310       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
320     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IC) GO TO 360
      DO 350 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        IBASE=0
        DO 340 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 330 IJK=1,LOT
            C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
            D(JA+J)=B(IA+I)+(B(IB+I)-B(IC+I))
            C(JB+J)=C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+    &
             I)+B(IC+I))))   -S1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60 &
             *(A(IB+I)-A(IC+I))))
            D(JB+J)=S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+    &
             I)+B(IC+I))))   +C1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60 &
             *(A(IB+I)-A(IC+I))))
            C(JC+J)=C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+    &
             I)+B(IC+I))))   -S2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60 &
             *(A(IB+I)-A(IC+I))))
            D(JC+J)=S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+    &
             I)+B(IC+I))))   +C2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60 &
             *(A(IB+I)-A(IC+I))))
            I=I+INC3
            J=J+INC4
330         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
340       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC-IINK
        JBASE=JBASE+JUMP
350     CONTINUE
      IF (IA.GT.IC) GO TO 900
360   CONTINUE
      IBASE=0
      DO 380 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 370 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
          C(JC+J)=-(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
          I=I+INC3
          J=J+INC4
370       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
380     CONTINUE
      GO TO 900
!
390   CONTINUE
      SSIN60=2.0*SIN60
      DO 394 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 392 IJK=1,LOT
          C(JA+J)=2.0*(A(IA+I)+A(IB+I))
          C(JB+J)=(2.0*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
          C(JC+J)=(2.0*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
          I=I+INC3
          J=J+INC4
392       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
394     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 4
!     -------------------
400   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
!
      IF (LA.EQ.M) GO TO 490
!
      DO 420 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 410 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(IC+I))+A(IB+I)
          C(JB+J)=(A(IA+I)-A(IC+I))-B(IB+I)
          C(JC+J)=(A(IA+I)+A(IC+I))-A(IB+I)
          C(JD+J)=(A(IA+I)-A(IC+I))+B(IB+I)
          I=I+INC3
          J=J+INC4
410       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
420     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      ID=ID-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.IC) GO TO 460
      DO 450 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        IBASE=0
        DO 440 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 430 IJK=1,LOT
            C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
            D(JA+J)=(B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))
            C(JC+J)=C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
                   -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
            D(JC+J)=S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
                   +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
            C(JB+J)=C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I))) &
                   -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
            D(JB+J)=S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I))) &
                   +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
            C(JD+J)=C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I))) &
                   -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
            D(JD+J)=S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I))) &
                   +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
            I=I+INC3
            J=J+INC4
430         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
440       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC-IINK
        ID=ID-IINK
        JBASE=JBASE+JUMP
450     CONTINUE
      IF (IB.GT.IC) GO TO 900
460   CONTINUE
      IBASE=0
      SIN45=SQRT(0.5)
      DO 480 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 470 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))
          C(JC+J)=B(IB+I)-B(IA+I)
          C(JD+J)=-SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))
          I=I+INC3
          J=J+INC4
470       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
480     CONTINUE
      GO TO 900
!
490   CONTINUE
      DO 494 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 492 IJK=1,LOT
          C(JA+J)=2.0*((A(IA+I)+A(IC+I))+A(IB+I))
          C(JB+J)=2.0*((A(IA+I)-A(IC+I))-B(IB+I))
          C(JC+J)=2.0*((A(IA+I)+A(IC+I))-A(IB+I))
          C(JD+J)=2.0*((A(IA+I)-A(IC+I))+B(IB+I))
          I=I+INC3
          J=J+INC4
492       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
494     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 5
!     -------------------
500   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC
      IE=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
!
      IF (LA.EQ.M) GO TO 590
!
      DO 520 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 510 IJK=1,LOT
          C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
          C(JB+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I &
           )))     -(SIN72*B(IB+I)+SIN36*B(IC+I))
          C(JC+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I &
           )))     -(SIN36*B(IB+I)-SIN72*B(IC+I))
          C(JD+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I &
           )))     +(SIN36*B(IB+I)-SIN72*B(IC+I))
          C(JE+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I &
           )))     +(SIN72*B(IB+I)+SIN36*B(IC+I))
          I=I+INC3
          J=J+INC4
510       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
520     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.ID) GO TO 560
      DO 550 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        IBASE=0
        DO 540 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 530 IJK=1,LOT
!
            A10(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))) &
             )         +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
            A20(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))) &
             )         -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
            B10(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))) &
             )         +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
            B20(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))) &
             )         -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
            A11(IJK)=SIN72*(B(IB+I)+B(IE+I))+SIN36*(B(IC+I)+B(ID+I))
            A21(IJK)=SIN36*(B(IB+I)+B(IE+I))-SIN72*(B(IC+I)+B(ID+I))
            B11(IJK)=SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))
            B21(IJK)=SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))
!
            C(JA+J)=A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))
            D(JA+J)=B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))
            C(JB+J)=C1*(A10(IJK)-A11(IJK))-S1*(B10(IJK)+B11(IJK))
            D(JB+J)=S1*(A10(IJK)-A11(IJK))+C1*(B10(IJK)+B11(IJK))
            C(JE+J)=C4*(A10(IJK)+A11(IJK))-S4*(B10(IJK)-B11(IJK))
            D(JE+J)=S4*(A10(IJK)+A11(IJK))+C4*(B10(IJK)-B11(IJK))
            C(JC+J)=C2*(A20(IJK)-A21(IJK))-S2*(B20(IJK)+B21(IJK))
            D(JC+J)=S2*(A20(IJK)-A21(IJK))+C2*(B20(IJK)+B21(IJK))
            C(JD+J)=C3*(A20(IJK)+A21(IJK))-S3*(B20(IJK)-B21(IJK))
            D(JD+J)=S3*(A20(IJK)+A21(IJK))+C3*(B20(IJK)-B21(IJK))
!
            I=I+INC3
            J=J+INC4
530         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
540       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC+IINK
        ID=ID-IINK
        IE=IE-IINK
        JBASE=JBASE+JUMP
550     CONTINUE
      IF (IB.GT.ID) GO TO 900
560   CONTINUE
      IBASE=0
      DO 580 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 570 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(IB+I))+A(IC+I)
          C(JB+J)=(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+I &
           )))     -(SIN36*B(IA+I)+SIN72*B(IB+I))
          C(JE+J)=-(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+ &
           I)))    -(SIN36*B(IA+I)+SIN72*B(IB+I))
          C(JC+J)=(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+I &
           )))     -(SIN72*B(IA+I)-SIN36*B(IB+I))
          C(JD+J)=-(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+ &
           I)))    -(SIN72*B(IA+I)-SIN36*B(IB+I))
          I=I+INC3
          J=J+INC4
570       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
580     CONTINUE
      GO TO 900
!
590   CONTINUE
      QQRT5=2.0*QRT5
      SSIN36=2.0*SIN36
      SSIN72=2.0*SIN72
      DO 594 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 592 IJK=1,LOT
          C(JA+J)=2.0*(A(IA+I)+(A(IB+I)+A(IC+I)))
          C(JB+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))
          C(JC+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))
          C(JD+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))
          C(JE+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))
          I=I+INC3
          J=J+INC4
592       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
594     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 6
!     -------------------
600   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC+2*M*INC1
      IE=IC
      IF=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
!
      IF (LA.EQ.M) GO TO 690
!
      DO 620 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 610 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))
          C(JD+J)=(A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))
          C(JB+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I))) &
                   -(SIN60*(B(IB+I)+B(IC+I)))
          C(JF+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I))) &
                   +(SIN60*(B(IB+I)+B(IC+I)))
          C(JC+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I))) &
                   -(SIN60*(B(IB+I)-B(IC+I)))
          C(JE+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I))) &
                   +(SIN60*(B(IB+I)-B(IC+I)))
          I=I+INC3
          J=J+INC4
610       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
620     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      IF=IF-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IC.EQ.ID) GO TO 660
      DO 650 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        KF=KE+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        C5=TRIGS(KF+1)
        S5=TRIGS(KF+2)
        IBASE=0
        DO 640 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 630 IJK=1,LOT
!
            A11(IJK)= (A(IE+I)+A(IB+I))+(A(IC+I)+A(IF+I))
            A20(IJK)=(A(IA+I)+A(ID+I))-0.5*A11(IJK)
            A21(IJK)=SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IF+I)))
            B11(IJK)= (B(IB+I)-B(IE+I))+(B(IC+I)-B(IF+I))
            B20(IJK)=(B(IA+I)-B(ID+I))-0.5*B11(IJK)
            B21(IJK)=SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IF+I)))
!
            C(JA+J)=(A(IA+I)+A(ID+I))+A11(IJK)
            D(JA+J)=(B(IA+I)-B(ID+I))+B11(IJK)
            C(JC+J)=C2*(A20(IJK)-B21(IJK))-S2*(B20(IJK)+A21(IJK))
            D(JC+J)=S2*(A20(IJK)-B21(IJK))+C2*(B20(IJK)+A21(IJK))
            C(JE+J)=C4*(A20(IJK)+B21(IJK))-S4*(B20(IJK)-A21(IJK))
            D(JE+J)=S4*(A20(IJK)+B21(IJK))+C4*(B20(IJK)-A21(IJK))
!
            A11(IJK)=(A(IE+I)-A(IB+I))+(A(IC+I)-A(IF+I))
            B11(IJK)=(B(IE+I)+B(IB+I))-(B(IC+I)+B(IF+I))
            A20(IJK)=(A(IA+I)-A(ID+I))-0.5*A11(IJK)
            A21(IJK)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
            B20(IJK)=(B(IA+I)+B(ID+I))+0.5*B11(IJK)
            B21(IJK)=SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IF+I)))
!
            C(JD+J)=C3*((A(IA+I)-A(ID+I))+A11(IJK))-S3*((B(IA+I)+B(ID+I &
             ))-B11(IJK))
            D(JD+J)=S3*((A(IA+I)-A(ID+I))+A11(IJK))+C3*((B(IA+I)+B(ID+I &
             ))-B11(IJK))
            C(JB+J)=C1*(A20(IJK)-B21(IJK))-S1*(B20(IJK)-A21(IJK))
            D(JB+J)=S1*(A20(IJK)-B21(IJK))+C1*(B20(IJK)-A21(IJK))
            C(JF+J)=C5*(A20(IJK)+B21(IJK))-S5*(B20(IJK)+A21(IJK))
            D(JF+J)=S5*(A20(IJK)+B21(IJK))+C5*(B20(IJK)+A21(IJK))
!
            I=I+INC3
            J=J+INC4
630         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
640       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC+IINK
        ID=ID-IINK
        IE=IE-IINK
        IF=IF-IINK
        JBASE=JBASE+JUMP
650     CONTINUE
      IF (IC.GT.ID) GO TO 900
660   CONTINUE
      IBASE=0
      DO 680 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 670 IJK=1,LOT
          C(JA+J)=A(IB+I)+(A(IA+I)+A(IC+I))
          C(JD+J)=B(IB+I)-(B(IA+I)+B(IC+I))
          C(JB+J)=(SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I))
          C(JF+J)=-(SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I))
          C(JC+J)=SIN60*(B(IC+I)-B(IA+I))+(0.5*(A(IA+I)+A(IC+I))-A(IB+I))
          C(JE+J)=SIN60*(B(IC+I)-B(IA+I))-(0.5*(A(IA+I)+A(IC+I))-A(IB+I))
          I=I+INC3
          J=J+INC4
670       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
680     CONTINUE
      GO TO 900
!
690   CONTINUE
      SSIN60=2.0*SIN60
      DO 694 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 692 IJK=1,LOT
          C(JA+J)=(2.0*(A(IA+I)+A(ID+I)))+(2.0*(A(IB+I)+A(IC+I)))
          C(JD+J)=(2.0*(A(IA+I)-A(ID+I)))-(2.0*(A(IB+I)-A(IC+I)))
          C(JB+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I))) &
                   -(SSIN60*(B(IB+I)+B(IC+I)))
          C(JF+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I))) &
                   +(SSIN60*(B(IB+I)+B(IC+I)))
          C(JC+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I))) &
                   -(SSIN60*(B(IB+I)-B(IC+I)))
          C(JE+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I))) &
                   +(SSIN60*(B(IB+I)-B(IC+I)))
          I=I+INC3
          J=J+INC4
692       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
694     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 8
!     -------------------
800   CONTINUE
      IBAD=3
      IF (LA.NE.M) GO TO 910
      IA=1
      IB=IA+LA*INC1
      IC=IB+2*LA*INC1
      ID=IC+2*LA*INC1
      IE=ID+2*LA*INC1
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
      JG=JF+JINK
      JH=JG+JINK
      SSIN45=SQRT(2.0)
!
      DO 820 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 810 IJK=1,LOT
          C(JA+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
          C(JE+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
          C(JC+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
          C(JG+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
          C(JB+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I)) &
                   +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JF+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I)) &
                   -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JD+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I)) &
                   -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          C(JH+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I)) &
                   +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          I=I+INC3
          J=J+INC4
810       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
820     CONTINUE
!
!     RETURN
!     ------
900   CONTINUE
      IBAD=0
910   CONTINUE
      IERR=IBAD
      RETURN
      END


