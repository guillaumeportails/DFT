C---------------------------------------------------------------------------
C
C  gfortran -std=f95 -Wall sorensen.f && ./a.out | tee tmp.m && octave tmp.m
C

      PROGRAM SORENSEN
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: M = 6
      INTEGER, PARAMETER :: N = 2**M
      REAL X(N),Y(N)
C
      WRITE (*,*) "% SORENSEN TEST"

      CALL FILL
      CALL FFTDIF(X,Y,N,M)
      CALL DISP('Y')
      WRITE (*,*) "F = fft(X); fprintf('ERR = %g\n',max(abs(Y-F)));"

      CALL FILL
      CALL FFTDIT(X,Y,N,M)
      CALL DISP('Y')
      WRITE (*,*) "F = fft(X); fprintf('ERR = %g\n',max(abs(Y-F)));"

C---------------------------------------------------------------------------

      CONTAINS

      REAL FUNCTION URAND()
      IMPLICIT NONE
      INTEGER, SAVE :: SEED = 1234
      SEED = 9821*SEED + 211327
      URAND = REAL(SEED)/2147483647.0
      END FUNCTION URAND

      SUBROUTINE FILL
      IMPLICIT NONE
      INTEGER I
      DO I = 1,N
            X(I) = 6*URAND()-3
            Y(I) = 6*URAND()-3
      END DO
      CALL DISP('X')
      END SUBROUTINE FILL

      SUBROUTINE DISP(C)
      IMPLICIT NONE
      CHARACTER C
      INTEGER I
      DO I = 1,N
            WRITE (*,*) "",C,"(",I,") = ",X(I)," + i*",Y(I),";"
      END DO
      END SUBROUTINE DISP

C---------------------------------------------------------------------------
C   A DUHAMEL-HOLLMAN SPLIT-RADIX DIF FFT
C   REF: ELECTRONICS LETTERS, JAN. 5, 1984
C   COMPLEX INPUT AND OUTPUT DATA IN ARRAYS X AND Y
C   LENGTH IS N = 2 ** M
C   C.S. BURRUS, RICE UNIV. DEC. 1984
C---------------------------------------------------------------------------
C
      SUBROUTINE FFTDIF (X,Y,N,M)
      IMPLICIT NONE
      INTEGER N,M
      REAL X(N),Y(N)
C
      INTEGER I,J,K, IS,ID, I0,I1,I2,I3, N1,N2,N4
      REAL R1,R2,S1,S2,S3, XT
      REAL CC1,CC3,SS1,SS3
      REAL E,A,A3
C
      N2 = 2*N
      DO 10 K = 1, M-1
        N2 = N2/2
        N4 = N2/4
        E = 6.283185307179586 / N2
        A = 0
        DO 20 J = 1, N4
            A3 = 3*A
            CC1 = COS(A)
            SS1 = SIN(A)
            CC3 = COS(A3)
            SS3 = SIN(A3)
            A = J*E
            IS = J
            ID = 2*N2
  40        DO 30 I0 = IS, N-1, ID
                I1 = I0 + N4
                I2 = I1 + N4
                I3 = I2 + N4 
                R1    = X(I0) - X(I2)
                X(I0) = X(I0) + X(I2)
                R2    = X(I1) - X(I3)
                X(I1) = X(I1) + X(I3) 
                S1    = Y(I0) - Y(I2)
                Y(I0) = Y(I0) + Y(I2) 
                S2    = Y(11) - Y(I3) 
                Y(I1) = Y(I1) + Y(I3)
C
                S3 = R1 - S2 
                R1 = R1 + S2
                S2 = R2 - S1
                R2 = R2 + S1
                X(I2) =  R1*CC1 - S2*SS1 
                Y(I2) = -S2*CC1 - R1*SS1 
                X(I3) =  S3*CC3 + R2*SS3 
                Y(I3) =  R2*CC3 - S3*SS3 
  30        CONTINUE 
            IS = 2*ID - N2 + J 
            ID = 4*ID 
            IF (IS.LT.N) GOTO 40 
  20    CONTINUE 
  10  CONTINUE 

C---------------LAST STAGE, LENGTH-2 BUTTERFLY-----------------------

      IS = 1
      ID = 4
  50  DO 60 I0 = IS, N, ID
        I1 = I0 + 1
        R1 = X(I0)
        X(I0) = R1 + X(I1) 
        X(I1) = R1 - X(I1)
        R1    = Y(I0)
        Y(I0) = R1 + Y(I1)
        Y(11) = R1 - Y(I1)
  60  CONTINUE
      IS = 2*ID - 1
      ID = 4*ID
      IF (IS.LT.N) GOTO 50
C
C------------------BIT REVERSE COUNTER---------------

      J = 1
      N1 = N - 1
      DO 104 I=1, N1
        IF (I.GE.J) GOTO 101
        XT   = X(J)
        X(J) = X(I)
        X(I) = XT
        XT   = Y(J)
        Y(J) = Y(I)
        Y(I) = XT
  101   K = N/2
  102   IF (K.GE.J) GOTO 103
        J = J - K
        K = K/2
        GOTO 102
  103   J = J + K
  104 CONTINUE
      RETURN
      END SUBROUTINE FFTDIF


C---------------------------------------------------------------------------
C   A DUHAMEL-HOLLMAN SPLIT-RADIX DIT FFT
C   REF: Electronics Letters, Jan. 5, 1984
C   Complex input and output data in arrays X and Y
C   Length is N = 2 ** M
C
C   H.V. Sorensen, Rice University, Jan. 4 1985 C
C 
C---------------------------------------------------------------------------
C
      SUBROUTINE FFTDIT (X,Y,N,M)
      IMPLICIT NONE
      INTEGER N,M
      REAL X(N),Y(N)
C
      INTEGER I,J,K, IS,ID, I0,I1,I2,I3, N1,N2,N4
      REAL R1,R2,R3,S1,S2, XT
      REAL CC1,CC3,SS1,SS3
      REAL E,A,A3
C
      J = 1
      N1 = N - 1
      DO 104 I=1, N1
        IF (I.GE.J) GOTO 101
        XT = X(J)
        X(J) = X(I)
        X(I) = XT
        XT = Y(J)
        Y(J) = Y(I)
        Y(I) = XT
  101   K = N/2
  102   IF (K.GE.J) GOTO 103
            J = J - K
            K = K/2
            GOTO 102
  103     J = J + K
  104 CONTINUE
C
C-------------------------Lebgth two transform---------------------
C
      IS = 1
      ID = 4
  70  DO 60 I0 = IS,N,ID
        I1 = I0 + 1
        R1 = X(I0)
        X(I0) = R1 + X(I1)
        X(I1) = R1 - X(I1)
        R1 = Y(I0)
        Y(I0) = R1 + Y(I1)
        Y(I1) = R1 - Y(I1)
  60  CONTINUE
      IS = 2*ID - 1
      ID = 4*ID
      IF (IS.LT.N) GOTO 70
C
C-------------------------L SHAPED BUTTERFLIES---------------------
C
      N2 = 2
      DO 10 K = 2, M
        N2 = N2*2
        N4 = N2/4
        E = 6.283185307179586 / N2
        A = 0
        DO 20 J = 1, N4 
            A3 = 3*A
            CC1 = COS(A)
            SS1 = SIN(A)
            CC3 = COS(A3)
            SS3 = SIN(A3)
            A = J*E
            IS = J
            ID = 2*N2
  40        DO 30 I0 = IS, N-1, ID
                I1 = I0 + N4
                I2 = I1 + N4
                I3 = I2 + N4
C
                R1 = X(I2)*CC1 + Y(I2)*SS1
                S1 = Y(I2)*CC1 - X(I2)*SS1
                R2 = X(I3)*CC3 + Y(I3)*SS3
                S2 = Y(I3)*CC3 - X(I3)*SS3
                R3 = R1 + R2
                R2 = R1 - R2
                R1 = S1 + S2
                S2 = S1 - S2
C
                X(I2) = X(I0) - R3
                X(I0) = X(I0) + R3
                X(I3) = X(I1) - S2
                X(I1) = X(I1) + S2
                Y(I2) = Y(I0) - R1
                Y(I0) = Y(I0) + R1
                Y(I3) = Y(I1) + R2
                Y(I1) = Y(I1) - R2
C
  30        CONTINUE
            IS = 2*ID - N2 + J
            ID = 4*ID
            IF (IS.LT.N) GOTO 40
  20    CONTINUE
  10  CONTINUE
      RETURN
      END SUBROUTINE FFTDIT

C---------------------------------------------------------------------------

       END PROGRAM SORENSEN
