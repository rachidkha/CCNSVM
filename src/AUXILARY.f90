!--------------------------------------------------------------------------------------
 
!----------------Derive of loss function

 
      SUBROUTINE DWDdrv (nobs, nvars, x, y, r, vl,delta)
         IMPLICIT NONE
         INTEGER :: nobs, nvars, i
         DOUBLE PRECISION :: x (nobs, nvars), y (nobs)
         DOUBLE PRECISION :: r (nobs), vl (nvars), dl (nobs)
            DOUBLE PRECISION :: delta
         vl = 0.0
         DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)- 1.0D0) / delta
                  END IF
               END DO
         vl = Matmul(dl * y, x) / nobs
      END SUBROUTINE DWDdrv
!----calculate nw et ws
 










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
        & ITRAN, LIVE, ITER, WSS, IFAULT)
      INTEGER ::ITER 
      INTEGER M,N,K ,IFAULT 
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      DOUBLE PRECISION A(M,N), D(M), C(K,N), AN1(K), AN2(K), WSS(K)

      DOUBLE PRECISION DT(2), ZERO, ONE
      INTEGER I,IL,J,L,INDX,IJ,II
      DOUBLE PRECISION BIG, DA, TEMP, DB, DC,AA

      DATA BIG /1.E30/, ZERO /0.0/, ONE /1.0/

      IFAULT = 3

      IF (K .LE. 1 .OR. K .GE. M) RETURN
      IFAULT = 0

      DO 50 I = 1, M
        IC1(I) = 1
        IC2(I) = 2
        DO 10 IL = 1, 2
          DT(IL) = ZERO
          DO 10 J = 1, N
            DA = A(I,J) - C(IL,J)
            DT(IL) = DT(IL) + DA*DA
   10   CONTINUE
        IF (DT(1) .GT. DT(2)) THEN
          IC1(I) = 2
          IC2(I) = 1
          TEMP = DT(1)
          DT(1) = DT(2)
          DT(2) = TEMP
        END IF
        DO 50 L = 3, K
          DB = ZERO
          DO 30 J = 1, N
            DC = A(I,J) - C(L,J)
            DB = DB + DC*DC
            IF (DB .GE. DT(2)) GO TO 50
   30     CONTINUE
          IF (DB .LT. DT(1)) GO TO 40
          DT(2) = DB
          IC2(I) = L
          GO TO 50
   40     DT(2) = DT(1)
          IC2(I) = IC1(I)
          DT(1) = DB
          IC1(I) = L
   50 CONTINUE

      DO 70 L = 1, K
        NC(L) = 0
        DO 60 J = 1, N
   60   C(L,J) = ZERO
   70 CONTINUE
      DO 90 I = 1, M
        L = IC1(I)
        NC(L) = NC(L) + 1
        DO 80 J = 1, N
   80   C(L,J) = C(L,J) + A(I,J)
   90 CONTINUE

      DO 120 L = 1, K
        IF (NC(L) .EQ. 0) THEN
          IFAULT = 1
          RETURN
        END IF
        AA = NC(L)
        DO 110 J = 1, N
  110   C(L,J) = C(L,J) / AA

        AN2(L) = AA / (AA + ONE)
        AN1(L) = BIG
        IF (AA .GT. ONE) AN1(L) = AA / (AA - ONE)
        ITRAN(L) = 1
        NCP(L) = -1
  120 CONTINUE
      INDX = 0
      DO 140 IJ = 1, ITER

        CALL OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &        ITRAN, LIVE, INDX)

        IF (INDX .EQ. M) GO TO 150

        CALL QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &       ITRAN, INDX)

        IF (K .EQ. 2) GO TO 150

        DO 130 L = 1, K
  130   NCP(L) = 0
  140 CONTINUE

      IFAULT = 2

  150 DO 160 L = 1, K
        WSS(L) = ZERO
        DO 160 J = 1, N
          C(L,J) = ZERO
  160 CONTINUE
      DO 170 I = 1, M
        II = IC1(I)
        DO 170 J = 1, N
          C(II,J) = C(II,J) + A(I,J)
  170 CONTINUE
      DO 190 J = 1, N
        DO 180 L = 1, K
  180   C(L,J) = C(L,J) / FLOAT(NC(L))
        DO 190 I = 1, M
          II = IC1(I)
          DA = A(I,J) - C(II,J)
          WSS(II) = WSS(II) + DA*DA
  190 CONTINUE

      RETURN
      END
!---------------------------------------------------------------------------------------
      SUBROUTINE OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &      ITRAN, LIVE, INDX)

      INTEGER M,N,K,INDX
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      DOUBLE PRECISION    A(M,N), D(M), C(K,N), AN1(K), AN2(K)

      INTEGER L,I,L1,L2,LL,J
      DOUBLE PRECISION ZERO, ONE
      DOUBLE PRECISION BIG,DE,DF,DA,DB,R2,RR,DC,DD,AL1,ALW,AL2,ALT

      DATA BIG /1.0E30/, ZERO /0.0/, ONE/1.0/

      DO 10 L = 1, K
        IF (ITRAN(L) .EQ. 1) LIVE(L) = M + 1
   10 CONTINUE
      DO 100 I = 1, M
        INDX = INDX + 1
        L1 = IC1(I)
        L2 = IC2(I)
        LL = L2

        IF (NC(L1) .EQ. 1) GO TO 90

        IF (NCP(L1) .EQ. 0) GO TO 30
        DE = ZERO
        DO 20 J = 1, N
          DF = A(I,J) - C(L1,J)
          DE = DE + DF*DF
   20   CONTINUE
        D(I) = DE * AN1(L1)

   30   DA = ZERO
        DO 40 J = 1, N
          DB = A(I,J) - C(L2,J)
          DA = DA + DB*DB
   40   CONTINUE
        R2 = DA * AN2(L2)
        DO 60 L = 1, K

          IF (I .GE. LIVE(L1) .AND. I .GE. LIVE(L) .OR. L .EQ. L1 .OR.&
          &   L .EQ. LL) GO TO 60
          RR = R2 / AN2(L)
          DC = ZERO
          DO 50 J = 1, N
            DD = A(I,J) - C(L,J)
            DC = DC + DD*DD
            IF (DC .GE. RR) GO TO 60
   50     CONTINUE
          R2 = DC * AN2(L)
          L2 = L
   60     CONTINUE
          IF (R2 .LT. D(I)) GO TO 70

          IC2(I) = L2
          GO TO 90

   70     INDX = 0
          LIVE(L1) = M + I
          LIVE(L2) = M + I
          NCP(L1) = I
          NCP(L2) = I
          AL1 = NC(L1)
          ALW = AL1 - ONE
          AL2 = NC(L2)
          ALT = AL2 + ONE
          DO 80 J = 1, N
            C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
            C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
   80     CONTINUE
          NC(L1) = NC(L1) - 1
          NC(L2) = NC(L2) + 1
          AN2(L1) = ALW / AL1
          AN1(L1) = BIG
          IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
          AN1(L2) = ALT / AL2
          AN2(L2) = ALT / (ALT + ONE)
          IC1(I) = L2
          IC2(I) = L1
   90   CONTINUE
        IF (INDX .EQ. M) RETURN
  100 CONTINUE
      DO 110 L = 1, K

        ITRAN(L) = 0
        LIVE(L) = LIVE(L) - M
  110 CONTINUE

      RETURN
      END
!--------------------------------------------------------------------------------------------
      SUBROUTINE QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
     &    ITRAN, INDX)

      INTEGER M,N,K,INDX
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K)
      DOUBLE PRECISION A(M,N), D(M), C(K,N), AN1(K), AN2(K)

      DOUBLE PRECISION ZERO, ONE
      INTEGER ICOUN,ISTEP,I,L1,L2,J
      DOUBLE PRECISION BIG,DA,DB,DD,AL1,ALW,AL2,ALT,R2,DE

      DATA BIG /1.0E30/, ZERO /0.0/, ONE /1.0/

      ICOUN = 0
      ISTEP = 0
   10 DO 70 I = 1, M
        ICOUN = ICOUN + 1
        ISTEP = ISTEP + 1
        L1 = IC1(I)
        L2 = IC2(I)

        IF (NC(L1) .EQ. 1) GO TO 60

        IF (ISTEP .GT. NCP(L1)) GO TO 30
        DA = ZERO
        DO 20 J = 1, N
          DB = A(I,J) - C(L1,J)
          DA = DA + DB*DB
   20   CONTINUE
        D(I) = DA * AN1(L1)

   30   IF (ISTEP .GE. NCP(L1) .AND. ISTEP .GE. NCP(L2)) GO TO 60
        R2 = D(I) / AN2(L2)
        DD = ZERO
        DO 40 J = 1, N
          DE = A(I,J) - C(L2,J)
          DD = DD + DE*DE
          IF (DD .GE. R2) GO TO 60
   40   CONTINUE

        ICOUN = 0
        INDX = 0
        ITRAN(L1) = 1
        ITRAN(L2) = 1
        NCP(L1) = ISTEP + M
        NCP(L2) = ISTEP + M
        AL1 = NC(L1)
        ALW = AL1 - ONE
        AL2 = NC(L2)
        ALT = AL2 + ONE
        DO 50 J = 1, N
          C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
          C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
   50   CONTINUE
        NC(L1) = NC(L1) - 1
        NC(L2) = NC(L2) + 1
        AN2(L1) = ALW / AL1
        AN1(L1) = BIG
        IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
        AN1(L2) = ALT / AL2
        AN2(L2) = ALT / (ALT + ONE)
        IC1(I) = L2
        IC2(I) = L1

   60   IF (ICOUN .EQ. M) RETURN
   70 CONTINUE
      GO TO 10
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE sample(m,k,b)
            IMPLICIT NONE
               integer::m
               integer::k
               integer::i,bo,j,h         
               DOUBLE PRECISION::a
               INTEGER  :: b(k)
               
               
                i = 1
                do while   (i<=k)
                    call random_number(a)
                    if (i == 1) then 
                           b(i) = floor(a*m) + 1 
                           i = i+1
                    else                         
                       bo = floor(a*m) + 1
                       h=1
                       do  j=1,(i-1)
                           if (bo==b(j)) then 
                                  h=0
                           end if
                       end do
                       if (h==1) then 
                             b(i) = bo
                             i = i+1 
                       end if
                    end if
                     
                end do
              
              RETURN
        END   SUBROUTINE sample
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE kmeans(A, M, N, CENTER ,  K, groupe,IC1, IC2, NC, AN1, AN2, NCP, D,&
       & ITRAN, LIVE, ITER, WSS, IFAULT)
            implicit none
            INTEGER :: M,N,K,ITER ,IFAULT,jer,l,ll
            INTEGER :: i,j
            
              INTEGER,PARAMETER::nstart2=5
            INTEGER :: IC(M), IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
            DOUBLE PRECISION :: A(M,N), D(M),  AN1(K), AN2(K), WSS(K),CENTER(K,N)
            DOUBLE PRECISION :: best,z  
            INTEGER :: groupe(M) 
            INTEGER ,DIMENSION(:),allocatable ::b
            
            DOUBLE PRECISION ,DIMENSION(:, :),allocatable ::C
            ALLOCATE(b(K),stat=jer)
            
            ALLOCATE (C(K,N), STAT=jer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
            CALL KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
                     & ITRAN, LIVE, ITER, WSS, IFAULT)
              
             best =  sum(WSS)
             
              ll = 0
              DO 
                ll = ll + 1
                IF(ll > nstart2) EXIT
                Call sample(M,K,b)
              DO  i=1,K
                  C(i,:)=A(b(i),:)
              END DO
                 CALL KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,&
                     & ITRAN, LIVE, ITER, WSS, IFAULT)
                  IC = IC1
                  z = sum(WSS)
                  
             if(z .LE. best) then 
               
               best = z
               IC1 =  IC
               
             END IF
                    
           END DO 
             
             groupe =IC1   
      DEALLOCATE (b , C) 
       return
  END subroutine kmeans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE standard(nobs,nvars,x,ju,isd,xmean,xnorm,maj)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER::  nobs
    INTEGER::nvars
    INTEGER::isd
    INTEGER::ju(nvars)
    DOUBLE PRECISION::  x(nobs,nvars)
    DOUBLE PRECISION::xmean(nvars)
    DOUBLE PRECISION::xnorm(nvars)
    DOUBLE PRECISION::maj(nvars)
    ! - - - local declarations - - -
    INTEGER:: j
! - - - begin - - -
    DO j=1,nvars
        IF(ju(j)==1) THEN
            xmean(j)=sum(x(:,j))/nobs     !mean
            x(:,j)=x(:,j)-xmean(j)
            maj(j)=dot_product(x(:,j),x(:,j))/nobs
              IF(isd==1) THEN
                xnorm(j)=sqrt(maj(j))    !standard deviation
                x(:,j)=x(:,j)/xnorm(j)
                maj(j)=1.0D0
            ENDIF
        ENDIF
    END DO
END SUBROUTINE standard

!----calculate nw et ws
SUBROUTINE COVAR(k,nvars,nobs,clust,adj,b,nw,ws)
IMPLICIT NONE
INTEGER:: k,h,nw,nobs,nvars
DOUBLE PRECISION :: ws
integer:: clust(nvars)
DOUBLE PRECISION:: adj(nobs,nvars)
DOUBLE PRECISION:: b(nvars)
                            nw = 1
                            ws = 0.0
                       DO h = 1, nvars


                        IF ((clust(h).EQ.clust(k)).AND.(h.ne. k))  THEN

                               nw = nw + 1
                              !w = dot_product(adj(:,k),adj(:,h))
                              ws = ws + adj(k,h) * b(h)

                        END IF

                      END DO
END SUBROUTINE
! --------------------------------------------------
SUBROUTINE chkvars (nobs, nvars, x, ju)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: ju (nvars)
      DOUBLE PRECISION :: x (nobs, nvars)
    ! - - - local declarations - - -
      INTEGER :: i
      INTEGER :: j
      DOUBLE PRECISION :: t
! - - - begin - - -
      DO j = 1, nvars
         ju (j) = 0
         t = x (1, j)
         DO i = 2, nobs
            IF (x(i, j) /= t) THEN
               ju (j) = 1
               EXIT
            END IF
         END DO
      END DO
END SUBROUTINE chkvars
! --------------------------------------------------
