

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE hsvmcluster (delta, lam2, nobs, nvars, x, y, tom,  KK, jd, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
& nbeta, alam, npass, jerr,cluster)
! --------------------------------------------------

 

      IMPLICIT NONE
    ! - - - arg types - - -

      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: KK
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd,i
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps 
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION,intent(in) :: x (nobs, nvars), tom(nvars,nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: cluster(nvars,nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
     
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
     
   
! - - - begin - - -
! - - - allocate variables - - -
    !  ALLOCATE(TOM(nvars,nvars))
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr

      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0
      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
       
      CALL hsvmclusterpathc (delta, lam2, maj, nobs, nvars, x,  y, tom,  KK , ju, &
     & pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, &
     & ibeta, nbeta, alam, npass, jerr, cluster)
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam

         nk = nbeta (l)
        IF (isd == 1) THEN
           DO j = 1, nk
               beta (j, l) = beta (j, l) / xnorm (ibeta(j))
            END DO
         END IF
         b0 (l) = b0 (l) - dot_product (beta(1:nk, l), &
        & xmean(ibeta(1:nk)))
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE hsvmcluster
! --------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE hsvmclusterpathc (delta, lam2, maj, nobs, nvars, x,  y, tom,  KK , ju, &
& pf, pf2, dfmax, pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, &
& nbeta, alam, npass, jerr, cluster)
! --------------------------------------------------

      IMPLICIT NONE
        ! - - - arg types - - -

      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: KK
      INTEGER :: Nsort
      INTEGER :: NiSort
      INTEGER :: NvSort
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr,clust(nvars)
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: delta
      DOUBLE PRECISION,intent(in) :: x (nobs, nvars), tom(nvars,nvars)
      DOUBLE PRECISION :: chol(nvars,nvars)
      DOUBLE PRECISION :: MAT (nvars, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: vl (nvars), ga(nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars),tlam
      DOUBLE PRECISION :: cluster(nvars,nlam) 

!!!!-------parameters for kmeans
      INTEGER, parameter ::ITER = 10
      INTEGER :: IFAULT 
      INTEGER :: IC2(nvars), NC(KK), NCP(KK), ITRAN(KK), LIVE(KK)
      DOUBLE PRECISION ::  ND(nvars), AN1(KK), AN2(KK), WSS(KK), CENTER(KK,nvars)
      DOUBLE PRECISION ::  weight
    ! - - - local declarations - - -
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: u, ul
      DOUBLE PRECISION :: w
      DOUBLE PRECISION :: ws
      DOUBLE PRECISION :: rhokk,omega
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al,al0
      DOUBLE PRECISION :: alf, TRA(nvars,nvars)
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: dl (nobs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta, oldcoeff
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r, newcoeff
      INTEGER :: jx, jxx (nvars)
     

      INTEGER :: i
      INTEGER :: nz
      INTEGER :: cI
      INTEGER :: k
      INTEGER :: j
      INTEGER :: jj
      INTEGER :: l
      INTEGER :: h
      INTEGER :: nw
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr;
      INTEGER :: ni
      INTEGER :: me, bn
      INTEGER :: COUNTER
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
  
       DOUBLE PRECISION::start,finish
! - - - begin - - -
! - - - allocate variables - - -
  
      ALLOCATE (b(0:nvars), STAT=jerr)
      ALLOCATE (newcoeff(0:nvars), STAT=jerr)
      ALLOCATE (oldcoeff(0:nvars), STAT=ierr)
       jerr = jerr + ierr
      ALLOCATE (oldbeta(0:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs), STAT=ierr)
      jerr = jerr + ierr
     

  


      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      
      r = 0.0D0
      b = 0.0D0
      oldbeta = 0.0D0
      oldcoeff = 0.0D0
      newcoeff=0.0D0
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      maj = 2.0 * maj / delta
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(DBLE(nlam) - 1.0D0))

      END IF
         !!!!!!!!!!!!!!!!!!!!!!initialiser les centers et recalculer la tom
                  
                    DO  i=1,KK
                      DO j = 1,nvars
                      CENTER(i,j)=x(i,j)
                      END DO
                   END DO
          !  call cpu_time(start)
            chol = matmul(TRANSPOSE(tom),tom)
        !   call dgemm("T","N",nvars,nvars,nvars,1.0,tom,nvars,tom,nvars,0.0,chol,nvars)
        !   call cpu_time(finish)
         !   print*, finish-start, 'seconds.'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
! --------- lambda loop ----------------------------


    DO l = 1, nlam

        al0 = al 
        
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
           
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al0 = 0.0D0
          CALL DWDdrv (nobs, nvars, x, y, r, vl,delta)
                ga = Abs(vl)
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN                      
                        al0 = Max (al0, ga(j) /pf(j))
                       
                     END IF
                  END IF
               END DO
               al = al0 * alf  
                
            END IF
         END IF
           
           counter = 0
               
!!!!!!!!!!!!!!!!!!!!!!!!! début boocle kmeans
         DO
              
                 oldcoeff(m(1:ni))  = newcoeff(m(1:ni))
                oldcoeff(0) = newcoeff(0)
                   

              !1111111111111111111111111111111 call kmeans
         IF (count( ABS(oldcoeff (m(1:ni))) /= 0.0D0) < KK) THEN     

            clust=1                
            
        ELSE
                DO i = 1,ni
                     j=m(i)     
                     MAT(:,j) = tom(:,j) * oldcoeff(j)  
                END DO  
                     
             CALL kmns(MAT,nvars, nvars, CENTER, KK, clust,IC2, NC, AN1, AN2, NCP, ND,&
                                      & ITRAN, LIVE, ITER, WSS, IFAULT)
                      
        END IF
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
               
               
             
    ! --------- outer loop ----------------------------
           ctr = 0
          
         DO
    
            oldbeta (0) = b (0)
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
        
        ! --middle loop-------------------------------------

          
            DO
            
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)-1.0D0) / delta
                  END IF
                        u = u + dl (i) * y (i) * x (i, k)

                END DO
                       ! --somme rho jl beta l-------------------------------------

                          nw = 1       
                          ws = 0.0d0    
                       DO h = 1, nvars
                          

                        IF ((clust(h).EQ.clust(k)).AND.(h.ne. k))  THEN 
                               nw = nw + 1 
                               ul = dble(chol(h,k)) * b(h) !dble(chol(h,k))
                               ws = ws + ul
                        END IF
                            
                            
                      END DO
                    
                       
  !-------------------------------------------------------
                     

                     u = maj (k) * b (k) - u / nobs + pf2(k)*lam2*ws/nw

                     v = al * pf (k)
                     v = Abs (u) - v
                     omega = maj(k) + pf2(k)* lam2*  (nw - 1)/nw
                     
                     IF (v > 0.0D0) THEN
                     	b (k) = sign (v, u) / omega

                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = b (k) - oldb
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r + y * x (:, k) * d
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
            END DO
               IF (ni > pmax) EXIT
               d = 0.0D0
               DO i = 1, nobs
               IF (r(i) > 1.0D0) THEN
                  dl (i) = 0.0D0
               ELSE IF (r(i) <= (1-delta)) THEN
                  dl (i) = - 1.0D0
               ELSE
                  dl (i) = (r(i)-1.0D0) / delta
               END IF
                  d = d + dl (i) * y (i)
               END DO !!!!!!!!!!!!!end update b0
              d = - 0.5D0 * delta * d / nobs
               IF (d /= 0.0D0) THEN
                  b (0) = b (0) +  d
                  r = r + y * d
                  dif = Max (dif, 2.0*d**2/delta)
               END IF
               IF (dif < eps) EXIT
      
        ! --inner loop----------------------
              
               DO
            

                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     u = 0.0D0
                     DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)-1.0D0) / delta
                  END IF
                        u = u + dl (i) * y (i) * x (i, k)

                     END DO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               
                           nw = 1
                           ws = 0.0d0
                      DO bn = 1, ni
                        h = m(bn)
                       IF ((clust(h).EQ.clust(k)).AND.(h.ne. k))  THEN                
                              nw = nw + 1 
                              ul = dble(chol(h,k)) * b(h)
                              ws = ws + ul
                       END IF
                      END DO
                      
               
  !-----------------------------------------------------------
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     u = maj (k) * b (k) - u / nobs + pf2(k)*lam2* ws / nw
                     v = al * pf (k)
                     v = Abs (u) - v
                     omega= maj(k) + pf2(k)*lam2* (nw - 1)/nw
                     IF (v > 0.0D0) THEN
                        b (k) = sign (v, u) / omega

                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = b (k) - oldb
                     IF (Abs(d) > 0.0D0) THEN
                        dif = Max (dif, 2.0*d**2/delta)
                        r = r + y * x (:, k) * d
                     END IF
                  END DO
                  d = 0.0D0
                  DO i = 1, nobs
                  IF (r(i) > 1.0D0) THEN
                     dl (i) = 0.0D0
                  ELSE IF (r(i) <= (1-delta)) THEN
                     dl (i) = - 1.0D0
                  ELSE
                     dl (i) = (r(i)-1.0D0) / delta
                  END IF
                     d = d + dl (i) * y (i)
                  END DO
                   d = - 0.5D0 * delta * d / nobs
                  IF (d /= 0.0D0) THEN
                     b (0) = b (0) + d
                     r = r + y * d
                     dif = Max (dif, 2.0*d**2/delta)
                  END IF
                  IF (dif < eps) EXIT
               IF(npass > maxit) EXIT
               END DO ! END FIRST LOOP
               IF(npass > maxit) EXIT
               END DO ! END SECOND LOOP
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
             vrg = 1
            IF ((b(0)-oldbeta(0))**2 >= eps) vrg = 0
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            
            IF (vrg == 1) EXIT
            ctr = ctr + 1
            IF (ctr > maxit) EXIT
             
           END DO !!!!!!!!!!!!!!!!!!!!!!!!fin coordinante descente
          IF (ni > 0) newcoeff (m(1:ni)) = b(m(1:ni))      
                      newcoeff(0) = b(0)
                
                
       !!!!!!!!!!!!!!!!!!!!The secon final check
        
          counter = counter + 1
         IF(sum((newcoeff-oldcoeff)**2)  / (mfl +sum( newcoeff**2)) < eps.or.(counter > maxit)) EXIT
        
        
        END DO ! END kmeans LOOP  
       !!!!!!!!!!!!!!!!!!!!!!save results
         IF (ni > pmax) THEN
          jerr = - 10000 - l
          EXIT
        END IF
              
         beta (1:ni, l) = newcoeff(m(1:ni))    
         b0 (l) = newcoeff(0)
         nbeta (l) = ni
         alam (l) = al
         nalam = l
         cluster(:, l) = clust
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (ABS(beta(1:ni, l)) > 0.0D0)
         IF (me > dfmax) EXIT
              
   END DO   !!!!!!!!!!!!111!!!!!!!!!!end  lambda loop
      DEALLOCATE ( b, oldbeta) 
 !!!!!!!!!!!!!!!!!!fin boocle lambda
      DEALLOCATE ( r, mm)
      DEALLOCATE (oldcoeff,newcoeff)
      RETURN
     

END SUBROUTINE hsvmclusterpathc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------
