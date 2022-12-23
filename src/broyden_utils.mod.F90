MODULE broyden_utils
  USE broy,                            ONLY: broy1
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE ovlap_utils,                     ONLY: dmatc
  USE parac,                           ONLY: parai,&
                                             paral
  USE spin,                            ONLY: clsd
  USE summat_utils,                    ONLY: summat
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: inversemat
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: broyden
  PUBLIC :: give_scr_broyden
  PUBLIC :: kerker

CONTAINS

  ! ==================================================================
  SUBROUTINE broyden(xmix,it,rinm1,rin0,rout0,rinp,nhg,nhgb,&
       ilsd,hg)
    ! ==--------------------------------------------------------------==
    ! == GENERATE NEXT ITERATION USING BROYDENS  METHOD               ==
    ! == (FROM DD JOHNSON, PHYS. REV. B, 38,12807,(1988)              ==
    ! == NHGB IS THE NUMBER OF PLANE-WAVES TO BE MIXED                ==
    ! == RIN0 IS INPUT DENSITY FROM CURRENT STEP (IT)                 ==
    ! == ROUT0 IS OUTPUT DENSITY FROM CURRENT STEP                    ==
    ! == RINM1 IS INPUT DENSITY FROM PREVIOUS STEP                    ==
    ! == RINP  IS NEW MIXED DENSITY                                   ==
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: xmix
    INTEGER                                  :: it
    COMPLEX(real_8)                          :: rinm1(*), rin0(*), rout0(*), &
                                                rinp(*)
    INTEGER                                  :: nhg, nhgb, ilsd
    REAL(real_8)                             :: hg(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'broyden'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: f(:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: deltaf(:,:,:), fm1(:,:), &
                                                u(:,:,:)
    INTEGER                                  :: i, ierr, ig, info, isub, j, &
                                                k, lbroyden
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: aux, sn
    REAL(real_8), ALLOCATABLE                :: a(:,:), b(:,:), c(:), GAMMA(:)
    REAL(real_8), ALLOCATABLE, SAVE          :: w(:,:)

    CALL tiset('   BROYDEN',isub)
    IF (it.GT.broy1%ibreset) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' WARNING| IT GT IBRESET IN BROYDEN. RESETTING IT=1'
       it=1
    ENDIF
    IF (ifirst.EQ.0) THEN

       ! Allocate FM1 (use to store F for the next generation).
       ALLOCATE(fm1(nhgb,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! Allocate U   (use to store Delta RHO).
       ALLOCATE(u(nhgb,broy1%ibreset,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! Allocate DELTAF (use to store F(:)-FM1(:))
       ALLOCATE(deltaf(nhgb,broy1%ibreset,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! Allocate W   (use to store weights)
       ALLOCATE(w(broy1%ibreset,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ifirst=1
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Test SCR
    CALL give_scr_broyden(lbroyden,tag,nhgb)

    ! ==--------------------------------------------------------------==
    ! Allocation inside SCR
    ! IF: QUESTION: why backward direction???      

    ALLOCATE(f(nhgb),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! TODO check stat
    ALLOCATE(a(broy1%ibreset,broy1%ibreset),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(b(broy1%ibreset,broy1%ibreset),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(c(broy1%ibreset),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(GAMMA(broy1%ibreset),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! For SUMMAT, A(:,:)=0 to avoid floating point exception
    CALL zeroing(a)!,broy1%ibreset*broy1%ibreset)
    ! ..F=nout-n
    !$omp parallel do private(IG)
    DO ig=1,nhgb
       f(ig)=rout0(ig)-rin0(ig)
    ENDDO
    ! ..For the first iteration use linear mixing
    IF (it.EQ.1) THEN
       IF (broy1%tadbroy) THEN
          ! ..Calculate W(IT)
          w(1,ilsd)=dotp(nhgb,f,f)
          CALL mp_sum(w(1,ilsd),parai%allgrp)
          w(1,ilsd)=1._real_8/SQRT(w(1,ilsd))
          w(1,ilsd)=MAX(1._real_8,w(1,ilsd))
       ELSE
          w(1,ilsd)=1._real_8
       ENDIF
       ! ..Mixed density on first iteration
       IF (broy1%kermix.EQ.0._real_8) THEN
          !$omp parallel do private(IG)
          DO ig=1,nhgb
             rinp(ig)=rin0(ig)+xmix*f(ig)! Linear mixing
          ENDDO
       ELSE
          !$omp parallel do private(IG)
          DO ig=1,nhgb
             rinp(ig)=rin0(ig)+xmix*f(ig)&
                  *hg(ig)/(hg(ig)+broy1%kermix*broy1%kermix) ! Kerker mixing
          ENDDO
       ENDIF
       ! ..Store F in FM1 for next iteration
       CALL dcopy(2*nhgb,f,1,fm1(1,ilsd),1)
    ENDIF
    ! ==============================================================
    IF (it.GT.1) THEN
       ! ..DELTA F
       !$omp parallel do private(IG)
       DO ig=1,nhgb
          deltaf(ig,it-1,ilsd)=f(ig)-fm1(ig,ilsd)
       ENDDO
       ! ..Store F in FM1 for next iteration
       CALL dcopy(2*nhgb,f,1,fm1(1,ilsd),1)
       sn=dotp(nhgb,deltaf(:,it-1,ilsd),deltaf(:,it-1,ilsd))
       CALL mp_sum(sn,parai%allgrp)
       sn=1._real_8/SQRT(sn)
       CALL dscal(nhgb*2,sn,deltaf(1,it-1,ilsd),1)
       ! ..DELTA u
       !$omp parallel do private(IG)
       DO ig=1,nhgb
          u(ig,it-1,ilsd)=rin0(ig)-rinm1(ig)
       ENDDO
       CALL dscal(nhgb*2,sn,u(1,it-1,ilsd),1)
       ! ..|U(n)>=G(1)|DELTAF(n)+|DELTA n(n)>
       CALL daxpy(nhgb*2,xmix,deltaf(1,it-1,ilsd),1,u(1,it-1,ilsd),1)
       IF (broy1%tadbroy) THEN
          ! ..Calculate W(IT-1)
          w(it-1,ilsd)=dotp(nhgb,f,f)
          CALL mp_sum(w(it-1,ilsd),parai%allgrp)
          w(it-1,ilsd)=1._real_8/SQRT(w(it-1,ilsd))
          w(it-1,ilsd)=MAX(1._real_8,w(it-1,ilsd))
       ELSE
          w(it-1,ilsd)=1._real_8
       ENDIF
       ! ..Calculate aij array
       CALL dsyrk('U','T',it-1,2*nhgb,2._real_8,deltaf(1,1,ilsd),2*nhgb,&
            0._real_8,a,broy1%ibreset)
       CALL dmatc('U',it-1,a,broy1%ibreset)
       IF (geq0)&
            CALL dger(it-1,it-1,-1._real_8,deltaf(1,1,ilsd),2*nhgb,&
            deltaf(1,1,ilsd),2*nhgb,a,broy1%ibreset)
       CALL summat(a,broy1%ibreset)
       DO i=1,it-1
          DO j=1,it-1
             a(j,i)=w(i,ilsd)*w(j,ilsd)*a(j,i)
          ENDDO
          a(i,i)=broy1%w02broy+a(i,i)
       ENDDO
       ! ..Invert A array. Output in A.
       CALL inversemat(it-1,a,broy1%ibreset,b,info)
       ! ..Compute C(k)=w(k)*<deltaf(k)|f(m)>        
       DO k=1,it-1
          c(k)=w(k,ilsd)*dotp(nhgb,deltaf(:,k,ilsd),f)
       ENDDO
       CALL mp_sum(c,it-1,parai%allgrp)
       ! ..Compute Gamma(L)=sum_k c(k) a(k,l)
       CALL dgemv('T',it-1,it-1,1._real_8,a,broy1%ibreset,c,1,0._real_8,gamma,1)
       ! ..Compute mixed density
       !$omp parallel do private(IG)
       DO ig=1,nhgb
          rinp(ig)=rin0(ig)+xmix*f(ig)
       ENDDO
       DO i=1,it-1
          aux=w(i,ilsd)*GAMMA(i)
          CALL daxpy(2*nhgb,-aux,u(1,i,ilsd),1,rinp,1)
       ENDDO
    ENDIF
    ! ..Copy the last nhg-nhgb elements of F to rinp
    IF (nhgb.LT.nhg)&
         CALL zcopy(nhg-nhgb,rout0(nhgb+1),1,rinp(nhgb+1),1)
    CALL tihalt('   BROYDEN',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(f,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__) ! TODO check stat
    DEALLOCATE(a,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(b,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(c,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gamma,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! FIXME when to deallocate this???
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE broyden
  ! ==================================================================
  SUBROUTINE give_scr_broyden(lbroyden,tag,nhgb)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lbroyden
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nhgb

! ==--------------------------------------------------------------==

    lbroyden=&
         2*nhgb+              & ! F
         broy1%ibreset*broy1%ibreset+     & ! A
         broy1%ibreset*broy1%ibreset+     & ! B
         broy1%ibreset+             & ! C
         broy1%ibreset              ! GAMMA
    tag='2*NHGB+IBRESET*(2*NHGB+...)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_broyden
  ! ==================================================================
  SUBROUTINE kerker(xmix,ymix,rin0,rout0,rinp,nhg,hg)
    ! ==--------------------------------------------------------------==
    ! == GENERATE NEXT ITERATION WITH KERKER DAMPING FACTOR           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: xmix, ymix
    COMPLEX(real_8)                          :: rin0(*), rout0(*), rinp(*)
    INTEGER                                  :: nhg
    REAL(real_8)                             :: hg(*)

    INTEGER                                  :: ig

    !$omp parallel do private(IG)
    DO ig=1,nhg
       rinp(ig)=rin0(ig)+xmix*(rout0(ig)-rin0(ig))&
            *hg(ig)/(hg(ig)+ymix*ymix)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kerker
  ! ==================================================================

END MODULE broyden_utils
