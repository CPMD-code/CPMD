MODULE davidson_utils
  USE cvan,                            ONLY: dvan,&
                                             qq
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE gsortho_utils,                   ONLY: gs_ortho
  USE hfx_drivers,                     ONLY: hfxpsi
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE soft,                            ONLY: soft_com
  USE sort_utils,                      ONLY: sort2
  USE spsi_utils,                      ONLY: give_scr_spsi,&
                                             spsi
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             fpar,&
                                             maxsys,&
                                             ncpw
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: dspevy
  USE vgsortho_utils,                  ONLY: vgsortho
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: davidson
  PUBLIC :: give_scr_davidson

CONTAINS

  SUBROUTINE check_added_vectors(c0, ncurr, nnew)
    COMPLEX(real_8)                          :: c0(:,:)
    INTEGER, INTENT(IN)                      :: ncurr, nnew

    CHARACTER(*), PARAMETER :: procedureN = 'check_added_vectors'
    REAL(real_8), PARAMETER                  :: droptol = 1.0E-6_real_8

    INTEGER                                  :: i, ierr, j, seed_size
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: seed
    REAL(real_8)                             :: norm, r1, r2
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: smat

    ALLOCATE(smat(1,ncurr+nnew),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(smat)

    CALL RANDOM_SEED(size=seed_size)
    ALLOCATE(seed(seed_size),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DO i=1,seed_size
       seed(i)=i*parai%me
    END DO
    CALL RANDOM_SEED(put=seed)

    DO i = 1, nnew

       CALL ovlap2(ncpw%ngw,1,ncurr+i-1,smat,c0(:,ncurr+i),c0(:,1:ncurr+i-1),.TRUE.)
       CALL mp_sum(smat,SIZE(smat),parai%allgrp)
       CALL dgemm('N','T',2*ncpw%ngw,1,ncurr+i-1,-1._real_8,c0,2*ncpw%ngw,smat,1,1._real_8,c0(1,ncurr+i),2*ncpw%ngw)

       ! check norm of ncurr+i
       norm=dotp(ncpw%ngw,c0(:,ncurr+i),c0(:,ncurr+i))
       CALL mp_sum(norm,parai%allgrp)
       norm=SQRT(norm)
       !TL IF (paral%io_parent) WRITE(*,*)"NORM of vector (",i,") is : ",norm
       DO WHILE (norm < droptol) 
          IF (paral%io_parent) WRITE(*,*)"NORM of vector (",i,") is lower than droptol! Resetting with random components..",norm
          DO j = 1, ncpw%ngw
             CALL RANDOM_NUMBER(r1)
             CALL RANDOM_NUMBER(r2)
             c0(j,ncurr+i)= CMPLX(1.0e-5_real_8 * ( r1 - 0.5_real_8 ), 1.0e-5_real_8 * ( r2 - 0.5_real_8 ), KIND=real_8)
          END DO
          ! renormalize vector
          norm=dotp(ncpw%ngw,c0(:,ncurr+i),c0(:,ncurr+i))
          CALL mp_sum(norm,parai%allgrp)
          norm=SQRT(norm)
          c0(:,ncurr+i) = c0(:,ncurr+i) / norm
          !
          CALL ovlap2(ncpw%ngw,1,ncurr+i-1,smat,c0(:,ncurr+i),c0(:,1:ncurr+i-1),.TRUE.)
          CALL mp_sum(smat,SIZE(smat),parai%allgrp)
          CALL dgemm('N','T',2*ncpw%ngw,1,ncurr+i-1,-1._real_8,c0,2*ncpw%ngw,smat,1,1._real_8,c0(1,ncurr+i),2*ncpw%ngw)
          norm=dotp(ncpw%ngw,c0(:,ncurr+i),c0(:,ncurr+i))
          CALL mp_sum(norm,parai%allgrp)
          norm=SQRT(norm)
       END DO

       ! renormalize new vector
       norm=dotp(ncpw%ngw,c0(:,ncurr+i),c0(:,ncurr+i))
       CALL mp_sum(norm,parai%allgrp)
       norm=SQRT(norm)
       c0(:,ncurr+i) = c0(:,ncurr+i) / norm

    END DO

    DEALLOCATE(seed,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

  END SUBROUTINE check_added_vectors

  ! ==================================================================
  SUBROUTINE davidson(ndiag,nadd,c0,cs,cr,sc0,cscr,vpot,&
       cgs,focc,nstate, psi,edav,vpp)
    ! ==--------------------------------------------------------------==
    ! ==  NDIAG      NUMBER OF STATES                                 ==
    ! ==  NADD       DIMENSION OF DAVIDSON MATRIX                     ==
    ! ==  C0         WAVEFUNCTIONS                                    ==
    ! ==  CS         SIGMA VECTORS                                    ==
    ! ==  CR         RESIDUALS                                        ==
    ! ==  SC0        S*C0 (only Vanderbilt)                           ==
    ! ==  CSCR       SCRATCH                                          ==
    ! ==  EDAV       EIGENVALUES OF ADAV                              ==
    ! ==--------------------------------------------------------------==
    ! ==  SPECIAL INPUT NEEDED FOR HARTREE-FOCK EXCHANGE              ==
    ! ==  CGS        GROUND STATE ORBITALS -> DEFINING VPOT           ==
    ! ==  FOCC       OCCUPATION NUMBERS OF CGS                        ==
    ! ==  NSTATE     NUMBER OF CGS ORBITALS                           ==
    ! ==--------------------------------------------------------------==
    ! ==  ADAV       DAVIDSON MATRIX                                  ==
    ! ==  UDAV       EIGENVECTORS OF ADAV                             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndiag, nadd
    COMPLEX(real_8) :: c0(ncpw%ngw,nadd), cs(ncpw%ngw,nadd), &
      cr(ncpw%ngw,nadd), sc0(ncpw%ngw,*), cscr(ncpw%ngw,nadd)
    REAL(real_8)                             :: vpot(:)
    COMPLEX(real_8)                          :: cgs(:,:)
    REAL(real_8)                             :: focc(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: edav(nadd), vpp(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'davidson'
    INTEGER, PARAMETER                       :: ispin = 1

    INTEGER                                  :: i, icycle, ierr, ig, iopt, &
                                                iprie, ir, is, isub, ixabs, &
                                                izamax, j, k, nconv, ncurr, &
                                                ndima, nnew, nrot2
    INTEGER, ALLOCATABLE                     :: INDEX(:)
    LOGICAL                                  :: tlsd2
    REAL(real_8)                             :: e1, e1n, e1z, eimax, ermax, &
                                                norm, sign, tcpu, time1
    REAL(real_8), ALLOCATABLE                :: adav(:), aux(:), udav(:,:)

    IF (ndiag.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('  DAVIDSON',isub)
    ALLOCATE(INDEX(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(adav(nadd*(nadd+1)/2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(udav(nadd, nadd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(aux(6*nadd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    sign=1._real_8
    IF (cntl%tlsd) sign=2._real_8
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    time1=m_walltime()
    ! ==--------------------------------------------------------------==
    ! ==  SHIFT POTENTIAL                                             ==
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       vpot(ir)=vpot(ir)+cntr%dshift
    ENDDO
    !$omp parallel do private(IG)
    DO ig=1,ncpw%ngw
       vpp(ig)=vpp(ig)+cntr%dshift
    ENDDO
    IF (pslo_com%tivan)&
         CALL daxpy(maxsys%nhxs*maxsys%nhxs*maxsys%nsx,cntr%dshift,qq(1,1,1),1,dvan(1,1,1),1)
    ! ==--------------------------------------------------------------==
    ! ==  INITIALIZE SIGMA VECTORS                                    ==
    ! ==--------------------------------------------------------------==
    CALL hpsi(c0,cs,sc0,vpot,psi,ndiag,1,ispin)
    CALL hfxpsi(cgs,c0,cs,focc,sign,psi,nstate,ndiag)
    ! ==--------------------------------------------------------------==
    ! ==  DAVIDSON MATRIX                                             ==
    ! ==--------------------------------------------------------------==
    iopt=1
    k=0
    DO j=1,ndiag
       DO i=j,ndiag
          k=k+1
          adav(k)=dotp(ncpw%ngw,c0(:,i),cs(:,j))
       ENDDO
    ENDDO
    CALL mp_sum(adav,(ndiag*(ndiag+1))/2,parai%allgrp)
    CALL dspevy(iopt,adav,edav,udav,nadd,ndiag,aux,SIZE(aux))
    ! ==--------------------------------------------------------------==
    ! ==  RESIDUALS                                                   ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(cr)!,ngw*ndiag)
    IF (pslo_com%tivan) THEN            ! cmb- avoid IF inside loop for good vectorization
       DO i=1,ndiag
          DO j=1,ndiag
             DO ig=1,ncpw%ngw
                cr(ig,i)=cr(ig,i)+(cs(ig,j)-sc0(ig,j)*edav(i))*udav(j,i)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO i=1,ndiag
          DO j=1,ndiag
             DO ig=1,ncpw%ngw
                cr(ig,i)=cr(ig,i)+(cs(ig,j)-c0(ig,j)*edav(i))*udav(j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  ROTATE ORBITALS AND SIGMA VECTORS                           ==
    ! ==--------------------------------------------------------------==
    CALL dgemm('N','N',2*ncpw%ngw,ndiag,ndiag,1.0_real_8,c0(1,1),2*ncpw%ngw,&
         udav(1,1),nadd,0.0_real_8,cscr(1,1),2*ncpw%ngw)
    CALL dcopy(2*ncpw%ngw*ndiag,cscr(1,1),1,c0(1,1),1)
    CALL dgemm('N','N',2*ncpw%ngw,ndiag,ndiag,1.0_real_8,cs(1,1),2*ncpw%ngw,&
         udav(1,1),nadd,0.0_real_8,cscr(1,1),2*ncpw%ngw)
    CALL dcopy(2*ncpw%ngw*ndiag,cscr(1,1),1,cs(1,1),1)
    IF (pslo_com%tivan) THEN
       CALL dgemm('N','N',2*ncpw%ngw,ndiag,ndiag,1.0_real_8,sc0(1,1),2*ncpw%ngw,&
            udav,nadd,0.0_real_8,cscr(1,1),2*ncpw%ngw)
       CALL dcopy(2*ncpw%ngw*ndiag,cscr(1,1),1,sc0(1,1),1)
    ENDIF
    ncurr=ndiag
    ndima=ndiag
    nconv=0
    ! ==--------------------------------------------------------------==
    ! ==          THE BASIC LOOP FOR DAVIDSON DIAGONALIZATION         ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/,A,A)') ' ITER  STATES  SUBSPACE    STATE     ENERGY ',&
         '    RESIDUAL     TCPU  '
    DO icycle=0,cnti%n_fries
       DO i=1,ndiag
          ixabs=izamax(ncpw%ngw,cr(1,i),1)
          eimax=ABS(cr(ixabs,i))
          CALL mp_max(eimax,parai%allgrp)
          IF (eimax.GT.cntr%epsdav) THEN
             nconv=i-1
             GOTO 200
          ENDIF
       ENDDO
       nconv=ndiag
200    CONTINUE
       ixabs=izamax(ncpw%ngw*ndiag,cr,1)
       j = ( ixabs - MOD( ixabs - 1, ncpw%ngw ) - 1 ) / ncpw%ngw + 1
       i = ixabs - ( j - 1 ) * ncpw%ngw
       ermax=ABS(cr(i,j))
       CALL mp_max(ermax,parai%allgrp)
       tcpu=(m_walltime()-time1)*0.001_real_8
       iprie=nconv+1
       IF (nconv.EQ.ndiag) iprie=ndiag
       IF (paral%io_parent)&
            WRITE(6,'(I5,3X,I4,4X,I4,7X,I4,F13.6,2X,E10.2,F9.2)')&
            icycle,nconv,ndima,iprie,edav(iprie)-cntr%dshift,ermax,tcpu
       time1=m_walltime()
       IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
       IF (nconv.EQ.ndiag .OR. soft_com%exsoft) GOTO 100
       ! ==--------------------------------------------------------------==
       ! ==  ADD NEW VECTORS TO THE SUBSPACE                             ==
       ! ==--------------------------------------------------------------==
       nnew=ndiag-nconv
       nnew=MIN(nnew,nadd-ncurr)
       ndima=ncurr+nnew-nconv
       DO i=ncurr+1,ncurr+nnew
          is=nconv+i-ncurr
          DO ig=1,ncpw%ngw
             cscr(ig,1)=1.0_real_8/(vpp(ig)-edav(is))
             c0(ig,i)=c0(ig,is) !*cscr(ig,1)
          ENDDO
          e1z=dotp(ncpw%ngw,c0(:,i),cr(:,is))
          e1n=dotp(ncpw%ngw,c0(:,i),c0(:,is))
          CALL mp_sum(e1z,parai%allgrp)
          CALL mp_sum(e1n,parai%allgrp)
          e1=e1z/e1n
          IF (pslo_com%tivan) THEN
             !$omp parallel do private(IG)
             DO ig=1,ncpw%ngw
                c0(ig,i)=(cr(ig,is)-e1*sc0(ig,is))*cscr(ig,1)
             ENDDO
          ELSE
             !$omp parallel do private(IG)
             DO ig=1,ncpw%ngw
                c0(ig,i)=(cr(ig,is)-e1*c0(ig,is))*cscr(ig,1)
             ENDDO
          ENDIF
       ENDDO

       ! renormalize new vector to avoid too small/large entries
       DO i=ncurr+1,ncurr+nnew
          norm=dotp(ncpw%ngw,c0(:,i),c0(:,i))
          CALL mp_sum(norm,parai%allgrp)
          norm=SQRT(norm)
          c0(:,i) = c0(:,i) / norm
       ENDDO

       CALL check_added_vectors(c0,ncurr,nnew)

       ! ==--------------------------------------------------------------==
       ! ==  ORTHOGONALIZE NEW VECTORS                                   ==
       ! ==--------------------------------------------------------------==
       IF (pslo_com%tivan) THEN
          CALL rnlsm(c0(:,ncurr+1:ncurr+nnew),nnew,1,1,.FALSE.)
          CALL dcopy(2*ncpw%ngw*nnew,c0(1,ncurr+1),1,sc0(1,ncurr+1),1)
          CALL spsi(nnew,sc0(1,ncurr+1))
          CALL vgsortho(c0,sc0,ncpw%ngw,ncurr+1,ncurr+nnew)
       ELSE
          CALL gs_ortho(c0,ncurr,c0(:,ncurr+1:ncurr+nnew),nnew)
       ENDIF

       ! ==--------------------------------------------------------------==
       ! ==  NEW SIGMA VECTORS                                           ==
       ! ==--------------------------------------------------------------==
       CALL hpsi(c0(:,ncurr+1:ncurr+nnew),cs(:,ncurr+1:ncurr+nnew),sc0(1,ncurr+1),&
            vpot,psi,nnew,1,ispin)
       CALL hfxpsi(cgs,c0(:,ncurr+1:ncurr+nnew),cs(:,ncurr+1:ncurr+nnew),focc,sign,&
            psi,nstate,nnew)
       ! ==--------------------------------------------------------------==
       ! ==  DAVIDSON MATRIX                                             ==
       ! ==--------------------------------------------------------------==
       iopt=1
       k=0
       DO j=nconv+1,ncurr+nnew
          DO i=j,ncurr+nnew
             k=k+1
             adav(k)=dotp(ncpw%ngw,c0(:,i),cs(:,j))
          ENDDO
       ENDDO
       CALL mp_sum(adav,k,parai%allgrp)
       CALL dspevy(iopt,adav,edav(nconv+1),udav,nadd,ndima,aux,SIZE(aux))
       ! ==--------------------------------------------------------------==
       ! ==  RESIDUALS                                                   ==
       ! ==--------------------------------------------------------------==
       CALL zeroing(cr)!,ngw*ndiag)
       IF (pslo_com%tivan) THEN        ! cmb- avoid IF inside loop for good vectorization
          DO i=nconv+1,ndiag
             DO j=nconv+1,ncurr+nnew
                DO ig=1,ncpw%ngw
                   cr(ig,i)=cr(ig,i)+(cs(ig,j)-sc0(ig,j)*edav(i))*&
                        udav(j-nconv,i-nconv)
                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO i=nconv+1,ndiag
             DO j=nconv+1,ncurr+nnew
                DO ig=1,ncpw%ngw
                   cr(ig,i)=cr(ig,i)+(cs(ig,j)-c0(ig,j)*edav(i))*&
                        udav(j-nconv,i-nconv)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       ! ==--------------------------------------------------------------==
       ! ==  ROTATE ORBITALS AND SIGMA VECTORS                           ==
       ! ==--------------------------------------------------------------==
       nrot2=ndima
       IF (ndima+ndiag.GT.nadd) nrot2=ndiag-nconv
       CALL dgemm('N','N',2*ncpw%ngw,nrot2,ndima,1.0_real_8,c0(1,nconv+1),&
            2*ncpw%ngw,udav(1,1),nadd,0.0_real_8,cscr(1,1),2*ncpw%ngw)
       CALL dcopy(2*ncpw%ngw*nrot2,cscr(1,1),1,c0(1,nconv+1),1)
       CALL dgemm('N','N',2*ncpw%ngw,nrot2,ndima,1.0_real_8,cs(1,nconv+1),&
            2*ncpw%ngw,udav(1,1),nadd,0.0_real_8,cscr(1,1),2*ncpw%ngw)
       CALL dcopy(2*ncpw%ngw*nrot2,cscr(1,1),1,cs(1,nconv+1),1)
       IF (pslo_com%tivan) THEN
          CALL dgemm('N','N',2*ncpw%ngw,nrot2,ndima,1.0_real_8,sc0(1,nconv+1),&
               2*ncpw%ngw,udav(1,1),nadd,0.0_real_8,cscr(1,1),2*ncpw%ngw)
          CALL dcopy(2*ncpw%ngw*nrot2,cscr(1,1),1,sc0(1,nconv+1),1)
       ENDIF

       !if (ndima+ndiag.gt.nadd) then
       !  ! when the buffer is full, add random noise.
       !  if(pslo_com%tivan) call stopgm(procedureN,'NYI', &
       !      __LINE__,__FILE__)
       !  call random_seed(size=seed_size)
       !  allocate(seed(seed_size))
       !  call random_seed(get=seed)
       !  do i=1,seed_size
       !    seed(i)=seed(i)+10**i*parai%me
       !  enddo
       !  call random_seed(put=seed)
       !  do i=nconv+1,nadd
       !     do ig=1,ngw
       !       call random_number(r)
       !       c0(ig,i) = c0(ig,i) + 1.0e-5_real_8 * ( r - 0.5_real_8 )
       !       !call random_number(r)
       !       !cs(ig,i) = cs(ig,i) + 1.0e-5_real_8 * ( r - 0.5_real_8 )
       !     enddo
       !  enddo
       !  do i=nconv+1,ndiag
       !     do ig=1,ngw
       !       call random_number(r)
       !       cr(ig,i) = cr(ig,i) + 1.0e-8_real_8 * ( r - 0.5_real_8 )
       !     enddo
       !  enddo
       !  call gs_ortho(c0,nconv,c0(1,nconv+1),nadd-nconv-1)
       !  deallocate(seed)
       !endif

       ncurr=nconv+nrot2
       ! ==--------------------------------------------------------------==
       ! ==                      END OF MAIN LOOP                        ==
       ! ==--------------------------------------------------------------==
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' ! !!!!!! NOT ALL ROOTS ARE CONVERGED !!!!!!! '
100 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==  ORDER THE STATES WITH RESPECT TO ENERGY                     ==
    ! ==--------------------------------------------------------------==
    CALL sort2(edav,ndiag,index)
    CALL dcopy(2*ncpw%ngw*ndiag,c0(1,1),1,cs(1,1),1)
    DO i=1,ndiag
       j=INDEX(i)
       CALL dcopy(2*ncpw%ngw,cs(1,j),1,c0(1,i),1)
    ENDDO
    cntl%tlsd=tlsd2
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       vpot(ir)=vpot(ir)-cntr%dshift
    ENDDO
    !$omp parallel do private(IG)
    DO ig=1,ncpw%ngw
       vpp(ig)=vpp(ig)-cntr%dshift
    ENDDO
    !$omp parallel do private(I)
    DO i=1,ndiag
       edav(i)=edav(i)-cntr%dshift
    ENDDO
    IF (pslo_com%tivan)&
         CALL daxpy(maxsys%nhxs*maxsys%nhxs*maxsys%nsx,-cntr%dshift,qq(1,1,1),1,dvan(1,1,1),1)
    ! 
    ! ==--------------------------------------------------------------==
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(adav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(udav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('  DAVIDSON',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE davidson
  ! ==================================================================
  SUBROUTINE give_scr_davidson(ldavidson,tag,ndiag,nadd)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldavidson
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ndiag, nadd

    INTEGER                                  :: ldspevy, lhpsi, lrnlsm, &
                                                lscr2, lspsi

    lscr2 =  ndiag+        & ! INDEX
         nadd*(nadd+1)/2+     & ! ADAV
         nadd*nadd            ! UDAV
    CALL give_scr_hpsi(lhpsi,tag,ndiag)
    ldspevy=3*ndiag
    IF (pslo_com%tivan) THEN
       CALL give_scr_rnlsm(lrnlsm,tag,ndiag,.FALSE.)
    ELSE
       lrnlsm=0
    ENDIF
    CALL give_scr_spsi(lspsi,tag)
    ldavidson=lscr2+MAX(lhpsi,ldspevy,lrnlsm,lspsi)+200
    tag='LSCR+MAX(HPSI,DIAG,RNLSM,SPSI)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_davidson
  ! ==================================================================

END MODULE davidson_utils
