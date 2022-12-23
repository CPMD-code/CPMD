MODULE friesner_c_utils
  USE adjmu_utils,                     ONLY: convfrie
  USE ehpsi_utils,                     ONLY: ehpsi_c
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fint,                            ONLY: fint1
  USE frsblk_c_utils,                  ONLY: reorder_c
  USE gsortho_utils,                   ONLY: gs_ortho_c
  USE hpsi_utils,                      ONLY: hpsi
  USE jacobi_utils,                    ONLY: jacobi
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE randtowf_utils,                  ONLY: randtowf
  USE sort_utils,                      ONLY: sort2
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: friesner_c
  PUBLIC :: krylov_ref_c

CONTAINS

  ! ==================================================================
  SUBROUTINE friesner_c(ndiag,nel,nconv,nhpsi,c0,cs,sc0,cscr,&
       vpot,psi,edav,&
       jspin,ikind,ikk,trefine,tinfo)
    ! ==--------------------------------------------------------------==
    ! ==  Kohn-Sham Matrix Diagonalization by Krylov-Space Method     ==
    ! ==  W.T. Pollard and R.A. Friesner                              ==
    ! ==  J. Chem. Phys. 99, 6742 (1993)                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndiag
    REAL(real_8)                             :: nel
    INTEGER                                  :: nconv, nhpsi
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                cs(nkpt%ngwk,*), &
                                                sc0(nkpt%ngwk,*), &
                                                cscr(nkpt%ngwk,*)
    REAL(real_8)                             :: vpot(*)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: edav(*)
    INTEGER                                  :: jspin, ikind, ikk
    LOGICAL                                  :: trefine, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'friesner_c'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)
    INTEGER, PARAMETER                       :: ispin = 1 
    REAL(real_8), PARAMETER :: big = 1.e33_real_8 , edavmax = 300._real_8, &
      edavmin = -300._real_8, eps = 1.e-14_real_8

    CHARACTER(len=100)                       :: formlanc
    CHARACTER(len=20)                        :: charspin
    COMPLEX(real_8), ALLOCATABLE             :: adav(:,:), work(:)
    INTEGER :: i, iaux, icycle, ierr, iformlanc, ig, isub, itest, j, lwork, &
      nconvold, ncurr, ngwk2, nhpsiold, nleft, ntest
    INTEGER, ALLOCATABLE                     :: INDEX(:)
    LOGICAL                                  :: tcheck, tconv, tdebug, tlsd2, &
                                                ttest
    REAL(real_8)                             :: amu, aux, b2, b2max, b2min, &
                                                fac, tim, tim1, tim2, &
                                                trotmax, trotmin
    REAL(real_8), ALLOCATABLE                :: alpha0(:), beta1(:), rwork(:)
    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==
! With LSD and 1 electron...

    IF (ndiag.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('FRIESNER_C',isub)
    ! TODO check stat
    ! TODO align for BG
    ALLOCATE(INDEX(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(adav(ndiag, ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(alpha0(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(beta1(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    tim1=m_walltime()
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    tconv=.FALSE.
    tdebug=.FALSE.
    ngwk2=2*nkpt%ngwk
    nconv=0
    nconvold=0
    nhpsiold=0
    IF (fint1%ttrot) THEN
       trotmin=EXP(-edavmax*fint1%betap)
       trotmax=EXP(-edavmin*fint1%betap)
    ENDIF
    CALL zeroing(index)!,ndiag)
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.tinfo) THEN
       IF (jspin.EQ.0) THEN
          charspin='<<<<<<<<<<<<<<<<<<<<'
       ELSEIF (jspin.EQ.1) THEN
          charspin='<<<<<<< ALPHA SPIN<<'
       ELSEIF (jspin.EQ.2) THEN
          charspin='<<<<<<<< BETA SPIN<<'
       ENDIF
       iformlanc=INT(LOG10(REAL(nkpt%nkpts,kind=real_8)))+1
       IF (paral%io_parent)&
            WRITE(formlanc,'(A,I1,A,I1,A,I2,A)')&
            '(/,1X,"<<",I',iformlanc,',":",I',iformlanc,",",&
            16-2*iformlanc,'("<")," LANCZOS DIAGONALIZATION ",A20)'
       IF (paral%io_parent)&
            WRITE(6,formlanc) ikk,nkpt%nkpts,charspin
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Subspace Matrix
    tcheck=.TRUE.
    ntest=0

    ! Workspace query call: determine optimal LWORK        
    ALLOCATE(rwork(MAX(1, 3*ndiag-2)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! for workspace query call
    lwork = -1
    CALL zheev('V','U',ndiag,adav,ndiag,edav,&
         work,lwork,rwork,ierr)
    lwork = work(1)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__) ! was allocated (1) for query call
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    DO WHILE(tcheck)
       IF (fint1%ttrot) THEN
          CALL ehpsi_c(c0,cs,vpot,psi,ndiag,ikind,maxfft)
       ELSE
          CALL hpsi(c0(:,1:ndiag),cs,sc0,vpot,psi,ndiag,ikind,ispin)
       ENDIF
       nhpsi=ndiag
       CALL ovlap_h(ndiag,adav,c0,cs)
       CALL sumhmat(adav,ndiag)
       CALL zheev('V','U',ndiag,adav,ndiag,edav,&
            work,lwork,rwork,ierr)
       ! 
       CALL reorder_c(ndiag,ndiag,edav,adav)
       ! Rotate Orbitals
       CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,c0(1,1),nkpt%ngwk,&
            adav(1,1),ndiag,zzero,cscr(1,1),nkpt%ngwk)
       CALL dcopy(ngwk2*ndiag,cscr(1,1),1,c0(1,1),1)
       tcheck=.FALSE.
       itest=0
       DO i=1,ndiag
          IF (fint1%ttrot) THEN
             ttest = (edav(i).GT.trotmax).OR.(edav(i).LT.trotmin)
          ELSE
             ttest = ABS(edav(i)).LT.eps
          ENDIF
          IF (ttest) THEN
             tcheck=.TRUE.
             IF (paral%io_parent)&
                  WRITE(6,*) 'FRIESNER_C| EIGENVECTOR',i,' IS VERY BAD! '
             ! Not useful: Now we use H|C0>.
             CALL randtowf(c0(:,i:i),1,ikind,ikk)
             itest=itest+1
             edav(i)=big
          ENDIF
       ENDDO
       IF (tcheck) THEN
          ntest=ntest+1
          IF (ntest.GT.10.AND.paral%parent)&
               CALL stopgm(' FRIESNER_C', 'CAN''T FIND GUEST VECTORS',& 
               __LINE__,__FILE__)
          CALL sort2(edav,ndiag,index)
          ! Put at the end randomized wavefunctions.
          CALL dcopy(ngwk2*ndiag,c0,1,cscr,1)
          DO i=1,ndiag
             j=INDEX(i)
             IF (edav(i).EQ.big) THEN
                CALL dcopy(ngwk2,cs(1,j),1,c0(1,i),1)
             ELSE
                CALL dcopy(ngwk2,cscr(1,j),1,c0(1,i),1)
             ENDIF
          ENDDO
          ! Orthogonalize.
          CALL gs_ortho_c(c0,ndiag-itest,&
               c0(1,ndiag-itest+1),itest,adav)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Rotate also H|C0> (IF TCHECK=.FALSE., ADAV unchanged).
    CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,cs(1,1),nkpt%ngwk,&
         adav(1,1),ndiag,zzero,cscr(1,1),nkpt%ngwk)
    CALL dcopy(ngwk2*ndiag,cscr(1,1),1,cs(1,1),1)
    tim2=m_walltime()
    tim=(tim2-tim1)*0.001_real_8
    IF (paral%io_parent.AND.tinfo) THEN
       WRITE(6,'(A,T58,F8.2)')&
            ' >> TIME FOR INITIAL SUBSPACE DIAGONALIZATION:  ',tiM
       WRITE(6,'(A,8X,A,8X,A)') ' >> CYCLE     NCONV',&
            'B2MAX','B2MIN     #HPSI      TIME'
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL zeroing(index)!,ndiag)
    DO icycle=1,cnti%n_fries
       tim1=m_walltime()
       b2max=0._real_8
       b2min=1.e30_real_8
       ! We calculate again the residue for all vectors
       nconv=0
       ncurr=nconv+1
       ! First Lanczos Step
       DO i=ncurr,ndiag
          DO ig=1,nkpt%ngwk
             cscr(ig,1)=cs(ig,i)-edav(i)*c0(ig,i)
          ENDDO
          b2=ddot(ngwk2,cscr,1,cscr,1)
          CALL mp_sum(b2,parai%allgrp)
          beta1(i)=SQRT(b2)
          alpha0(i)=edav(i)
          IF (b2.LT.cntr%b2limit) THEN
             nconv=nconv+1
             INDEX(i)=nconv
             IF (b2.GT.b2max) b2max=b2
             IF (b2.LT.b2min) b2min=b2
          ELSE
             CALL dcopy(ngwk2,cscr,1,cs(1,i),1)
             fac=1._real_8/beta1(i)
             CALL dscal(ngwk2,fac,cs(1,i),1)
          ENDIF
       ENDDO
       ! Order states: converged first
       DO i=ncurr,ndiag
          j=INDEX(i)
          IF (j.NE.0.AND.j.NE.i) THEN
             CALL dswap(ngwk2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngwk2,cs(1,i),1,cs(1,j),1)
             INDEX(j)=j
             INDEX(i)=0
             aux = beta1(i)
             beta1(i) = beta1(j)
             beta1(j) = aux
             aux = alpha0(i)
             alpha0(i) = alpha0(j)
             alpha0(j) = aux
          ENDIF
       ENDDO
       ncurr=nconv+1
       ! Lanczos Refinement
       IF (icycle.EQ.1.AND.ncurr.GT.ndiag) THEN
          trefine=.FALSE.
       ELSE
          trefine=.TRUE.
       ENDIF
       CALL krylov_ref_c(ncurr,ndiag,nconv,nhpsi,c0,cs,sc0,cscr,&
            vpot,psi,index,alpha0,beta1,&
            b2min,b2max,ikind)
       ! Order states : converged first
       DO i=ncurr,ndiag
          j=INDEX(i)
          IF (j.NE.0.AND.j.NE.i) THEN
             CALL dswap(ngwk2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngwk2,cs(1,i),1,cs(1,j),1)
             INDEX(j)=j
             INDEX(i)=0
             ! Swap BETA1
             aux=beta1(j)
             beta1(j)=beta1(i)
             beta1(i)=aux
          ENDIF
       ENDDO
       ! Reorthogonalize
       IF (ncurr.LE.ndiag) THEN
          CALL gs_ortho_c(c0,ncurr-1,c0(1,ncurr),ndiag-ncurr+1,adav)
       ENDIF
       ! Calculate new forces; only for states that entered Lanczos
       nleft=ndiag-ncurr+1
       IF (ncurr.LE.ndiag) THEN
          IF (fint1%ttrot) THEN
             CALL ehpsi_c(c0(1,ncurr),cs(1,ncurr),vpot,psi,&
                  nleft,ikind,maxfft)
          ELSE
             CALL hpsi(c0(:,ncurr:ncurr+nleft-1),cs(1,ncurr),sc0(1,ncurr),&
                  vpot,psi,nleft,ikind,ispin)
          ENDIF
          nhpsi=nhpsi+nleft
       ENDIF
       ! Diagonalize Kohn-Sham Matrix
       CALL ovlap_h(ndiag,adav,c0,cs)
       CALL sumhmat(adav,ndiag)
       ! We inverse the order
       ! (not the same one for ZHEEV than for this routine).
       CALL zheev('V','U',ndiag,adav,ndiag,edav,&
            work,lwork,rwork,ierr)
       CALL reorder_c(ndiag,ndiag,edav,adav)
       CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,c0(1,1),nkpt%ngwk,&
            adav(1,1),ndiag,zzero,cscr(1,1),nkpt%ngwk)
       CALL dcopy(ngwk2*ndiag,cscr(1,1),1,c0(1,1),1)
       CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,cs(1,1),nkpt%ngwk,&
            adav(1,1),ndiag,zzero,cscr(1,1),nkpt%ngwk)
       CALL dcopy(ngwk2*ndiag,cscr(1,1),1,cs(1,1),1)
       ! If one eigenvalue is equal to zero (very rare: bug ??).
       cnti%ntrans=nconv
       DO i=1,nconv
          IF (fint1%ttrot) THEN
             ttest = (edav(i).GT.trotmax).OR.(edav(i).LT.trotmin)
          ELSE
             ttest = ABS(edav(i)).LT.eps
          ENDIF
          IF (ttest) THEN
             nconv=nconv-1
             cnti%ntrans=cnti%ntrans+1
             IF (paral%io_parent)&
                  WRITE(6,*) ' FRIESNER_C| ','VERY BAD EIGENVALUE !!'
             IF (cnti%ntrans.GT.ndiag) THEN
                CALL randtowf(c0(:,i:i),1,ikind,ikk)
             ELSE
                CALL dswap(ngwk2,c0(1,i),1,c0(1,cnti%ntrans),1)
                CALL dswap(ngwk2,cs(1,i),1,cs(1,cnti%ntrans),1)
             ENDIF
          ENDIF
       ENDDO
       tim2=m_walltime()
       tim=(tim2-tim1)*0.001_real_8
       IF (paral%io_parent.AND.tinfo) THEN
          WRITE(6,'(4X,I5,5X,I5,2(1PE13.3),0PF10.2,0PF10.2)')&
               icycle,nconv,b2max,b2min,&
               REAL(nhpsi-nhpsiold,kind=real_8)/REAL(ndiag,kind=real_8),tiM
       ENDIF
       nhpsiold=nhpsi
       ! Test if the number of eigenvalues is enough to stop.
       IF ((.NOT.fint1%tfral).AND.cntl%tfint.AND.&
            (nconv.NE.nconvold).AND.&
            (nconv.NE.ndiag)) THEN
          ! Spin or not.
          IF (tlsd2.AND.(nconv.GE.nel).OR.&
               (.NOT.tlsd2.AND.(nconv.GE.(nel+1)/2)) ) THEN
             tconv=.TRUE.
             nconvold=nconv
             CALL convfrie(ndiag,nconv,nel,edav,cntl%tlsd,amu,tconv)
          ENDIF
       ENDIF
       IF (tconv.OR.nconv.EQ.ndiag) THEN
          GOTO 200
       ENDIF
       ! End of Krylov Space Cycles
    ENDDO
200 CONTINUE
    IF (paral%io_parent.AND.((.NOT.tconv).AND.nconv.NE.ndiag).AND.tdebug)&
         THEN
       WRITE(6,'(1X,64("!"))')
       WRITE(6,'(" !!",A,T64,"!!")')&
            ' FRIESNER_C| NOT ALL ROOTS ARE CONVERGED'
       WRITE(6,'(1X,64("!"))')
    ENDIF
    ! Order Roots with respect to energy
    CALL sort2(edav,ndiag,index)
    CALL dcopy(ngwk2*ndiag,c0(1,1),1,cs(1,1),1)
    IF (.NOT.fint1%ttrot) CALL dscal(ndiag,-1._real_8,edav(1),1)
    DO i = 1 , ndiag/2
       aux = edav(ndiag-i+1)
       edav(ndiag-i+1) = edav(i)
       edav(i) = aux
       iaux = INDEX(ndiag-i+1)
       INDEX(ndiag-i+1) = INDEX(i)
       INDEX(i) = iaux
    ENDDO
    DO i=1,ndiag
       j=INDEX(i)
       CALL dcopy(ngwk2,cs(1,j),1,c0(1,i),1)
    ENDDO
    cntl%tlsd=tlsd2
    CALL tihalt('FRIESNER_C',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(adav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(alpha0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(beta1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE friesner_c
  ! ==================================================================
  SUBROUTINE krylov_ref_c(ncurr,ndiag,nconv,nhpsi,c0,cs,sc0,cscr,&
       vpot,psi,index,alpha0,beta1,&
       b2min,b2max,ikind)
    ! ==--------------------------------------------------------------==
    ! ==  LANCZOS REFINEMENT LOOPS                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ncurr, ndiag, nconv, nhpsi
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*), &
                                                cs(nkpt%ngwk,*), &
                                                sc0(nkpt%ngwk,*), &
                                                cscr(nkpt%ngwk,*)
    REAL(real_8)                             :: vpot(*)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: INDEX(*)
    REAL(real_8)                             :: alpha0(*), beta1(*), b2min, &
                                                b2max
    INTEGER                                  :: ikind

    CHARACTER(*), PARAMETER                  :: procedureN = 'krylov_ref_c'
    INTEGER, PARAMETER                       :: ispin = 1 

    COMPLEX(real_8), EXTERNAL                :: zdotc
    INTEGER                                  :: i, ib, ibn0, ibn1, ibn2, &
                                                ierr, ii, ik, ikk, iovlm, j, &
                                                jkk, jl, nkb0, nkb1, nkbeff, &
                                                nkbeffnew, nkry
    INTEGER, ALLOCATABLE                     :: ind_k(:), jconv(:)
    REAL(real_8)                             :: b22, fac, zmax
    REAL(real_8), ALLOCATABLE                :: alpha(:,:), b2(:), beta(:,:), &
                                                tmat(:,:), ymat(:,:), &
                                                yomat(:,:), zeta(:)
    REAL(real_8), EXTERNAL                   :: ddot

    ALLOCATE(jconv(cnti%nkry_block),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ind_k(cnti%nkry_block),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(tmat(cnti%nkry_max, cnti%nkry_max),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(beta(cnti%nkry_max, cnti%nkry_block),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(alpha(cnti%nkry_max, cnti%nkry_block),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ymat(cnti%nkry_max, cnti%nkry_max),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(zeta(cnti%nkry_max),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(b2(cnti%nkry_block),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(yomat(cnti%nkry_max, cnti%nkry_block),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DO i=ncurr,ndiag,cnti%nkry_block
       nkbeff=MIN(cnti%nkry_block,ndiag-i+1)
       CALL dcopy(nkbeff*nkpt%ngwk*2,c0(1,i),1,cscr(1,1),1)
       CALL dcopy(nkbeff*nkpt%ngwk*2,cs(1,i),1,cscr(1,nkbeff+1),1)
       DO ib=1,nkbeff
          jconv(ib)=0
          ind_k(ib)=i+ib-1
          beta(1,ib)=0.0_real_8
          beta(2,ib)=beta1(i+ib-1)
          alpha(1,ib)=alpha0(i+ib-1)
       ENDDO
       DO nkry=3,cnti%nkry_max
          nkb0=(nkry-1)*nkbeff+1
          nkb1=(nkry-2)*nkbeff+1
          IF (fint1%ttrot) THEN
             CALL ehpsi_c(cscr(1,nkb1),cscr(1,nkb0),vpot,psi,&
                  nkbeff,ikind,maxfft)
          ELSE
             CALL hpsi(cscr(:,nkb1:nkb1+nkbeff-1),cscr(1,nkb0),sc0,&
                  vpot,psi,nkbeff,ikind,ispin)
          ENDIF
          nhpsi=nhpsi+nkbeff
          DO ib=1,nkbeff
             ibn0=(nkry-1)*nkbeff+ib
             ibn1=(nkry-2)*nkbeff+ib
             ibn2=(nkry-3)*nkbeff+ib
             CALL daxpy(nkpt%ngwk*2,-beta(nkry-1,ib),cscr(1,ibn2),1,&
                  cscr(1,ibn0),1)

             alpha(nkry-1,ib)=REAL(&
                  zdotc(nkpt%ngwk,cscr(1,ibn1),1,cscr(1, ibn0),1))
             CALL mp_sum(alpha(nkry-1,ib),parai%allgrp)
             CALL daxpy(nkpt%ngwk*2,-alpha(nkry-1,ib),cscr(1,ibn1),1,&
                  cscr(1,ibn0),1)
             beta(nkry,ib)=ddot(nkpt%ngwk*2,cscr(1,ibn0),1,&
                  cscr(1,ibn0),1)
             CALL mp_sum(beta(nkry,ib),parai%allgrp)
             beta(nkry,ib)=SQRT(beta(nkry,ib))
             fac=1._real_8/beta(nkry,ib)
             CALL dscal(nkpt%ngwk*2,fac,cscr(1,ibn0),1)
             ! Diagonalize the Lanczos Matrix
             CALL zeroing(tmat)!,cnti%nkry_max*cnti%nkry_max)
             DO j=1,nkry-1
                tmat(j,j)=alpha(j,ib)
                tmat(j+1,j)=beta(j+1,ib)
                tmat(j,j+1)=beta(j+1,ib)
             ENDDO
             CALL jacobi(cnti%nkry_max,nkry-1,tmat,zeta,ymat,ierr)
             ! Calculate beta^2 for the eigenvector
             ! with biggest overlap with initial vector
             zmax=-1.e30_real_8
             DO ii=1,nkry-1
                IF (zeta(ii).GE.zmax) THEN
                   zmax=zeta(ii)
                   iovlm=ii
                   CALL dcopy(nkry-1,ymat(1,iovlm),1,yomat(1,ib),1)
                ENDIF
             ENDDO
             b2(ib)=beta(nkry,ib)*ymat(nkry-1,iovlm)
             b2(ib)=b2(ib)*b2(ib)
             IF (b2(ib).LT.cntr%b2limit) THEN
                jconv(ib)=1
                nconv=nconv+1
                INDEX(ind_k(ib))=nconv
                IF (nkbeff.NE.1) THEN
                   CALL zeroing(c0(:,ind_k(ib)))!,nkpt%ngwk)
                   DO ik=1,nkry-1
                      ikk=(ik-1)*nkbeff+ib
                      CALL daxpy(nkpt%ngwk*2,ymat(ik,iovlm),cscr(1,ikk),1,&
                           c0(1,ind_k(ib)),1)
                   ENDDO
                ELSE
                   CALL dgemv('N',nkpt%ngwk*2,nkry-1,1._real_8,cscr(1,1),nkpt%ngwk*2,&
                        ymat(1,iovlm),1,0._real_8,c0(1,ind_k(1)),1)
                ENDIF
             ENDIF
          ENDDO
          jl=0
          DO ib=1,nkbeff
             IF (jconv(ib).EQ.0) THEN
                jl=jl+1
                IF (jl.NE.ib) THEN
                   CALL dcopy(cnti%nkry_max,beta(1,ib),1,beta(1,jl),1)
                   CALL dcopy(cnti%nkry_max,alpha(1,ib),1,alpha(1,jl),1)
                   CALL dcopy(cnti%nkry_max,yomat(1,ib),1,yomat(1,jl),1)
                   b22=b2(ib)
                   b2(ib)=b2(jl)
                   b2(jl)=b22
                   ii=ind_k(ib)
                   ind_k(ib)=ind_k(jl)
                   ind_k(jl)=ii
                ENDIF
             ENDIF
          ENDDO
          nkbeffnew=jl
          IF (nkbeff.NE.nkbeffnew.AND.nkbeffnew.NE.0) THEN
             DO ik=1,nkry
                jl=0
                DO ib=1,nkbeff
                   IF (jconv(ib).EQ.0) THEN
                      jl=jl+1
                      ikk=(ik-1)*nkbeff+ib
                      jkk=(ik-1)*nkbeffnew+jl
                      CALL dcopy(nkpt%ngwk*2,cscr(1,ikk),1,cscr(1,jkk),1)
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          nkbeff=nkbeffnew
          DO ib=1,nkbeff
             jconv(ib)=0
          ENDDO
          IF (nkbeff.EQ.0) GOTO 100
       ENDDO
       IF (nkbeff.NE.1) THEN
          DO ib=1,nkbeff
             CALL zeroing(c0(:,ind_k(ib)))!,nkpt%ngwk)
             DO ik=1,nkry-2
                ikk=(ik-1)*nkbeff+ib
                CALL daxpy(nkpt%ngwk*2,yomat(ik,ib),cscr(1,ikk),1,&
                     c0(1,ind_k(ib)),1)
             ENDDO
          ENDDO
       ELSE
          CALL dgemv('N',nkpt%ngwk*2,nkry-2,1._real_8,cscr(1,1),nkpt%ngwk*2,&
               yomat(1,1),1,0._real_8,c0(1,ind_k(1)),1)
       ENDIF
100    CONTINUE
       DO ib=1,MIN(cnti%nkry_block,ndiag-i+1)
          IF (b2(ib).GT.b2max) b2max=b2(ib)
          IF (b2(ib).LT.b2min) b2min=b2(ib)
          beta1(ind_k(ib))=SQRT(b2(ib))
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(jconv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ind_k,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(tmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(beta,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(alpha,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ymat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(zeta,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(b2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(yomat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE krylov_ref_c
  ! ==================================================================

END MODULE friesner_c_utils
