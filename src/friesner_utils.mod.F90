MODULE friesner_utils
  USE adjmu_utils,                     ONLY: convfrie
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE frsblk_utils,                    ONLY: reorder
  USE gsortho_utils,                   ONLY: gs_ortho
  USE hfx_drivers,                     ONLY: hfxpsi
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE jacobi_utils,                    ONLY: jacobi
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE randtowf_utils,                  ONLY: randtowf
  USE sort_utils,                      ONLY: sort2
  USE summat_utils,                    ONLY: give_scr_summat,&
                                             summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             fpar,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: friesner
  PUBLIC :: krylov_ref
  PUBLIC :: give_scr_friesner
  PUBLIC :: give_scr_krylov_ref

CONTAINS

  ! ==================================================================
  SUBROUTINE friesner(ndiag,nel,nconv,nhpsi,nstate,&
       c0,cs,sc0,cscr,cgs,focc,&
       vpot,psi,edav,jspin,trefine,tinfo)
    ! ==--------------------------------------------------------------==
    ! ==  Kohn-Sham Matrix Diagonalization by Krylov-Space Method     ==
    ! ==  W.T. Pollard and R.A. Friesner                              ==
    ! ==  J. Chem. Phys. 99, 6742 (1993)                              ==
    ! ==--------------------------------------------------------------==


    INTEGER                                  :: ndiag
    REAL(real_8)                             :: nel
    INTEGER                                  :: nconv, nhpsi, nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), &
                                                cs(ncpw%ngw,*), &
                                                sc0(ncpw%ngw,*), &
                                                cscr(ncpw%ngw,*), &
                                                cgs(ncpw%ngw,*)
    REAL(real_8)                             :: focc(:), vpot(*)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: edav(*)
    INTEGER                                  :: jspin
    LOGICAL                                  :: trefine, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'friesner'
    INTEGER, PARAMETER                       :: ispin = 1 
    REAL(real_8), PARAMETER :: big = 1.e33_real_8 , edavmax = 300._real_8, &
      edavmin = -300._real_8, eps = 1.e-14_real_8

    CHARACTER(len=20)                        :: charspin
    INTEGER                                  :: i, iaux, icycle, ierr, ig, &
                                                isub, itest, j, lwork, &
                                                nconvold, ncurr, ngw2, &
                                                nhpsiold, nleft, ntest
    INTEGER, ALLOCATABLE                     :: INDEX(:)
    LOGICAL                                  :: tcheck, tconv, tdebug, tlsd2, &
                                                ttest
    REAL(real_8)                             :: amu, aux, b2, b2max, b2min, &
                                                fac, sign, tim, tim1, tim2, &
                                                trotmax, trotmin
    REAL(real_8), ALLOCATABLE                :: adav(:,:), alpha0(:), &
                                                beta1(:), work(:)

! ==--------------------------------------------------------------==
! With LSD and 1 electron...

    IF (ndiag.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('  FRIESNER',isub)
    ! ==--------------------------------------------------------------==
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
    sign=-1._real_8
    IF (cntl%tlsd) sign=-2._real_8
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    tconv=.FALSE.
    tdebug=.FALSE.
    ngw2=ncpw%ngw*2
    nconv=0
    nconvold=0
    nhpsiold=0
    IF (fint1%ttrot) THEN
       trotmin=EXP(-edavmax*fint1%betap)
       trotmax=EXP(-edavmin*fint1%betap)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%parent.AND.tinfo) THEN
       IF (jspin.EQ.0) THEN
          charspin='<<<<<<<<<<<<<<<<<<<<'
       ELSEIF (jspin.EQ.1) THEN
          charspin='<<<<<<< ALPHA SPIN<<'
       ELSEIF (jspin.EQ.2) THEN
          charspin='<<<<<<<< BETA SPIN<<'
       ENDIF
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,19("<"),A,A20)')&
            ' LANCZOS DIAGONALIZATION ',charspin
    ENDIF
    ! ==--------------------------------------------------------------==
    ! LWORK query ! TODO verify that this is correct
    lwork = -1
    ALLOCATE(work(1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL dsyev('V','U',ndiag,adav,ndiag,edav,work,lwork,ierr)
    lwork = work(1)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ! Subspace Matrix
    tcheck=.TRUE.
    ntest=0
    DO WHILE(tcheck)
       IF (fint1%ttrot) THEN
          CALL ehpsi(c0,cs,vpot,psi,ndiag)
       ELSE
          CALL hpsi(c0(:,1:ndiag),cs,sc0,vpot,psi,ndiag,1,ispin)
          CALL hfxpsi(cgs(:,1:nstate),c0(:,1:ndiag),cs(:,1:ndiag),focc(1:nstate),sign,psi,nstate,ndiag)
       ENDIF
       nhpsi=ndiag
       CALL ovlap2(ncpw%ngw,ndiag,ndiag,adav,c0,cs,.TRUE.)
       CALL summat(adav,ndiag)
       CALL dsyev('V','U',ndiag,adav,ndiag,edav,work,lwork,ierr)
       ! IF (IO_PARENT.AND.(LWORK < WORK(1))) THEN
       ! WRITE(6,*) ' WARNING| LWORK SUBOPTIMAL FOR DSYEV'
       ! WRITE(6,*) '    LWORK:',LWORK
       ! WRITE(6,*) ' OPTIMAL:',INT(WORK(1))
       ! ENDIF
       CALL reorder(ndiag,ndiag,edav,adav)
       ! Rotate Orbitals
       CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,c0(1,1),ngw2,&
            adav(1,1),ndiag,0._real_8,cscr(1,1),ngw2)
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,c0(1,1),1)
       ! Check if an eigenvalue is not equal to 0 (dim < NDIAG)
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
                  WRITE(6,*) 'FRIESNER| EIGENVECTOR',i,' IS VERY BAD! '
             ! Not useful: Now we use H|C0>.
             ! CALL RANDTOWF(C0(1,I),1,0,0)
             itest=itest+1
             edav(i)=1.e33_real_8
          ENDIF
       ENDDO
       IF (tcheck) THEN
          ntest=ntest+1
          IF (ntest.GT.10.AND.paral%parent)&
               CALL stopgm(' FRIESNER', 'CAN''T FIND GUEST VECTORS',& 
               __LINE__,__FILE__)
          CALL sort2(edav,ndiag,index)
          ! Put at the end randomized wavefunctions.
          CALL dcopy(ngw2*ndiag,c0,1,cs,1)
          DO i=1,ndiag
             j=INDEX(i)
             IF (edav(i).EQ.big) THEN
                CALL dcopy(ngw2,cs(1,j),1,c0(1,i),1)
             ELSE
                CALL dcopy(ngw2,cscr(1,j),1,c0(1,i),1)
             ENDIF
          ENDDO
          ! Orthogonalize.
          CALL gs_ortho(c0,ndiag-itest,c0(1,ndiag-itest+1),itest)
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    ! Rotate also H|C0> (IF TCHECK=.FALSE., ADAV unchanged)
    CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,cs(1,1),ngw2,&
         adav(1,1),ndiag,0._real_8,cscr(1,1),ngw2)
    CALL dcopy(ngw2*ndiag,cscr(1,1),1,cs(1,1),1)
    tim2=m_walltime()
    tim=(tim2-tim1)*0.001_real_8
    IF (paral%parent.AND.tinfo) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T58,F8.2)')&
            ' >> TIME FOR INITIAL SUBSPACE DIAGONALIZATION:  ',tim
       IF (paral%io_parent)&
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
          DO ig=1,ncpw%ngw
             cscr(ig,1)=cs(ig,i)-edav(i)*c0(ig,i)
          ENDDO
          b2=dotp(ncpw%ngw,cscr(:,1),cscr(:,1))
          CALL mp_sum(b2,parai%allgrp)
          beta1(i)=SQRT(b2)
          alpha0(i)=edav(i)
          IF (b2.LT.cntr%b2limit) THEN
             nconv=nconv+1
             INDEX(i)=nconv
             IF (b2.GT.b2max) b2max=b2
             IF (b2.LT.b2min) b2min=b2
          ELSE
             CALL dcopy(ngw2,cscr(1,1),1,cs(1,i),1)
             fac=1._real_8/beta1(i)
             CALL dscal(ngw2,fac,cs(1,i),1)
          ENDIF
       ENDDO
       ! Order states: converged first
       DO i=ncurr,ndiag
          j=INDEX(i)
          IF (j.NE.0.AND.j.NE.i) THEN
             CALL dswap(ngw2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngw2,cs(1,i),1,cs(1,j),1)
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
       CALL krylov_ref(ncurr,ndiag,nconv,nhpsi,nstate,&
            c0,cs,sc0,cscr,cgs,focc,sign,&
            vpot,psi,index,alpha0,beta1,b2min,b2max)
       ! Order states : converged first
       DO i=ncurr,ndiag
          j=INDEX(i)
          IF (j.NE.0.AND.j.NE.i) THEN
             CALL dswap(ngw2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngw2,cs(1,i),1,cs(1,j),1)
             INDEX(j)=j
             INDEX(i)=0
             aux = beta1(i)
             beta1(i) = beta1(j)
             beta1(j) = aux
          ENDIF
       ENDDO
       ! Reorthogonalize
       IF (ncurr.LE.ndiag) THEN
          CALL gs_ortho(c0,ncurr-1,c0(1,ncurr),ndiag-ncurr+1)
       ENDIF
       ! Calculate new forces; only for states that entered Lanczos
       nleft=ndiag-ncurr+1
       IF (ncurr.LE.ndiag) THEN
          IF (fint1%ttrot) THEN
             CALL ehpsi(c0(1,ncurr),cs(1,ncurr),vpot,psi,nleft)
          ELSE
             CALL hpsi(c0(:,ncurr:ncurr+nleft),cs(1,ncurr),sc0(1,ncurr),&
                  vpot,psi,nleft,1,ispin)
             CALL hfxpsi(cgs(:,1:nstate),c0(:,ncurr:ncurr+nleft-1),cs(:,ncurr:ncurr+nleft-1),focc(1:nstate),sign,&
                  psi,nstate,nleft)
          ENDIF
          nhpsi=nhpsi+nleft
       ENDIF
       ! Diagonalize Kohn-Sham Matrix
       CALL ovlap2(ncpw%ngw,ndiag,ndiag,adav,c0,cs,.TRUE.)
       CALL summat(adav,ndiag)
       ! This function does not preserve the order of eigenfunctions.
       CALL dsyev('V','U',ndiag,adav,ndiag,edav,work,lwork,ierr)
       ! We inverse the order
       ! (not the same one for DSYEV than for this routine).
       CALL reorder(ndiag,ndiag,edav,adav)
       CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,c0(1,1),ngw2,&
            adav(1,1),ndiag,0._real_8,cscr(1,1),ngw2)
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,c0(1,1),1)
       CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,cs(1,1),ngw2,&
            adav(1,1),ndiag,0._real_8,cscr(1,1),ngw2)
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,cs(1,1),1)
       ! If one eigenvalue is equal to zero (like friesner_c).
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
                  WRITE(6,*) ' FRIESNER| ','VERY BAD EIGENVALUE ! !'
             IF (cnti%ntrans.GT.ndiag) THEN
                CALL randtowf(c0(:,i:i),1,0,0)
             ELSE
                CALL dswap(ngw2,c0(1,i),1,c0(1,cnti%ntrans),1)
                CALL dswap(ngw2,cs(1,i),1,cs(1,cnti%ntrans),1)
             ENDIF
          ENDIF
       ENDDO
       tim2=m_walltime()
       tim=(tim2-tim1)*0.001_real_8
       IF (paral%parent.AND.tinfo) THEN
          IF (paral%io_parent)&
               WRITE(6,'(4X,I5,5X,I5,2(1PE13.3),0PF10.2,0PF10.2)')&
               icycle,nconv,b2max,b2min,&
               REAL(nhpsi-nhpsiold,kind=real_8)/REAL(ndiag,kind=real_8),tim
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
             CALL convfrie(ndiag,nconv,nel,edav,cntl%tlsd,&
                  amu,tconv)
          ENDIF
       ENDIF
       IF (tconv.OR.nconv.EQ.ndiag) THEN
          GOTO 200
       ENDIF
       ! End of Krylov Space Cycles
    ENDDO
200 CONTINUE
    IF (paral%parent.AND.((.NOT.tconv).AND.nconv.NE.ndiag).AND.tdebug) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("! "))')
       IF (paral%io_parent)&
            WRITE(6,'(" ! !",A,T64,"!!")')&
            ' FRIESNER| NOT ALL ROOTS ARE CONVERGED'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("! "))')
    ENDIF
    ! Order Roots with respect to energy
    CALL sort2(edav,ndiag,index)
    CALL dcopy(ngw2*ndiag,c0(1,1),1,cs(1,1),1)
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
       CALL dcopy(2*ncpw%ngw,cs(1,j),1,c0(1,i),1)
    ENDDO
    cntl%tlsd=tlsd2
    ! ==--------------------------------------------------------------==
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
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
    ! ==--------------------------------------------------------------==
    CALL tihalt('  FRIESNER',isub)
    RETURN
  END SUBROUTINE friesner
  ! ==================================================================
  SUBROUTINE krylov_ref(ncurr,ndiag,nconv,nhpsi,nstate,&
       c0,cs,sc0,cscr,cgs,focc,sign,&
       vpot,psi,index,alpha0,beta1,b2min,b2max)
    ! ==--------------------------------------------------------------==
    ! ==  LANCZOS REFINEMENT LOOPS                                    ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ncurr, ndiag, nconv, nhpsi, &
                                                nstate
    COMPLEX(real_8) :: c0(ncpw%ngw,ndiag), cs(ncpw%ngw,ndiag), &
      sc0(ncpw%ngw,*), cscr(ncpw%ngw,cnti%nkry_block*cnti%nkry_max), &
      cgs(ncpw%ngw,*)
    REAL(real_8)                             :: focc(*), sign, vpot(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    INTEGER                                  :: INDEX(ndiag)
    REAL(real_8)                             :: alpha0(ndiag), beta1(ndiag), &
                                                b2min, b2max

    CHARACTER(*), PARAMETER                  :: procedureN = 'krylov_ref'
    INTEGER, PARAMETER                       :: ispin = 1 

    INTEGER                                  :: i, ib, ibn0, ibn1, ibn2, &
                                                ierr, ii, ik, ikk, iovlm, j, &
                                                jkk, jl, nkb0, nkb1, nkbeff, &
                                                nkbeffnew, nkry
    INTEGER, ALLOCATABLE                     :: ind_k(:), jconv(:)
    REAL(real_8)                             :: b22, fac, zmax
    REAL(real_8), ALLOCATABLE                :: alpha(:,:), b2(:), beta(:,:), &
                                                tmat(:,:), ymat(:,:), &
                                                yomat(:,:), zeta(:)

! ==--------------------------------------------------------------==
! TODO check stat
! TODO align for BG

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
       CALL dcopy(nkbeff*ncpw%ngw*2,c0(1,i),1,cscr(1,1),1)
       CALL dcopy(nkbeff*ncpw%ngw*2,cs(1,i),1,cscr(1,nkbeff+1),1)
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
             CALL ehpsi(cscr(1,nkb1),cscr(1,nkb0),vpot,psi,nkbeff)
          ELSE
             CALL hpsi(cscr(:,nkb1:nkb1+nkbeff-1),cscr(1,nkb0),sc0,&
                  vpot,psi,nkbeff,1,ispin)
             CALL hfxpsi(cgs(:,1:nstate),cscr(:,nkb1:nkb1+nkbeff-1),cscr(:,nkb0:nkb0+nkbeff-1),&
                  focc(1:nstate),sign,psi,nstate,nkbeff)
          ENDIF
          nhpsi=nhpsi+nkbeff
          DO ib=1,nkbeff
             ibn0=(nkry-1)*nkbeff+ib
             ibn1=(nkry-2)*nkbeff+ib
             ibn2=(nkry-3)*nkbeff+ib
             CALL daxpy(ncpw%ngw*2,-beta(nkry-1,ib),cscr(1,ibn2),1,&
                  cscr(1,ibn0),1)
             alpha(nkry-1,ib)=dotp(ncpw%ngw,cscr(:,ibn1),cscr(:,ibn0))
             CALL mp_sum(alpha(nkry-1,ib),parai%allgrp)
             CALL daxpy(ncpw%ngw*2,-alpha(nkry-1,ib),cscr(1,ibn1),1,&
                  cscr(1,ibn0),1)
             beta(nkry,ib)=dotp(ncpw%ngw,cscr(:,ibn0),cscr(:,ibn0))
             CALL mp_sum(beta(nkry,ib),parai%allgrp)
             beta(nkry,ib)=SQRT(beta(nkry,ib))
             fac=1._real_8/beta(nkry,ib)
             CALL dscal(ncpw%ngw*2,fac,cscr(1,ibn0),1)
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
                   CALL zeroing(c0(:,ind_k(ib)))!,ngw)
                   DO ik=1,nkry-1
                      ikk=(ik-1)*nkbeff+ib
                      CALL daxpy(ncpw%ngw*2,ymat(ik,iovlm),cscr(1,ikk),1,&
                           c0(1,ind_k(ib)),1)
                   ENDDO
                ELSE
                   CALL dgemv('N',ncpw%ngw*2,nkry-1,1._real_8,cscr(1,1),ncpw%ngw*2,&
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
                      CALL dcopy(ncpw%ngw*2,cscr(1,ikk),1,cscr(1,jkk),1)
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
             CALL zeroing(c0(:,ind_k(ib)))!,ngw)
             DO ik=1,nkry-2
                ikk=(ik-1)*nkbeff+ib
                CALL daxpy(ncpw%ngw*2,yomat(ik,ib),cscr(1,ikk),1,&
                     c0(1,ind_k(ib)),1)
             ENDDO
          ENDDO
       ELSE
          CALL dgemv('N',ncpw%ngw*2,nkry-2,1._real_8,cscr(1,1),ncpw%ngw*2,&
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
  END SUBROUTINE krylov_ref
  ! ==================================================================
  SUBROUTINE give_scr_friesner(lfriesner,tag,ndiag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lfriesner
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ndiag

    INTEGER                                  :: ldiag, lhpsi, lkrylov_ref, &
                                                lscr, lsummat

    lscr=0
    ldiag=0
    lkrylov_ref=0
    lsummat=0
    lhpsi=0
    IF (.NOT.fint1%ttrot) THEN
       CALL give_scr_hpsi(lhpsi,tag,ndiag)
    ENDIF
    CALL give_scr_summat(lsummat,tag,ndiag)
    CALL give_scr_krylov_ref(lkrylov_ref,tag,ndiag)
    IF (tkpts%tkpnt) THEN
       ldiag=((3*ndiag-2) + 5*(2*ndiag))! ZHEEV
       lscr=ndiag/2+1+2*ndiag*ndiag+2*ndiag! INDEX ADAV ALPHA0 BETA1
    ELSE
       ldiag=5*ndiag         ! DSYEV
       lscr=ndiag/2+1+ndiag*ndiag+2*ndiag ! INDEX ADAV ALPHA0 BETA1
    ENDIF
    ldiag=MAX(ldiag,2*ndiag)  ! CONVFRIE
    lfriesner=lscr+MAX(lhpsi,lsummat,lkrylov_ref,ldiag)+10
    tag='LSCR+MAX(HPSI,SUMHMAT,...)' 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_friesner
  ! ==================================================================
  SUBROUTINE give_scr_krylov_ref(lkrylov_ref,tag,ndiag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lkrylov_ref
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ndiag

    INTEGER                                  :: lhpsi

    lkrylov_ref=&
         (cnti%nkry_block/2+1)+(cnti%nkry_block/2+1)+ & ! JCONV IND_K
         cnti%nkry_max*cnti%nkry_max+cnti%nkry_max*cnti%nkry_block+   & ! TMAT BETA
         cnti%nkry_max*cnti%nkry_block+cnti%nkry_max*cnti%nkry_max+   & ! ALPHA YMAT
         cnti%nkry_max+cnti%nkry_block+cnti%nkry_max*cnti%nkry_block+10 ! ZETA B2 YOMAT
    IF (fint1%ttrot) THEN
       lhpsi=0
    ELSE
       CALL give_scr_hpsi(lhpsi,tag,ndiag)
    ENDIF
    lkrylov_ref=lkrylov_ref+lhpsi
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_krylov_ref
  ! ==================================================================

END MODULE friesner_utils
