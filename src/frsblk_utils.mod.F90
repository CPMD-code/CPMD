MODULE frsblk_utils
  USE adjmu_utils,                     ONLY: convfrie
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE geq0mod,                         ONLY: geq0
  USE gsortho_utils,                   ONLY: gs_ortho
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu_vec_cmplx
  USE rgs_utils,                       ONLY: rgs
  USE summat_utils,                    ONLY: give_scr_summat,&
                                             summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: frsblk
  PUBLIC :: give_scr_frsblk
  PUBLIC :: puttmat
  PUBLIC :: reorder
  PUBLIC :: rsdblk
  PUBLIC :: fmdiag
  PUBLIC :: prjcnv

CONTAINS

  SUBROUTINE frsblk(ndiag,nel,nconv,nhpsi,c0,cs,sc0,cscr,vpot,psi,&
       edav,jspin,trefine,tinfo)
    ! ==--------------------------------------------------------------==
    ! ==  Kohn-Sham Matrix Diagonalization by Block Krylov-Space      == 
    ! ==  Method                                                      ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ndiag
    REAL(real_8)                             :: nel
    INTEGER                                  :: nconv, nhpsi
    COMPLEX(real_8)                          :: c0(ncpw%ngw,ndiag), &
                                                cs(ncpw%ngw,*), &
                                                sc0(ncpw%ngw,*), &
                                                cscr(ncpw%ngw,*)
    REAL(real_8)                             :: vpot(*)
    COMPLEX(real_8)                          :: psi(:)
    REAL(real_8)                             :: edav(*)
    INTEGER                                  :: jspin
    LOGICAL                                  :: trefine, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'frsblk'
    INTEGER, PARAMETER                       :: ispin = 1
    REAL(real_8), PARAMETER                  :: edavmax = 300._real_8, &
                                                edavmin = -300._real_8, &
                                                eps = 1.e-14_real_8

    CHARACTER(len=20)                        :: charspin
    INTEGER                                  :: i, icycle, ierr, ig, isub, &
                                                itest, j, lwork, nbleff, &
                                                nblk, nconvold, ncurr, ngw2, &
                                                nhpsiold, nk, nleft, ntest
    INTEGER, ALLOCATABLE                     :: INDEX(:)
    LOGICAL                                  :: tcheck, tconv, tlsd2, ttest
    REAL(real_8)                             :: amu, aux, b2, b2max, b2min, &
                                                tim, tim1, tim2, trotmax, &
                                                trotmin
    REAL(real_8), ALLOCATABLE                :: adav(:,:), alpha(:), beta(:), &
                                                tmat(:), work(:), wt(:)

! Variables
! ==--------------------------------------------------------------==
! With LSD and 1 electron...

    IF (ndiag.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('    FRSBLK',isub)
    nk=cnti%nkry_block*cnti%nkry_max
    ! TODO check stat
    ! TODO align for BG
    ALLOCATE(INDEX(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(tmat(3 * nk * nk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(alpha(cnti%nkry_block * cnti%nkry_block * (cnti%nkry_max + 1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(beta(cnti%nkry_block * cnti%nkry_block * cnti%nkry_max),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(wt(nk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(adav(ndiag, ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    lwork=MAX(1,3*ndiag-1)
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    tim1=m_walltime()
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    tconv=.FALSE.
    ngw2=ncpw%ngw*2
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
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,13("<"),A,A20)')&
            ' BLOCK LANCZOS DIAGONALIZATION ',charspin
    ENDIF
    ! Subspace Matrix
    tcheck=.TRUE.
    ntest=0
    DO WHILE(tcheck)
       IF (fint1%ttrot) THEN
          CALL ehpsi(c0,cs,vpot,psi,ndiag)
       ELSE
          CALL stopgm(procedureN,'fix me: wrong dummy args + provide a regtest', &
               __LINE__,__FILE__)
          !fix me call hpsi(c0,cs,sc0,vpot,psi,ndiag,1,ispin)
       ENDIF
       nhpsi=ndiag
       CALL ovlap2(ncpw%ngw,ndiag,ndiag,adav,c0,cs,.TRUE.)
       CALL dsyev('V','U',ndiag,adav,ndiag,edav,work,lwork,ierr)
       CALL reorder(ndiag,ndiag,edav,adav)
       ! Rotate Orbitals
       CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,c0(1,1),ngw2,&
            adav(1,1),ndiag,0._real_8,cscr(1,1),ngw2)
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,c0(1,1),1)
       CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,cs(1,1),ngw2,&
            adav(1,1),ndiag,0._real_8,cscr(1,1),ngw2)
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,cs(1,1),1)
       ! Check if an eigenvalue is not equal to 0 (dim < NDIAG)
       tcheck=.FALSE.
       DO i=1,ndiag
          IF (fint1%ttrot) THEN
             ttest = (edav(i).GT.trotmax).OR.(edav(i).LT.trotmin)
          ELSE
             ttest = ABS(edav(i)).LT.eps
          ENDIF
          IF (ttest) THEN
             tcheck=.TRUE.
             IF (paral%io_parent)&
                  WRITE(6,*) 'FRSBLK| EIGENVECTOR',i,' IS VERY BAD! '
             CALL repprngu_vec_cmplx(ncpw%ngw,c0(:,i))
             IF (geq0) c0(1,i)=CMPLX(REAL(c0(1,i)),0._real_8,kind=real_8)
             itest=i
          ENDIF
       ENDDO
       IF (tcheck) THEN
          ntest=ntest+1
          IF (ntest.GT.10.AND.paral%parent)&
               CALL stopgm(' FRSBLK', 'CAN''T FIND GUEST VECTORS',& 
               __LINE__,__FILE__)
          IF (itest.EQ.1.AND.ndiag.NE.1)&
               CALL dswap(ngw2,c0(1,1),1,c0(1, ndiag),1)
          CALL rgs(c0,ndiag,adav)
       ENDIF
    ENDDO
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
          IF (b2.LT.cntr%b2limit) THEN
             nconv=nconv+1
             INDEX(i)=nconv
          ENDIF
          IF (b2.GT.b2max) b2max=b2
          IF (b2.LT.b2min) b2min=b2
       ENDDO
       ! Order states: converged first
       DO i=ncurr,ndiag
          j=INDEX(i)
          IF (j.NE.0.AND.j.NE.i) THEN
             CALL dswap(ngw2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngw2,cs(1,i),1,cs(1,j),1)
             INDEX(j)=j
             INDEX(i)=0
             aux = edav(i)
             edav(i)=edav(j)
             edav(j)=aux
          ENDIF
       ENDDO
       ncurr=nconv+1
       nleft=ndiag-nconv
       ! Lanczos Refinement
       IF (icycle.EQ.1.AND.ncurr.GT.ndiag) THEN
          trefine=.FALSE.
       ELSE
          trefine=.TRUE.
       ENDIF
       ! ..New
       DO i=ncurr,ndiag,cnti%nkry_block
          nbleff=MIN(cnti%nkry_block,ndiag-i+1)
          ! ==----------------------------------------------------------==
          ! Preparation for refinement
          ! ==----------------------------------------------------------==
          ! NOTE: The type of the actual argument differs from the type of
          ! the dummy argument.  [C0, CS]
          CALL prpkrv(ngw2,nbleff,cnti%nkry_max,nconv,c0,c0(:,i),&
               cs(1,i),cscr,alpha,beta,edav(i),&
               nhpsi,vpot,psi,sc0,fint1%ttrot)
          ! ==----------------------------------------------------------==
          ! Refinement Loop
          ! ==----------------------------------------------------------==
          nblk=nbleff*cnti%nkry_max
          ! NOTE: The type of the actual argument differs from the type of
          ! the dummy argument.  [C0, CSCR]
          CALL kryref(ngw2,nbleff,cnti%nkry_max,nconv,c0,edav,cscr,alpha,&
               beta,tmat,nblk*nblk,wt,nblk,nhpsi,&
               vpot, psi,sc0,fint1%ttrot)
          ! ==----------------------------------------------------------==
          ! ..V=[V_1 V_2.. V_L] Y
          CALL dgemm('N','N',ngw2,nbleff,nblk,1._real_8,cscr,ngw2,tmat,&
               nblk,0._real_8,c0(1,i),ngw2)
          ! ==----------------------------------------------------------==
          ! == End of refinement over states                            ==
          ! ==----------------------------------------------------------==
       ENDDO
       ! ..New
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
             CALL stopgm(procedureN,'fix me: wrong dummy args + provide a regtest', &
                  __LINE__,__FILE__)
             !fix me call hpsi(c0(1,ncurr),cs(1,ncurr),sc0(1,ncurr),&
             !     vpot,psi,nleft,1,ispin)
          ENDIF
          nhpsi=nhpsi+nleft
       ENDIF
       ! Diagonalize Kohn-Sham Matrix
       CALL ovlap2(ncpw%ngw,ndiag,ndiag,adav,c0,cs,.TRUE.)
       CALL dsyev('V','U',ndiag,adav,ndiag,edav,work,lwork,ierr)
       CALL reorder(ndiag,ndiag,edav,adav)
       CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,c0(1,1),ngw2,adav(1,1),&
            ndiag,0._real_8,cscr(1,1),ngw2)
       CALL dcopy(ngw2*ndiag,cscr(1,1),1,c0(1,1),1)
       CALL dgemm('N','N',ngw2,ndiag,ndiag,1._real_8,cs(1,1),ngw2,adav(1,1),&
            ndiag,0._real_8,cscr(1,1),ngw2)
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
                  WRITE(6,*) ' FRSBLK| ','VERY BAD EIGENVALUE !!'
             IF (cnti%ntrans.GT.ndiag) THEN
                CALL repprngu_vec_cmplx(ncpw%ngw,c0(:,i))
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
    IF (paral%parent.AND.((.NOT.tconv).AND.nconv.NE.ndiag)) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
       IF (paral%io_parent)&
            WRITE(6,'(" !!",A,T64,"!!")')&
            ' FRSBLK| NOT ALL ROOTS ARE CONVERGED'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("!"))')
    ENDIF
    cntl%tlsd=tlsd2
    CALL tihalt('    FRSBLK',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    DEALLOCATE(tmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(alpha,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(beta,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(wt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(adav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE frsblk
  ! ==================================================================
  SUBROUTINE give_scr_frsblk(lfrsblk,tag,ndiag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lfrsblk
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ndiag

    INTEGER                                  :: lhpsi, lscr, lscrdiag, &
                                                lsummat, nk

    lhpsi=0
    lsummat=0
    lscr=0
    lscrdiag=0
    CALL give_scr_hpsi(lhpsi,tag,ndiag)
    CALL give_scr_summat(lsummat,tag,ndiag)
    nk=cnti%nkry_block*cnti%nkry_max
    IF (tkpts%tkpnt) THEN
       lscrdiag=(3*ndiag-2) + 2*(2*ndiag)
       lscr=(6*nk/2+1)+(ndiag/2+1)+         & ! ISCR INDEX
            6*nk*nk+nk+                             & ! TMAT WT
            2*(cnti%nkry_block*cnti%nkry_block*(cnti%nkry_max+1))+ & ! ALPHA
            2*(cnti%nkry_block*cnti%nkry_block*cnti%nkry_max)+     & ! BETA
            2*(ndiag*ndiag)                         ! ADAV
    ELSE
       lscrdiag=8*ndiag*cnti%nkry_block
       lscr=(6*nk/2+1)+(ndiag/2+1)+         & ! ISCR INDEX
            3*nk*nk+nk+                             & ! TMAT WT
            (cnti%nkry_block*cnti%nkry_block*(cnti%nkry_max+1))+   & ! ALPHA
            (cnti%nkry_block*cnti%nkry_block*cnti%nkry_max)+       & ! BETA
            (ndiag*ndiag)                           ! ADAV
    ENDIF
    lfrsblk=lscr+MAX(lhpsi,lsummat,lscrdiag)+10
    tag='LSCR+MAX(HPSI,SUMMAT,DIAG)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_frsblk
  ! ==================================================================
  SUBROUTINE puttmat(char,t,m,n,a,ma,na,ibeg,jbeg)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: char
    INTEGER                                  :: m, n
    REAL(real_8)                             :: t(m,n)
    INTEGER                                  :: ma, na
    REAL(real_8)                             :: a(ma,na)
    INTEGER                                  :: ibeg, jbeg

    INTEGER                                  :: i, j

    IF (ma*ibeg.GT.m) CALL stopgm('PUTTMAT',' MA*IBEG.GT.M',& 
         __LINE__,__FILE__)
    IF (na*jbeg.GT.n) CALL stopgm('PUTTMAT',' NA*JBEG.GT.N',& 
         __LINE__,__FILE__)
    IF (char.EQ.'N') THEN
       DO j=1,na
          DO i=1,ma
             t(i+ma*(ibeg-1),j+na*(jbeg-1))=a(i,j)
          ENDDO
       ENDDO
    ELSEIF (char.EQ.'T') THEN
       DO j=1,na
          DO i=1,ma
             t(i+ma*(ibeg-1),j+na*(jbeg-1))=a(j,i)
          ENDDO
       ENDDO
    ELSE
       CALL stopgm('PUTTMAT','ILLEGAL CHAR',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE puttmat
  ! ================================================================
  SUBROUTINE getmat(char,t,m,n,a,ma,na,ibeg,jbeg)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: char
    INTEGER                                  :: m, n
    REAL(real_8)                             :: t(m,n)
    INTEGER                                  :: ma, na
    REAL(real_8)                             :: a(ma,na)
    INTEGER                                  :: ibeg, jbeg

    INTEGER                                  :: i, j

    IF (ma*ibeg.GT.m) CALL stopgm('GETMAT',' MA*IBEG.GT.M',& 
         __LINE__,__FILE__)
    IF (na*jbeg.GT.n) CALL stopgm('GETMAT',' NA*JBEG.GT.N',& 
         __LINE__,__FILE__)
    IF (char.EQ.'N') THEN
       DO j=1,na
          DO i=1,ma
             a(i,j)=t(i+ma*(ibeg-1),j+na*(jbeg-1))
          ENDDO
       ENDDO
    ELSEIF (char.EQ.'T') THEN
       DO j=1,na
          DO i=1,ma
             a(j,i)=t(i+ma*(ibeg-1),j+na*(jbeg-1))
          ENDDO
       ENDDO
    ELSE
       CALL stopgm('GETMAT','ILLEGAL CHAR',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getmat
  ! =================================================================
  SUBROUTINE reorder(m,n,w,a)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n
    REAL(real_8)                             :: w(n), a(m,n)

    INTEGER                                  :: j
    REAL(real_8)                             :: aux

    DO j=1,n/2
       CALL dswap(m,a(1,j),1,a(1,n-j+1),1)
       aux=w(n-j+1)
       w(n-j+1)=w(j)
       w(j)=aux
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE reorder
  ! ==================================================================
  SUBROUTINE puttab(n,kd,a)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: n, kd
    REAL(real_8)                             :: a(n,n)

    CHARACTER(*), PARAMETER                  :: procedureN = 'puttab'

    INTEGER                                  :: i, ibeg, ierr, j
    REAL(real_8), ALLOCATABLE                :: aux(:)

    ALLOCATE(aux(n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(aux)!,n)
    DO j=1,n
       CALL dcopy(n,a(1,j),1,aux,1)
       CALL zeroing(a(:,j))!,n)
       ibeg=MAX(j-kd,1)
       DO i=ibeg,j
          a(1+kd+i-j,j)=aux(i)
       ENDDO
    ENDDO
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE puttab
  ! ==================================================================
  SUBROUTINE rsdblk(char,m,n,v1,a,v2)
    ! ==--------------------------------------------------------------==
    ! == V_I = V_I - V_(I-1) A_(I-1)                                  ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: char
    INTEGER                                  :: m, n
    REAL(real_8)                             :: v1(m,n), a(n,n), v2(m,n)

    INTEGER                                  :: isub

    CALL tiset('    RSDBLK',isub)
    IF (char.EQ.'N') THEN
       CALL dgemm('N','N',m,n,n,-1._real_8,v1,m,a,n,1._real_8,v2,m)
    ELSEIF (char.EQ.'T') THEN
       CALL dgemm('N','T',m,n,n,-1._real_8,v1,m,a,n,1._real_8,v2,m)
    ELSE
       CALL stopgm('RSDBLK','ILLEGAL CHAR',& 
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt('    RSDBLK',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rsdblk
  ! ==================================================================
  SUBROUTINE fmdiag(m,n,t,wt)
    ! ==--------------------------------------------------------------==
    ! Returns the n largest eigenvalues of MxM symtric matrix T       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n
    REAL(real_8)                             :: t(m,m,*), wt(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fmdiag'

    INTEGER                                  :: ierr, info, isub, lwork, meval
    INTEGER, ALLOCATABLE                     :: ifail(:), iwork(:)
    REAL(real_8), ALLOCATABLE                :: work(:)

    CALL tiset('    FMDIAG',isub)
    IF (m.EQ.n) THEN
       lwork=MAX(1,3*n-1)! TODO optimal LWORK
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL dsyev('V','U',n,t,n,wt,work,lwork,info)
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSE
       lwork=8*n! TODO optimal LWORK
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(iwork(5*n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(ifail(n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL dsyevx('V','I','U',m,t,m,1._real_8,1._real_8,m-n+1,m,0._real_8,meval,wt,&
            t(1,1,2),m,work,lwork,iwork,ifail,info)
       DEALLOCATE(ifail,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(iwork,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       IF (meval.LT.n) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' WARNING| DSYEVX RETURNED MEVAL < N',meval,n
       ENDIF
       CALL dcopy(m*n,t(1,1,2),1,t(1,1,1),1)
    ENDIF
    IF (info.NE.0) CALL stopgm('FMDIAG','FMDIAG FAILED',& 
         __LINE__,__FILE__)
    CALL reorder(m,n,wt,t)
    CALL tihalt('    FMDIAG',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fmdiag
  ! ==================================================================
  SUBROUTINE bandm(in,n,t,wt)
    ! ==--------------------------------------------------------------==
    ! == Returns the n largest eigenvalues                            ==
    ! == of MxM symetric banded matrix T                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: in, n
    REAL(real_8)                             :: t(in,in,3), wt(in)

    CHARACTER(*), PARAMETER                  :: procedureN = 'bandm'

    INTEGER                                  :: ierr, info, isub, meval
    INTEGER, ALLOCATABLE                     :: ifail(:), iwork(:)
    REAL(real_8), ALLOCATABLE                :: work(:)

    CALL tiset('     BANDM',isub)
    ! ..Put T is form suitable for banded matrix diagonalisation 
    CALL puttab(in,n,t)
    ALLOCATE(work(7*in),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(iwork(5*in),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ifail(n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! CALL DSBEV('V','U',IN,N,T,IN,WT,T(1,2),IN,SCR,INFO)
    CALL dsbevx('V','I','U',in,n,t,in,t(1,1,2),in,1._real_8,1._real_8,&
         in-n+1,in,1.e-10_real_8,meval,wt,t(1,1,3),in,work,iwork,ifail,info)
    DEALLOCATE(ifail,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(iwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (meval.LT.n) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' WARNING| DSBEVX RETURNED MEVAL < N',meval,n
    ENDIF
    CALL dcopy(in*in,t(1,1,3),1,t(1,1,1),1)
    CALL reorder(in,n,wt,t)
    CALL tihalt('     BANDM',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE bandm
  ! ==================================================================
  SUBROUTINE prjcnv(m,np,n0,v0,w,v)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, np, n0
    REAL(real_8)                             :: v0(m,n0), w(n0), v(m,np)

    CHARACTER(*), PARAMETER                  :: procedureN = 'prjcnv'

    INTEGER                                  :: ierr, isub
    REAL(real_8), ALLOCATABLE                :: smat(:,:)

    IF (n0.EQ.0) THEN
       RETURN
    ELSE
       RETURN
    ENDIF
    CALL tiset('    PRJCNV',isub)
    ALLOCATE(smat(n0,np),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ..SMAT=V0^T.V 
    CALL ovlap2(m/2,n0,np,smat,v0,v,.TRUE.)
    CALL mp_sum(smat,np*n0,parai%allgrp)
    ! deb      DO J=1,N0
    ! deb        CALL DSCAL(NP,W(J),SMAT(J,1),N0)
    ! deb      ENDDO
    CALL dgemm('N','N',m,np,n0,-1._real_8,v0,m,smat,n0,1._real_8,v,m)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('    PRJCNV',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prjcnv
  ! ==================================================================
  SUBROUTINE write_matrix(char,m,n,a)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=10)                        :: char
    INTEGER                                  :: m, n
    REAL(real_8)                             :: a(m,n)

    INTEGER                                  :: i, j

    IF (paral%io_parent)&
         WRITE(6,*) char
    DO i=1,m
       IF (paral%io_parent)&
            WRITE(6,'(12(1PE15.6))') (a(i,j),j=1,n)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE write_matrix
  ! ==================================================================


END MODULE frsblk_utils


SUBROUTINE prpkrv(m,n,nkry,nconv,vconv,v0,vs,v,am,bm,w,&
     nhpsi,vpot,psi,sc0,ttrot)
  ! ==---------------------------------------------------------------==
  ! == Returns A_1,A_2,B_2 and V in a form suitable for Krylov       ==
  ! == refinement                                                    ==
  ! ==---------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE frsblk_utils, ONLY : prjcnv, rsdblk, puttmat
  USE hpsi_utils, ONLY : hpsi
  USE fft_maxfft,                      ONLY: maxfftn
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: m, n, nkry, nconv
  REAL(real_8)                               :: vconv(m,*), v0(m,n), vs(m,n), &
                                                v(m,n,nkry+1), &
                                                am(n,n,nkry+1), bm(n,n,nkry), &
                                                w(n)
  INTEGER                                    :: nhpsi
  REAL(real_8)                               :: vpot(*)
  COMPLEX(real_8)                            :: psi(maxfftn), sc0(*)
  LOGICAL                                    :: ttrot

  CHARACTER(*), PARAMETER                    :: procedureN = 'prpkrv'
  INTEGER, PARAMETER                         :: ispin = 1

  INTEGER                                    :: i, isub

  CALL tiset('    PRPKRV',isub)
  ! ==---------------------------------------------------------------==
  CALL dcopy(m*n,v0,1,v(1,1,1),1)
  CALL dcopy(m*n,vs,1,v(1,1,2),1)
  ! ..Compute residual 
  DO i=1,n
     CALL daxpy(m,-w(i),v(1,i,1),1,v(1,i,2),1)
  ENDDO
  ! ..V2 B2=R  
  CALL zeroing(bm)!,n*n*nkry)
  CALL mgs(v(1,1,2),m,n,bm(1,1,2))
  ! ..   Setup A_1=diag(w)
  CALL zeroing(am)!,n*n)
  DO i=1,n
     am(i,i,1)=w(i)
  ENDDO
  ! ..V_3=H.V2
  IF (ttrot) THEN
     CALL ehpsi(v(:,:,2),v(:,:,3),vpot,psi,n)
  ELSE
     CALL stopgm(procedureN,'fix me: wrong dummy args + provide a regtest', &
          __LINE__,__FILE__)
     !fix me call hpsi(v(:,:,2),v(:,:,3),sc0,vpot,psi,n,1,ispin)
  ENDIF
  nhpsi=nhpsi+n
  ! .. project out converged states
  CALL prjcnv(m,n,nconv,vconv,w,v(1,1,3))
  ! ..V_3=V_3-V_1 B_2^T
  CALL rsdblk('T',m,n,v(1,1,1),bm(1,1,2),v(1,1,3))
  ! ..Ovlap 
  CALL ovlap2(m/2,n,n,am(1,1,2),v(1,1,2),v(1,1,3),.TRUE.)
  ! ==--------------------------------------------------------------==
  CALL tihalt('    PRPKRV',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE prpkrv
! ==================================================================
SUBROUTINE kryref(m,n,nkry,nconv,v0,w,v,am,bm,t,ldt,wt,ldwt,&
     nhpsi,vpot,psi,sc0,ttrot)
  ! ==================================================================
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE frsblk_utils, ONLY :  rsdblk, prjcnv, puttmat, fmdiag
  USE hpsi_utils, ONLY : hpsi
  USE summat_utils, ONLY : summat
  USE fft_maxfft,                      ONLY: maxfftn
  USE zeroing_utils,                   ONLY: zeroing
  IMPLICIT NONE
  INTEGER                                    :: m, n, nkry, nconv
  REAL(real_8)                               :: v0(m,*), w(n), v(m,n,nkry+1), &
                                                am(n,n,nkry+1), bm(n,n,nkry), &
                                                t(*)
  INTEGER                                    :: ldt, ldwt
  REAL(real_8)                               :: wt(ldwt)
  INTEGER                                    :: nhpsi
  REAL(real_8)                               :: vpot(*)
  COMPLEX(real_8)                            :: psi(maxfftn), sc0(*)
  LOGICAL                                    :: ttrot

  CHARACTER(*), PARAMETER                    :: procedureN = 'kryref'
  INTEGER, PARAMETER                         :: ispin = 1

  INTEGER                                    :: i, isub, j

  CALL tiset('    KRYREF',isub)
  DO i=3,nkry
     ! ..V_I = V_I - V_(I-1) A_(I-1) 
     CALL rsdblk('N',m,n,v(1,1,i-1),am(1,1,i-1),v(1,1,i))
     ! ..V_I B_I = V_I
     CALL mgs(v(1,1,i),m,n,bm(1,1,i))
     ! ..HPSI:  H V_I
     IF (ttrot) THEN
        CALL ehpsi(v(1,1,i),v(1,1,i+1),vpot,psi,n)
     ELSE
        CALL stopgm(procedureN,'fix me: wrong dummy args + provide a regtest', &
             __LINE__,__FILE__)
        !fix me call hpsi(v(1,1,i),v(1,1,i+1),sc0,vpot,psi,n,1,ispin)
     ENDIF
     nhpsi=nhpsi+n
     ! .. project out converged states
     CALL prjcnv(m,n,nconv,v0,w,v(1,1,i+1))
     ! ..V_(I+1) = H V_I - V_(I-1)B_I 
     CALL rsdblk('T',m,n,v(1,1,i-1),bm(1,1,i),v(1,1,i+1))
     ! ..Ovlap:  A_I=V_I^T V_(I+1)
     CALL ovlap2(m/2,n,n,am(1,1,i),v(1,1,i),v(1,1,i+1),.TRUE.)
     CALL summat(am(1,1,i),n)
  ENDDO
  ! ..Setup T-matrix. First diagonal terms      
  i=nkry
  CALL zeroing(t(1:i*n*i*n))!,i*n*i*n)
  CALL puttmat('N',t,i*n,i*n,am(1,1,1),n,n,1,1)
  DO j=2,i
     CALL puttmat('N',t,i*n,i*n,am(1,1,j),n,n,j,j)
     CALL puttmat('T',t,i*n,i*n,bm(1,1,j),n,n,j-1,j)
     CALL puttmat('N',t,i*n,i*n,bm(1,1,j),n,n,j,j-1)
  ENDDO
  ! ..Full matrix diag. 
  CALL fmdiag(i*n,n,t,wt)
  ! ..banded matrix diagonalisation. This seems to be slower than 
  ! ..full diag.  
  ! CALL BANDM(I*N,N,T,WT,SCR,LSCR,ISCR,5*I*N,ISCR(5*I*N+1))
  ! ==--------------------------------------------------------------==
  CALL tihalt('    KRYREF',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE kryref
! ==================================================================
! ==================================================================
SUBROUTINE mgs(a,m,n,r)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE rgs_utils, ONLY : rgs
  IMPLICIT NONE
  INTEGER                                    :: m, n
  COMPLEX(real_8)                            :: a(m,n)
  REAL(real_8)                               :: r(n,*)

  INTEGER                                    :: info, isub

  CALL tiset('        MGS',isub)
  CALL rgs(a,n,r)
  CALL dtrtri('U','N',n,r,n,info)
  CALL tihalt('        MGS',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE mgs
