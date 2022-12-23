MODULE frsblk_c_utils
  USE adjmu_utils,                     ONLY: convfrie
  USE ehpsi_utils,                     ONLY: ehpsi_c
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fint,                            ONLY: fint1
  USE geq0mod,                         ONLY: geq0
  USE gsortho_utils,                   ONLY: gs_ortho_c
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE prng_utils,                      ONLY: repprngu_vec_cmplx
  USE rgs_utils,                       ONLY: rgs,&
                                             rgs_c
  USE summat_utils,                    ONLY: give_scr_summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: frsblk_c
  PUBLIC :: give_scr_frsblk_c
  !public :: mgs_c
  !public :: puttmat_c
  !public :: getmat_c
  PUBLIC :: reorder_c
  !public :: kryref_c
  !public :: prpkrv_c
  !public :: puttab_c
  !public :: rsdblk_c
  !public :: fmdiag_c
  !public :: bandm_c
  !public :: prjcnv_c
  !public :: write_matrix_c

CONTAINS

  ! ==================================================================
  SUBROUTINE frsblk_c(ndiag,nel,nconv,nhpsi,c0,cs,sc0,cscr,vpot,psi,&
       edav,jspin,ikind,ikk,&
       trefine,tinfo)
    ! ==--------------------------------------------------------------==
    ! ==  Kohn-Sham Matrix Diagonalization by Block Krylov-Space      == 
    ! ==  Method                                                      ==
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: ndiag
    REAL(real_8)                             :: nel
    INTEGER                                  :: nconv, nhpsi
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,ndiag), &
                                                cs(nkpt%ngwk,*), &
                                                sc0(nkpt%ngwk,*), &
                                                cscr(nkpt%ngwk,*)
    REAL(real_8)                             :: vpot(*)
    COMPLEX(real_8)                          :: psi(maxfftn)
    REAL(real_8)                             :: edav(*)
    INTEGER                                  :: jspin, ikind, ikk
    LOGICAL                                  :: trefine, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'frsblk_c'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)
    INTEGER, PARAMETER                       :: ispin = 1 
    REAL(real_8), PARAMETER                  :: edavmax = 300._real_8, &
                                                edavmin = -300._real_8 , &
                                                eps = 1.e-14_real_8

    CHARACTER(len=100)                       :: formlanc
    CHARACTER(len=20)                        :: charspin
    COMPLEX(real_8), ALLOCATABLE             :: adav(:,:), alpha(:), beta(:), &
                                                tmat(:), work(:)
    INTEGER :: i, icycle, ierr, iformlanc, ig, isub, itest, j, liscr, lwork, &
      nbleff, nblk, nconvold, ncurr, ngwk2, nhpsiold, nk, nleft, ntest
    INTEGER, ALLOCATABLE                     :: INDEX(:), iscr(:)
    LOGICAL                                  :: tcheck, tconv, tlsd2, ttest
    REAL(real_8)                             :: amu, aux, b2, b2max, b2min, &
                                                tim, tim1, tim2, trotmax, &
                                                trotmin
    REAL(real_8), ALLOCATABLE                :: radav(:,:), rwork(:), wt(:)
    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==
! With LSD and 1 electron...

    IF (ndiag.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('  FRSBLK_C',isub)
    nk=cnti%nkry_block*cnti%nkry_max
    liscr=6*nk
    ! TODO check stat
    ! TODO align for BG
    ALLOCATE(iscr(liscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(INDEX(ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(tmat(3*nk*nk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(wt(nk),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(alpha(cnti%nkry_block*cnti%nkry_block*(cnti%nkry_max + 1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(beta(cnti%nkry_block*cnti%nkry_block*cnti%nkry_max),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(adav(ndiag, ndiag),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    tim1=m_walltime()
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    tconv=.FALSE.
    ngwk2=nkpt%ngwk*2
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
            11-2*iformlanc,&
            '("<")," BLOCK LANCZOS DIAGONALIZATION ",A19)'
       IF (paral%io_parent)&
            WRITE(6,formlanc) ikk,nkpt%nkpts,charspin
    ENDIF

    ALLOCATE(rwork(3*ndiag-2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(2*ndiag-1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    lwork = -1
    CALL zheev('V','U',ndiag,adav,ndiag,edav,work,lwork,&
         rwork,ierr)
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
          CALL ehpsi_c(c0,cs,vpot,psi,ndiag,ikind,maxfftn)
       ELSE
          CALL hpsi(c0,cs,sc0,vpot,psi,ndiag,ikind,ispin)
       ENDIF
       nhpsi=ndiag
       CALL ovlap2_c(nkpt%ngwk,ndiag,ndiag,adav,c0,cs)
       CALL sumhmat(adav,ndiag)
       CALL zheev('V','U',ndiag,adav,ndiag,edav,work,lwork,&
            rwork,ierr)
       CALL reorder_c(ndiag,ndiag,edav,adav)
       ! Rotate Orbitals
       CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,c0(1,1),nkpt%ngwk,adav(1,1),&
            ndiag,zzero,cscr(1,1),nkpt%ngwk)
       CALL dcopy(ngwk2*ndiag,cscr(1,1),1,c0(1,1),1)
       CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,cs(1,1),nkpt%ngwk,adav(1,1),&
            ndiag,zzero,cscr(1,1),nkpt%ngwk)
       CALL dcopy(ngwk2*ndiag,cscr(1,1),1,cs(1,1),1)
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
                  WRITE(6,*) 'FRSBLK_C| EIGENVECTOR',i,' IS VERY BAD! '
             CALL repprngu_vec_cmplx(nkpt%ngwk,c0(:,i))
             IF (geq0) c0(1,i)=CMPLX(REAL(c0(1,i)),0._real_8,kind=real_8)
             itest=i
          ENDIF
       ENDDO
       IF (tcheck) THEN
          ntest=ntest+1
          IF (ntest.GT.10.AND.paral%parent)&
               CALL stopgm(' FRSBLK_C', 'CAN''T FIND GUEST VECTORS',& 
               __LINE__,__FILE__)

          IF (itest.EQ.1.AND.ndiag.NE.1)&
               CALL dswap(ngwk2,c0(1,1),1,c0(1,ndiag),1)
          ALLOCATE(radav(ndiag, ndiag),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL rgs(c0,ndiag,radav)
          DEALLOCATE(radav,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDDO



    tim2=m_walltime()
    tim=(tim2-tim1)*0.001_real_8
    IF (paral%io_parent.AND.tinfo) THEN
       WRITE(6,'(A,T58,F8.2)')&
            ' >> TIME FOR INITIAL SUBSPACE DIAGONALIZATION:  ',tiM
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
          DO ig=1,nkpt%ngwk
             cscr(ig,1)=cs(ig,i)-edav(i)*c0(ig,i)
          ENDDO
          b2=ddot(ngwk2,cscr,1,cscr,1)
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
             CALL dswap(ngwk2,c0(1,i),1,c0(1,j),1)
             CALL dswap(ngwk2,cs(1,i),1,cs(1,j),1)
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
          CALL prpkrv_c(nkpt%ngwk,nbleff,cnti%nkry_max,nconv,c0(1,1),c0(1,i),&
               cs(1,i),cscr,alpha,beta,edav(i),&
               nhpsi,vpot,psi,sc0,ikind,fint1%ttrot)
          ! ==----------------------------------------------------------==
          ! Refinement Loop
          ! ==----------------------------------------------------------==
          nblk=nbleff*cnti%nkry_max
          CALL kryref_c(nkpt%ngwk,nbleff,cnti%nkry_max,nconv,c0,edav,cscr,alpha,&
               beta,tmat,nblk*nblk,wt,nblk,nhpsi,&
               vpot, psi,sc0,ikind,fint1%ttrot)
          ! ==----------------------------------------------------------==
          ! ..V=[V_1 V_2.. V_L] Y
          CALL zgemm('N','N',nkpt%ngwk,nbleff,nblk,zone,cscr,nkpt%ngwk,tmat,&
               nblk,zzero,c0(1,i),nkpt%ngwk)
          ! ==----------------------------------------------------------==
          ! ==  End of refinement over states                           ==
          ! ==----------------------------------------------------------==
       ENDDO
       ! ..New
       ! Reorthogonalize
       IF (ncurr.LE.ndiag) THEN
          CALL gs_ortho_c(c0,ncurr-1,c0(1,ncurr),ndiag-ncurr+1,adav)
       ENDIF
       ! Calculate new forces; only for states that entered Lanczos
       nleft=ndiag-ncurr+1
       IF (ncurr.LE.ndiag) THEN
          IF (fint1%ttrot) THEN
             CALL ehpsi_c(c0(1,ncurr),cs(1,ncurr),vpot,psi,nleft,ikind,maxfftn)
          ELSE
             CALL hpsi(c0(:,ncurr:ncurr+nleft-1),cs(1,ncurr),sc0(1,ncurr),&
                  vpot,psi,nleft,ikind,ispin)
          ENDIF
          nhpsi=nhpsi+nleft
       ENDIF
       ! Diagonalize Kohn-Sham Matrix
       CALL ovlap2_c(nkpt%ngwk,ndiag,ndiag,adav,c0,cs)
       CALL sumhmat(adav,ndiag)
       CALL zheev('V','U',ndiag,adav,ndiag,edav,work,lwork,&
            rwork,ierr)
       ! We inverse the order
       ! (not the same one for DSYEV than for this routine).
       CALL reorder_c(ndiag,ndiag,edav,adav)
       CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,c0(1,1),nkpt%ngwk,adav(1,1),&
            ndiag,zzero,cscr(1,1),nkpt%ngwk)
       CALL dcopy(ngwk2*ndiag,cscr(1,1),1,c0(1,1),1)
       CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,cs(1,1),nkpt%ngwk,adav(1,1),&
            ndiag,zzero,cscr(1,1),nkpt%ngwk)
       CALL dcopy(ngwk2*ndiag,cscr(1,1),1,cs(1,1),1)
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
                  WRITE(6,*) ' FRSBLK| ','VERY BAD EIGENVALUE ! !'
             IF (cnti%ntrans.GT.ndiag) THEN
                CALL repprngu_vec_cmplx(nkpt%ngwk,c0(:,i))
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
    IF (paral%io_parent.AND.((.NOT.tconv).AND.nconv.NE.ndiag)) THEN
       WRITE(6,'(1X,64("! "))')
       WRITE(6,'(" ! !",A,T64,"!!")')&
            ' FRSBLK| NOT ALL ROOTS ARE CONVERGED'
       WRITE(6,'(1X,64("! "))')
    ENDIF
    cntl%tlsd=tlsd2
    CALL tihalt('  FRSBLK_C',isub)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(iscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(index,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(tmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(wt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(alpha,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(beta,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(adav,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE frsblk_c
  ! ==================================================================
  SUBROUTINE give_scr_frsblk_c(lfrsblk,tag,ndiag)
    ! ==--------------------------------------------------------------==

    INTEGER                                  :: lfrsblk
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: ndiag

    INTEGER                                  :: lhpsi, lsummat

    CALL give_scr_hpsi(lhpsi,tag,ndiag)
    CALL give_scr_summat(lsummat,tag,ndiag)
    lfrsblk=MAX(lhpsi,lsummat,23*ndiag*cnti%nkry_block)+10
    tag='MAX(LHPSI,LSUMHMAT,23*NDIAG...' 
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_frsblk_c
  ! ==================================================================
  SUBROUTINE mgs_c(a,n,r)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: a(*)
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: r(n,*)

    INTEGER                                  :: info, isub

! Variables
! ==================================================================

    CALL tiset('     MGS_C',isub)
    CALL rgs_c(a,n,r)
    CALL ztrtri('U','N',n,r,n,info)
    CALL tihalt('     MGS_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mgs_c
  ! ===================================================================
  SUBROUTINE puttmat_c(char,t,m,n,a,ma,na,ibeg,jbeg)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: char
    INTEGER                                  :: m, n
    COMPLEX(real_8)                          :: t(m,n)
    INTEGER                                  :: ma, na
    COMPLEX(real_8)                          :: a(ma,na)
    INTEGER                                  :: ibeg, jbeg

    INTEGER                                  :: i, j

    IF (ma*ibeg.GT.m) CALL stopgm('PUTTMAT_C',' MA*IBEG.GT.M',& 
         __LINE__,__FILE__)
    IF (na*jbeg.GT.n) CALL stopgm('PUTTMAT_C',' NA*JBEG.GT.N',& 
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
             t(i+ma*(ibeg-1),j+na*(jbeg-1))=CONJG(a(j,i))
          ENDDO
       ENDDO
    ELSE
       CALL stopgm('PUTTMAT_C','ILLEGAL CHAR',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE puttmat_c
  ! ================================================================
  SUBROUTINE getmat_c(char,t,m,n,a,ma,na,ibeg,jbeg)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: char
    INTEGER                                  :: m, n
    COMPLEX(real_8)                          :: t(m,n)
    INTEGER                                  :: ma, na
    COMPLEX(real_8)                          :: a(ma,na)
    INTEGER                                  :: ibeg, jbeg

    INTEGER                                  :: i, j

    IF (ma*ibeg.GT.m) CALL stopgm('GETMAT_C',' MA*IBEG.GT.M',& 
         __LINE__,__FILE__)
    IF (na*jbeg.GT.n) CALL stopgm('GETMAT_C',' NA*JBEG.GT.N',& 
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
             a(j,i)=CONJG(t(i+ma*(ibeg-1),j+na*(jbeg-1)))
          ENDDO
       ENDDO
    ELSE
       CALL stopgm('GETMAT_C','ILLEGAL CHAR',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getmat_c
  ! ==================================================================
  SUBROUTINE reorder_c(m,n,w,a)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n
    REAL(real_8)                             :: w(n)
    COMPLEX(real_8)                          :: a(m,n)

    INTEGER                                  :: j
    REAL(real_8)                             :: aux

    DO j=1,n/2
       CALL dswap(2*m,a(1,j),1,a(1,n-j+1),1)
       aux=w(n-j+1)
       w(n-j+1)=w(j)
       w(j)=aux
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE reorder_c
  ! ==================================================================
  SUBROUTINE kryref_c(m,n,nkry,nconv,v0,w,v,am,bm,t,ldt,wt,ldwt,&
       nhpsi,vpot,psi,sc0,ikind,ttrot)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n, nkry, nconv
    COMPLEX(real_8)                          :: v0(m,*)
    REAL(real_8)                             :: w(n)
    COMPLEX(real_8)                          :: v(m,n,nkry+1), &
                                                am(n,n,nkry+1), bm(n,n,nkry), &
                                                t(*)
    INTEGER                                  :: ldt, ldwt
    REAL(real_8)                             :: wt(ldwt)
    INTEGER                                  :: nhpsi
    REAL(real_8)                             :: vpot(*)
    COMPLEX(real_8)                          :: psi(:), sc0(*)
    INTEGER                                  :: ikind
    LOGICAL                                  :: ttrot

    INTEGER, PARAMETER                       :: ispin = 1 

    INTEGER                                  :: i, isub, j

! ==--------------------------------------------------------------==

    CALL tiset('  KRYREF_C',isub)
    DO i=3,nkry
       ! ..V_I = V_I - V_(I-1) A_(I-1) 
       CALL rsdblk_c('N',m,n,v(1,1,i-1),am(1,1,i-1),v(1,1,i))
       ! ..V_I B_I = V_I
       CALL mgs_c(v(1,1,i),n,bm(1,1,i))
       ! ..HPSI:  H V_I
       IF (ttrot) THEN
          CALL ehpsi_c(v(1,1,i),v(1,1,i+1),vpot,psi,n,ikind,maxfftn)
       ELSE
          CALL hpsi(v(:,:,i),v(1,1,i+1),sc0,vpot,psi,n,ikind,&
               ispin)
       ENDIF
       nhpsi=nhpsi+n
       ! .. project out converged states
       CALL prjcnv_c(m,n,nconv,v0,w,v(1,1,i+1))
       ! ..V_(I+1) = H V_I - V_(I-1)B_I 
       CALL rsdblk_c('T',m,n,v(1,1,i-1),bm(1,1,i),v(1,1,i+1))
       ! ..Ovlap:  A_I=V_I^T V_(I+1)
       CALL ovlap2_c(m,n,n,am(1,1,i),v(1,1,i),v(1,1,i+1))
       CALL sumhmat(am(1,1,i),n)
    ENDDO
    ! ..Setup T-matrix. First diagonal terms      
    i=nkry
    CALL zeroing(t(1:i*n*i*n))!,i*n*i*n)
    CALL puttmat_c('N',t,i*n,i*n,am(1,1,1),n,n,1,1)
    DO j=2,i
       CALL puttmat_c('N',t,i*n,i*n,am(1,1,j),n,n,j,j)
       CALL puttmat_c('T',t,i*n,i*n,bm(1,1,j),n,n,j-1,j)
       CALL puttmat_c('N',t,i*n,i*n,bm(1,1,j),n,n,j,j-1)
    ENDDO
    ! ..Full matrix diag. 
    CALL fmdiag_c(i*n,n,t,wt)
    ! ..banded matrix diagonalisation. This seems to be slower than 
    ! ..full diag.  
    ! CALL BANDM(I*N,N,T,WT,SCR,LSCR,ISCR,5*I*N,ISCR(5*I*N+1))
    ! ..
    ! ==--------------------------------------------------------------==
    CALL tihalt('  KRYREF_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kryref_c
  ! ==================================================================
  SUBROUTINE prpkrv_c(m,n,nkry,nconv,vconv,v0,vs,v,am,bm,w,&
       nhpsi,vpot,psi,sc0,ikind,ttrot)
    ! ==--------------------------------------------------------------==
    ! == Returns A_1,A_2,B_2 and V in a form suitable for Krylov      ==
    ! == refinement                                                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n, nkry, nconv
    COMPLEX(real_8)                          :: vconv(m,*), v0(m,n), vs(m,n), &
                                                v(m,n,nkry+1), &
                                                am(n,n,nkry+1), bm(n,n,nkry)
    REAL(real_8)                             :: w(n)
    INTEGER                                  :: nhpsi
    REAL(real_8)                             :: vpot(*)
    COMPLEX(real_8)                          :: psi(:), sc0(*)
    INTEGER                                  :: ikind
    LOGICAL                                  :: ttrot

    INTEGER, PARAMETER                       :: ispin = 1 

    INTEGER                                  :: i, isub

! ==---------------------------------------------------------------==

    CALL tiset('  PRPKRV_C',isub)
    ! ==---------------------------------------------------------------==
    CALL dcopy(2*m*n,v0,1,v(1,1,1),1)
    CALL dcopy(2*m*n,vs,1,v(1,1,2),1)
    ! ..Compute residual 
    DO i=1,n
       CALL daxpy(2*m,-w(i),v(1,i,1),1,v(1,i,2),1)
    ENDDO
    ! ..V2 B2=R  
    CALL zeroing(bm)!,n*n*nkry)
    CALL mgs_c(v(1,1,2),n,bm(1,1,2))
    ! ..   Setup A_1=diag(w)
    CALL zeroing(am)!,n*n)
    DO i=1,n
       am(i,i,1)=CMPLX(w(i),0._real_8,kind=real_8)
    ENDDO
    ! ..V_3=H.V2
    IF (ttrot) THEN
       CALL ehpsi_c(v(1,1,2),v(1,1,3),vpot,psi,n,ikind,maxfftn)
    ELSE
       CALL hpsi(v(:,:,2),v(1,1,3),sc0,vpot,psi,&
            n,ikind,ispin)
    ENDIF
    nhpsi=nhpsi+n
    ! .. project out converged states
    CALL prjcnv_c(m,n,nconv,vconv,w,v(1,1,3))
    ! ..V_3=V_3-V_1 B_2^T
    CALL rsdblk_c('T',m,n,v(1,1,1),bm(1,1,2),v(1,1,3))
    ! ..Ovlap 
    CALL ovlap2_c(m,n,n,am(1,1,2),v(1,1,2),v(1,1,3))
    CALL sumhmat(am(1,1,2),n)
    ! ==--------------------------------------------------------------==
    CALL tihalt('  PRPKRV_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prpkrv_c
  ! ==================================================================
  SUBROUTINE puttab_c(n,kd,a)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, kd
    COMPLEX(real_8)                          :: a(n,n)

    CHARACTER(*), PARAMETER                  :: procedureN = 'puttab_c'

    COMPLEX(real_8), ALLOCATABLE             :: aux(:)
    INTEGER                                  :: i, ibeg, ierr, j

    ALLOCATE(aux(n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(aux)!,n)
    DO j=1,n
       CALL zcopy(n,a(1,j),1,aux,1)
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
  END SUBROUTINE puttab_c
  ! ==================================================================
  SUBROUTINE rsdblk_c(char,m,n,v1,a,v2)
    ! ==--------------------------------------------------------------==
    ! == V_I = V_I - V_(I-1) A_(I-1)                                  ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: char
    INTEGER                                  :: m, n
    COMPLEX(real_8)                          :: v1(m,n), a(n,n), v2(m,n)

    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) 

    INTEGER                                  :: isub

! ==--------------------------------------------------------------==

    CALL tiset('  RSDBLK_C',isub)
    IF (char.EQ.'N') THEN
       CALL zgemm('N','N',m,n,n,-zone,v1,m,a,n,zone,v2,m)
    ELSEIF (char.EQ.'T') THEN
       CALL zgemm('N','C',m,n,n,-zone,v1,m,a,n,zone,v2,m)
    ELSE
       CALL stopgm('RSDBLK_C','ILLEGAL CHAR',& 
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt('  RSDBLK_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rsdblk_c
  ! ==================================================================
  SUBROUTINE fmdiag_c(m,n,t,wt)
    ! ==--------------------------------------------------------------==
    ! == Returns the n largest eigenvalues of MxM symetric matrix T   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n
    COMPLEX(real_8)                          :: t(m,m,*)
    REAL(real_8)                             :: wt(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fmdiag_c'
    REAL(real_8), PARAMETER                  :: small = 1.e-14_real_8

    COMPLEX(real_8), ALLOCATABLE             :: work(:)
    INTEGER                                  :: ierr, info, isub, lwork, meval
    INTEGER, ALLOCATABLE                     :: ifail(:), iwork(:)
    REAL(real_8), ALLOCATABLE                :: rwork(:)

    CALL tiset('  FMDIAG_C',isub)
    IF (m.EQ.n) THEN
       ALLOCATE(rwork(MAX(1,3*n-2)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       lwork = MAX(1,2*n-1)! TODO determine optimal LWORK
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL zheev('V','U',n,t,n,wt,work,lwork,rwork,info)
       DEALLOCATE(rwork,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(rwork(7*m),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       lwork = 2*n! TODO determine optimal LWORK
       ALLOCATE(work(lwork),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(iwork(5*n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(ifail(n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL zheevx('V','I','U',m,t,m,1._real_8,1._real_8,m-n+1,m,small,&
            meval,wt,t(1,1,2),m,work,lwork,rwork,iwork,ifail,info)
       DEALLOCATE(work,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(rwork,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(iwork,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(ifail,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       IF (meval.LT.n) THEN
          IF (paral%io_parent)&
               WRITE(6,*) ' WARNING| ZHEEVX RETURNED MEVAL < N',meval,n
       ENDIF
       CALL dcopy(2*m*n,t(1,1,2),1,t(1,1,1),1)
    ENDIF
    IF (info.NE.0) CALL stopgm('FMDIAG_C','FMDIAG FAILED',& 
         __LINE__,__FILE__)
    CALL reorder_c(m,n,wt,t)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(iwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('  FMDIAG_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fmdiag_c
  ! ==================================================================
  SUBROUTINE bandm_c(in,n,t,wt,ifail)
    ! ==--------------------------------------------------------------==
    ! == Returns the n largest eigenvalues                            ==
    ! == of MxM symetric banded matrix T                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: in, n
    COMPLEX(real_8)                          :: t(in,in,3)
    REAL(real_8)                             :: wt(in)
    INTEGER                                  :: ifail(in)

    CHARACTER(*), PARAMETER                  :: procedureN = 'bandm_c'

    COMPLEX(real_8), ALLOCATABLE             :: work(:)
    INTEGER                                  :: ierr, info, isub, meval
    INTEGER, ALLOCATABLE                     :: iwork(:)
    REAL(real_8), ALLOCATABLE                :: rwork(:)

    CALL tiset('   BANDM_C',isub)
    ! ..Put T is form suitable for banded matrix diagonalisation 
    CALL puttab_c(in,n,t)

    ALLOCATE(work(in),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rwork(7*in),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(iwork(5*in),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zhbevx('V','I','U',in,n,t,in,t(1,1,2),in,1._real_8,1._real_8,&
         in-n+1,in,1.e-10_real_8,meval,wt,t(1,1,3),in,work,rwork,iwork,&
         ifail,info)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(iwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    IF (meval.LT.n) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' WARNING| ZHBEVX RETURNED MEVAL < N',meval,n
    ENDIF
    CALL dcopy(2*in*in,t(1,1,3),1,t(1,1,1),1)
    CALL reorder_c(in,n,wt,t)
    CALL tihalt('   BANDM_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE bandm_c
  ! ==================================================================
  SUBROUTINE prjcnv_c(m,np,n0,v0,w,v)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, np, n0
    COMPLEX(real_8)                          :: v0(m,n0)
    REAL(real_8)                             :: w(n0)
    COMPLEX(real_8)                          :: v(m,np)

    CHARACTER(*), PARAMETER                  :: procedureN = 'prjcnv_c'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8)

    COMPLEX(real_8), ALLOCATABLE             :: smat(:,:)
    INTEGER                                  :: ierr, isub

    IF (n0.EQ.0) THEN
       RETURN
    ELSE
       RETURN
    ENDIF
    CALL tiset('  PRJCNV_C',isub)
    ALLOCATE(smat(n0,np),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ..SMAT=V0^T.V 
    CALL ovlap2_c(m,n0,np,smat,v0,v)
    CALL mp_sum(smat,2*np*n0,parai%allgrp)
    ! deb      DO J=1,N0
    ! deb        CALL DSCAL(NP,W(J),SMAT(J,1),N0)
    ! deb      ENDDO
    CALL zgemm('N','N',m,np,n0,-zone,v0,m,smat,n0,zone,v,m)
    DEALLOCATE(smat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('  PRJCNV_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE prjcnv_c
  ! ==================================================================
  SUBROUTINE write_matrix_c(char,m,n,a)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=10)                        :: char
    INTEGER                                  :: m, n
    COMPLEX(real_8)                          :: a(m,n)

    INTEGER                                  :: i, j

    IF (.NOT.paral%io_parent) RETURN
    WRITE(6,*) char
    DO i=1,m
       WRITE(6,'(12(1PE10.3,1PE10.3))') (a(i,j),j=1,n)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE write_matrix_c
  ! ==================================================================

END MODULE frsblk_c_utils
