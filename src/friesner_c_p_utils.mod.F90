MODULE friesner_c_p_utils
  USE ehpsi_utils,                     ONLY: ehpsi_c
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fint,                            ONLY: fint1
  USE frsblk_c_utils,                  ONLY: reorder_c
  USE gsortho_utils,                   ONLY: gs_ortho_c
  USE hpsi_utils,                      ONLY: hpsi
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE sort_utils,                      ONLY: sort2
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: frie_c_p

CONTAINS

  ! ==================================================================
  SUBROUTINE frie_c_p(ndiag,nel,nconv,nhpsi,c0,cs,sc0,cscr,vpot,psi,&
       edav,jspin,ikind,ikk,trefine,tinfo)
    ! ==--------------------------------------------------------------==
    ! ==  Kohn-Sham Matrix Diagonalization by Krylov-Space Method     ==
    ! ==  W.T. Pollard and R.A. Friesner                              ==
    ! ==  J. Chem. Phys. 99, 6742 (1993)                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ndiag
    REAL(real_8)                             :: nel
    INTEGER                                  :: nconv, nhpsi
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,ndiag), &
                                                cs(nkpt%ngwk,ndiag), &
                                                sc0(nkpt%ngwk,*), &
                                                cscr(nkpt%ngwk,*)
    REAL(real_8)                             :: vpot(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    REAL(real_8)                             :: edav(ndiag)
    INTEGER                                  :: jspin, ikind, ikk
    LOGICAL                                  :: trefine, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'frie_c_p'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)
    INTEGER, PARAMETER                       :: ispin = 1 

    CHARACTER(len=100)                       :: formlanc
    CHARACTER(len=20)                        :: charspin
    COMPLEX(real_8), ALLOCATABLE             :: adav(:,:), work(:)
    INTEGER                                  :: i, iaux, ierr, iformlanc, &
                                                istate, isub, j, lwork, &
                                                ncurr, ngwk2
    INTEGER, ALLOCATABLE                     :: INDEX(:)
    LOGICAL                                  :: tlsd2
    REAL(real_8)                             :: aux
    REAL(real_8), ALLOCATABLE                :: alpha0(:), beta1(:), rwork(:)

! ==--------------------------------------------------------------==
! With LSD and 1 electron...

    IF (ndiag.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('FRIE_C_mia',isub)
    ! TODO check stat
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
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    ngwk2=2*nkpt%ngwk
    nconv=0
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
    IF (fint1%ttrot) THEN
       CALL ehpsi_c(c0,cs,vpot,psi,ndiag,ikind,maxfft)
    ELSE
       CALL hpsi(c0,cs,sc0,vpot,psi,ndiag,ikind,ispin)
    ENDIF
    nhpsi=ndiag
    CALL ovlap_h(ndiag,adav,c0,cs)
    CALL sumhmat(adav,ndiag)

    ! ==----------------------------------------------------------------==
    ! ==                 TEST  HAMILTONIAN                              ==
    ! 
    IF (paral%io_parent) THEN
       WRITE(6,'("HAMILTONIAN LEFT :"/)')
       j = 1
       DO istate = 1,ndiag
          WRITE(6,'("HAM ",I5,I5,5(f11.6,f11.6))')&
               istate,j,adav(istate,j),adav(istate,j+1),&
               adav(istate,j+2),adav(istate,j+3),&
               adav(istate,j+4)
       ENDDO
       WRITE(6,'(/"HAMILTONIAN RIGHT :"/)')
       j = ndiag/2+1
       DO istate = 1,ndiag
          WRITE(6,'("HAM ",I5,I5,5(f11.6,f11.6))')&
               istate,j,adav(istate,j),adav(istate,j+1),&
               adav(istate,j+2),adav(istate,j+3),&
               adav(istate,j+4)
       ENDDO

    ENDIF

    ALLOCATE(work(2*ndiag-1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rwork(3*ndiag-2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    lwork = -1
    CALL zheev('V','U',ndiag,adav,ndiag,edav,&
         work,lwork,rwork,ierr)
    lwork = work(1)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(work(lwork),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zheev('V','U',ndiag,adav,ndiag,edav,&
         work,lwork,rwork,ierr)
    DEALLOCATE(work,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rwork,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL reorder_c(ndiag,ndiag,edav,adav)
    ! Rotate Orbitals
    CALL zgemm('N','N',nkpt%ngwk,ndiag,ndiag,zone,c0(1,1),nkpt%ngwk,&
         adav(1,1),ndiag,zzero,cscr(1,1),nkpt%ngwk)
    CALL dcopy(ngwk2*ndiag,cscr(1,1),1,c0(1,1),1)
    ncurr = 1
    CALL gs_ortho_c(c0,ncurr-1,c0(1,ncurr),ndiag-ncurr+1,adav)

    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
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
    CALL tihalt('FRIE_C_mia',isub)
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
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE frie_c_p
  ! ==================================================================


END MODULE friesner_c_p_utils
