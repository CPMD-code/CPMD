MODULE vibana_utils
  USE adat,                            ONLY: elem
  USE coor,                            ONLY: tau0,&
                                             taup
  USE cotr,                            ONLY: cotc0,&
                                             hess
  USE detdof_utils,                    ONLY: detdof
  USE error_handling,                  ONLY: stopgm
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE hessin_utils,                    ONLY: hessin
  USE hessout_utils,                   ONLY: hessout
  USE initrun_driver,                  ONLY: initrun
  USE initrun_utils,                   ONLY: give_scr_initrun
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE lscal,                           ONLY: mode
  USE parac,                           ONLY: paral
  USE rmas,                            ONLY: rmass
  USE secder_utils,                    ONLY: molvib,&
                                             purged,&
                                             vibeig
  USE store_types,                     ONLY: restart1
  USE symm,                            ONLY: symmi
  USE symtrz_utils,                    ONLY: give_scr_symmat,&
                                             symmat
  USE system,                          ONLY: maxsys
  USE utils,                           ONLY: symma
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vibana
  PUBLIC :: hesfile
  PUBLIC :: give_scr_vibana
  PUBLIC :: give_scr_hesfile

CONTAINS

  ! ==================================================================
  SUBROUTINE vibana(sder)
    ! ==--------------------------------------------------------------==
    ! ==  VIBRATIONAL ANALYSIS                                        ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: sder(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'vibana'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: i, ia, ierr, info, is, ix, j, &
                                                k, lscr, naux
    REAL(real_8)                             :: xmass
    REAL(real_8), ALLOCATABLE                :: scr(:), vibe(:), xma(:)

    IF (.NOT.paral%parent) RETURN
    CALL give_scr_vibana(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vibe(3*ions1%nat),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(xma(3*maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! Symmetrize Hessian
    IF (symmi%indpg.NE.0.AND.(symmi%nrot.NE.1.OR.symmi%ntvec.GT.1)) THEN
       CALL symmat(sder,2)
    ENDIF
    CALL symma(sder,3*ions1%nat)
    ! Write Hessian to file HESSIAN
    IF (cotc0%nodim.EQ.3*ions1%nat) THEN
       CALL dcopy(cotc0%nodim*cotc0%nodim,sder(1,1),1,hess(1,1),1)
       CALL hessout
    ELSE
       CALL dcopy(9*ions1%nat*ions1%nat,sder(1,1),1,hess(1,1),1)
    ENDIF
    ! Write output file for MOLVIB program
    CALL molvib(tau0,sder)
    ! Mass weighted force constants
    ix=0
    DO is=1,ions1%nsp
       xmass=SQRT(1._real_8/rmass%pma0(is))
       DO ia=1,ions0%na(is)
          DO k=1,3
             ix=ix+1
             xma(ix)=xmass
          ENDDO
       ENDDO
    ENDDO
    !$omp parallel do private(I,J)
    DO i=1,3*ions1%nat
       DO j=1,3*ions1%nat
          sder(i,j)=sder(i,j)*xma(i)*xma(j)
       ENDDO
    ENDDO
    ! Diagonalization
    naux=9*ions1%nat*ions1%nat
    info=0
    CALL dsyev('V','U',3*ions1%nat,sder,3*ions1%nat,vibe,scr,naux,info)
    IF (info.NE.0) CALL stopgm('SECDER','DSYEV INFO',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I)
    DO i=1,3*ions1%nat
       vibe(i)=SIGN(5140.487_real_8*SQRT(ABS(vibe(i))),vibe(i))
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(/," ",64("*"),/,A,/)')&
         ' HARMONIC FREQUENCIES [cm**-1]:'
    IF (paral%io_parent)&
         WRITE(6,'(4(F16.4))') (vibe(i),i=1,3*ions1%nat)
    ! Write output file for Eigenvectors   
    CALL vibeig(vibe,sder,3*ions1%nat,.TRUE.)
    ! Write a selected mode to standard output
    IF (mode.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",64("*"),/,A,I3,A,/)')&
            ' POSITIONS AND HESSIAN EIGENMODE ', mode,&
            ' [atomic units]:'
       ix=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ix=ix+1
             IF (paral%io_parent)&
                  WRITE(6,'(1X,I3,1X,A2,3F8.4,1X,3(1PE11.3))')&
                  ix,elem%el(ions0%iatyp(is)),&
                  (tau0(k,ia,is),k=1,3),&
                  (sder(k+(ix-1)*3,mode),k=1,3)
          ENDDO
       ENDDO
    ENDIF
    ! Purification
    IF (paral%io_parent)&
         WRITE(6,'(/,A)') ' PURIFICATION OF DYNAMICAL MATRIX'
    CALL dcopy(9*ions1%nat*ions1%nat,hess,1,sder,1)
    ! tu   purification only AFTER mass weighting of the Hessian
    ! tu   CALL PURGED(SDER,HESS,TAU0)
    !$omp parallel do private(I,J)
    DO i=1,3*ions1%nat
       DO j=1,3*ions1%nat
          sder(i,j)=sder(i,j)*xma(i)*xma(j)
       ENDDO
    ENDDO
    CALL purged(sder,hess,tau0)
    ! Diagonalization
    naux=9*ions1%nat*ions1%nat
    info=0
    CALL dsyev('V','U',3*ions1%nat,sder,3*ions1%nat,vibe,scr,naux,info)
    IF (info.NE.0) CALL stopgm('SECDER','DSYEV INFO',& 
         __LINE__,__FILE__)
    !$omp parallel do private(I)
    DO i=1,3*ions1%nat
       vibe(i)=SIGN(5140.487_real_8*SQRT(ABS(vibe(i))),vibe(i))
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(/," ",64("*"),/,A,/)')&
         ' HARMONIC FREQUENCIES [cm**-1]:'
    IF (paral%io_parent)&
         WRITE(6,'(4(F16.4))') (vibe(i),i=1,3*ions1%nat)
    IF (paral%io_parent)&
         WRITE(6,'(A,E16.8)') 'ChkSum(FREQ) = ',SUM(ABS(vibe(1:3*ions1%nat)))
    ! Write output file for Eigenvectors   
    CALL vibeig(vibe,sder,3*ions1%nat,.FALSE.)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(xma,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vibe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vibana
  ! ==================================================================
  SUBROUTINE hesfile(sder,eigv,c0,c2,sc0)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: sder(cotc0%nodim,cotc0%nodim),&
                                                eigv(*)
    COMPLEX(real_8)                          :: c0(:,:,:), c2(*), sc0(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hesfile'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: psi(:,:)
    INTEGER                                  :: ierr, il_psi_1d, il_psi_2d, &
                                                il_rhoe, il_rhoe_1d, &
                                                il_rhoe_2d, lscr
    INTEGER, DIMENSION(100) :: irec = (/1,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,&
      1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
      0,0,0,0,0,0,0,0,0/)
    REAL(real_8), ALLOCATABLE                :: rhoe(:,:), scr(:), tscr(:,:,:)

    CALL rhoe_psi_size(il_rhoe,il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d,&
         il_rhoe_1d=il_rhoe_1d,il_rhoe_2d=il_rhoe_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_hesfile(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       ALLOCATE(tscr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(tscr)!,3*maxsys%nax*maxsys%nsx)
    ELSE
       ! Avoid 'not allocated' runtime error
       ALLOCATE(taup(1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! 
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    IF (paral%parent) THEN
       CALL detdof(tau0,tscr)
       CALL hessin(tau0,restart1%rhe)
       CALL dcopy(9*ions1%nat*ions1%nat,hess,1,sder,1)
    ENDIF
    ! 
    IF (paral%parent) THEN
       DEALLOCATE(taup,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(tscr,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hesfile
  ! ==================================================================
  SUBROUTINE give_scr_vibana(lvibana,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvibana
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lsymmat

! ==--------------------------------------------------------------==

    lsymmat=0
    CALL give_scr_symmat(lsymmat,tag)
    lvibana=MAX(9*ions1%nat*ions1%nat+9*ions1%nat,lsymmat)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vibana
  ! ==================================================================
  SUBROUTINE give_scr_hesfile(lhesfile,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lhesfile
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: linitrun

    CALL give_scr_initrun(linitrun,tag)
    lhesfile=linitrun
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_hesfile
  ! ==================================================================

END MODULE vibana_utils
