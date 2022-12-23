MODULE hpsi_utils
  USE coor,                            ONLY: tau0,fion
  USE error_handling,                  ONLY: stopgm
  USE hubbardu,                        ONLY: hubbu
  USE hubbardu_utils,                  ONLY: hubbardUcorrection,&
                                             give_scr_hubbardu,&
                                             add_hubbardu
  USE fft_maxfft,                      ONLY: maxfft
  USE fnonloc_utils,                   ONLY: fnonloc,&
                                             give_scr_fnonloc
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE pslo,                            ONLY: pslo_com
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE spin,                            ONLY: lspin2
  USE spsi_utils,                      ONLY: spsi
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE vpsi_utils,                      ONLY: vpsi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hpsi
  PUBLIC :: give_scr_hpsi

CONTAINS

  ! ==================================================================
  SUBROUTINE hpsi(c0,c2,sc0,vpot,psi,nstate,ikind,ispin)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE GENERALIZED FORCE C2=CMPLX(DFR,DFI) ACTING ON ALL THE   ==
    ! ==  ELECTRONIC STATES REPRESENTED BY THE VECTOR C0=CMPLX(CR,CI) ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:                                                       ==
    ! ==   C0(NGWK,NSTATE) Input Wavefunctions                        ==
    ! ==   VPOT    Local Potential (in real space)                    ==
    ! ==   LSCR    Dimension of scratch array SCR                     ==
    ! ==   NSTATE  Number of states                                   ==
    ! ==   IKIND   Index of k point                                   ==
    ! ==   ISPIN   Need with LSD option for diagonalization scheme    ==
    ! ==         dimension of VPOT(NNR1,ISPIN)                        ==
    ! == OUTPUT:                                                      ==
    ! ==   C2 = H|C0>                                                 ==
    ! ==   PSI              Used for FFT                              ==
    ! ==   SCR(LSCR)        Scratch array                             ==
    ! ==   SC0(NGWK,NSTATE) Used only by TIVAN                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(nkpt%ngwk,nstate), &
                                                c2(nkpt%ngwk,nstate)
    INTEGER                                  :: ikind, ispin
    REAL(real_8)                             :: vpot(fpar%nnr1,ispin)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hpsi'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: auxc(:), pab(:)
    INTEGER                                  :: i, ierr, il_auxc, il_ddia, &
                                                il_pab, isub, lrnlsm
    REAL(real_8), ALLOCATABLE                :: foc(:)
    COMPLEX(real_8), ALLOCATABLE             :: C2U(:,:)
    INTEGER                                  :: ISTATE, IG
! VARIABLES
! ==--------------------------------------------------------------==

    CALL tiset('      HPSI',isub)
    IF (lspin2%tlse) THEN
       il_pab=2*maxfft
    ELSE
       il_pab = 2
    ENDIF
    ALLOCATE(foc(nstate),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == INITIALIZE ELEC-FORCE ARRAYS                                 ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(c2)!,nkpt%ngwk*nstate)
    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    lrnlsm=lrnlsm+100
    ALLOCATE(pab(il_pab/2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    ALLOCATE(auxc(lrnlsm),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL rnlsm(c0,nstate,1,ikind,.FALSE.)
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DO i=1,nstate
       foc(i)=1.0_real_8
    ENDDO
    IF (pslo_com%tivan) THEN
       CALL dcopy(2*nkpt%ngwk*nstate,c0(1,1),1,sc0(1,1),1)
       CALL spsi(nstate,sc0)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == Compute the force on the electronic degrees of freedom due   ==
    ! == to the local potential (stored in VPOT)                      ==
    ! ==--------------------------------------------------------------==
    CALL vpsi(c0,c2,foc,vpot,psi,nstate,ikind,ispin,.TRUE.)
    IF(cntl%thubb)THEN
      ALLOCATE(c2u(ncpw%ngw,nstate),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
      CALL hubbardUcorrection(C0,C2U,TAU0,FION,NSTATE,PSI,.FALSE.,ISPIN)
      CALL add_hubbardu(c2,c2u,nstate)
      DEALLOCATE(c2u,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    END IF 
    ! ==--------------------------------------------------------------==
    ! == Compute the force on the electronic degrees of freedom due   ==
    ! == to the non-local part of the potential, and add it to the    ==
    ! == other piece, coming from the local contribution.             ==
    ! ==--------------------------------------------------------------==
    CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
    CALL fnonloc(c2,foc,nstate,ikind,ispin,.TRUE.)
    ! kpt  Seems to be unnecessary.
    IF (geq0) THEN
       IF (tkpts%tkpnt) THEN
          CALL zclean_k(c2,nstate,ncpw%ngw)
       ELSE
          CALL zclean(c2,nstate,ncpw%ngw)
       ENDIF
    ENDIF
    IF (cntl%tdavi) CALL dscal(2*nkpt%ngwk*nstate,-1.0_real_8,c2(1,1),1)
    DEALLOCATE(pab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(foc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('      HPSI',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hpsi
  ! ==================================================================
  SUBROUTINE give_scr_hpsi(lhpsi,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lhpsi
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: il_auxc, il_ddia, lrnlsm
    INTEGER                                  :: lhubbu

    CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_fnonloc(il_auxc,il_ddia,nstate)
    lhpsi=nstate+             & ! FOC
         il_ddia+             & ! DDIA
         MAX(il_auxc,lrnlsm)+  & ! AUXC and SCR
         100                  ! SCRATCH tool
    IF (lspin2%tlse) lhpsi=lhpsi+2*maxfft
  
    
    CALL give_scr_hubbardu(nstate,lhubbu,tag)
    lhpsi=MAX(lhpsi,lhubbu)
    tag ='NSTATE+IL_DDIA+MAX(IL_AUXC...)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_hpsi
  ! ==================================================================

END MODULE hpsi_utils
