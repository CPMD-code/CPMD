MODULE rinitwf_driver
  USE atomwf_utils,                    ONLY: atomwf
  USE atwf,                            ONLY: atwp,&
                                             loadc_foc_array_size
  USE copot_utils,                     ONLY: copot
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE forcedr_driver,                  ONLY: forcedr
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_sum
  USE newcell_utils,                   ONLY: newcell
  USE nlcc,                            ONLY: corel
  USE ortho_utils,                     ONLY: ortho
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE pslo,                            ONLY: pslo_com
  USE randtowf_utils,                  ONLY: randtowf
  USE reshaper,                        ONLY: reshape_inplace
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: bsnfi,&
                                             iteropt,&
                                             ropt_mod
  USE setbasis_utils,                  ONLY: loadc
  USE sphe,                            ONLY: maskgw,&
                                             tsphere
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2bye
  USE utils,                           ONLY: setkwf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rinitwf

CONTAINS

  ! ==================================================================
  SUBROUTINE rinitwf(c0,c2,sc0,nstate,tau0,fion,rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION OF WAVEFUNCTION                               ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:,:), c2(*), sc0(*)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rinitwf'

    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: isub
    LOGICAL                                  :: qmmm_s

!(nkpt%ngwk,nstate,*)
! RHOE(NNR1,*)
! Variables
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    tag=''
    IF (cntl%tprcp.OR.cntl%tpres) CALL newcell
    IF (tkpts%tkblock) CALL rkpt_swap(c0,nstate,1,'HGKP HGKM MASKGW TWNL')
    IF (cnti%inwfun.EQ.1) THEN
       qmmm_s=lqmmm%qmmm
       lqmmm%qmmm=.FALSE.
       CALL randwf(c0,c2,sc0,nstate,tau0,fion,rhoe,psi)
       lqmmm%qmmm=qmmm_s
    ELSEIF (cnti%inwfun.EQ.2) THEN
       CALL atomwf(c0,nstate,tau0,fion,rhoe,psi)
    ELSEIF (cnti%inwfun.EQ.3) THEN
       CALL simplewf(c0,c2,nstate,tau0)
    ELSE
       IF (paral%io_parent) WRITE(6,*) ' RINITWF| UNKNOWN OPTION'
       CALL stopgm('RINITWF',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rinitwf
  ! ==================================================================
  SUBROUTINE simplewf(c0,c2,nstate,tau0)
    ! ==--------------------------------------------------------------==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate,*)
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'simplewf'

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: CATOM_loc
    INTEGER                                  :: ia, iaorb, iat, ierr, is, &
                                                isub, ixx, natst
    REAL(real_8)                             :: foc(loadc_foc_array_size), sfc
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: overlap

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! ==  PHASE FACTORS                                               ==
    ! ==--------------------------------------------------------------==
    CALL phfac(tau0)
    ! ==--------------------------------------------------------------==

    ALLOCATE(CATOM_loc(nkpt%ngwk,atwp%nattot),overlap(atwp%nattot,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Allocation problem' ,& 
         __LINE__,__FILE__)
    overlap(:,:)=0.0_real_8
    CATOM_loc(:,:)=0.0_real_8

    ! load atomic guess to PW basis
    iaorb=1
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          CALL loadc(CATOM_loc(1,iaorb),foc,ncpw%ngw,ncpw%ngw,atwp%nattot-iaorb+1,&
               SIZE(foc),is,iat,natst)
          DO ixx=iaorb,iaorb+natst-1
             sfc=dotp(ncpw%ngw,CATOM_loc(:,ixx),catom_loc(:,ixx))
             CALL mp_sum(sfc,parai%allgrp)
             IF (sfc.EQ.0._real_8) THEN
                IF (paral%io_parent) WRITE(6,'(A,A,I5,A)') ' ATRHO|',&
                     ' THE NORM OF ATOMIC ORBITAL (',ixx,') IS NULL'
                CALL stopgm(procedureN,'WRONG ATOMIC ORBITAL',& 
                     __LINE__,__FILE__)
             ELSE
                sfc=1._real_8/SQRT(sfc)
             ENDIF
             CALL dscal(2*ncpw%ngw,sfc,CATOM_loc(1,ixx),1)
          ENDDO
          iaorb=iaorb+natst
       ENDDO
    ENDDO

    ! Init random C0
    CALL randtowf(c2(:,:,1),nstate,0,1)

    ! Overlap matrix
    CALL ovlap2(ncpw%ngw,atwp%nattot,nstate,overlap,CATOM_loc,c2,.FALSE.)
    CALL mp_sum(overlap,atwp%nattot*nstate,parai%allgrp)

    ! Project
    IF (ncpw%ngw>0) THEN
       CALL dgemm('N','N',2*ncpw%ngw,nstate,atwp%nattot,1.0_real_8,CATOM_loc,2*ncpw%ngw,&
            overlap,atwp%nattot,0.0_real_8,c0,2*ncpw%ngw)
    ENDIF

    ! Ortho
    CALL ortho(nstate,c0(:,:,1),c2)

    ! The other kpoints are the same atomic wf (C0 Gamma).
    IF (tkpts%tkpnt) THEN
       CALL stopgm( procedureN, 'Dont know if that is working' ,& 
            __LINE__,__FILE__)
       CALL setkwf(ncpw%ngw,nstate,c0)
    ENDIF

    DEALLOCATE(CATOM_loc,overlap,stat=ierr)
    IF (ierr.NE.0) CALL stopgm( procedureN, 'Deallocation problem' ,& 
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE simplewf
  ! ==================================================================
  SUBROUTINE randwf(c0,c2,sc0,nstate,tau0,fion,rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! == PERFORMS AN INITIAL STEP IN SELF CONSISTENCY STARTING WITH   ==
    ! == THE RANDOM 'C0'                                              ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(nkpt%ngwk,nstate,*)
    COMPLEX(real_8), TARGET                  :: c2(nkpt%ngwk,nstate,*)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'randwf'

    COMPLEX(real_8), DIMENSION(:, :), &
      POINTER                                :: c2_ptr
    INTEGER                                  :: ierr, ik, ikk, ikpt, isub, &
                                                kbeg, kend, kinc, nkpoint
    LOGICAL                                  :: tfor
    REAL(real_8)                             :: cfac
    REAL(real_8), ALLOCATABLE                :: eigv(:)

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    tfor=.FALSE.
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    iteropt%nfi=0
    bsnfi=0
    ! ==--------------------------------------------------------------==
    ! ==  PHASE FACTORS                                               ==
    ! ==--------------------------------------------------------------==
    CALL phfac(tau0)
    ! ==--------------------------------------------------------------==
    ! ==  RANDOM INITIALIZATION COEFFICIENTS FOR THE WAVEFUNCTIONS    ==
    ! ==--------------------------------------------------------------==
    cfac=dt2bye
    ! Loop over k points
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       nkpoint=nkpbl(ikpt)
       IF (tkpts%tkblock.AND.tsphere) CALL rkpt_swap(maskgw,1,ikpt,'MASKGW')
       DO ik=1,nkpoint
          ikk=kpbeg(ikpt)+ik
          CALL randtowf(c0(:,:,ik),nstate,ik,ikk)
          ! Orthogonalization
          IF (pslo_com%tivan) CALL rnlsm(c0(:,:,ik),nstate,&
               ikpt,ik,.FALSE.)
          CALL ortho(nstate,c0(:,:,ik),c2(:,:,ik))
       ENDDO
       IF (tkpts%tkblock) THEN
          CALL wkpt_swap(c0,nstate,ikpt,'C0')
          IF (tkpts%tknoswap) GOTO 100
       ENDIF
    ENDDO
100 CONTINUE
    IF (corel%tinlc) CALL copot(rhoe,psi,.FALSE.)
    ! ==--------------------------------------------------------------==
    ! ==  GRADIENT WITH RESPECT TO WAVEFUNCTIONS                      ==
    ! ==--------------------------------------------------------------==
    IF (.NOT.cntl%tmdbo) THEN   ! AK IF(.NOT.cntl%tdiag) THEN
       ALLOCATE(eigv(nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (tkpts%tkpnt) CALL stopgm('RANDWF',&
            'K POINTS AND NONORT NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       !vw avoiding segfaulting (assumed shape)
       CALL reshape_inplace(c2, (/SIZE(c2,1),SIZE(c2,2)/), c2_ptr)
       !CALL forcedr(c0(:,:,1),c2(:,:,1),sc0,rhoe,psi,tau0,fion,eigv,&
       CALL forcedr(c0(:,:,1),c2_ptr(:,:),sc0,rhoe,psi,tau0,fion,eigv,&
            nstate,1,.TRUE.,tfor)
       ! ==------------------------------------------------------------==
       ! ==  STEEPEST DESCENT STEP                                     ==
       ! ==------------------------------------------------------------==
       CALL daxpy(2*nstate*nkpt%ngwk*nkpoint,cfac,c2,1,c0,1)
       ! ==------------------------------------------------------------==
       ! ==  ORTHOGONALIZATION                                         ==
       ! ==------------------------------------------------------------==
       IF (.NOT.cntl%nonort) THEN
          DO ik=1,nkpoint
             IF (pslo_com%tivan) CALL rnlsm(c0(:,:,ik),nstate,&
                  1,ik,.FALSE.)
             CALL ortho(nstate,c0(:,:,ik),c2(:,:,ik))
          ENDDO
       ENDIF
       DEALLOCATE(eigv,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF                     ! End of .NOT.cntl%tmdbo
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE randwf
END MODULE rinitwf_driver
