MODULE vofrhot_utils
  USE eextern_utils,                   ONLY: eextern
  USE efld,                            ONLY: textfld
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE potfor_utils,                    ONLY: potfor
  USE ppener_utils,                    ONLY: ppener
  USE simulmod,                        ONLY: vploc
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE td_input,                        ONLY: td_prop
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vofrhot
  PUBLIC :: give_scr_vofrhot

CONTAINS

  ! ==================================================================
  SUBROUTINE vofrhot(tau0,fion,rhoe,v,vtemp,tfor)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ONE-PARTICLE POTENTIAL V IN REAL SPACE                  ==
    ! ==  THE TOTAL ENERGY ETOT                                       ==
    ! ==  THE FORCES FION ACTING ON THE IONS                          ==
    ! ==--------------------------------------------------------------==
    ! ==  RHOE: in  electronic density in real space                  ==
    ! ==  V:    out electronic density in G space                     ==
    ! ==--------------------------------------------------------------==
    ! HER[
    ! EHR]
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(fpar%nnr1)
    COMPLEX(real_8)                          :: v(:), vtemp(ncpw%nhg)
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrhot'

    COMPLEX(real_8)                          :: ee, eh, ei, eps
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:)
    INTEGER                                  :: ierr, ir, isub, nnrs

! Variables
! TODO refactor these arrays 
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! SCR partition
    ! ==--------------------------------------------------------------==
    ! Arrays inside SCR.
    ALLOCATE(eivps(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(eirop(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

    nnrs  = spar%nr1s*spar%nr2s*spar%nr3s
    CALL eicalc(eivps,eirop)
    ! ..External Field or forces of external potential
    ener_com%eext=0._real_8
    ! EHR[
    IF ((textfld.OR.cntl%texadd).AND.(.NOT.td_prop%td_extpot)) THEN
       ! EHR]
       CALL eextern(rhoe,v,eirop,tau0,fion,ener_com%eext,tfor)
       CALL mp_sum(ener_com%eext,parai%allgrp)
    ENDIF
    ! TRANSFORM THE DENSITY TO G SPACE
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir)=CMPLX(rhoe(ir),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! ==        MARTYNA & TUCKERMAN CLUSTER METHOD                    ==
    ! ==--------------------------------------------------------------==
    ! COMPUTE POTENTIAL ENERGY CONTRIBUTION TO FORCES ON IONS
    IF (tfor) CALL potfor(fion,v,eirop)
    ! COMPUTE ELECTROSTATIC AND PSEUDOPOTENTIAL ENERGIES.
    CALL ppener(eh,ei,ee,eps,v,vtemp,eivps,eirop)
    ! Hartree term for all charges
    ener_com%ehep = REAL(eh)*parm%omega
    ener_com%ehii = REAL(ei)*parm%omega + crge%nel*vploc
    ener_com%ehee = REAL(ee)*parm%omega
    IF (cntl%tdiag) THEN
       ! Ion-Ion interaction - e-e (double counting)
       ener_com%eht =  ener_com%ehii-ener_com%ehee&
            + ener_com%esr&
            - ener_com%eself/REAL(parai%nproc,kind=real_8)
    ELSE
       ener_com%eht =  ener_com%ehep&
            + ener_com%esr&
            - ener_com%eself/REAL(parai%nproc,kind=real_8)
    ENDIF
    ener_com%epseu=2._real_8*REAL(eps)*parm%omega
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vofrhot
  ! ==================================================================
  SUBROUTINE give_scr_vofrhot(lvofrhot,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrhot
    CHARACTER(len=30)                        :: tag

    lvofrhot = 2*ncpw%nhg*2
    tag      ='2*NHG*2'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrhot
  ! ==================================================================

END MODULE vofrhot_utils
