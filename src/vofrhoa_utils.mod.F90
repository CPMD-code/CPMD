#include "cpmd_global.h"

MODULE vofrhoa_utils
  USE eextern_utils,                   ONLY: eextern
  USE efld,                            ONLY: textfld
  USE eicalc_utils,                    ONLY: eicalc
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: fwfftn
  USE htrstr_utils,                    ONLY: htrstr
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nvtx_utils
  USE parac,                           ONLY: parai
  USE potfor_utils,                    ONLY: potfor
  USE ppener_utils,                    ONLY: ppener
  USE simulmod,                        ONLY: vploc
  USE state_utils,                     ONLY: copy_to_re
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE td_input,                        ONLY: td_prop
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vlocst_utils,                    ONLY: vlocst
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vofrhoa
  PUBLIC :: give_scr_vofrhoa

CONTAINS

  ! ==================================================================
  SUBROUTINE vofrhoa(tau0,fion,rhoe,v,vtemp,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ONE-PARTICLE POTENTIAL V IN REAL SPACE                  ==
    ! ==  THE TOTAL ENERGY ETOT                                       ==
    ! ==  THE FORCES FION ACTING ON THE IONS                          ==
    ! ==--------------------------------------------------------------==
    ! ==  RHOE: in  electronic density in real space                  ==
    ! ==  V:    out electronic density in G space                     ==
    ! ==--------------------------------------------------------------==
    ! EHR[
    ! EHR]
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(:)
    COMPLEX(real_8)                          :: v(:), vtemp(:)
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrhoa'

    COMPLEX(real_8)                          :: ee, eh, ei, eps
    COMPLEX(real_8), ALLOCATABLE             :: eirop(:), eivps(:)
    INTEGER                                  :: ierr, isub, nnrs
    REAL(real_8)                             :: ehs

! Variables
! these arrays are in modules now 
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

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
    CALL zeroing( v )
    CALL copy_to_re ( fpar%nnr1 ,rhoe , v )
!!$omp parallel do private(IR)
    !DO ir=1,fpar%nnr1
    !   v(ir)=CMPLX(rhoe(ir),0._real_8,kind=real_8)
    !ENDDO
    CALL  fwfftn(v,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! ==                       PERIODIC SYSTEMS                       ==
    ! ==--------------------------------------------------------------==
    ! COMPUTE POTENTIAL ENERGY CONTRIBUTION TO FORCES ON IONS
    IF (tfor) CALL potfor(fion,v,eirop)
    ! COMPUTE ELECTROSTATIC AND PSEUDOPOTENTIAL ENERGIES.
    CALL ppener(eh,ei,ee,eps,v(:),vtemp,eivps,eirop)
    ! Hartree term for all charges
    ener_com%ehep = REAL(eh)*parm%omega
    ! Ion-Ion interaction (T.D. add VPLOC)
    ! NEL*VPLOC is the G=0 part of the Hartree energy, for compatibility
    ! reasons we add it to the band energy.
    ! (Thierry Deutsch) to have independent term with Raggio
    ! add NEL*VPLOC with Hartree energy(IGEQ0 processor has VPLOC!=0)
    ener_com%ehii = REAL(ei)*parm%omega + crge%nel*vploc
    ! Electron-Electron interaction
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
    ! Local pseudopotential energies
    ener_com%epseu=2._real_8*REAL(eps)*parm%omega
    IF (tstress) THEN
       ! Hartree contribution to stress
       ehs=REAL(eh)*parm%omega
       CALL htrstr(ehs,v,eirop)
       ! Local pseudopotential contribution to stress
       CALL vlocst(ener_com%epseu,v,eivps)
    ENDIF

    DEALLOCATE(eivps,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(eirop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE vofrhoa
  ! ==================================================================
  SUBROUTINE give_scr_vofrhoa(lvofrhoa,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrhoa
    CHARACTER(len=30)                        :: tag

    lvofrhoa = 2*ncpw%nhg*2 + 200
    tag      ='2*NHG*2'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrhoa
  ! ==================================================================

END MODULE vofrhoa_utils
