MODULE dynit_utils
  USE clas,                            ONLY: clas3,&
                                             clas4,&
                                             tclas
  USE cnst,                            ONLY: scmass
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntr
  USE tpar,                            ONLY: &
       dt2_elec, dt2_ions, dt2bye, dt2bym, dt2hbe, dt_elec, dt_ions, dtb2me, &
       dtb2mi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dynit

CONTAINS

  ! ==================================================================
  SUBROUTINE dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ! ==--------------------------------------------------------------==
    ! == INITIALIZATION OF STEPS FOR GEOMETRY OPTIMIZATION OR         ==
    ! == MOLECULAR DYNAMICS (quantities inside tpar.inc)              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekincp, ekin1, ekin2, temp1, &
                                                temp2, ekinh1, ekinh2

    INTEGER                                  :: is

! Variables
! ==================================================================
! ==  BOUNDS FOR ELECTRONIC AND IONIC DYNAMICS                    ==
! ==--------------------------------------------------------------==

    ekincp=0._real_8
    ekin1=cntr%ekinw+cntr%toll
    ekin2=cntr%ekinw-cntr%toll
    temp1=cntr%tempw+cntr%tolp
    temp2=cntr%tempw-cntr%tolp
    ekinh1=cntr%ekinhr+cntr%tolc
    ekinh2=cntr%ekinhr-cntr%tolc
    ! ==================================================================
    ! ==  TIME STEP FUNCTIONS                                         ==
    ! ==--------------------------------------------------------------==
    ! For electrons
    dt_elec=cntr%delt_elec
    dt2_elec=cntr%delt_elec*cntr%delt_elec
    dt2bye=dt2_elec/cntr%emass
    dt2hbe=0.5_real_8*dt2_elec/cntr%emass
    dtb2me=cntr%delt_elec/(2.0_real_8*cntr%emass)
    ! For ions
    dt_ions=cntr%delt_ions
    dt2_ions=cntr%delt_ions*cntr%delt_ions
    DO is=1,ions1%nsp
       dt2bym(is)=dt2_ions/rmass%pma(is)
    ENDDO
    DO is=1,ions1%nsp
       dtb2mi(is)=cntr%delt_ions/(2.0_real_8*rmass%pma(is))
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (tclas) THEN
       DO is=1,clas3%ncltyp
          clas4%cstep(is)=cntr%delt_ions/(2.0_real_8*clas4%clmas(is)*scmass)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dynit
  ! ==================================================================

END MODULE dynit_utils
