MODULE setbsstate_utils
  USE bsym,                            ONLY: bsclcs,&
                                             hdoel,&
                                             hsdown,&
                                             hspin,&
                                             hsup,&
                                             hupel
  USE elct,                            ONLY: crge
  USE hubbardu,                        ONLY: hubbu
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE ropt,                            ONLY: ropt_mod
  USE spin,                            ONLY: spin_mod,&
                                             tdsp1
  USE system,                          ONLY: cntl
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setbsstate

!TODO:  reimplement Hubbard U for Broken Symmetry 
!       This was implemented in cpmd 3.x (r1887 Bochum) 
!       contact: niklas.siemer@theochem.rub.de

CONTAINS

  ! ==================================================================
  SUBROUTINE setbsstate
    ! CB: Sets occupation in BS calc, dependening on BSCLCS (BS:1, HS:2)
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i

! ==--------------------------------------------------------------==
! ==--------------------------------------------------------------==

    IF (paral%parent) THEN
       IF (bsclcs.EQ.1) THEN
          spin_mod%nspin=0
          tdsp1%nupel=crge%n/2
          tdsp1%ndoel=crge%n/2
          spin_mod%nsup=crge%n/2
          spin_mod%nsdown=crge%n/2
       ELSE
          spin_mod%nspin=hspin
          tdsp1%nupel=hupel
          tdsp1%ndoel=hdoel
          spin_mod%nsup=hsup
          spin_mod%nsdown=hsdown
       ENDIF
       CALL zeroing(crge%f)!,n)
       DO i=1,INT(tdsp1%nupel)
          crge%f(i,1)=1._real_8
       ENDDO
       DO i=spin_mod%nsup+1,spin_mod%nsup+INT(tdsp1%ndoel)
          crge%f(i,1) = 1._real_8
       ENDDO
    ENDIF
    CALL mp_bcast(spin_mod%nspin,parai%source,parai%allgrp)
    CALL mp_bcast(tdsp1%nupel,parai%source,parai%allgrp)
    CALL mp_bcast(tdsp1%ndoel,parai%source,parai%allgrp)
    CALL mp_bcast(spin_mod%nsup,parai%source,parai%allgrp)
    CALL mp_bcast(spin_mod%nsdown,parai%source,parai%allgrp)
    CALL mp_bcast(crge%f,SIZE(crge%f),parai%source,parai%allgrp)
    IF(cntl%thubb) CALL mp_bcast(hubbu%u,hubbu%nuatm,parai%source,parai%allgrp)
    ropt_mod%sdiis=.TRUE.
    ropt_mod%spcg=.TRUE.
  END SUBROUTINE setbsstate

END MODULE setbsstate_utils
