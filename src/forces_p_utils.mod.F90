MODULE forces_p_utils
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE h0psi1_p_utils,                  ONLY: give_scr_h0psi1,&
                                             h0psi1_p
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE norm,                            ONLY: gnmax,&
                                             gnorm
  USE parac,                           ONLY: parai
  USE perturbation_p_utils,            ONLY: energy
  USE response_pmod,                   ONLY: dmbi,&
                                             ener1,&
                                             nmr_options,&
                                             response1
  USE rhoofr_p_utils,                  ONLY: give_scr_rhoofr_p
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE v1ofrho_p_utils,                 ONLY: give_scr_v1ofrho,&
                                             v1ofrho_p
  USE vpsi_p_utils,                    ONLY: v1psi0_p
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: forces_p
  PUBLIC :: give_scr_forces_p

CONTAINS

  ! ==================================================================
  SUBROUTINE forces_p(c0,c1,c12,h1nl,z11,drhoe,v1_loc,eirop1,&
       tau0,fion,rhoe,psi,nstate,firstcall)
    ! =================================================================
    ! This routine calculates the gradient with respect to the first
    ! order perturbation wavefunctions psi_1 (C1) of the second order
    ! energy functional to be optimized:
    ! 
    ! etot = sum_a<c1_a|h0|c1_a>-
    ! sum_a,b <c1_a|z11_b,a|c1_b>+
    ! sum_a <c0|h1|c1>+<c1|h1|c0>+
    ! 1/2 int rr' {n1(r) n1(r')  [ 1/|r-r'|  +  dmu_xc/dn d(r-r')]}
    ! 
    ! xavier gonze  PRA 52 p1096 (1995) formula (91)
    ! =================================================================
    ! The resulting force is NOT orthogonalized with respect to
    ! the ground state orbitals. The maximum component is NOT calculated.
    ! 
    ! =================================================================
    REAL(real_8)                             :: drhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: v1_loc(ncpw%nhg,clsd%nlsd), &
                                                eirop1(ncpw%nhg)
    REAL(real_8)                             :: tau0(:,:,:), fion(:,:,:), &
                                                rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: z11(nstate,nstate)
    COMPLEX(real_8) :: h1nl(nkpt%ngwk,nstate), c12(nkpt%ngwk,nstate), &
      c1(nkpt%ngwk,nstate), c0(nkpt%ngwk,nstate)
    LOGICAL                                  :: firstcall

    INTEGER                                  :: is1, isub

! variables
! TESTING, in progress...
! ==--------------------------------------------------------------==

    CALL tiset('  forces_p',isub)
    gnmax  = 0.0_real_8
    gnorm  = 0.0_real_8
    ener1%eh0    = 0.0_real_8
    ener1%elag1  = 0.0_real_8
    ener1%exc1   = 0.0_real_8
    ener1%enl1   = 0.0_real_8
    ! ==--  |C12>  =  H0  |C1>   -------------------------------------==
    IF (.NOT. firstcall) THEN ! because IF FIRSTCALL, then c1=0 anyway
       CALL h0psi1_p(c1,c0,drhoe,crge%f,tau0,fion,psi,&
            nstate,c12)
       ener1%eh0 = - energy(c1,c12,nstate)
    ELSE
       IF (response1%trho1) CALL zeroing(drhoe)!,nnr1*clsd%nlsd)
    ENDIF                     ! if not firstcall
    ! This routine also computes the FIRST ORDER DENSITY (drhoe)

    ! ==--  Lagrange multipliers  ------------------------------------==
    ! ==--  |C12>  +=  eps_ij  |C1>  ---------------------------------==
    ! elag1  =  sum_a,b  <c1_a|  z11_b,a  |c1_b>
    IF (.NOT. firstcall) THEN ! because IF FIRSTCALL, then c1=0 anyway
       IF (((response1%tnmr.OR.response1%tepr) .AND. nmr_options%tlocalize) .OR. (response1%tinteraction)) THEN
          ! call DSYMM ( 'R','U', 2*ngw, nstate, 1._real_8,
          ! &        z11, nstate, c1, 2*ngw, 1._real_8, c12, 2*ngw)
          CALL rotate(1._real_8, c1, 1._real_8, c12, z11, nstate, 2*ncpw%ngw,&
               cntl%tlsd, spin_mod%nsup, spin_mod%nsdown)
       ELSE
          DO is1=1,nstate
             CALL daxpy(2*ncpw%ngw, z11(is1,is1), c1(1,is1),1, c12(1,is1),1)
          ENDDO
       ENDIF
       ener1%elag1 = - energy(c1,c12,nstate) - ener1%eh0
    ENDIF                     ! if not firstcall

    ! ==--  The perturbation potential  ------------------------------==
    ! ==--  external (v1_loc) and internal  v[n^1]  ------------------==
    ! ==--  |C12>  +=  v^1 |C0>  -------------------------------------==
    IF (response1%trho1) THEN
       IF (dmbi%fastdexc .AND. cntl%tgc) THEN
          ! calculate d2Exc/dn1dn1 using LDA...
          cntl%tgc=.FALSE.
          CALL v1ofrho_p(v1_loc,eirop1,drhoe,rhoe,psi,&
               firstcall)
          cntl%tgc=.TRUE.
       ELSE
          CALL v1ofrho_p(v1_loc,eirop1,drhoe,rhoe,psi,&
               firstcall)
       ENDIF

       ! \ returns the total potential in rhoe
       CALL v1psi0_p(c0,c12,rhoe,psi,nstate)
       ! Calculates also the following energies:
       ! \  EHT1  =    <1|  v_Coul[n^(1)]  |0>  the implicit Coulomb potential 
       ! \                                      of rho1
       ! \  ELOC1 =    <1|  v_ext |0>           the explicit external potential
       ! \  EXC1  =    <1|  v_xc [n^(1)]  |0>   the implicit xc-potential of rho1
    ENDIF

    ! ==--  The nonlocal part   --------------------------------------==
    ! ==--  |C12>  +=  |H1NL>   --------------------------------------==
    IF (response1%tnonlocal) THEN
       CALL daxpy(2*nstate*ncpw%ngw,1._real_8,h1nl(1,1),1,c12(1,1),1)
       ener1%enl1 = - 2._real_8 * energy(c1,h1nl,nstate)
    ENDIF
    ! ==--  ENERGIES  ------------------------------------------------==
    ! gathering energies from all procs
    CALL mp_sum(ener1%eloc1,parai%allgrp)
    CALL mp_sum(ener1%enl1,parai%allgrp)
    CALL mp_sum(ener1%eht1,parai%allgrp)
    CALL mp_sum(ener1%exc1,parai%allgrp)
    CALL mp_sum(ener1%elag1,parai%allgrp)
    CALL mp_sum(ener1%eh0,parai%allgrp)

    ener_com%etot =  ener1%eh0 + ener1%eloc1 + ener1%eht1 + ener1%enl1 + ener1%exc1 + ener1%elag1
    ! ==--------------------------------------------------------------==
    CALL tihalt('  forces_p',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE forces_p
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE give_scr_forces_p(l_forces,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: l_forces
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l_h0psi1, l_rhoofr_p, &
                                                l_v1ofrho

! variables
! ==--------------------------------------------------------------==

    CALL give_scr_h0psi1(l_h0psi1,tag,nstate)

    IF (response1%phonon .OR. response1%tlanphon .OR. response1%teigens .OR. cntl%tinr&
         .OR. response1%teigensystem .OR. response1%thardness .OR. response1%traman .OR. response1%tfukui&
         .OR. response1%tdummyatom .OR. response1%tinteraction .OR. response1%tvoa) THEN
       CALL give_scr_rhoofr_p(l_rhoofr_p,tag)
       CALL give_scr_v1ofrho(l_v1ofrho,tag)
    ELSE
       l_rhoofr_p=0
       l_v1ofrho=0
    ENDIF


    l_forces = MAX(l_h0psi1,l_rhoofr_p,l_v1ofrho)
    tag     = 'perturbation forces'
    RETURN
  END SUBROUTINE give_scr_forces_p
  ! ==================================================================

END MODULE forces_p_utils
