MODULE simple_model_p_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE system,                          ONLY: ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: simple_ortho_p

CONTAINS


  ! ==================================================================
  SUBROUTINE simple_ortho_p(nstate,c0,c1,S_star)
    ! ==================================================================
    ! omS_matrix is (1-S) matrix
    ! omS2_matrix is (1-S)**2 matrix
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate), &
                                                c1(nkpt%ngwk,nstate)
    REAL(real_8)                             :: S_star(nstate,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'simple_ortho_p'

    INTEGER                                  :: ierr, isub
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8), ALLOCATABLE, SAVE          :: overlap(:,:), star_ovl(:,:)

    CALL tiset('simp_orth_p',isub)

    IF (ifirst .EQ. 0) THEN
       ifirst = 1
       ALLOCATE(overlap(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(star_ovl(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(   overlap)!,nstate*nstate)
       CALL zeroing(star_ovl)!,nstate*nstate)
    ENDIF

    ! first calculate <c1|c0>
    CALL ovlap(nstate,overlap,c0,c1)
    CALL mp_sum(overlap,nstate*nstate,parai%allgrp)

    ! then state by state orthog... 
    ! zero order term
    CALL dgemm('N','N',2*ncpw%ngw,nstate,nstate,-1._real_8,c0,&
         2*ncpw%ngw,overlap,nstate,+1._real_8,c1,2*ncpw%ngw)

    ! first+second order correction in S matrix
    ! calculating (S_star)(overlap) matrix result in star_ovl
    CALL dgemm('N','N',nstate,nstate,nstate,1._real_8,S_star,&
         NSTATE,overlap,NSTATE,0._real_8,star_ovl,NSTATE)
    ! Now applying star_ovl to C0 to get the final projector
    ! and subtract that to C1
    CALL dgemm('N','N',2*ncpw%ngw,nstate,nstate,-1._real_8,c0,&
         2*ncpw%ngw,star_ovl,nstate,+1._real_8,c1,2*ncpw%ngw)

    CALL tihalt('simp_orth_p',isub)
    RETURN
  END SUBROUTINE simple_ortho_p

END MODULE simple_model_p_utils
