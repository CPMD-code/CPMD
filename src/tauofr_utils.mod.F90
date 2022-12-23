MODULE tauofr_utils
  USE cp_grp_utils,                    ONLY: cp_grp_redist
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: gk,&
                                             indzs,&
                                             nzhs
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: invfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai
  USE part_1d,                         ONLY: part_1d_get_el_in_blk,&
                                             part_1d_nbr_el_in_blk
  USE pslo,                            ONLY: pslo_com
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm
  USE tauf,                            ONLY: itaur,&
                                             tau,&
                                             vtau
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tauofr
  PUBLIC :: dpsisc
  !public :: tauadd

CONTAINS

  ! ==================================================================
  SUBROUTINE tauofr(c0,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==          COMPUTES   TAU = SUM_i |nabla PSI_i|^2              ==
    ! ==--------------------------------------------------------------==
    ! == WARNING: ALL WAVEFUNCTIONS C0 HAVE TO BE ORTHOGONAL          ==
    !
    ! Added support for CP groups
    !                         17.10.2019 M.P. Bircher @ Comp-Phys/UniVie
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'tauofr'

    INTEGER                                  :: ierr, is, is1, is2, isub, isub2
    INTEGER, SAVE                            :: ifirst = 0

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF (pslo_com%tivan) CALL stopgm(procedureN,"USPP NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,"LSE NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (ifirst.EQ.0) THEN
       ifirst=1
       ALLOCATE(tau(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vtau(fpar%nnr1,clsd%nlsd),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    itaur=1
    ! ==--------------------------------------------------------------==
    ! Initialize
    CALL zeroing(tau)!,clsd%nlsd*nnr1)
    CALL zeroing(vtau) ! Since vtau is not deallocated between iterations, we zero it here
    ! Loop over the electronic states
    DO is = 1,part_1d_nbr_el_in_blk(nstate,parai%cp_inter_me,parai%cp_nogrp),2
       is1 = part_1d_get_el_in_blk(is,nstate,parai%cp_inter_me,parai%cp_nogrp)
       is2 = nstate+1
       IF (is+1.LE.part_1d_nbr_el_in_blk(nstate,parai%cp_inter_me,parai%cp_nogrp))&
            is2 = part_1d_get_el_in_blk(is+1,nstate,parai%cp_inter_me,parai%cp_nogrp)
       ! x direction
       CALL zeroing(psi)!,maxfft)
       CALL dpsisc(c0,psi,1,is1,is2,nstate)
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       CALL tauadd(psi,tau,is1,is2,nstate)
       ! y direction
       CALL zeroing(psi)!,maxfft)
       CALL dpsisc(c0,psi,2,is1,is2,nstate)
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       CALL tauadd(psi,tau,is1,is2,nstate)
       ! z direction
       CALL zeroing(psi)!,maxfft)
       CALL dpsisc(c0,psi,3,is1,is2,nstate)
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       CALL tauadd(psi,tau,is1,is2,nstate)
    ENDDO
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_grps_b',isub2)
       CALL cp_grp_redist(tau,fpar%nnr1,clsd%nlsd)
       CALL tihalt(procedureN//'_grps_b',isub2)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE tauofr
  ! ==================================================================
  SUBROUTINE dpsisc(c0,psi,k,is1,is2,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(maxfft)
    INTEGER                                  :: k, is1, is2, nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    INTEGER                                  :: ig

    IF (is2.GT.nstate) THEN
       !$omp parallel do private(IG) shared(k,is1)
       DO ig=1,ncpw%ngw
          psi(nzhs(ig))=gk(k,ig)*c0(ig,is1)
          psi(indzs(ig))=-gk(k,ig)*CONJG(c0(ig,is1))
       ENDDO
    ELSE
       !$omp parallel do private(IG) shared(k,is1,is2)
       DO ig=1,ncpw%ngw
          psi(nzhs(ig))=gk(k,ig)*(c0(ig,is1)+uimag*c0(ig,is2))
          psi(indzs(ig))=-gk(k,ig)*(CONJG(c0(ig,is1))+&
               uimag*CONJG(c0(ig,is2)))
       ENDDO
    ENDIF
    IF (geq0) psi(nzhs(1))=CMPLX(0._real_8,0._real_8,kind=real_8)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE dpsisc
  ! ==================================================================
  SUBROUTINE tauadd(psi,tau,is1,is2,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psi(maxfft)
    REAL(real_8)                             :: tau(fpar%nnr1,clsd%nlsd)
    INTEGER                                  :: is1, is2, nstate

    INTEGER                                  :: ir, ispin1, ispin2
    REAL(real_8)                             :: coef1, coef2, r1, r2

    coef1=0.5_real_8*parm%tpiba2*crge%f(is1,1)/parm%omega
    IF (is2.GT.nstate) THEN
       coef2=0._real_8
    ELSE
       coef2=0.5_real_8*parm%tpiba2*crge%f(is2,1)/parm%omega
    ENDIF
    IF (cntl%tlsd) THEN
       ispin1=1
       ispin2=1
       IF (is1.GT.spin_mod%nsup) ispin1=2
       IF (is2.GT.spin_mod%nsup) ispin2=2
       !$omp parallel do private(IR,R1,R2)
       DO ir=1,fpar%nnr1
          r1=AIMAG(psi(ir))
          r2=REAL(psi(ir))
          tau(ir,ispin1)=tau(ir,ispin1)+coef1*r1*r1
          tau(ir,ispin2)=tau(ir,ispin2)+coef2*r2*r2
       ENDDO
    ELSE
       !$omp parallel do private(IR,R1,R2)
       DO ir=1,fpar%nnr1
          r1=AIMAG(psi(ir))
          r2=REAL(psi(ir))
          tau(ir,1)=tau(ir,1)+coef1*r1*r1+coef2*r2*r2
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE tauadd
  ! ==================================================================

END MODULE tauofr_utils
